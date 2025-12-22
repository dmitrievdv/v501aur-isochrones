using Makie
using CairoMakie, GLMakie
using DataFrames, CSV
using Polynomials
using LaTeXStrings
using Printf
using FITSIO
using Dierckx
using Optim
using Interpolations

include("kurucz-int.jl")

star_masses = [1:0.15:7; 3.0:0.05:5.0;]
sort!(star_masses)
mesa_dir = "mesa-data"
mesa_processed_dir = "mesa-data-processed"
mesa_processed_temp = "mesa-data-processed-temp"
mesa_processed_tess_dir = "mesa-data-processed-tess"

mkpath(mesa_dir)
mkpath(mesa_processed_dir)
mkpath(mesa_processed_temp)
mkpath(mesa_processed_tess_dir)

function process_mesa_history_data(dir)
    tracks_dir = "$dir/mesa_tracks"
    models = readdir(tracks_dir)
    for model in models
        println(model)
        mass, z, rot = [parse(Float64, m.match) for m in eachmatch(r"[0-9]+.[0-9]+", model)]
        output_file = @sprintf "ROT%.2fM%.2fZ%.4f.csv" rot mass z

        model_dir = "$tracks_dir/$model"
        suffixes = [split(dir, '_')[2] for dir in readdir(model_dir)]
        model_df = foldr(vcat, [get_mesa_df("$model_dir/LOGS_$suffix/history.data") for suffix in suffixes])
        sort!(model_df, [:model_number])
        CSV.write("$mesa_processed_temp/$output_file", model_df)
    end
end

function get_mesa_df(mesa_file :: AbstractString)
    lines = readlines(mesa_file)
    splitted_lines = [split(line) for line in lines]
    col_names = splitted_lines[6]
    n_cols = length(col_names)
    cols = [[parse(Float64, splitted_line[i_col]) for splitted_line in splitted_lines[7:end]] for i_col in 1:n_cols]
    mesa_df = DataFrame([Pair(name_col_tuple...) for name_col_tuple in zip(col_names, cols)])
end

function get_mesa_df(star_mass)
    mesa_file_name = "track_M$(star_mass)_Z0.0147_ROT0.00_history.data"

    lines = readlines("$mesa_dir/$mesa_file_name")
    splitted_lines = [split(line) for line in lines]
    col_names = splitted_lines[6]
    n_cols = length(col_names)
    cols = [[parse(Float64, splitted_line[i_col]) for splitted_line in splitted_lines[7:end]] for i_col in 1:n_cols]
    mesa_df = DataFrame([Pair(name_col_tuple...) for name_col_tuple in zip(col_names, cols)])
end

function get_mesa_df_processed(star_mass :: Real)
    star_mass_str = @sprintf "%.2f" star_mass
    mesa_file_name = "ROT0.00M$(star_mass_str)Z0.0147.csv"

    get_mesa_df_processed(mesa_file_name)
end

# function get_mesa_df_processed(mesa_file_name)
#     CSV.read("$mesa_processed_dir/$mesa_file_name", DataFrame)
# end

function get_mesa_df_processed(mesa_file_name :: AbstractString, mesa_dir = mesa_processed_dir :: AbstractString)
    CSV.read("$mesa_dir/$mesa_file_name", DataFrame)
end

function get_mesa_df_processed(mesa_file_name, tess_response_spl, Av, mesa_dir = mesa_processed_dir :: AbstractString)
    mesa_df = get_mesa_df_processed(mesa_file_name, mesa_dir)

    n_mesa = nrow(mesa_df)

    log_tess_luminosity = zeros(n_mesa)
    log_tess_luminosity_noext = zeros(n_mesa)

    for i_mesa = 1:n_mesa
        log_tess_luminosity[i_mesa] = try
            log10(mesa_df.radius[i_mesa]^2 * calc_tess_kurucz_flux_with_extinction(mesa_df.Teff[i_mesa], mesa_df.log_g[i_mesa], Av, tess_response_spl))
        catch 
            -1e2
        end

        log_tess_luminosity_noext[i_mesa] = try
            log10(mesa_df.radius[i_mesa]^2 * calc_tess_kurucz_flux_no_extinction(mesa_df.Teff[i_mesa], mesa_df.log_g[i_mesa], tess_response_spl))
        catch 
            -1e2
        end
    end
    # println(log_tess_luminosity)

    mesa_df[!, :log_TESS] = log_tess_luminosity
    mesa_df[!, :log_TESS_noext] = log_tess_luminosity_noext

    mkpath("$mesa_processed_dir-tess")
    CSV.write("$mesa_processed_dir-tess/$mesa_file_name", mesa_df)

    mesa_df
end

function get_mesa_df_processed(star_mass :: Real, tess_response_spl, Av)
    star_mass_str = @sprintf "%.2f" star_mass
    mesa_file_name = "ROT0.00M$(star_mass_str)Z0.0147.csv"

    get_mesa_df_processed(mesa_file_name, tess_response_spl, Av)
end

function get_age_df(mesa_df, lg_age)
    lg_ages = log10.(mesa_df.star_age)
    i_age = findfirst(x -> x > lg_age, lg_ages)
    # println(i_age)

    if isnothing(i_age)
        df = DataFrame([Pair(name_val_tuple...) for name_val_tuple in zip(names(mesa_df), collect(mesa_df[end,:]))])
        df.model_number .= -1.0
        df.star_age .= 10^lg_age
        return df
    end

    vals_prev_age = collect(mesa_df[i_age-1,:])
    vals_next_age = collect(mesa_df[i_age,:])

    prev_age = mesa_df.star_age[i_age-1]
    next_age = mesa_df.star_age[i_age]

    vals_age = vals_prev_age + (10^lg_age - prev_age)/(next_age - prev_age)*(vals_next_age - vals_prev_age)
    return DataFrame([Pair(name_val_tuple...) for name_val_tuple in zip(names(mesa_df), vals_age)])
end

function get_isochrone_df(mesa_dfs, lg_age)
    foldl(vcat, get_age_df.(mesa_dfs, lg_age))
end

function interpolate_isochrone_df(isochrone_df, M)
    isochrone_masses = isochrone_df.star_mass[isochrone_df.model_number .> 0]
    min_mass = minimum(isochrone_masses)
    max_mass = maximum(isochrone_masses)
    if ((M > max_mass) | (M < min_mass))
        df = DataFrame([Pair(name_val_tuple...) for name_val_tuple in zip(names(isochrone_df), collect(isochrone_df[end,:]))])
        df[:,:] .= NaN
        df.model_number .= -1.0
        return df
    end

    i_mass = findfirst(x -> x > M, isochrone_df.star_mass)
    vals_prev = collect(isochrone_df[i_mass-1,:])
    vals_next = collect(isochrone_df[i_mass,:])

    mass_prev = isochrone_df.star_mass[i_mass-1]
    mass_next = isochrone_df.star_mass[i_mass]

    vals = vals_prev + (M - mass_prev)/(mass_next - mass_prev)*(vals_next - vals_prev)
    return DataFrame([Pair(name_val_tuple...) for name_val_tuple in zip(names(isochrone_df), vals)])
end

function interpolate_isochrone_df(isochrone_df, M)
    isochrone_masses = isochrone_df.star_mass[isochrone_df.model_number .> 0]
    min_mass = minimum(isochrone_masses)
    max_mass = maximum(isochrone_masses)
    if ((M > max_mass) | (M < min_mass))
        df = DataFrame([Pair(name_val_tuple...) for name_val_tuple in zip(names(isochrone_df), collect(isochrone_df[end,:]))])
        df[:,:] .= NaN
        df.model_number .= -1.0
        return df
    end

    i_mass = findfirst(x -> x > M, isochrone_df.star_mass)
    vals_prev = collect(isochrone_df[i_mass-1,:])
    vals_next = collect(isochrone_df[i_mass,:])

    mass_prev = isochrone_df.star_mass[i_mass-1]
    mass_next = isochrone_df.star_mass[i_mass]

    vals = vals_prev + (M - mass_prev)/(mass_next - mass_prev)*(vals_next - vals_prev)
    return DataFrame([Pair(name_val_tuple...) for name_val_tuple in zip(names(isochrone_df), vals)])
end

function interpolate_isochrone_lg_L(isochrone_df, M)
    isochrone_masses = isochrone_df.star_mass[isochrone_df.model_number .> 0]
    min_mass = minimum(isochrone_masses)
    max_mass = maximum(isochrone_masses)
    if ((M > max_mass) | (M < min_mass))
        return NaN
    end

    i_mass = findfirst(x -> x > M, isochrone_df.star_mass)
    lg_L_prev = isochrone_df.log_L[i_mass-1]
    lg_L_next = isochrone_df.log_L[i_mass]

    mass_prev = isochrone_df.star_mass[i_mass-1]
    mass_next = isochrone_df.star_mass[i_mass]

    return lg_L_prev + (M - mass_prev)/(mass_next - mass_prev)*(lg_L_next - lg_L_prev)
end

function find_secondary_mass(mass_function, primary_mass)
    mass_poly = Polynomial([-mass_function*primary_mass^2, -2mass_function*primary_mass, -mass_function, 1])
    mass_poly_roots = Polynomials.roots(mass_poly)
    return real.(mass_poly_roots[abs.(imag.(mass_poly_roots)) .≈ 0.0])[1]
end

function fit_max_mass(max_age, max_mass, poly_deg)

    lg_age_to_fit = log10.(max_age[1e7 .< max_age .< 4.1e8])
    mass_to_fit = max_mass[1e7 .< max_age .< 4.1e8]

    function f_to_fit(x)
        f_vals = Polynomial(x).(lg_age_to_fit)
        sum((log10.(mass_to_fit) - f_vals) .^ 2)
    end


    init_val = zeros(poly_deg+1)

    res = optimize(f_to_fit, init_val, LBFGS())
    res.minimizer
end


# def __init__(self):
#         self.mass_c = 0.079
#         self.sigma_lm = 0.69
#         self.k_kroupa = 4.53
#         self.a_coef = 0.2791

#     def prob(self, mass):
#         if 0.1 <= mass <= 1.0:
#             return self.k_kroupa * 1.0 / mass * np.exp(
#                 -np.power(np.log10(mass) - np.log10(self.mass_c), 2.0) / 2.0 / np.power(self.sigma_lm, 2.0))
#         elif 1.0 < mass <= 150:
#             return self.k_kroupa * self.a_coef * np.power(mass, -2.35)
#         else:
#             return None

function IMF(m)
    mass_c = 0.079
    σ_lm = 0.69
    k_kroupa = 4.53/9.954211189
    a_coef = 0.2791
    if 0.1 ≤ m ≤ 1.0
        k_kroupa / m * exp(-(m - mass_c)^2 / 2 / σ_lm^2)
    elseif 1.0 < m ≤ 150
        k_kroupa * a_coef / m^(2.35)
    else
        0.0
    end
end

age_1 = 1e8
age_2 = 4e8

function prior_μ(lg_age, μ_1, μ_2, poly_max_fit)
    m = calc_max_mass(lg_age, poly_max_fit)
    m^2/(age_2 - age_1)*IMF(m*μ_1)*IMF(m*μ_2)*log(10)*10^lg_age
end

function extract_important_data(mesa_df, poly_max_fit)
    return hcat(log10.(mesa_df.star_age), mesa_df.star_mass ./ calc_max_mass.(log10.(mesa_df.star_age), Ref(poly_max_fit)), mesa_df.log_TESS)
end

function extract_important_data_isochrone(log_age, mesa_dfs, poly_max_fit)
    isochrone = fill(1e4, length(mesa_dfs)+1, 2)
    last_points = zeros(length(mesa_dfs)+1, 3)
    i_last = 0
    for (i_mesa, mesa_df) in enumerate(mesa_dfs)
        data = extract_important_data(mesa_df, poly_max_fit)
        last_points[end-i_mesa+1, :] = data[end, :]
        i_first = findfirst(x -> x > log_age, data[:,1])
        if isnothing(i_first)
            isochrone[i_mesa, 1] = 1e4
            isochrone[i_mesa, 2] = 1e4
            continue
        end
        i_last += 1

        n_data = size(data)[1]

        itp_indeces = if i_first < 6
            1:10
        elseif i_first > size(data)[1] - 4
            n_data-9:n_data
        else 
            i_first-5:i_first+4
        end

        itp_mass = interpolate(data[itp_indeces, 1], data[itp_indeces, 2], AkimaMonotonicInterpolation())
        itp_flux = interpolate(data[itp_indeces, 1], data[itp_indeces, 3], AkimaMonotonicInterpolation())

        isochrone[i_mesa,1] = itp_mass(log_age)
        isochrone[i_mesa,2] = itp_flux(log_age)
    end

    last_points_sorted = sortslices(last_points, dims = 1)

    i_first = findfirst(x -> x > log_age, last_points_sorted[:,1])
    itp_indeces = if i_first < 6
        1:10
    elseif i_first > length(mesa_dfs) - 4
        n_data-9:n_data
    else 
        i_first-5:i_first+4
    end

    # sorted_itp = sortperm(last_points[itp_indeces, 1])

    itp_ages_points = last_points_sorted[itp_indeces,1]
    itp_mass_points = last_points_sorted[itp_indeces,2]
    itp_flux_points = last_points_sorted[itp_indeces,3]


    itp_mass = interpolate(itp_ages_points, itp_mass_points, AkimaMonotonicInterpolation())
    itp_flux = interpolate(itp_ages_points, itp_flux_points, AkimaMonotonicInterpolation())

    isochrone[i_last+1,1] = itp_mass(log_age)
    isochrone[i_last+1,2] = itp_flux(log_age)

    
    return isochrone
end

function find_tess_flux(lg_age, iso_mass :: Real, mesa_dfs, poly_max_fit)
    isochrone = extract_important_data_isochrone(lg_age, mesa_dfs, poly_max_fit)
    # sort_is = sortperm(isochrone[:,1])
    itp_flux = interpolate(isochrone[:,1], isochrone[:,2], AkimaMonotonicInterpolation())
    itp_flux(iso_mass)
end

function find_tess_flux(lg_age, iso_mass, mesa_dfs, poly_max_fit)
    isochrone = extract_important_data_isochrone(lg_age, mesa_dfs, poly_max_fit)
    sort_is = sortperm(isochrone[:,1])
    itp_flux = interpolate(isochrone[sort_is,1], isochrone[sort_is,2], AkimaMonotonicInterpolation())
    min_mass = isochrone[1,1]
    # println(min_mass * calc_max_mass(lg_age, poly_max_fit))
    map(iso_mass) do μ
        if (μ < min_mass) | (μ > 0.995)
            return 1e4
        else
            itp_flux(μ)
        end
    end
end

function calc_iso_mass_function(lg_age, μ_1, μ_2, poly_max_fit)
    max_mass = calc_max_mass(lg_age, poly_max_fit)
    m_1 = max_mass*μ_1; m_2 = max_mass*μ_2
    calc_mass_function(m_1, m_2)
end

calc_mass_function(m_1, m_2) = m_2^3/(m_1 + m_2)^2

function calc_likelihood(lg_age, μs_1, μs_2, lg_flux_rel, lg_flux_err, mass_function, mass_function_err, mesa_dfs, poly_max_fit) # assuming same grid for μ_1, μ_2
    fluxes = find_tess_flux(lg_age, [μs_1; μs_2], mesa_dfs, poly_max_fit)
    n_μ_1 = length(μs_1)
    fluxes_1 = fluxes[1:n_μ_1]
    n_μ_2 = length(μs_2)
    fluxes_2 = fluxes[n_μ_1+1:end]
    # println(fluxes)
    # n_μ = length(μs)
    likelihood = Matrix{Float64}(undef, n_μ_1, n_μ_2)
    max_mass = calc_max_mass(lg_age, poly_max_fit)
    for (i_μ_1, μ_1) in enumerate(μs_1)
        if fluxes[i_μ_1] > 9e3
            likelihood[i_μ_1, :] .= 0.0
            continue
        end
        m_1 = max_mass*μ_1
        for (i_μ_2, μ_2) in enumerate(μs_2)
            if fluxes[i_μ_2] > 9e3
                likelihood[i_μ_1, i_μ_2] = 0.0
                continue
            end
            model_lg_flux_rel = fluxes_1[i_μ_1] - fluxes_2[i_μ_2]
            m_2 = max_mass*μ_2
            model_mass_function = calc_mass_function(m_1, m_2)
            # println("$m_1 $m_2 $model_mass_function $model_lg_flux_rel")
            # println(prior_μ(lg_age, μ_1, μ_2, poly_max_fit))
            # println(exp(-(model_lg_flux_rel - lg_flux_rel)^2/2lg_flux_err^2))
            # println(exp(-(model_mass_function - mass_function)^2/2mass_function_err^2))
            likelihood[i_μ_1, i_μ_2] = prior_μ(lg_age, μ_1, μ_2, poly_max_fit)
            likelihood[i_μ_1, i_μ_2] *= exp(-(model_lg_flux_rel - lg_flux_rel)^2/2lg_flux_err^2)
            likelihood[i_μ_1, i_μ_2] *= exp(-(model_mass_function - mass_function)^2/2mass_function_err^2)
            if isnan(likelihood[i_μ_1, i_μ_2])
                likelihood[i_μ_1, i_μ_2] = 0.0
            end
        end
    end

    return likelihood
end


function calc_likelihood_mass_function(lg_age, μs_1, mass_functions, mass_matrix, lg_flux_rel, lg_flux_err, mass_function, mass_function_err, mesa_dfs, poly_max_fit) # assuming same grid for μ_1, μ_2
    n_μ_1 = length(μs_1)
    n_f = length(mass_functions)
    
    max_mass = calc_max_mass(lg_age, poly_max_fit)
    mass_matrix /= max_mass
    # mass_matrix = zeros(n_μ_1, n_f+1)
    # mass_matrix[:,1] .= μs_1

    # for (i_μ_1, μ_1) in enumerate(μs_1)
    #     for (i_f, f) in enumerate(mass_functions) 
    #         m_1 = max_mass*μ_1
    #         m_2 = find_secondary_mass(f, m_1)
    #         # println("$f $m_1 $m_2")
    #         mass_matrix[i_μ_1, i_f+1] = m_2/max_mass
    #     end
    # end
    

    fluxes = find_tess_flux(lg_age, mass_matrix, mesa_dfs, poly_max_fit)
    
    fluxes_1 = fluxes[:,1]
    # n_μ_2 = length(μs_2)
    fluxes_2 = fluxes[:, 2:end]
    # println(fluxes)
    # n_μ = length(μs)
    likelihood = Matrix{Float64}(undef, n_μ_1, n_f)
    
    for (i_μ_1, μ_1) in enumerate(μs_1)
        if fluxes_1[i_μ_1,1] > 9e3
            likelihood[i_μ_1, :] .= 0.0
            continue
        end
        for (i_f, f) in enumerate(mass_functions)
            if fluxes_2[i_μ_1, i_f] > 9e3
                likelihood[i_μ_1, i_f] = 0.0
                continue
            end
            μ_2 = mass_matrix[i_μ_1, i_f + 1]
            model_lg_flux_rel = fluxes_1[i_μ_1] - fluxes_2[i_μ_1, i_f]
            model_mass_function = f
            # println("$m_1 $m_2 $model_mass_function $model_lg_flux_rel")
            # println(prior_μ(lg_age, μ_1, μ_2, poly_max_fit))
            # println(exp(-(model_lg_flux_rel - lg_flux_rel)^2/2lg_flux_err^2))
            # println(exp(-(model_mass_function - mass_function)^2/2mass_function_err^2))
            likelihood[i_μ_1, i_f] = prior_μ(lg_age, μ_1, μ_2, poly_max_fit)
            likelihood[i_μ_1, i_f] *= exp(-(model_lg_flux_rel - lg_flux_rel)^2/2lg_flux_err^2)
            likelihood[i_μ_1, i_f] *= exp(-(model_mass_function - mass_function)^2/2mass_function_err^2)
            likelihood[i_μ_1, i_f] /= max_mass*(μ_2/(μ_1 + μ_2)^2)*(3 - 2μ_2/(μ_1+μ_2))
            if isnan(likelihood[i_μ_1, i_f])
                likelihood[i_μ_1, i_f] = 0.0
            end
        end
    end

    return likelihood
end


mesa_model_dir = "/home/coloboquito/work/progs/mesa/v501aur"
# process_mesa_history_data(mesa_model_dir)

tess_response_spl = read_tess_response_spl()
Av = 0.54*3.1
# mesa_dfs = get_mesa_df_processed.(readdir(mesa_processed_temp), Ref(tess_response_spl), Av, mesa_processed_temp)
mesa_dfs = get_mesa_df_processed.(readdir(mesa_processed_tess_dir), Ref("mesa-data-processed-tess"))
star_masses = [mesa_df.star_mass[1] for mesa_df in mesa_dfs]

begin

max_age = [mesa_df.star_age[end] for mesa_df in mesa_dfs]
max_mass = [mesa_df.star_mass[end] for mesa_df in mesa_dfs]

max_mass_fit = fit_max_mass(max_age, max_mass, 3)
max_mass_fit[1] -= 8e-4

poly_max_fit = Polynomial(max_mass_fit)

calc_max_mass(lg_age, poly_max_fit) = 10^poly_max_fit(lg_age)

lg_flux_rel = 1.91; lg_flux_err = 0.05
mass_function = 0.1373; mass_function_err = 0.0002

n_μ_1 = 500; n_μ_2 = 500; n_age = 1600
n_m_1 = 1600; n_m_2 = 10; n_f = 10

μ_1_start = 0.8; μ_1_end = 0.99; μ_1_step = (μ_1_end-μ_1_start)/n_μ_1
μ_2_start = 0.3; μ_2_end = 0.5; μ_2_step = (μ_2_end-μ_2_start)/n_μ_2
lg_age_step = (log10(age_2) - log10(age_1))/n_age
f_step = 10*mass_function_err/n_f

m_1_start = 2.5; m_1_end = 5.0; m_1_step = (m_1_end-m_1_start)/n_m_1
m_2_start = 1.3; m_2_end = 1.9; m_2_step = (m_2_end-m_2_start)/n_m_2

μs_1 = [μ_1_start + (i_μ - 0.5)*μ_1_step for i_μ = 1:n_μ_1]
μs_2 = [μ_2_start + (i_μ - 0.5)*μ_2_step for i_μ = 1:n_μ_2]

ms_1 = [m_1_start + (i_m - 0.5)*m_1_step for i_m = 1:n_m_1]
ms_2 = [m_2_start + (i_m - 0.5)*m_2_step for i_m = 1:n_m_2]
lg_ages = [log10(age_1) + (i - 0.5)*lg_age_step for i = 1:n_age]
max_masses = calc_max_mass.(lg_ages, Ref(poly_max_fit))

mass_functions = [mass_function - 5*mass_function_err + (i_f-0.5)*f_step for i_f = 1:n_f]

mass_matrix = zeros(n_m_1, n_f+1)
mass_matrix[:,1] .= ms_1
for (i_m_1, m_1) in enumerate(ms_1)
    for (i_f, f) in enumerate(mass_functions) 
        m_2 = find_secondary_mass(f, m_1)
        mass_matrix[i_m_1, i_f+1] = m_2
    end
end

likelihood = stack([calc_likelihood_mass_function(lg_age, ms_1/max_mass, mass_functions, mass_matrix, lg_flux_rel, lg_flux_err, 
                            mass_function, mass_function_err, mesa_dfs, poly_max_fit)/max_mass
                                     for (lg_age, max_mass) in zip(lg_ages, max_masses)])
likelihood = likelihood / sum(likelihood)/ lg_age_step/m_1_step/f_step

smooth_n = 4
smoothed_likelihood = zeros(n_m_1 ÷ smooth_n, n_f, n_age ÷ smooth_n)
smoothed_ms_1 = zeros(n_m_1 ÷ smooth_n)
smoothed_lg_ages = zeros(n_age ÷ smooth_n) 

for i_m_1 = 1:n_m_1÷smooth_n
    m_1_is = ((i_m_1-1)*smooth_n+1):(i_m_1*smooth_n)
    smoothed_ms_1[i_m_1] = sum(ms_1[m_1_is])/smooth_n
    for i_age = 1:n_age÷smooth_n
        age_is = ((i_age-1)*smooth_n+1):(i_age*smooth_n)
        smoothed_likelihood[i_m_1,:,i_age] .= sum(likelihood[m_1_is,:,age_is], dims = (1,3))[1,:,1]/smooth_n^2
        smoothed_lg_ages[i_age] = sum(lg_ages[age_is])/smooth_n
    end
end
end

begin 
fig = Figure()
grid = fig[1,1] = GridLayout()
ax_m1 = Axis(grid[1,1])
hidexdecorations!(ax_m1, grid = false)

ax_lgt = Axis(grid[2,2])
hideydecorations!(ax_lgt, grid = false)

ax_hm = Axis(grid[2,1], xlabel = "Primary mass", ylabel = "Age")

linkyaxes!(ax_hm, ax_lgt)
linkxaxes!(ax_hm, ax_m1)

colgap!(grid, 10)
rowgap!(grid, 10)

rowsize!(grid, 2, Auto(2))
colsize!(grid, 1, Auto(2))

lines!(ax_lgt, sum(smoothed_likelihood, dims = (1,2))[1,1,:]*m_1_step*f_step*smooth_n, smoothed_lg_ages)
lines!(ax_m1, smoothed_ms_1, sum(smoothed_likelihood, dims = (2,3))[:,1,1]*f_step*lg_age_step*smooth_n)

hm = heatmap!(ax_hm, smoothed_ms_1, smoothed_lg_ages, sum(smoothed_likelihood, dims = (2))[:,1,:]*f_step, colormap = :binary)
max_masses = calc_max_mass.(smoothed_lg_ages, Ref(poly_max_fit))

lines!(ax_hm, max_masses[max_masses .< 5], smoothed_lg_ages[max_masses .< 5], label = "Death line", linestyle = :dash, color = :red)

Colorbar(fig[2,1], hm, vertical = false)
Legend(grid[1,2], ax_hm, tellwidth = false)


fig
end

# lg_ages = [8.0:0.1:8.7;]
# isochrone_dfs = get_isochrone_df.(Ref(mesa_dfs), lg_ages)

# GLMakie.activate!()

# fig_binary = Figure(size = (900, 900))
# ax_mass = Axis(fig_binary[1,1], xlabel = L"M_1", ylabel = L"\lg F_1 - \lg F_2", 
#                 xticks = WilkinsonTicks(6;k_min = 5,k_max = 7), yticks = WilkinsonTicks(6;k_min = 5,k_max = 7))
# ax_2mass = Axis(fig_binary[1,2], xlabel = L"M_2", ylabel = L"\lg F_1 - \lg F_2", 
#                 xticks = WilkinsonTicks(6;k_min = 5,k_max = 7), yticks = WilkinsonTicks(6;k_min = 5,k_max = 7))
# ax_logg = Axis(fig_binary[2,1], xlabel = L"\lg g", ylabel = L"\lg F_1 - \lg F_2", 
#                 xticks = WilkinsonTicks(6;k_min = 5,k_max = 7), yticks = WilkinsonTicks(6;k_min = 5,k_max = 7))
# ylims!(ax_logg, (1.4, 2.4))
# xlims!(ax_logg, (1.4, 3.0))
# ax_teff = Axis(fig_binary[2,2], xlabel = L"T_\mathrm{eff}", ylabel = L"\lg F_1 - \lg F_2", 
#                     xticks = WilkinsonTicks(6,k_min = 5,k_max = 7), yticks = WilkinsonTicks(6,k_min = 5,k_max = 7))
# ylims!(ax_teff, (1.4, 2.4))
# xlims!(ax_teff, (4000, 6000))



# for (i_age, lg_age) in enumerate(lg_ages)
#     isochrone_df = isochrone_dfs[i_age]
#     isochrone_df.log_L[isochrone_df.model_number .< 0] .= NaN
#     primary_masses = [2.25:0.01:7.0;]
#     secondary_masses = find_secondary_mass.(mass_function, primary_masses)
#     primary_dfs = foldl(vcat, interpolate_isochrone_df.(Ref(isochrone_df), primary_masses))
#     secondary_dfs = foldl(vcat, interpolate_isochrone_df.(Ref(isochrone_df), secondary_masses))

#     lg_TESS_primary = primary_dfs.log_TESS
#     lg_TESS_secondary = secondary_dfs.log_TESS
#     lg_L_primary = primary_dfs.log_L
#     lg_L_secondary = secondary_dfs.log_L

#     Δlg_L = lg_L_primary - lg_L_secondary
#     Δlg_TESS = lg_TESS_primary - lg_TESS_secondary
#     if isempty(Δlg_L[@. !isnan(Δlg_L)])
#         continue
#     end
#     lines!(ax_mass, primary_dfs.star_mass, Δlg_TESS, label = L"\lg t = %$lg_age")
#     lines!(ax_2mass, secondary_dfs.star_mass, Δlg_TESS, label = L"\lg t = %$lg_age")
#     lines!(ax_logg, primary_dfs.log_g, Δlg_TESS)
#     lines!(ax_teff, primary_dfs.Teff, Δlg_TESS)

#     scatter!(ax_mass, primary_dfs.star_mass, Δlg_TESS, label = L"\lg t = %$lg_age")
#     scatter!(ax_2mass, secondary_dfs.star_mass, Δlg_TESS, label = L"\lg t = %$lg_age")
#     scatter!(ax_logg, primary_dfs.log_g, Δlg_TESS)
#     scatter!(ax_teff, primary_dfs.Teff, Δlg_TESS)
# end
# axislegend(ax_mass, position = :lt, merge = true)



# fig = Figure()
# ax_tracks = Axis(fig[1,1], xreversed = true)
# ax_age = Axis(fig[1,3])
# ax_iso = Axis(fig[1,2])

# for (i_mass, star_mass) in enumerate(star_masses)
#     lines!(ax_tracks, mesa_dfs[i_mass].Teff, mesa_dfs[i_mass].log_L, label = "M = $star_mass")
#     lines!(ax_age, mesa_dfs[i_mass].star_mass, log10.(mesa_dfs[i_mass].star_age))
# end

# for (i_age, lg_age) in enumerate(lg_ages)
#     isochrone_df = isochrone_dfs[i_age]
#     isochrone_df.log_L[isochrone_df.model_number .< 0] .= NaN
#     masses = star_masses
#     lg_L = interpolate_isochrone_lg_L.(Ref(isochrone_df), masses)
#     lines!(ax_iso, masses, lg_L, label = "lg t = $lg_age")
# end

# axislegend(ax_tracks)
# axislegend(ax_iso)
# fig_binary

