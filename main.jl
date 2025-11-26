using Makie
using CairoMakie, GLMakie
using DataFrames, CSV
using Polynomials
using LaTeXStrings
using Printf
using FITSIO
using Dierckx


include("kurucz-int.jl")

star_masses = [1:0.15:7; 3.0:0.05:5.0;]
sort!(star_masses)
mesa_dir = "mesa-data"
mesa_processed_dir = "mesa-data-processed"
mesa_processed_temp = "mesa-data-processed-temp"

mkpath(mesa_dir)
mkpath(mesa_processed_dir)
mkpath(mesa_processed_temp)

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

mesa_model_dir = "/home/coloboquito/work/progs/mesa/v501aur"
# process_mesa_history_data(mesa_model_dir)

tess_response_spl = read_tess_response_spl()
Av = 0.54*3.1
# mesa_dfs = get_mesa_df_processed.(readdir(mesa_processed_temp), Ref(tess_response_spl), Av, mesa_processed_temp)
mesa_dfs = get_mesa_df_processed.(readdir(mesa_processed_dir), Ref("mesa-data-processed-tess"))
star_masses = [mesa_df.star_mass[1] for mesa_df in mesa_dfs]
lg_ages = [8.0:0.1:8.6;]

luminosity_rel = 100
mass_function = 0.1373


isochrone_dfs = get_isochrone_df.(Ref(mesa_dfs), lg_ages)

CairoMakie.activate!()

fig_binary = Figure(size = (900, 900))
ax_mass = Axis(fig_binary[1,1], xlabel = L"M_1", ylabel = L"\lg F_1 - \lg F_2", 
                xticks = WilkinsonTicks(6;k_min = 5,k_max = 7), yticks = WilkinsonTicks(6;k_min = 5,k_max = 7))
ax_2mass = Axis(fig_binary[1,2], xlabel = L"M_2", ylabel = L"\lg F_1 - \lg F_2", 
                xticks = WilkinsonTicks(6;k_min = 5,k_max = 7), yticks = WilkinsonTicks(6;k_min = 5,k_max = 7))
ax_logg = Axis(fig_binary[2,1], xlabel = L"\lg g", ylabel = L"\lg F_1 - \lg F_2", 
                xticks = WilkinsonTicks(6;k_min = 5,k_max = 7), yticks = WilkinsonTicks(6;k_min = 5,k_max = 7))
ylims!(ax_logg, (1.4, 2.4))
xlims!(ax_logg, (1.4, 3.0))
ax_teff = Axis(fig_binary[2,2], xlabel = L"T_\mathrm{eff}", ylabel = L"\lg F_1 - \lg F_2", 
                    xticks = WilkinsonTicks(6,k_min = 5,k_max = 7), yticks = WilkinsonTicks(6,k_min = 5,k_max = 7))
ylims!(ax_teff, (1.4, 2.4))
xlims!(ax_teff, (4000, 6000))



for (i_age, lg_age) in enumerate(lg_ages)
    isochrone_df = isochrone_dfs[i_age]
    isochrone_df.log_L[isochrone_df.model_number .< 0] .= NaN
    primary_masses = [2.25:0.01:7.0;]
    secondary_masses = find_secondary_mass.(mass_function, primary_masses)
    primary_dfs = foldl(vcat, interpolate_isochrone_df.(Ref(isochrone_df), primary_masses))
    secondary_dfs = foldl(vcat, interpolate_isochrone_df.(Ref(isochrone_df), secondary_masses))

    lg_TESS_primary = primary_dfs.log_TESS
    lg_TESS_secondary = secondary_dfs.log_TESS
    lg_L_primary = primary_dfs.log_L
    lg_L_secondary = secondary_dfs.log_L

    Δlg_L = lg_L_primary - lg_L_secondary
    Δlg_TESS = lg_TESS_primary - lg_TESS_secondary
    if isempty(Δlg_L[@. !isnan(Δlg_L)])
        continue
    end
    lines!(ax_mass, primary_dfs.star_mass, Δlg_TESS, label = L"\lg t = %$lg_age")
    lines!(ax_2mass, secondary_dfs.star_mass, Δlg_TESS, label = L"\lg t = %$lg_age")
    lines!(ax_logg, primary_dfs.log_g, Δlg_TESS)
    lines!(ax_teff, primary_dfs.Teff, Δlg_TESS)
end
axislegend(ax_mass, position = :lt)



fig = Figure()
ax_tracks = Axis(fig[1,1], xreversed = true)
ax_age = Axis(fig[1,3])
ax_iso = Axis(fig[1,2])

for (i_mass, star_mass) in enumerate(star_masses)
    lines!(ax_tracks, mesa_dfs[i_mass].Teff, mesa_dfs[i_mass].log_L, label = "M = $star_mass")
    lines!(ax_age, mesa_dfs[i_mass].star_mass, log10.(mesa_dfs[i_mass].star_age))
end

for (i_age, lg_age) in enumerate(lg_ages)
    isochrone_df = isochrone_dfs[i_age]
    isochrone_df.log_L[isochrone_df.model_number .< 0] .= NaN
    masses = star_masses
    lg_L = interpolate_isochrone_lg_L.(Ref(isochrone_df), masses)
    lines!(ax_iso, masses, lg_L, label = "lg t = $lg_age")
end

axislegend(ax_tracks)
axislegend(ax_iso)
fig_binary