using Makie
using GLMakie
using DataFrames
using Polynomials
using LaTeXStrings

star_masses = [2.5:0.25:7.0;]
mesa_dir = "mesa-data"

function get_mesa_df(star_mass)
    mesa_file_name = "track_M$(star_mass)_Z0.0147_ROT0.00_history.data"

    lines = readlines("$mesa_dir/$mesa_file_name")
    splitted_lines = [split(line) for line in lines]
    col_names = splitted_lines[6]
    n_cols = length(col_names)
    cols = [[parse(Float64, splitted_line[i_col]) for splitted_line in splitted_lines[7:end]] for i_col in 1:n_cols]
    mesa_df = DataFrame([Pair(name_col_tuple...) for name_col_tuple in zip(col_names, cols)])
end

function get_age_df(mesa_df, lg_age)
    lg_ages = log10.(mesa_df.star_age)
    i_age = findfirst(x -> x > lg_age, lg_ages)
    println(i_age)

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
    mass_poly_roots = roots(mass_poly)
    return real.(mass_poly_roots[abs.(imag.(mass_poly_roots)) .≈ 0.0])[1]
end

mesa_dfs = get_mesa_df.(star_masses)
lg_ages = [7.5:0.1:8.5;]

luminosity_rel = 100
mass_function = 0.4

isochrone_dfs = get_isochrone_df.(Ref(mesa_dfs), lg_ages)

fig_binary = Figure()
ax_mass = Axis(fig_binary[1,1], xlabel = L"M_1", ylabel = L"\lg L_1 - \lg L_2")
ax_logg = Axis(fig_binary[1,2], xlabel = L"\lg g", ylabel = L"\lg L_1 - \lg L_2")
ax_teff = Axis(fig_binary[1,3], xlabel = L"T_\mathrm{eff}", ylabel = L"\lg L_1 - \lg L_2")

for (i_age, lg_age) in enumerate(lg_ages)
    isochrone_df = isochrone_dfs[i_age]
    isochrone_df.log_L[isochrone_df.model_number .< 0] .= NaN
    primary_masses = [2.25:0.25:7.0;]
    secondary_masses = find_secondary_mass.(mass_function, primary_masses)
    primary_dfs = foldl(vcat, interpolate_isochrone_df.(Ref(isochrone_df), primary_masses))
    secondary_dfs = foldl(vcat, interpolate_isochrone_df.(Ref(isochrone_df), secondary_masses))
    lg_L_primary = primary_dfs.log_L
    lg_L_secondary = secondary_dfs.log_L
    Δlg_L = lg_L_primary - lg_L_secondary
    if isempty(Δlg_L[@. !isnan(Δlg_L)])
        continue
    end
    lines!(ax_mass, primary_dfs.star_mass, Δlg_L, label = "lg t = $lg_age")
    lines!(ax_logg, primary_dfs.log_g, Δlg_L)
    lines!(ax_teff, primary_dfs.Teff, Δlg_L)
end
axislegend(ax_mass)



fig = Figure()
ax_tracks = Axis(fig[1,1], xreversed = true)
ax_iso = Axis(fig[1,2])

for (i_mass, star_mass) in enumerate(star_masses)
    lines!(ax_tracks, mesa_dfs[i_mass].Teff, mesa_dfs[i_mass].log_L, label = "M = $star_mass")
end

for (i_age, lg_age) in enumerate(lg_ages)
    isochrone_df = isochrone_dfs[i_age]
    isochrone_df.log_L[isochrone_df.model_number .< 0] .= NaN
    masses = [2.25:0.01:7.0;]
    lg_L = interpolate_isochrone_lg_L.(Ref(isochrone_df), masses)
    lines!(ax_iso, masses, lg_L, label = "lg t = $lg_age")
end

axislegend(ax_tracks)
axislegend(ax_iso)
fig_binary