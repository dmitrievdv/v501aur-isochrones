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
    mesa_df = DataFrame([(name_col_tuple[1] => name_col_tuple[2]) for name_col_tuple in zip(col_names, cols)])
end

function get_age_df(mesa_df, lg_age)
    lg_ages = log10.(mesa_df.star_age)
    i_age = findfirst(x -> x > lg_age, lg_ages)
    println(i_age)

    if isnothing(i_age)
        df = DataFrame([(name_val_tuple[1] => name_val_tuple[2]) for name_val_tuple in zip(names(mesa_df), collect(mesa_df[end,:]))])
        df.model_number .= -1.0
        df.star_age .= 10^lg_age
        return df
    end

    vals_prev_age = collect(mesa_df[i_age-1,:])
    vals_next_age = collect(mesa_df[i_age,:])

    prev_age = mesa_df.star_age[i_age-1]
    next_age = mesa_df.star_age[i_age]

    vals_age = vals_prev_age + (10^lg_age - prev_age)/(next_age - prev_age)*(vals_next_age - vals_prev_age)
    return DataFrame([(name_val_tuple[1] => name_val_tuple[2]) for name_val_tuple in zip(names(mesa_df), vals_age)])
end

function get_isochrone_df(mesa_dfs, lg_age)
    foldl(vcat, get_age_df.(mesa_dfs, lg_age))
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

fig_transit = Figure()
ax_transit = Axis(fig_transit[1,1], xlabel = L"M_1", ylabel = L"\lg L_1 - \lg L_2")
for (i_age, lg_age) in enumerate(lg_ages)
    isochrone_df = isochrone_dfs[i_age]
    isochrone_df.log_L[isochrone_df.model_number .< 0] .= NaN
    primary_masses = [2.25:0.25:7.0;]
    secondary_masses = find_secondary_mass.(mass_function, primary_masses)
    lg_L_primary = interpolate_isochrone_lg_L.(Ref(isochrone_df), primary_masses)
    lg_L_secondary = interpolate_isochrone_lg_L.(Ref(isochrone_df), secondary_masses)
    Δlg_L = lg_L_primary - lg_L_secondary
    if isempty(Δlg_L[@. !isnan(Δlg_L)])
        continue
    end
    lines!(ax_transit, primary_masses, Δlg_L, label = "lg t = $lg_age")
end
axislegend(ax_transit)



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
fig