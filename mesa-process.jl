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