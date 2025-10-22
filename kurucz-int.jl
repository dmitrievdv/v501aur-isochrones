function find_near_nodes(x, nodes)
    x_2 = findfirst(y -> y ≥ 0.0, nodes .- x)
    if x_2 == 1
        return 2, 1
    end
    if (isnothing(x_2)) #& ((nodes[end] - x) ≥ -1e-3)
        return length(nodes), length(nodes) - 1
    end 
    return x_2 - 1, x_2
end

function linear_interpolation(x, x_1, x_2, y_1, y_2)
    return y_1 + (x - x_1)*(y_2 - y_1)/(x_2 - x_1)
end

function interpolate_kurucz(T_eff, logg)
    logg_nodes = [0.0:0.5:5.0;]
    T_eff_nodes = [3500:250.0:9750; 10000:5e2:12500; 13000:1e3:35000]

    i_T_eff_1, i_T_eff_2 = find_near_nodes(T_eff, T_eff_nodes)
    i_logg_1, i_logg_2 = find_near_nodes(logg, logg_nodes)

    kurucz_fits_1 = FITS("kurucz-p00/kp00_$(round(Int, T_eff_nodes[i_T_eff_1])).fits")
    kurucz_fits_2 = FITS("kurucz-p00/kp00_$(round(Int, T_eff_nodes[i_T_eff_2])).fits")

    λs_1 = read(kurucz_fits_1[2], "WAVELENGTH") 
    λs_2 = read(kurucz_fits_2[2], "WAVELENGTH")

    int_g1 = round(Int, logg_nodes[i_logg_1]*10)
    int_g2 = round(Int, logg_nodes[i_logg_2]*10)

    int_g1_str = @sprintf "%02i" int_g1
    int_g2_str = @sprintf "%02i" int_g2

    kurucz_spec_11 = read(kurucz_fits_1[2], "g$int_g1_str")
    kurucz_spec_12 = read(kurucz_fits_1[2], "g$int_g2_str")

    kurucz_spec_21 = read(kurucz_fits_2[2], "g$int_g1_str")
    kurucz_spec_22 = read(kurucz_fits_2[2], "g$int_g2_str")

    kurucz_spec_1 = linear_interpolation(logg, logg_nodes[[i_logg_1, i_logg_2]]..., kurucz_spec_11, kurucz_spec_12)
    kurucz_spec_2 = linear_interpolation(logg, logg_nodes[[i_logg_1, i_logg_2]]..., kurucz_spec_21, kurucz_spec_22)

    if !isnothing(findfirst(x -> abs(x) > 1e-3, λs_1 .- λs_2))
        kurucz_spec_2_spl = Spline1D(λs_2, kurucz_spec_2, k = 1)
        kurucz_spec_2 .= kurucz_spec_2_spl.(λs_1)
    end
        
    kurucz_spec = linear_interpolation(T_eff, T_eff_nodes[[i_T_eff_1, i_T_eff_2]]..., kurucz_spec_1, kurucz_spec_2)
    return λs_1, kurucz_spec
end

function read_tess_response_spl()
    tess_response_df = CSV.read("tess-response-function-v2.0.csv", DataFrame, comment = "#")
    tess_response_spl = Spline1D(tess_response_df[:,1] * 10, tess_response_df[:,2], k = 1)
end

function calc_tess_kurucz_flux_no_extinction(T_eff, logg, tess_response_spl)
    λs, kurucz_spec = interpolate_kurucz(T_eff, logg)
# wls, flux, cntm = Korg.synth(Teff = T_eff,
#                              logg = logg, 
#                              m_H = 0.0, 
#                              wavelengths = Korg.Wavelengths(5000, 12000, 1), 
#                              linelist = Korg.read_linelist("gfallvac08oct17.dat", format = "kurucz_vac"), 
#                              R = 100)

    Δλs_kurucz = λs[2:end] - λs[1:end-1]
    λs_mid_kurucz = (λs[2:end] + λs[1:end-1])/2
    kurucz_spec_mid = (kurucz_spec[2:end] + kurucz_spec[1:end-1])/2
    return sum(@. kurucz_spec_mid * tess_response_spl(λs_mid_kurucz) * Δλs_kurucz)
end


function calc_tess_kurucz_flux_no_extinction(T_eff, logg)
    tess_response_spl = read_tess_response_spl()
    calc_tess_kurucz_flux_no_extinction(T_eff, logg, tess_response_spl)
end

extinct_law(λ, Av) = Av/(2.5*log10(ℯ))/(λ/5e3)

function calc_tess_kurucz_flux_with_extinction(T_eff, logg, Av, tess_response_spl)
    λs, kurucz_spec = interpolate_kurucz(T_eff, logg)
# wls, flux, cntm = Korg.synth(Teff = T_eff,
#                              logg = logg, 
#                              m_H = 0.0, 
#                              wavelengths = Korg.Wavelengths(5000, 12000, 1), 
#                              linelist = Korg.read_linelist("gfallvac08oct17.dat", format = "kurucz_vac"), 
#                              R = 100)

    Δλs_kurucz = λs[2:end] - λs[1:end-1]
    λs_mid_kurucz = (λs[2:end] + λs[1:end-1])/2
    kurucz_spec_mid = (kurucz_spec[2:end] + kurucz_spec[1:end-1])/2
    return sum(@. kurucz_spec_mid * tess_response_spl(λs_mid_kurucz) * exp(-extinct_law(λs_mid_kurucz, Av)) * Δλs_kurucz)
end

function calc_tess_kurucz_flux_with_extinction(T_eff, logg, Av)
    tess_response_spl = read_tess_response_spl()
    calc_tess_kurucz_flux_with_extinction(T_eff, logg, Av, tess_response_spl)
end

# Av = 0.54 * 3.1

# tess_rel_no_Av = calc_tess_kurucz_flux_no_extinction(7800, 4.2)/calc_tess_kurucz_flux_no_extinction(4700, 2.2)*1.5^2/26.6^2
# tess_rel_Av = calc_tess_kurucz_flux_with_extinction(7800, 4.2, Av)/calc_tess_kurucz_flux_with_extinction(4700, 2.2, Av)*1.5^2/26.6^2

# logg_grid = [1.5:0.1:5.0;]
# T_eff_grid = [3500:1e2:10000;]

# n_logg = length(logg_grid)
# n_T_eff = length(T_eff_grid)

# tess_rel_Av_arr = zeros(n_logg, n_T_eff)
# tess_rel_arr = zeros(n_logg, n_T_eff)

# for (i_logg, logg) in enumerate(logg_grid)
#     for (i_T_eff, T_eff) in enumerate(T_eff_grid)
#         if (T_eff > 9250) & (logg < 2.0)
#             tess_rel_arr[i_logg, i_T_eff] = -1.0
#             tess_rel_Av_arr[i_logg, i_T_eff] = -1.0
#         else
#             tess_rel_arr[i_logg, i_T_eff] = calc_tess_kurucz_flux_no_extinction(T_eff, logg)
#             tess_rel_Av_arr[i_logg, i_T_eff] = calc_tess_kurucz_flux_with_extinction(T_eff, logg, Av)
#         end
#     end
# end

# mkpath("tess-tables/")
# for (i_T_eff, T_eff) in enumerate(T_eff_grid)
#     open("tess-tables/tess_flux_$(round(Int, T_eff)).dat", "w") do io
#         println(io, "# TESS fluxes for T_eff = $(round(Int, T_eff)).")
#         println(io, "# Fluxes are in arbitrary units, meant to only be used in relation to each other.")
#         println(io, "# Negative numbers mean that this (T_eff, logg) combination is invalid for Kurucz, 1993 models.")
#         println(io, "# Columns: logg, flux, flux with extinction. \n")
#         for i_logg in 1:n_logg
#             @printf io "%5f %10e %10e \n" logg_grid[i_logg] tess_rel_arr[i_logg, i_T_eff] tess_rel_Av_arr[i_logg, i_T_eff]
#         end
#     end
# end