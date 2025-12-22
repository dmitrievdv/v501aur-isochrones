


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
"""
функция масс
"""
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

# во многих функциях таскается poly_max_fit --- фит максимальных масс

"""
априор в параметрах lg t, m_1/max_mass, m_2/max_mass, где max_mass -- максимальная масса для возраста lg_age
"""
function prior_μ(lg_age, μ_1, μ_2, poly_max_fit)
    m = calc_max_mass(lg_age, poly_max_fit)
    m^2/(age_2 - age_1)*IMF(m*μ_1)*IMF(m*μ_2)*log(10)*10^lg_age
end

"""
получение важных для статистики параметров из даты (массы, )
"""
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

        itp_mass = AkimaInterpolation(data[itp_indeces, 2], data[itp_indeces, 1])
        itp_flux = AkimaInterpolation(data[itp_indeces, 3], data[itp_indeces, 1])

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


    itp_mass = AkimaInterpolation(itp_mass_points, itp_ages_points)
    itp_flux = AkimaInterpolation(itp_flux_points, itp_ages_points)

    isochrone[i_last+1,1] = itp_mass(log_age)
    isochrone[i_last+1,2] = itp_flux(log_age)

    
    return isochrone
end

function find_tess_flux(lg_age, iso_mass :: Real, mesa_dfs, poly_max_fit)
    isochrone = extract_important_data_isochrone(lg_age, mesa_dfs, poly_max_fit)
    # sort_is = sortperm(isochrone[:,1])
    itp_flux = AkimaInterpolation(isochrone[:,2], isochrone[:,1])
    itp_flux(iso_mass)
end

function find_tess_flux(lg_age, iso_mass, mesa_dfs, poly_max_fit)
    isochrone = extract_important_data_isochrone(lg_age, mesa_dfs, poly_max_fit)
    sort_is = sortperm(isochrone[:,1])
    itp_flux = AkimaInterpolation(isochrone[sort_is,2], isochrone[sort_is,1])
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

function calc_posterior(lg_age, μs_1, μs_2, lg_flux_rel, lg_flux_err, mass_function, mass_function_err, mesa_dfs, poly_max_fit) # assuming same grid for μ_1, μ_2
    fluxes = find_tess_flux(lg_age, [μs_1; μs_2], mesa_dfs, poly_max_fit)
    n_μ_1 = length(μs_1)
    fluxes_1 = fluxes[1:n_μ_1]
    n_μ_2 = length(μs_2)
    fluxes_2 = fluxes[n_μ_1+1:end]
    # println(fluxes)
    # n_μ = length(μs)
    posterior = Matrix{Float64}(undef, n_μ_1, n_μ_2)
    max_mass = calc_max_mass(lg_age, poly_max_fit)
    for (i_μ_1, μ_1) in enumerate(μs_1)
        if fluxes[i_μ_1] > 9e3
            posterior[i_μ_1, :] .= 0.0
            continue
        end
        m_1 = max_mass*μ_1
        for (i_μ_2, μ_2) in enumerate(μs_2)
            if fluxes[i_μ_2] > 9e3
                posterior[i_μ_1, i_μ_2] = 0.0
                continue
            end
            model_lg_flux_rel = fluxes_1[i_μ_1] - fluxes_2[i_μ_2]
            m_2 = max_mass*μ_2
            model_mass_function = calc_mass_function(m_1, m_2)
            # println("$m_1 $m_2 $model_mass_function $model_lg_flux_rel")
            # println(prior_μ(lg_age, μ_1, μ_2, poly_max_fit))
            # println(exp(-(model_lg_flux_rel - lg_flux_rel)^2/2lg_flux_err^2))
            # println(exp(-(model_mass_function - mass_function)^2/2mass_function_err^2))
            posterior[i_μ_1, i_μ_2] = prior_μ(lg_age, μ_1, μ_2, poly_max_fit)
            posterior[i_μ_1, i_μ_2] *= exp(-(model_lg_flux_rel - lg_flux_rel)^2/2lg_flux_err^2)
            posterior[i_μ_1, i_μ_2] *= exp(-(model_mass_function - mass_function)^2/2mass_function_err^2)
            if isnan(posterior[i_μ_1, i_μ_2])
                posterior[i_μ_1, i_μ_2] = 0.0
            end
        end
    end

    return posterior
end


function calc_posterior_mass_function(lg_age, μs_1, mass_functions, mass_matrix, lg_flux_rel, lg_flux_err, mass_function, mass_function_err, mesa_dfs, poly_max_fit) # assuming same grid for μ_1, μ_2
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
    posterior = Matrix{Float64}(undef, n_μ_1, n_f)
    
    for (i_μ_1, μ_1) in enumerate(μs_1)
        if fluxes_1[i_μ_1,1] > 9e3
            posterior[i_μ_1, :] .= 0.0
            continue
        end
        for (i_f, f) in enumerate(mass_functions)
            if fluxes_2[i_μ_1, i_f] > 9e3
                posterior[i_μ_1, i_f] = 0.0
                continue
            end
            μ_2 = mass_matrix[i_μ_1, i_f + 1]
            model_lg_flux_rel = fluxes_1[i_μ_1] - fluxes_2[i_μ_1, i_f]
            model_mass_function = f
            # println("$m_1 $m_2 $model_mass_function $model_lg_flux_rel")
            # println(prior_μ(lg_age, μ_1, μ_2, poly_max_fit))
            # println(exp(-(model_lg_flux_rel - lg_flux_rel)^2/2lg_flux_err^2))
            # println(exp(-(model_mass_function - mass_function)^2/2mass_function_err^2))
            posterior[i_μ_1, i_f] = prior_μ(lg_age, μ_1, μ_2, poly_max_fit)
            posterior[i_μ_1, i_f] *= exp(-(model_lg_flux_rel - lg_flux_rel)^2/2lg_flux_err^2)
            posterior[i_μ_1, i_f] *= exp(-(model_mass_function - mass_function)^2/2mass_function_err^2)
            posterior[i_μ_1, i_f] /= max_mass*(μ_2/(μ_1 + μ_2)^2)*(3 - 2μ_2/(μ_1+μ_2))
            if isnan(posterior[i_μ_1, i_f])
                posterior[i_μ_1, i_f] = 0.0
            end
        end
    end

    return posterior
end