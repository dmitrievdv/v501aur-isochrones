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

include("mesa-process.jl") # Обработка мезы

"""
find secondary mass from mass_function and primary mass
"""
function find_secondary_mass(mass_function, primary_mass)
    mass_poly = Polynomial([-mass_function*primary_mass^2, -2mass_function*primary_mass, -mass_function, 1])
    mass_poly_roots = Polynomials.roots(mass_poly)
    return real.(mass_poly_roots[abs.(imag.(mass_poly_roots)) .≈ 0.0])[1]
end


"""
фитирование максимальной массы для возраста полиномом степени poly_deg по дате max_age, max_mass
"""
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

# границы сетки по возрастам
age_1 = 1e8 
age_2 = 4e8

include("bayes.jl") # вся статистика тут

# mesa_model_dir = "/home/coloboquito/work/progs/mesa/v501aur" # папка с сырыми эволюционными моделями (актуально для компа в пулково)
# process_mesa_history_data(mesa_model_dir) # обработка сырых данных мезы, сохраняет в mesa-data-processed-temp

# расчет потоков тесс и считывание данных мезы

tess_response_spl = read_tess_response_spl() 
Av = 0.54*3.1
# mesa_dfs = get_mesa_df_processed.(readdir(mesa_processed_temp), Ref(tess_response_spl), Av, mesa_processed_temp) # добавление тесс потоков и сохранение в "mesa-data-processed-tess"
mesa_dfs = get_mesa_df_processed.(readdir(mesa_processed_tess_dir), Ref(mesa_processed_tess_dir)) # считывает из тесс папки
star_masses = [mesa_df.star_mass[1] for mesa_df in mesa_dfs] # начальные массы

begin # Байесовский блок -- вычисления статистики 
max_age = [mesa_df.star_age[end] for mesa_df in mesa_dfs] # конечные возраста
max_mass = [mesa_df.star_mass[end] for mesa_df in mesa_dfs] # конечные массы

max_mass_fit = fit_max_mass(max_age, max_mass, 3) # фит максимальной массы для возраста
max_mass_fit[1] -= 8e-4 # чтобы последняя точка изохроны всегда была больше 1 чуть опускаем весь фит (костыльно)
                        # т.к. в расчетах есть небольшой шум, из-за этого и фитируем 

poly_max_fit = Polynomial(max_mass_fit) # создаем объект-полином из фита

calc_max_mass(lg_age, poly_max_fit) = 10^poly_max_fit(lg_age)

# наблюдательные данные
lg_flux_rel = 1.91; lg_flux_err = 0.05 
mass_function = 0.1373; mass_function_err = 0.0002

# сетки по различным параметрам (μ = m/max_mass)
n_μ_1 = 500; n_μ_2 = 500; n_age = 400
n_m_1 = 400; n_m_2 = 10; n_f = 10

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

# расчет матрицы масс -- так как расчет массы из функции масс довольно затратен он выносится за расчет апостериора

mass_matrix = zeros(n_m_1, n_f+1)
mass_matrix[:,1] .= ms_1
for (i_m_1, m_1) in enumerate(ms_1)
    for (i_f, f) in enumerate(mass_functions) 
        m_2 = find_secondary_mass(f, m_1)
        mass_matrix[i_m_1, i_f+1] = m_2
    end
end

# расчет постериора по параметрам μ, f, lg t
posterior= stack([calc_posterior_mass_function(lg_age, ms_1/max_mass, mass_functions, mass_matrix, lg_flux_rel, lg_flux_err, 
                            mass_function, mass_function_err, mesa_dfs, poly_max_fit)/max_mass # деление на массу для перехода от μ к m
                                     for (lg_age, max_mass) in zip(lg_ages, max_masses)])
posterior = posterior / sum(posterior)/ lg_age_step/m_1_step/f_step

# сглаживание постериора расчетом вероятности в квадрате smooth_n*smooth_n точек сетки по m,lg t
smooth_n = 4
smoothed_posterior = zeros(n_m_1 ÷ smooth_n, n_f, n_age ÷ smooth_n)
smoothed_ms_1 = zeros(n_m_1 ÷ smooth_n)
smoothed_lg_ages = zeros(n_age ÷ smooth_n) 

for i_m_1 = 1:n_m_1÷smooth_n
    m_1_is = ((i_m_1-1)*smooth_n+1):(i_m_1*smooth_n)
    smoothed_ms_1[i_m_1] = sum(ms_1[m_1_is])/smooth_n
    for i_age = 1:n_age÷smooth_n
        age_is = ((i_age-1)*smooth_n+1):(i_age*smooth_n)
        smoothed_posterior[i_m_1,:,i_age] .= sum(posterior[m_1_is,:,age_is], dims = (1,3))[1,:,1]/smooth_n^2
        smoothed_lg_ages[i_age] = sum(lg_ages[age_is])/smooth_n
    end
end
end

begin # графика
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

lines!(ax_lgt, sum(smoothed_posterior, dims = (1,2))[1,1,:]*m_1_step*f_step*smooth_n, smoothed_lg_ages)
lines!(ax_m1, smoothed_ms_1, sum(smoothed_posterior, dims = (2,3))[:,1,1]*f_step*lg_age_step*smooth_n)

hm = heatmap!(ax_hm, smoothed_ms_1, smoothed_lg_ages, sum(smoothed_posterior, dims = (2))[:,1,:]*f_step, colormap = :binary)
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

