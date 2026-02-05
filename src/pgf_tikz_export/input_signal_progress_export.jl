#=
- Compute sequence μ_i and ratio μ_i / max μ_j for j ∈ {1,...,i} 
- Compute input signal u(t) for aluminum and steel 38Si7 for variable final iteration N
=#

L = 0.2; # Length of 1D rod

# Aluminium
λa = 237;  # Thermal conductivity
ρa = 2700; # Density
ca = 900;  # Specific heat capacity
α_a = λa / (ρa * ca) # Diffusivity
γ_a = L^2 / α_a

# Steel 38Si7
λs = 40;   # Thermal conductivity
ρs = 7800; # Density
cs = 460;  # Specific heat capacity
α_s = λs / (ρs * cs) # Diffusivity
γ_s = L^2 / α_s

η(L,α,i) = BigFloat(L)^(2i+1) / (BigFloat(α)^(i+1) * factorial(big(2i+1)))

idx_grid = 0:40;
eta_al = zeros(BigFloat, length(idx_grid))
eta_st = zeros(BigFloat, length(idx_grid))

for (n, iter) in enumerate(idx_grid)
    eta_al[n] = η(L,α_a,iter)
    eta_st[n] = η(L,α_s,iter)
end

using FastGaussQuadrature
w = 2;
T = 1000;
bump(t) = exp(-1 / (t/T - (t/T)^2)^w)
t_gq, weights_gq = FastGaussQuadrature.gausslegendre(1000)
p = T/2;
Ω_int = p *FastGaussQuadrature.dot( weights_gq ,bump.(p*t_gq .+ p))

diff_ref = 100; # (y_f - y_0) = 100 Kelvin

using DelimitedFiles
T = 1000; 
path_2_file = string("results/h_results/h_results_T_", T, ".txt")
h_data = readdlm(path_2_file, '\t', BigFloat, '\n')


u_al_all = λa*(diff_ref/Ω_int) * (h_data' .* eta_al)'
u_st_all = λs*(diff_ref/Ω_int) * (h_data' .* eta_st)'

u_al_sum = similar(u_al_all)
u_st_sum = similar(u_st_all)

u_al_norm_2 = zeros(size(u_al_all)[2])
u_st_norm_2 = zeros(size(u_st_all)[2])

using LinearAlgebra
for i in axes(u_al_all)[2]
    # Compute sum for input signal
    u_al_sum[:,i] = sum(u_al_all[:,1:i], dims=2)
    u_st_sum[:,i] = sum(u_st_all[:,1:i], dims=2)
    
    # Compute norm of input signal
    u_al_norm_2[i] = norm(u_al_all[:,i],2)
    u_st_norm_2[i] = norm(u_st_all[:,i],2)
end

u_al_norm_2_log10 = log10.(u_al_norm_2)
u_st_norm_2_log10 = log10.(u_st_norm_2)

# Ratio μ_i / max μ_j for j ∈ {1,...,i} 
# -> to find maximum iteration number
norm_rel_al = similar(u_al_norm_2)
norm_rel_st = similar(u_st_norm_2)

for (idx, el) in enumerate(u_al_norm_2)
    norm_rel_al[idx] = el / maximum(u_al_norm_2[begin:idx])
end

for (idx, el) in enumerate(u_st_norm_2)
    norm_rel_st[idx] = el / maximum(u_st_norm_2[begin:idx])
end

idx_grid = 0:size(u_al_all)[2]-1;
data_norm_abs = hcat(idx_grid, u_al_norm_2_log10, u_st_norm_2_log10)
data_norm_rel = hcat(idx_grid, norm_rel_al, norm_rel_st)

using DelimitedFiles
path_norm_abs = "results/latex_export/norm_abs_export.dat"
path_norm_rel = "results/latex_export/norm_rel_export.dat"

open(path_norm_abs, "w") do io
    writedlm(io,  data_norm_abs, ' ')
end

open(path_norm_rel, "w") do io
    writedlm(io,  data_norm_rel, ' ')
end

J = 1000;   
tgrid = T/J : T/J : T-T/J;
data_u_al = hcat(tgrid, u_al_sum[:,1], u_al_sum[:,3], u_al_sum[:,7])
data_u_st = hcat(tgrid, u_st_sum[:,5], u_st_sum[:,10], u_st_sum[:,15])

data_u_al =  round.(Float64.(data_u_al),digits=6)
data_u_st =  round.(Float64.(data_u_st),digits=6)

path_u_al = "results/latex_export/u_input_aluminum_progress.dat"
path_u_st = "results/latex_export/u_input_steel_progress.dat"

open(path_u_al, "w") do io
    writedlm(io,  data_u_al, ' ')
end

open(path_u_st, "w") do io
    writedlm(io,  data_u_st, ' ')
end

