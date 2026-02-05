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

log_eta_al = log10.(eta_al)
log_eta_st = log10.(eta_st)

# Ratio η_{i+1} / η_{i}
ratio_al = eta_al[begin+1:end] ./  eta_al[begin:end-1]
ratio_st = eta_st[begin+1:end] ./  eta_st[begin:end-1]

ratio_al_log = log10.(ratio_al)
ratio_st_log = log10.(ratio_st)

# η reaches its maximum
eta_max_al = η(L,α_a,9)
eta_max_st = η(L,α_s,29)

# η drops below 1
eta_0_al = η(L,α_a,28)
eta_0_st = η(L,α_s,83)

log_eta_al = round.(Float64.(log_eta_al),digits=6)
log_eta_st = round.(Float64.(log_eta_st),digits=6)
data_eta = hcat(idx_grid,log_eta_al, log_eta_st)

using DelimitedFiles
path_eta = "results/latex_export/eta_sequence_export.dat"

open(path_eta, "w") do io
    writedlm(io,  data_eta, ' ')
end

ratio_al_log = round.(Float64.(ratio_al_log),digits=6)
ratio_st_log = round.(Float64.(ratio_st_log),digits=6)
data_ratio = hcat(idx_grid[begin:end-1],ratio_al_log, ratio_st_log)

path_ratio = "results/latex_export/eta_ratio_export.dat"

open(path_ratio, "w") do io
    writedlm(io,  data_ratio, ' ')
end
