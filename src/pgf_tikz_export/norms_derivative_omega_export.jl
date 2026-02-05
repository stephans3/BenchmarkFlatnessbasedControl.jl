# Create Figure: Logarithmic norms of d/dt Ω(t) with respect to T and ω
using DelimitedFiles, LinearAlgebra, FastGaussQuadrature
t_gq, weights_gq = FastGaussQuadrature.gausslegendre(10000)

function bump_fun(t; T=1, w=2)

    if (t<=0) || (t>=T)
        return 0;
    else
        return exp(-1 / (t/T - (t/T)^2)^w)
    end
end

# Fixed ω=2 and variable T
T_set = [10, 100, 1000];
h_norm_T = zeros(41,length(T_set))
w = 2.0

for (idx, T) in enumerate(T_set)
    p = T/2;
    Ω_den = p *FastGaussQuadrature.dot( weights_gq ,bump_fun.(p*t_gq .+ p,T=T, w=w))
    
    path_2_file = string("results/h_results/h_results_T_", T, ".txt")
    h_data = readdlm(path_2_file, '\t', BigFloat, '\n')
    h_data = h_data[:,1:41]/Ω_den
    h_norm_T[:,idx] = log10.( mapslices(x->norm(x,2), h_data[2:end-1,:], dims=1))
end


# Fixed T=1000 and variable ω
w_set = [11, 15, 20, 25, 30]
h_norm_w = zeros(41,length(w_set))
T = 1000;
p = T/2;

for (idx, w) in enumerate(w_set)
    Ω_den = p *FastGaussQuadrature.dot( weights_gq ,bump_fun.(p*t_gq .+ p,T=T, w=0.1*w))

    path_2_file = string("results/h_results/h_results_T_1000_w_", w, ".txt")
    h_data = readdlm(path_2_file, '\t', BigFloat, '\n')
    h_data = h_data[:,1:41]/Ω_den
    h_norm_w[:,idx] = log10.( mapslices(x->norm(x,2), h_data[2:end-1,:], dims=1))
end



# log_eta_al = round.(Float64.(log_eta_al),digits=6)
# log_eta_st = round.(Float64.(log_eta_st),digits=6)

idx_grid = 0:size(h_norm_T)[1]-1;
data_norm_T = hcat(idx_grid, h_norm_T)

idx_grid = 0:size(h_norm_w)[1]-1;   
data_norm_w = hcat(idx_grid, h_norm_w)

replace!(data_norm_T, NaN => 0.0)
replace!(data_norm_w, NaN => 0.0)

path_w = "results/latex_export/norm_omega_der_fix_w_2.dat"
path_T = "results/latex_export/norm_omega_der_fix_T_1000.dat"

open(path_w, "w") do io
    writedlm(io,  data_norm_T, ' ')
end

open(path_T, "w") do io
    writedlm(io, data_norm_w, ' ')
end

