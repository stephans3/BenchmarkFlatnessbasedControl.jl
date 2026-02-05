using FastGaussQuadrature

t_gq, weights_gq = FastGaussQuadrature.gausslegendre(10000)

function bump_fun(t; T=1, w=2)

    if (t<=0) || (t>=T)
        return 0;
    else
        return exp(-1 / (t/T - (t/T)^2)^w)
    end
end

function transition_fun(t; T=1, w=2)
    if t<= 0
        return 0;
    elseif t >= T
        return 1;
    else
        q = t/2;
        p = T/2;
        Ω_num = q *FastGaussQuadrature.dot( weights_gq ,bump_fun.(q*t_gq .+ q,T=T, w=w))
        Ω_den = p *FastGaussQuadrature.dot( weights_gq ,bump_fun.(p*t_gq .+ p,T=T, w=w))
        return Ω_num / Ω_den 
    end

end

Tf = 100;
J = 1000;
tgrid = 0.0 : Tf/J : Tf;
p = Tf/2;
w_vec = [11,15,20,25,30]
# Transition
Φ = zeros(length(tgrid), length(w_vec))
d1Φ = zeros(length(tgrid), length(w_vec))

for (idx, w) in enumerate(w_vec)
    # Original function
    Φ[:,idx] = transition_fun.(tgrid,T=Tf,w=0.1*w) 

    # First derivative
    Ω_den = p *FastGaussQuadrature.dot( weights_gq ,bump_fun.(p*t_gq .+ p,T=Tf, w=0.1*w))
    Ω_num = bump_fun.(tgrid,T=Tf, w=0.1*w)
    d1Φ[:,idx] = Ω_num / Ω_den
end

data_Φ = hcat(tgrid, round.(Float64.(Φ),digits=6))
data_d1Φ = hcat(tgrid, round.(Float64.(d1Φ),digits=6))

using DelimitedFiles
path_Φ = "results/latex_export/transition_export.dat"

open(path_Φ, "w") do io
    writedlm(io,  data_Φ, ' ')
end

path_d1Φ = "results/latex_export/transition_bump_export.dat"

open(path_d1Φ, "w") do io
    writedlm(io,  data_d1Φ, ' ')
end

