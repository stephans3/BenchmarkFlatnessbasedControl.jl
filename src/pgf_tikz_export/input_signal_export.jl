using DelimitedFiles

N = 999;
ws = [1.5, 2.0, 2.5, 3.0]
data_u_al = zeros(BigFloat,N,0);
data_u_st = zeros(BigFloat,N,0);

for (i,w) in enumerate(ws)
    path_al = string("u_aluminum_w_", round(Int,w*10))
    path_st = string("u_steel_w_", round(Int,w*10))
    path_u_al = "results/u_input/"*path_al*".txt"
    path_u_st = "results/u_input/"*path_st*".txt"

    u_al = readdlm(path_u_al, '\t', BigFloat, '\n')
    u_st = readdlm(path_u_st, '\t', BigFloat, '\n')

    data_u_al = hcat(data_u_al, u_al)
    data_u_st = hcat(data_u_st, u_st)
end

replace!(data_u_al, NaN => 0.0)
replace!(data_u_st, NaN => 0.0)

data_u_al = vcat(zeros(1,length(ws)), data_u_al,zeros(1,length(ws)))
data_u_al = hcat(collect(0:1:N+1),data_u_al)

data_u_st = vcat(zeros(1,length(ws)), data_u_st,zeros(1,length(ws)))
data_u_st = hcat(collect(0:1:N+1),data_u_st)

data_u_al = round.(Float64.(data_u_al),digits=6)
data_u_st = round.(Float64.(data_u_st),digits=6)

path_u_al = "results/latex_export/u_aluminum_export.dat"
path_u_st = "results/latex_export/u_steel_export.dat"

open(path_u_al, "w") do io
    writedlm(io,  data_u_al, ' ')
end

open(path_u_st, "w") do io
    writedlm(io, data_u_st, ' ')
end