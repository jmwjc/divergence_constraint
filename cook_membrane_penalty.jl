using Revise
using TimerOutputs 
using SparseArrays, Pardiso
using CairoMakie
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵈᵢⱼσᵈᵢⱼdxdy, ∫∫qpdxdy, ∫∫p∇udxdy, ∫vᵢgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric

include("import_cook.jl")

const to = TimerOutput()
ps = MKLPardisoSolver()

ndiv = 2
@timeit to "import data" begin
n = 2
elements, nodes, nodes_p = import_elasticity_linear_mix("./msh/plate_with_hole_tri3_"*string(ndiv)*".msh","./msh/plate_with_hole_tri3_"*string(n)*".msh",n)
# nx = 7;ny = 3
# elements, nodes, nodes_p = import_elasticity_linear_mix("./msh/plate_with_hole_tri3_"*string(ndiv)*".msh","./msh/plate_with_hole_tri3_"*string(ny)*"_"*string(nx)*".msh",ny)
# elements, nodes, nodes_p = import_elasticity_quadratic_mix("./msh/plate_with_hole_tri6_"*string(ndiv)*".msh","./msh/plate_with_hole_tri6_"*string(n)*".msh",n)
# nx = 68;ny = 32
# elements, nodes, nodes_p = import_elasticity_quadratic_mix("./msh/plate_with_hole_tri6_"*string(ndiv)*".msh","./msh/plate_with_hole_tri3_"*string(ny)*"_"*string(nx)*".msh",ny)

nₚ = length(nodes_p)
end

nₑ = length(elements["Ωᵘ"])
# nₛ = 3
nᵤ = length(nodes)

# T = 1.0e3
# E = 3.0e6
P = 6.25
E = 70
# ν = 0.3
ν = 0.5-1e-8
Ē = E/(1.0-ν^2)
ν̄ = ν/(1.0-ν)
Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
Cᵢⱼᵢⱼ = E/(1+ν)/2

prescribe!(elements["Ωᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵘ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->P) 
prescribe!(elements["Γᵍᵘ"],:α=>(x,y,z)->1e12)
prescribe!(elements["Γᵍᵘ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍᵘ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₁₂=>(x,y,z)->0.0)

𝑎ᵘ = ∫∫εᵈᵢⱼσᵈᵢⱼdxdy=>elements["Ωᵘ"]
𝑎ᵖ = ∫∫qpdxdy=>elements["Ωᵖ"]
𝑏ᵖ = ∫∫p∇udxdy=>(elements["Ωᵖ"],elements["Ωᵘ"])
𝑎ᵘᵅ = ∫vᵢgᵢds=>elements["Γᵍᵘ"]
𝑓 = ∫vᵢtᵢds=>elements["Γᵗ"]
# 𝑓 = [
#     ∫∫vᵢbᵢdxdy=>elements["Ωᵘ"],
#     ∫vᵢtᵢds=>elements["Γᵗ"],
# ]

kᵘᵘ = zeros(2*nᵤ,2*nᵤ)
kᵖᵖ = zeros(nₚ,nₚ)
kᵖᵘ = zeros(nₚ,2*nᵤ)
fᵖ = zeros(nₚ)
fᵘ = zeros(2*nᵤ)

@timeit to "assembly" begin
𝑎ᵘ(kᵘᵘ)
𝑎ᵖ(kᵖᵖ)
𝑏ᵖ(kᵖᵘ)
𝑎ᵘᵅ(kᵘᵘ,fᵘ)
𝑓(fᵘ)
end
k =sparse([-kᵘᵘ kᵖᵘ';kᵖᵘ kᵖᵖ])
f = [-fᵘ;fᵖ]
d = zeros(2*nᵤ+nₚ)

set_matrixtype!(ps, -2)
k = get_matrix(ps,k,:N)
@timeit to "solve" pardiso(ps,d,k,f)

𝑢₁ = d[1:2:2*nᵤ]
𝑢₂ = d[2:2:2*nᵤ]
𝑝 = d[2*nᵤ+1:2*nᵤ+nₚ]

# @timeit to "plot figure" begin
# fig = Figure()
# ind = 100
# ax = Axis(fig[1,1], 
#     aspect = DataAspect(), 
#     xticksvisible = false,
#     xticklabelsvisible=false, 
#     yticksvisible = false, 
#     yticklabelsvisible=false,
# )
# hidespines!(ax)
# hidedecorations!(ax)
# xs = LinRange(0, 48, 4*ind)
# ys = LinRange(-6, 6, ind)
# zs = zeros(4*ind,ind)
# 𝗠 = zeros(21)
# for (i,x) in enumerate(xs)
#     for (j,y) in enumerate(ys)
#         indices = sp(x,y,0.0)
#         ni = length(indices)
#         𝓒 = [nodes_p[i] for i in indices]
#         data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
#         ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
#         𝓖 = [ξ]
#         a = type(𝓒,𝓖)
#         set𝝭!(a)
#         p = 0.0
#         N = ξ[:𝝭]
#         for (k,xₖ) in enumerate(𝓒)
#             p += N[k]*xₖ.p
#         end
#         zs[i,j] = p
#     end
# end
# surface!(xs,ys,zeros(4*ind,ind),color=zs,shading=NoShading,colormap=:lightrainbow)
# contour!(xs,ys,zs,levels=-1e3:200:1e3,color=:azure)
# Colorbar(fig[1,2], limits=(-900,900), colormap=:lightrainbow)
# save("./png/cantilever_mix_"*poly*"_"*string(ndiv)*"_"*string(nₚ)*".png",fig, px_per_unit = 10.0)
# end

show(to)
# fig