using Revise
using TimerOutputs 
using SparseArrays, Pardiso
using CairoMakie
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵈᵢⱼσᵈᵢⱼdxdy, ∫∫qpdxdy, ∫∫p∇udxdy, ∫vᵢgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric

include("import_plate_with_hole.jl")

const to = TimerOutput()
ps = MKLPardisoSolver()

ndiv = 2
@timeit to "import data" begin
n = 2
# elements, nodes, nodes_p = import_elasticity_linear_mix("./msh/plate_with_hole_tri3_"*string(ndiv)*".msh","./msh/plate_with_hole_tri3_"*string(n)*".msh",n)
# nx = 7;ny = 3
# elements, nodes, nodes_p = import_elasticity_linear_mix("./msh/plate_with_hole_tri3_"*string(ndiv)*".msh","./msh/plate_with_hole_tri3_"*string(ny)*"_"*string(nx)*".msh",ny)
elements, nodes, nodes_p = import_elasticity_quadratic_mix("./msh/plate_with_hole_tri6_"*string(ndiv)*".msh","./msh/plate_with_hole_tri6_"*string(n)*".msh",n)
# nx = 68;ny = 32
# elements, nodes, nodes_p = import_elasticity_quadratic_mix("./msh/plate_with_hole_tri6_"*string(ndiv)*".msh","./msh/plate_with_hole_tri3_"*string(ny)*"_"*string(nx)*".msh",ny)

nₚ = length(nodes_p)
end

nₑ = length(elements["Ωᵘ"])
# nₛ = 3
nᵤ = length(nodes)

# T = 1.0e3
# E = 3.0e6
T = 1.0
E = 1.0e4
# ν = 0.3
ν = 0.5-1e-8
Ē = E/(1.0-ν^2)
ν̄ = ν/(1.0-ν)
Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
Cᵢⱼᵢⱼ = E/(1+ν)/2
a = 1
b = 5

# n = 2
# u(x,y) = (1+2*x+3*y)^n
# v(x,y) = (4+5*x+6*y)^n
# ∂u∂x(x,y) = 2*n*(1+2*x+3*y)^abs(n-1)
# ∂u∂y(x,y) = 3*n*(1+2*x+3*y)^abs(n-1)
# ∂v∂x(x,y) = 5*n*(4+5*x+6*y)^abs(n-1)
# ∂v∂y(x,y) = 6*n*(4+5*x+6*y)^abs(n-1)
# ∂²u∂x²(x,y)  = 4*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²u∂x∂y(x,y) = 6*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²u∂y²(x,y)  = 9*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²v∂x²(x,y)  = 25*n*(n-1)*(4+5*x+6*y)^abs(n-2)
# ∂²v∂x∂y(x,y) = 30*n*(n-1)*(4+5*x+6*y)^abs(n-2)
# ∂²v∂y²(x,y)  = 36*n*(n-1)*(4+5*x+6*y)^abs(n-2)

r(x,y) = (x^2+y^2)^0.5
θ(x,y) = atan(y/x)
u(x,y) = T*a*(1+ν̄)/2/Ē*(r(x,y)/a*2/(1+ν̄)*cos(θ(x,y)) + a/r(x,y)*(4/(1+ν̄)*cos(θ(x,y))+cos(3*θ(x,y))) - a^3/r(x,y)^3*cos(3*θ(x,y)))
v(x,y) = T*a*(1+ν̄)/2/Ē*( -r(x,y)/a*2*ν̄/(1+ν̄)*sin(θ(x,y)) - a/r(x,y)*(2*(1-ν̄)/(1+ν̄)*sin(θ(x,y))-sin(3*θ(x,y))) - a^3/r(x,y)^3*sin(3*θ(x,y)) )
∂u∂x(x,y) = T/Ē*(1 + a^2/2/r(x,y)^2*((ν̄-3)*cos(2*θ(x,y))-2*(1+ν̄)*cos(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν̄)*cos(4*θ(x,y)))
∂u∂y(x,y) = T/Ē*(-a^2/r(x,y)^2*((ν̄+5)/2*sin(2*θ(x,y))+(1+ν̄)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν̄)*sin(4*θ(x,y)))
∂v∂x(x,y) = T/Ē*(-a^2/r(x,y)^2*((ν̄-3)/2*sin(2*θ(x,y))+(1+ν̄)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν̄)*sin(4*θ(x,y)))
∂v∂y(x,y) = T/Ē*(-ν̄ - a^2/2/r(x,y)^2*((1-3*ν̄)*cos(2*θ(x,y))-2*(1+ν̄)*cos(4*θ(x,y))) - 3*a^4/2/r(x,y)^4*(1+ν̄)*cos(4*θ(x,y)))

ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = 0.5*(∂u∂y(x,y) + ∂v∂x(x,y))
σ₁₁(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*ε₁₁(x,y) + ν*ε₂₂(x,y))
σ₂₂(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + (1-ν)*ε₂₂(x,y))
σ₃₃(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + ν*ε₂₂(x,y))
σ₁₂(x,y) = E/(1+ν)*ε₁₂(x,y)

∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))

∂σ₁₁∂x(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
∂σ₁₁∂y(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
∂σ₂₂∂x(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y) + (1-ν)*∂ε₂₂∂x(x,y))
∂σ₂₂∂y(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y) + (1-ν)*∂ε₂₂∂y(x,y))
∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)
p(x,y) = (σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3

prescribe!(elements["Ωᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵘ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵘ"],:b₁=>(x,y,z)->b₁(x,y))
prescribe!(elements["Ωᵘ"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂) 
# prescribe!(elements["Γᵗ"],:α=>(x,y,z)->1e12)
# prescribe!(elements["Γᵗ"],:g₁=>(x,y,z)->u(x,y)) 
# prescribe!(elements["Γᵗ"],:g₂=>(x,y,z)->v(x,y)) 
# prescribe!(elements["Γᵗ"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
# prescribe!(elements["Γᵗ"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)
# prescribe!(elements["Γᵗ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍᵘ"],:α=>(x,y,z)->1e12)
prescribe!(elements["Γᵍᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍᵘ"],:g₂=>(x,y,z)->v(x,y))
# prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z)->1.0)
# prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z,n₁,n₂)->abs(n₁))
prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z,n₁,n₂)->abs(n₂))
prescribe!(elements["Γᵍᵘ"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍᵘ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
prescribe!(elements["Ωᵍᵖ"],:p=>(x,y,z)->p(x,y))

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
# d = k\f

set_matrixtype!(ps, -2)
k = get_matrix(ps,k,:N)
@timeit to "solve" pardiso(ps,d,k,f)

𝑢₁ = d[1:2:2*nᵤ]
𝑢₂ = d[2:2:2*nᵤ]
𝑝 = d[2*nᵤ+1:2*nᵤ+nₚ]
push!(nodes,:d₁=>𝑢₁)
push!(nodes,:d₂=>𝑢₂)
push!(nodes_p,:p=>𝑝)

@timeit to "compute error" begin
Hₑ_𝒖, L₂_𝒖 = Hₑ_PlaneStress(elements["Ωᵍᵘ"])
Hₑ_dev = Hₑ_PlaneStrain_Deviatoric(elements["Ωᵍᵘ"])
L₂_𝑝 = L₂𝑝(elements["Ωᵍᵖ"])
end

println(log10(L₂_𝒖))
println(log10(Hₑ_𝒖))
println(log10(Hₑ_dev))
println(log10(L₂_𝑝))

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