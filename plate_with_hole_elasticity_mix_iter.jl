using Revise
using TimerOutputs 
using SparseArrays, Pardiso, XLSX
using CairoMakie
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫qpdxdy, ∫∫sᵢⱼsᵢⱼdxdy, ∫∫p∇udxdy, ∫∫sᵢⱼεᵢⱼdxdy, ∫pnᵢgᵢds, ∫sᵢⱼnⱼgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric

include("import_plate_with_hole.jl")

ndiv = 8

indices = 1:16
nₜ = length(indices)
L₂_𝒖   = zeros(nₜ)
Hₑ_𝒖   = zeros(nₜ)
Hₑ_dev = zeros(nₜ)
L₂_𝑝   = zeros(nₜ)

for (i,n) in enumerate(indices)
ps = MKLPardisoSolver()

elements, nodes, nodes_p = import_elasticity_linear_mix("./msh/plate_with_hole_tri3_"*string(ndiv)*".msh","./msh/plate_with_hole_tri3_"*string(n)*".msh",n)
# elements, nodes, nodes_p = import_elasticity_quadratic_mix("./msh/plate_with_hole_tri6_"*string(ndiv)*".msh","./msh/plate_with_hole_tri3_"*string(n)*".msh",n)
# nx = 131;ny = 32
# elements, nodes, nodes_p = import_linear_mix("./msh/cantilever_"*poly*"_"*string(ndiv)*".msh","./msh/cantilever_"*string(ny)*"_"*string(nx)*".msh",nx,ny)

nₚ = length(nodes_p)

nₑ = length(elements["Ωᵘ"])
nₛ = 1
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

# n = 1
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

# ∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
# ∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
# ∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
# ∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
# ∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
# ∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))

# ∂σ₁₁∂x(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
# ∂σ₁₁∂y(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
# ∂σ₂₂∂x(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y) + (1-ν)*∂ε₂₂∂x(x,y))
# ∂σ₂₂∂y(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y) + (1-ν)*∂ε₂₂∂y(x,y))
# ∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
# ∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
# b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
# b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)
p(x,y) = (σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3

prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν)
# prescribe!(elements["Ωᵘ"],:b₁=>(x,y,z)->b₁(x,y))
# prescribe!(elements["Ωᵘ"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂) 
prescribe!(elements["Γᵍᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍᵘ"],:g₂=>(x,y,z)->v(x,y))
# prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
# prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)
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

𝑎ˢ = ∫∫sᵢⱼsᵢⱼdxdy=>elements["Ωˢ"]
𝑎ᵖ = ∫∫qpdxdy=>elements["Ωᵖ"]
𝑏ˢ = ∫∫sᵢⱼεᵢⱼdxdy=>(elements["Ωˢ"],elements["Ωᵘ"])
𝑏ᵖ = ∫∫p∇udxdy=>(elements["Ωᵖ"],elements["Ωᵘ"])
𝑏ˢᵅ = ∫sᵢⱼnⱼgᵢds=>(elements["Γᵍˢ"],elements["Γᵍᵘ"])
𝑏ᵖᵅ = ∫pnᵢgᵢds=>(elements["Γᵍᵖ"],elements["Γᵍᵘ"])
𝑓 = ∫vᵢtᵢds=>elements["Γᵗ"]
# 𝑓 = [
#     ∫vᵢtᵢds=>elements["Γᵗ"],
#     ∫∫vᵢbᵢdxdy=>elements["Ωᵘ"]
# ]

kˢˢ = zeros(4*nₛ*nₑ,4*nₛ*nₑ)
kᵖᵖ = zeros(nₚ,nₚ)
kˢᵘ = zeros(4*nₛ*nₑ,2*nᵤ)
kᵖᵘ = zeros(nₚ,2*nᵤ)
fˢ = zeros(4*nₛ*nₑ)
fᵖ = zeros(nₚ)
fᵘ = zeros(2*nᵤ)

𝑎ˢ(kˢˢ)
𝑎ᵖ(kᵖᵖ)
𝑏ˢ(kˢᵘ)
𝑏ᵖ(kᵖᵘ)
𝑏ˢᵅ(kˢᵘ,fˢ)
𝑏ᵖᵅ(kᵖᵘ,fᵖ)
𝑓(fᵘ)
# k = [zeros(2*nᵤ,2*nᵤ) kᵖᵘ' kˢᵘ';kᵖᵘ kᵖᵖ zeros(nₚ,4*nₛ*nₑ);kˢᵘ zeros(4*nₛ*nₑ,nₚ) kˢˢ]
k = sparse([zeros(2*nᵤ,2*nᵤ) kᵖᵘ' kˢᵘ';kᵖᵘ kᵖᵖ zeros(nₚ,4*nₛ*nₑ);kˢᵘ zeros(4*nₛ*nₑ,nₚ) kˢˢ])
f = [-fᵘ;fᵖ;fˢ]
d = zeros(2*nᵤ+nₚ+4*nₛ*nₑ)
# d = k\f

set_matrixtype!(ps, -2)
k = get_matrix(ps,k,:N)
pardiso(ps,d,k,f)

𝑢₁ = d[1:2:2*nᵤ]
𝑢₂ = d[2:2:2*nᵤ]
𝑝 = d[2*nᵤ+1:2*nᵤ+nₚ]
push!(nodes,:d₁=>𝑢₁)
push!(nodes,:d₂=>𝑢₂)
push!(nodes_p,:p=>𝑝)

Hₑ_𝒖_, L₂_𝒖_ = Hₑ_PlaneStress(elements["Ωᵍᵘ"])
Hₑ_dev_ = Hₑ_PlaneStrain_Deviatoric(elements["Ωᵍᵘ"])
L₂_𝑝_ = L₂𝑝(elements["Ωᵍᵖ"])

L₂_𝒖[i] = log10(L₂_𝒖_)
Hₑ_𝒖[i] = log10(Hₑ_𝒖_)
Hₑ_dev[i] = log10(Hₑ_dev_)
L₂_𝑝[i] = log10(L₂_𝑝_)

println("n = $n, L₂_𝒖 = $L₂_𝒖_, Hₑ_𝒖 = $Hₑ_𝒖_, Hₑ_dev = $Hₑ_dev_, L₂_𝑝 = $L₂_𝑝_")

end


XLSX.openxlsx("./xlsx/plate_with_hole_linear_mix.xlsx", mode = "rw") do xf
    sheet = xf[1]
    row = "A"
    row_L₂_𝒖 = "B"
    row_Hₑ_𝒖 = "C"
    row_Hₑ_dev = "D"
    row_L₂_𝑝 = "E"
    for (n,L₂_𝒖_,Hₑ_𝒖_,Hₑ_dev_,L₂_𝑝_) in zip(indices,L₂_𝒖,Hₑ_𝒖,Hₑ_dev,L₂_𝑝)
        sheet[row*string(n)] = n
        sheet[row_L₂_𝒖*string(n)] = L₂_𝒖_
        sheet[row_Hₑ_𝒖*string(n)] = Hₑ_𝒖_
        sheet[row_Hₑ_dev*string(n)] = Hₑ_dev_
        sheet[row_L₂_𝑝*string(n)] = L₂_𝑝_
    end
end