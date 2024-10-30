
using ApproxOperator, LinearAlgebra, XLSX, SparseArrays
using ApproxOperator.Elasticity: ∫∫εᵈᵢⱼσᵈᵢⱼdxdy, ∫∫qpdxdy, ∫∫p∇udxdy, ∫vᵢgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝

include("import_patchtest.jl")

ndiv = 16

indices = 2:64

n_eig_nonzeros = zeros(Int,length(indices))
n_eig_real = zeros(Int,length(indices))
min_eig_nonzeros = zeros(length(indices))
min_eig_real = zeros(length(indices))

for (i,n) in enumerate(indices)

# elements, nodes, nodes_p = import_infsup_linear_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(n)*".msh",n,n)
# elements, nodes, nodes_p = import_infsup_linear_mix("./msh/patchtest_quad_"*string(ndiv)*".msh","./msh/patchtest_"*string(n)*".msh",n,n)
# elements, nodes, nodes_p = import_infsup_quadratic_mix("./msh/patchtest_tri6_"*string(ndiv)*".msh","./msh/patchtest_"*string(n)*".msh",n,n)
elements, nodes, nodes_p = import_infsup_quadratic_mix("./msh/patchtest_quad8_"*string(ndiv)*".msh","./msh/patchtest_"*string(n)*".msh",n,n)
# nx = 2;ny = n;
# elements, nodes, nodes_p = import_infsup_linear_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(nx)*"_"*string(ny)*".msh",nx,ny)

nₑ = length(elements["Ωᵘ"])
nᵤ = length(nodes)
nₚ = length(nodes_p)

E = 1.0
ν = 0.3
# ν = 0.5-1e3

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

# ε₁₁(x,y) = ∂u∂x(x,y)
# ε₂₂(x,y) = ∂v∂y(x,y)
# ε₁₂(x,y) = 0.5*(∂u∂y(x,y) + ∂v∂x(x,y))
# σ₁₁(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*ε₁₁(x,y) + ν*ε₂₂(x,y))
# σ₂₂(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + (1-ν)*ε₂₂(x,y))
# σ₃₃(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + ν*ε₂₂(x,y))
# σ₁₂(x,y) = E/(1+ν)*ε₁₂(x,y)
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
# p(x,y) = (σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3

prescribe!(elements["Ωᵘ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵘ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν, index=:𝑔)
# prescribe!(elements["Ωᵘ"],:b₁=>(x,y,z)->b₁(x,y))
# prescribe!(elements["Ωᵘ"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Γ⁴ᵘ"],:α=>(x,y,z)->1e10)
# prescribe!(elements["Γ¹ᵘ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
# prescribe!(elements["Γ¹ᵘ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
# prescribe!(elements["Γ²ᵘ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
# prescribe!(elements["Γ²ᵘ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
# prescribe!(elements["Γ³ᵘ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
# prescribe!(elements["Γ³ᵘ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γ⁴ᵘ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ⁴ᵘ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ⁴ᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴ᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴ᵘ"],:n₁₂=>(x,y,z)->0.0)
# prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->u(x,y))
# prescribe!(elements["Ωᵍᵘ"],:v=>(x,y,z)->v(x,y))
# prescribe!(elements["Ωᵍᵖ"],:p=>(x,y,z)->p(x,y))
# elements["Γᵗ"] = elements["Γ¹ᵘ"]∪elements["Γ²ᵘ"]∪elements["Γ³ᵘ"]

𝑎ᵘ = ∫∫εᵈᵢⱼσᵈᵢⱼdxdy=>elements["Ωᵘ"]
𝑎ᵖ = ∫∫qpdxdy=>elements["Ωᵖ"]
𝑏ᵖ = ∫∫p∇udxdy=>(elements["Ωᵖ"],elements["Ωᵘ"])
𝑎ᵘᵅ = ∫vᵢgᵢds=>elements["Γ⁴ᵘ"]
# 𝑓 = [
#     ∫∫vᵢbᵢdxdy=>elements["Ωᵘ"],
#     ∫vᵢtᵢds=>elements["Γᵗ"]
# ]

kᵘᵘ = zeros(2*nᵤ,2*nᵤ)
kᵖᵖ = zeros(nₚ,nₚ)
kᵖᵘ = zeros(nₚ,2*nᵤ)
fᵖ = zeros(nₚ)
fᵘ = zeros(2*nᵤ)

𝑎ᵘ(kᵘᵘ)
𝑎ᵖ(kᵖᵖ)
𝑏ᵖ(kᵖᵘ)
𝑎ᵘᵅ(kᵘᵘ,fᵘ)

# d = [-kᵘᵘ kᵖᵘ';kᵖᵘ kᵖᵖ]\[-fᵘ;fᵖ]

# 𝑢₁ = d[1:2:2*nᵤ]
# 𝑢₂ = d[2:2:2*nᵤ]
# 𝑝 = d[2*nᵤ+1:2*nᵤ+nₚ]
# push!(nodes,:d₁=>𝑢₁)
# push!(nodes,:d₂=>𝑢₂)
# push!(nodes_p,:p=>𝑝)

# L₂u = L₂(elements["Ωᵍᵘ"])
# L₂p = L₂𝑝(elements["Ωᵍᵖ"])

# println(L₂u)
# println(L₂p)

# kᵈ = kᵘᵘ
# kᵛ = kᵖᵘ\kᵖᵖ*kᵖᵘ
# val = eigvals(kᵛ,kᵈ)
val = eigvals(kᵖᵘ'*(kᵖᵖ\kᵖᵘ),kᵘᵘ)

# println(2*nᵤ-nₚ+1)
# println("Unsorted Eigenvalue")
# println.(val[2*nᵤ-nₚ.+(-2:4)]);

val_sign = zeros(2*nᵤ)
for (ii,v) in enumerate(val)
    if v isa Real
        val_sign[ii] = sign(v)
    else
        val_sign[ii] = sign(v.re) < -1e-8 ? -1.0 : 1.0
    end
end
val_real = val_sign .* abs.(val)
val_abs = abs.(val)
# println("Sorted Eigenvalue")
val_sort = sort(val_abs)
# println.(val_sort[2*nᵤ-nₚ.+(-2:4)]);

n_eig_real[i] = count(x-> abs(x)>1e-8, val_real)
n_eig_nonzeros[i] = count(x-> x > 1e-8,val_sort)
min_eig_real[i] = min(val_real[val_real.>1e-8]...)
min_eig_nonzeros[i] = val_sort[2*nᵤ - n_eig_nonzeros[i] + 1]

end

XLSX.openxlsx("./xlsx/infsup.xlsx", mode = "rw") do xf
    sheet = xf[1]
    for (n,n_eig_r,min_eig_r,n_eig_n,min_eig_n) in zip(indices,n_eig_real,min_eig_real,n_eig_nonzeros,min_eig_nonzeros)
        sheet["A"*string(n)] = n
        sheet["B"*string(n)] = n_eig_r
        sheet["C"*string(n)] = min_eig_r^0.5
        sheet["D"*string(n)] = n_eig_n
        sheet["E"*string(n)] = min_eig_n^0.5
    end
end