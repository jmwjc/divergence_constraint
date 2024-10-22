
using ApproxOperator, LinearAlgebra, XLSX
using ApproxOperator.Elasticity: ∫∫εᵈᵢⱼσᵈᵢⱼdxdy, ∫∫qpdxdy, ∫∫p∇udxdy, ∫vᵢgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝

include("import_patchtest.jl")

ndiv = 16

n = 36

# elements, nodes, nodes_p = import_patchtest_elasticity_penalty("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_c_"*string(nₚ)*".msh")
# elements, nodes, nodes_p = import_infsup_linear_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(n)*".msh",n,n)
elements, nodes, nodes_p = import_infsup_quadratic_mix("./msh/patchtest_tri6_"*string(ndiv)*".msh","./msh/patchtest_"*string(n)*".msh",n,n)

nₑ = length(elements["Ωᵘ"])
nᵤ = length(nodes)
nₚ = length(nodes_p)

E = 1.0e0
ν = 0.3
# ν = 0.5-1e1

prescribe!(elements["Ωᵘ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵘ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Γ⁴ᵘ"],:α=>(x,y,z)->1e10)
prescribe!(elements["Γ⁴ᵘ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ⁴ᵘ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ⁴ᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴ᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴ᵘ"],:n₁₂=>(x,y,z)->0.0)

𝑎ᵘ = ∫∫εᵈᵢⱼσᵈᵢⱼdxdy=>elements["Ωᵘ"]
𝑎ᵖ = ∫∫qpdxdy=>elements["Ωᵖ"]
𝑏ᵖ = ∫∫p∇udxdy=>(elements["Ωᵖ"],elements["Ωᵘ"])
𝑎ᵘᵅ = ∫vᵢgᵢds=>elements["Γ⁴ᵘ"]

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

val = eigvals(kᵖᵘ'*(kᵖᵖ\kᵖᵘ),kᵘᵘ)
# val = eigvals(kᵖᵘ'*inv(kᵖᵖ)*kᵖᵘ,kᵘᵘ)

# println(2*nᵤ-nₚ+1)
# println("Unsorted Eigenvalue")
# println.(val[2*nᵤ-nₚ.+(-2:4)]);

val_sign = zeros(2*nᵤ)
for (ii,v) in enumerate(val)
    if v isa Real
        val_sign[ii] = sign(v)
    else
        val_sign[ii] = sign(v.re) < -1e-10 ? -1.0 : 1.0
    end
end
val_real = val_sign .* abs.(val)
val_abs = abs.(val)
# println("Sorted Eigenvalue")
val_sort = sort(val_abs)
# println.(val_sort[2*nᵤ-nₚ.+(-2:4)]);

n_eig_positive = count(x-> isa(x,Real) ? x > 1e-10 : x.re > 1e-10,val)
n_eig_nonzeros = count(x-> x > 1e-10,val_sort)
min_eig_positive = isa(val[2*nᵤ - n_eig_positive + 1],Real) ? val[2*nᵤ - n_eig_positive + 1] : val[2*nᵤ - n_eig_positive + 1].re
min_eig_nonzeros = val_sort[2*nᵤ - n_eig_nonzeros + 1]
min_eig_real = min(val_real[abs.(val_real).>1e-10]...)