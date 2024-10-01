
using BenchmarkTools
using SparseArrays, Pardiso
using ApproxOperator
using ApproxOperator.Heat: ∫∫qᵢpᵢdxdy, ∫pᵢnᵢuds, ∫∫∇𝒑udxdy, ∫pᵢnᵢgⱼds, ∫vtdΓ, ∫vgdΓ, L₂, L₂𝒑, H₁

include("import_plate_with_hole.jl")

ps = MKLPardisoSolver()

ndiv = 32
poly = "tri3"
# n = 8
L₂_𝒑 = zeros(34)
H₁_𝑢 = zeros(34)
L₂_𝑢 = zeros(34)

# for n in 9:42
# elements, nodes, nodes_u = import_linear_mix("./msh/plate_with_hole_"*poly*"_"*string(ndiv)*".msh","./msh/plate_with_hole_"*poly*"_"*string(n)*".msh",n)
# # n₂ = 32
# # n₁ = 68
# # elements, nodes, nodes_u = import_linear_mix("./msh/plate_with_hole_"*poly*"_"*string(ndiv)*".msh","./msh/plate_with_hole_"*poly*"_"*string(n₂)*"_"*string(n₁)*".msh",n₂)
# # elements, nodes, nodes_u = import_linear_mix("./msh/plate_with_hole_"*poly*"_"*string(ndiv)*".msh","./msh/plate_with_hole_"*poly*"_"*string(ndiv)*".msh")

# nₚ = length(nodes)
# nᵤ = length(nodes_u)

# # n = 1
# # u(x,y) = (x+y)^n
# # ∂u∂x(x,y) = n*(x+y)^abs(n-1)
# # ∂u∂y(x,y) = n*(x+y)^abs(n-1)
# # ∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# # ∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
# # ∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# # b(x,y,z) = -∂²u∂x²(x,y)-∂²u∂y²(x,y)

# r(x,y) = (x^2+y^2)^0.5
# θ(x,y) = atan(y/x)
# u(x,y) = (r(x,y) + 1/r(x,y))cos(θ(x,y))
# ∂u∂x(x,y) = 1 - 1/r(x,y)^2 + 2*1/r(x,y)^2 * sin(θ(x,y))^2
# ∂u∂y(x,y) = - 2/r(x,y)^2 * sin(θ(x,y))*cos(θ(x,y))

# prescribe!(elements["Γᵍᵘ"],:α=>(x,y,z)->1e9)
# prescribe!(elements["Γᵍᵘ"],:g=>(x,y,z)->u(x,y))
# prescribe!(elements["Γᵗ"],:t=>(x,y,z,n₁,n₂)->∂u∂x(x,y)*n₁ + ∂u∂y(x,y)*n₂)
# prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->u(x,y))
# prescribe!(elements["Ωᵍᵘ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
# prescribe!(elements["Ωᵍᵘ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
# prescribe!(elements["Ωᵍᵘ"],:∂u∂z=>(x,y,z)->0.0)
# prescribe!(elements["Ωᵍᵖ"],:𝑝₁=>(x,y,z)->∂u∂x(x,y))
# prescribe!(elements["Ωᵍᵖ"],:𝑝₂=>(x,y,z)->∂u∂y(x,y))
# prescribe!(elements["Ωᵍᵖ"],:𝑝₃=>(x,y,z)->0.0)

# 𝑎 = ∫∫qᵢpᵢdxdy=>elements["Ωᵖ"]
# 𝑏 = [
#     ∫pᵢnᵢuds=>(elements["∂Ωᵖ"],elements["∂Ωᵘ"]),
#     ∫∫∇𝒑udxdy=>(elements["Ωᵖ"],elements["Ωᵘ"]),
# ]
# 𝑏ᵅ = ∫pᵢnᵢgⱼds=>(elements["Γᵍᵖ"],elements["Γᵍᵘ"])
# 𝑓 = ∫vtdΓ=>elements["Γᵗ"]

# kᵖᵖ = zeros(2*nₚ,2*nₚ)
# fᵖ = zeros(2*nₚ)
# kᵖᵘ = zeros(2*nₚ,nᵤ)
# fᵘ = zeros(nᵤ)

# 𝑎(kᵖᵖ)
# 𝑏(kᵖᵘ)
# 𝑏ᵅ(kᵖᵘ,fᵖ)
# 𝑓(fᵘ)

# k = sparse([kᵖᵖ kᵖᵘ;kᵖᵘ' zeros(nᵤ,nᵤ)])
# f = [fᵖ;-fᵘ]
# d = zeros(nᵤ+2*nₚ)
# set_matrixtype!(ps, -2)
# k = get_matrix(ps,k,:N)
# pardiso(ps,d,k,f)

# 𝑝₁ = d[1:2:2*nₚ]
# 𝑝₂ = d[2:2:2*nₚ]
# 𝑢 = d[2*nₚ+1:end]
# push!(nodes,:p₁=>𝑝₁)
# push!(nodes,:p₂=>𝑝₂)
# push!(nodes,:p₃=>zeros(nₚ))
# push!(nodes_u,:d=>𝑢)

# L₂_𝒑_ = L₂𝒑(elements["Ωᵍᵖ"])
# H₁_𝑢_, L₂_𝑢_ = H₁(elements["Ωᵍᵘ"])

# L₂_𝒑_ = log10(L₂_𝒑_)
# H₁_𝑢_ = log10(H₁_𝑢_)
# L₂_𝑢_ = log10(L₂_𝑢_)
# println("n = $n, L₂_𝒑 = $L₂_𝒑_, H₁_𝑢 = $H₁_𝑢_,L₂_𝑢 = $L₂_𝑢_")

# end

XLSX.openxlsx("./xlsx/plate_with_hole_heat.xlsx", mode = "rw") do xf
    sheet = xf[1]
    row = "A"
    row_L₂_𝑢 = "B"
    row_H₁_𝑢 = "C"
    row_L₂_𝒑 = "D"
    for (n,L₂_𝑢_,H₁_𝑢_,L₂_𝒑_) in zip(9:42,L₂_𝑢,H₁_𝑢,L₂_𝒑)
        sheet[row*string(n)] = n
        sheet[row_L₂_𝑢*string(n)] = L₂_𝑢_
        sheet[row_H₁_𝑢*string(n)] = H₁_𝑢_
        sheet[row_L₂_𝒑*string(n)] = L₂_𝒑_
    end
end