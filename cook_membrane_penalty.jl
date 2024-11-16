using Revise
using TimerOutputs 
using SparseArrays, Pardiso
using CairoMakie, XLSX, WriteVTK
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵈᵢⱼσᵈᵢⱼdxdy, ∫∫qpdxdy, ∫∫p∇udxdy, ∫vᵢgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric

include("import_cook.jl")

const to = TimerOutput()
ps = MKLPardisoSolver()

ndiv = 4
poly = "quad"
@timeit to "import data" begin
n = 4
elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cook_"*poly*"_"*string(ndiv)*".msh","./msh/cook_tri3_"*string(n)*".msh",n)
# elements, nodes, nodes_p, sp, type = import_quadratic_mix("./msh/cook_"*poly*"_"*string(ndiv)*".msh","./msh/cook_quad_"*string(n)*".msh",n)
nₚ = length(nodes_p)
end

nₑ = length(elements["Ωᵘ"])
# nₛ = 3
nᵤ = length(nodes)

# T = 1.0e3
# E = 3.0e6
P = 6.25
E = 70.0
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
push!(nodes,:u₁=>𝑢₁,:u₂=>𝑢₂)
push!(nodes_p,:p=>𝑝)

# @timeit to "plot figure" begin
# fig = Figure(figure_padding = 1,size = (400,600))
# ind = 100
# ax = Axis(fig[1,1], 
#     aspect = DataAspect(), 
#     xticksvisible = true,
#     xticklabelsvisible=false, 
#     yticksvisible = false, 
#     yticklabelsvisible=false,
#     backgroundcolor = :transparent,
# )
# hidespines!(ax)
# hidedecorations!(ax)
# index = [1,2,3,1]
# α = 1.0
# for elm in elements["Ωᵘ"]
#     x = [node.x+α*node.u₁ for node in elm.𝓒[index]]
#     y = [node.y+α*node.u₂ for node in elm.𝓒[index]]
#     lines!(ax,x,y,color=:black,linewidth = 3)
# end
# vertices = [[node.x+α*node.u₁ for node in nodes] [node.y+α*node.u₂ for node in nodes]]
colors = zeros(nᵤ)
𝗠 = zeros(21)
for (i,node) in enumerate(nodes)
    x = node.x
    y = node.y
    indices = sp(x,y,0.0)
    ni = length(indices)
    𝓒 = [nodes_p[i] for i in indices]
    data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
    ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
    𝓖 = [ξ]
    a = type(𝓒,𝓖)
    set𝝭!(a)
    p = 0.0
    N = ξ[:𝝭]
    for (k,xₖ) in enumerate(𝓒)
        p += N[k]*xₖ.p
    end
    colors[i] = p
end
# faces = zeros(Int,nₑ,3)
# for (e,elm) in enumerate(elements["Ωᵘ"])
#     faces[e,:] .= [xᵢ.𝐼 for xᵢ in elm.𝓒[1:3]]
# end
# mesh!(vertices,faces,color=colors,shading = NoShading,colormap=:haline,colorrange = (-50,15))

# coord = [[node.x+α*node.u₁ for node in nodes] [node.y+α*node.u₂ for node in nodes]]

# x = [node.x+α*node.u₁ for node in nodes]
# y = [node.y+α*node.u₂ for node in nodes]
# tricontourf!(ax,x,y,colors,levels=collect(-60:5:20), colormap=Reverse(:deep))
# surface!(xs,ys,zeros(4*ind,ind),color=zs,shading=NoShading,colormap=:lightrainbow)
# contour!(xs,ys,zs,levels=-1e3:200:1e3,color=:azure)
# Colorbar(fig[1,2], limits=(-50,15), colormap=:haline)
# Colorbar(fig[1,2], colormap=:haline)
# save("./png/cook_mix_"*poly*"_"*string(ndiv)*"_"*string(nₚ)*".png",fig, px_per_unit = 10.0)
# end

# XLSX.openxlsx("./xlsx/contour.xlsx", mode = "rw") do xf
#     sheet = xf[1]
#     for (i,node) in enumerate(nodes)
#         x = node.x
#         y = node.y
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
#         sheet["A"*string(i)] = round(x + α*node.u₁, digits=2)
#         sheet["B"*string(i)] = round(y + α*node.u₂, digits=2)
#         sheet["C"*string(i)] = round(p, digits=2)
#     end
# end

# points = [[node.x+α*node.u₁ for node in nodes]';[node.y+α*node.u₂ for node in nodes]';zeros(1,nᵤ)]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TRIANGLE,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# cells = [MeshCell(VTKCellTypes.VTK_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# vtk_grid("./vtk/cook_"*poly*"_"*string(ndiv)*"_"*string(nₚ),points,cells) do vtk
#     vtk["𝑝"] = colors
# end
show(to)
println(𝑢₂[3])
# fig