
using Gmsh, Statistics

function import_fem(filename::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    elements["Ω"] = getElements(nodes,entities["Ω"])
    elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"],normal=true)
    elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"],normal=true)

    push!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    push!(elements["Γᵗ"],:𝝭)
    push!(elements["Γᵍ"],:𝝭)

    gmsh.finalize()

    set∇𝝭!(elements["Ω"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵍ"])

    return elements, nodes
end

function import_linear_mix(filename1::String,filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_u = get𝑿ᵢ()
    xᵘ = nodes_u.x
    yᵘ = nodes_u.y
    zᵘ = nodes_u.z
    Ω = getElements(nodes_u, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 1.5*s*ones(length(nodes_u))
    push!(nodes_u,:s₁=>s,:s₂=>s,:s₃=>s)

    integrationOrder_Ω = 2
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 2

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ωᵖ"] = getElements(nodes,entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍᵖ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["∂Ωᵖ"] = getElements(nodes, entities["Γ"],   integrationOrder_Γ, normal = true)
    elements["Γ¹ᵖ"] = getElements(nodes,entities["Γ¹"], integrationOrder_Γ, normal = true)
    elements["Γ²ᵖ"] = getElements(nodes,entities["Γ²"], integrationOrder_Γ, normal = true)
    elements["Γ³ᵖ"] = getElements(nodes,entities["Γ³"], integrationOrder_Γ, normal = true)
    elements["Γ⁴ᵖ"] = getElements(nodes,entities["Γ⁴"], integrationOrder_Γ, normal = true)
    elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]

    push!(elements["Ωᵖ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["Ωᵍᵖ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
    push!(elements["∂Ωᵖ"],:𝝭)
    push!(elements["Γ¹ᵖ"],:𝝭)
    push!(elements["Γ²ᵖ"],:𝝭)
    push!(elements["Γ³ᵖ"],:𝝭)
    push!(elements["Γ⁴ᵖ"],:𝝭)

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵘ,yᵘ,zᵘ,n = 3,γ = 5)
    elements["Ωᵘ"] = getElements(nodes_u, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["∂Ωᵘ"] = getElements(nodes_u, entities["Γ"], type, integrationOrder_Γ, sp, normal = true)
    elements["Ωᵍᵘ"] = getElements(nodes_u, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γ¹ᵘ"] = getElements(nodes_u, entities["Γ¹"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵘ"] = getElements(nodes_u, entities["Γ²"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ³ᵘ"] = getElements(nodes_u, entities["Γ³"],type,  integrationOrder_Γ, sp, normal = true)
    elements["Γ⁴ᵘ"] = getElements(nodes_u, entities["Γ⁴"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γᵘ"] = elements["Γ¹ᵘ"]∪elements["Γ²ᵘ"]∪elements["Γ³ᵘ"]∪elements["Γ⁴ᵘ"]

    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵘ"], :𝝭)
    push!(elements["∂Ωᵘ"], :𝝭)
    push!(elements["Γ¹ᵘ"], :𝝭)
    push!(elements["Γ²ᵘ"], :𝝭)
    push!(elements["Γ³ᵘ"], :𝝭)
    push!(elements["Γ⁴ᵘ"], :𝝭)
    push!(elements["Ωᵘ"],  :𝗠=>𝗠)
    push!(elements["∂Ωᵘ"], :𝗠=>𝗠)
    push!(elements["Γ¹ᵘ"], :𝗠=>𝗠)
    push!(elements["Γ²ᵘ"], :𝗠=>𝗠)
    push!(elements["Γ³ᵘ"], :𝗠=>𝗠)
    push!(elements["Γ⁴ᵘ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵘ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y, :∂𝝭∂z)
    push!(elements["Ωᵍᵘ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    set∇𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["∂Ωᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵖ"])
    set𝝭!(elements["Ωᵘ"])
    set𝝭!(elements["∂Ωᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])
    set𝝭!(elements["Γᵘ"])

    gmsh.finalize()

    return elements, nodes, nodes_u
end

function cal_area_support(elms::Vector{ApproxOperator.AbstractElement})
    𝐴s = zeros(length(elms))
    for (i,elm) in enumerate(elms)
        x₁ = elm.𝓒[1].x
        y₁ = elm.𝓒[1].y
        x₂ = elm.𝓒[2].x
        y₂ = elm.𝓒[2].y
        x₃ = elm.𝓒[3].x
        y₃ = elm.𝓒[3].y
        𝐴s[i] = 0.5*(x₁*y₂ + x₂*y₃ + x₃*y₁ - x₂*y₁ - x₃*y₂ - x₁*y₃)
    end
    avg𝐴 = mean(𝐴s)
    var𝐴 = var(𝐴s)
    s = (4/3^0.5*avg𝐴)^0.5
    return s, var𝐴
end
