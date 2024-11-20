
using Gmsh 

function import_fem(filename::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    integrationOrder_Ω = 2
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 2

    elements["Ω"] = getElements(nodes,entities["Ω"],integrationOrder_Ω)
    elements["Ωᵍ"] = getElements(nodes,entities["Ω"],integrationOrder_Ωᵍ)
    elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"],integrationOrder_Γ,normal=true)
    elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"],integrationOrder_Γ,normal=true)
    elements["Γʳ"] = getElements(nodes,entities["Γʳ"],integrationOrder_Γ,normal=true)

    push!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    push!(elements["Ωᵍ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    push!(elements["Γᵗ"],:𝝭)
    push!(elements["Γᵍ"],:𝝭)
    push!(elements["Γʳ"],:𝝭)

    gmsh.finalize()

    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γʳ"])

    return elements, nodes
end

function import_linear_mix(filename1::String,filename2::String,n)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    s = 1.5
    s₁ = s/n*ones(length(nodes_p))
    s₂ = s/n*ones(length(nodes_p))
    s₃ = s/n*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)

    integrationOrder_Ω = 2
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 2

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"], integrationOrder_Γ, normal = true)
    elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"], integrationOrder_Γ, normal = true)
    elements["Γʳ"] = getElements(nodes,entities["Γʳ"],integrationOrder_Γ,normal=true)

    push!(elements["Ωᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    push!(elements["Ωᵍᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    push!(elements["Γᵗ"],:𝝭)
    push!(elements["Γᵍ"],:𝝭)
    push!(elements["Γʳ"],:𝝭)

    set∇𝝭!(elements["Ωᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γʳ"])

    type = ReproducingKernel{:Linear3D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"],type,  integrationOrder_Γ, sp, normal = true)

    nₘ = 10
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    ∂𝗠∂z = zeros(nₘ)
    push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y, :∂𝝭∂z)
    push!(elements["Γᵍᵖ"], :𝝭)
    push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂𝗠∂z=>∂𝗠∂z)
    push!(elements["Γᵍᵖ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵖ"], :𝝭)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠)

    set∇𝝭!(elements["Ωᵖ"])
    # set𝝭!(elements["Ωᵍᵖ"])
    # set𝝭!(elements["Γᵍᵖ"])

    gmsh.finalize()

    return elements, nodes, nodes_p, sp, type
end

function import_quadratic_mix(filename1::String,filename2::String,n)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    Ω = getElements(nodes_p, entities["Ω"])
    s = 2.5
    s₁ = s/n*ones(length(nodes_p))
    s₂ = s/n*ones(length(nodes_p))
    s₃ = s/n*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)

    integrationOrder_Ω = 4
    integrationOrder_Ωᵍ = 8
    integrationOrder_Γ = 4

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements["Ωᵘ"] = getElements(nodes,entities["Ω"], integrationOrder_Ω)
    elements["Ωᵍᵘ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"], integrationOrder_Γ, normal = true)
    elements["Γᵍᵘ"] = getElements(nodes,entities["Γᵍ"], integrationOrder_Γ, normal = true)

    push!(elements["Ωᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    push!(elements["Ωᵍᵘ"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    push!(elements["Γᵗ"],:𝝭)
    push!(elements["Γᵍᵘ"],:𝝭)

    set∇𝝭!(elements["Ωᵘ"])
    set∇𝝭!(elements["Ωᵍᵘ"])
    set𝝭!(elements["Γᵗ"])
    set𝝭!(elements["Γᵍᵘ"])

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    type = ReproducingKernel{:Quadratic3D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["Ωᵍᵖ"] = getElements(nodes_p, entities["Ω"], type,  integrationOrder_Ωᵍ, sp)
    elements["Γᵍᵖ"] = getElements(nodes_p, entities["Γᵍ"],type,  integrationOrder_Γ, sp, normal = true)

    nₘ = 55
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ωᵖ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y, :∂𝝭∂z)
    push!(elements["Γᵍᵖ"], :𝝭)
    push!(elements["Ωᵖ"],  :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂𝗠∂z=>∂𝗠∂z)
    push!(elements["Γᵍᵖ"], :𝗠=>𝗠)
    push!(elements["Ωᵍᵖ"], :𝝭)
    push!(elements["Ωᵍᵖ"], :𝗠=>𝗠)

    set∇𝝭!(elements["Ωᵖ"])
    set𝝭!(elements["Ωᵍᵖ"])
    set𝝭!(elements["Γᵍᵖ"])

    gmsh.finalize()

    return elements, nodes, nodes_p, sp, type
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
