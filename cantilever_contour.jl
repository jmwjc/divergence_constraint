
using CairoMakie

L = 48.0
D = 12.0
P = 1000
E = 3e6
ν = 0.4999999
Ē = E/(1.0-ν^2)
ν̄ = ν/(1.0-ν)
I = D^3/12
EI = Ē*I

Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
Cᵢⱼᵢⱼ = E/(1+ν)/2

u(x,y) = -P*y/6/EI*((6*L-3*x)*x + (2+ν̄)*(y^2-D^2/4))
v(x,y) = P/6/EI*(3*ν̄*y^2*(L-x) + (4+5*ν̄)*D^2*x/4 + (3*L-x)*x^2)
∂u∂x(x,y) = -P/EI*(L-x)*y
∂u∂y(x,y) = -P/6/EI*((6*L-3*x)*x + (2+ν̄)*(3*y^2-D^2/4))
∂v∂x(x,y) = P/6/EI*((6*L-3*x)*x - 3*ν̄*y^2 + (4+5*ν̄)*D^2/4)
∂v∂y(x,y) = P/EI*(L-x)*y*ν̄

ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = ∂u∂y(x,y) + ∂v∂x(x,y)
σ₁₁(x,y) = -P*(L-x)*y/I
σ₂₂(x,y) = 0.0
σ₁₂(x,y) = P/2/I*(D^2/4-y^2)
σ̄₁₁(x,y) = Cᵢᵢᵢᵢ*ε₁₁(x,y) + Cᵢᵢⱼⱼ*ε₂₂(x,y)
σ̄₂₂(x,y) = Cᵢᵢⱼⱼ*ε₁₁(x,y) + Cᵢᵢᵢᵢ*ε₂₂(x,y)
σ̄₃₃(x,y) = Cᵢᵢⱼⱼ*ε₁₁(x,y) + Cᵢᵢⱼⱼ*ε₂₂(x,y)
σ̄₁₂(x,y) = Cᵢⱼᵢⱼ*ε₁₂(x,y)
𝑝(x,y) = (σ̄₁₁(x,y) + σ̄₂₂(x,y) + σ̄₃₃(x,y))/3

fig = Figure()
ind = 100
ax = Axis(fig[1,1], 
    aspect = DataAspect(), 
    xticksvisible = false,
    xticklabelsvisible=false, 
    yticksvisible = false, 
    yticklabelsvisible=false,
)
hidespines!(ax)
hidedecorations!(ax)
xs = LinRange(0, 48, 4*ind)
ys = LinRange(-6, 6, ind)
zs = [𝑝(x,y) for x in xs, y in ys]
surface!(xs,ys,zeros(4*ind,ind),color=zs,colorrange=(-1000,1000),shading=NoShading,colormap=:lightrainbow)
contour!(xs[1:end-1],ys,zs[1:end-1,:],levels=-1e3:200:1e3,color=:azure)
# Colorbar(fig[1,2], limits=(-900,900), colormap=:lightrainbow)
save("./png/cantilever_exact_solution.png",fig, px_per_unit = 10.0)
fig