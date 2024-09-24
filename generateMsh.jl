
using BenchmarkExample
import Gmsh: gmsh
using BubbleMsh

# n = 1

# filename = "patchtest_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1)

# filename = "patchtest_quad_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1, quad=true)

# filename = "cantilever_tri3_"
# BenchmarkExample.CantileverBeam.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1, quad=false, order=1)

# filename = "cantilever_quad_"
# BenchmarkExample.CantileverBeam.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1, quad=true, order=1)

# for n in 49:56
#     filename = "cantilever_"
#     BenchmarkExample.CantileverBeam.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1, quad=false, order=1)
# end

# filename = "plate_with_hole_tri3_"
# BenchmarkExample.PlateWithHole.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1)

# ndiv = 2
# filename = "./msh/plate_with_hole_b_20.msh"
# bubblemsh(filename,[2.5,2.5,0.0],[1.5,1.5,0.0],70-20,0.4,0.4, maxiter=1000)
# nidv = 4
# filename = "./msh/plate_with_hole_b_58.msh"
# bubblemsh(filename,[2.5,2.5,0.0],[1.5,1.5,0.0],195-58,0.24,0.2, maxiter=1000)
# nidv = 8
# filename = "./msh/plate_with_hole_b_96.msh"
# bubblemsh(filename,[2.5,2.5,0.0],[1.5,1.5,0.0],709-115,0.15,0.15, maxiter=1000)
# nidv = 16
# filename = "./msh/plate_with_hole_b_304.msh"
# bubblemsh(filename,[2.5,2.5,0.0],[1.5,1.5,0.0],2698-304,0.085,0.07, maxiter=2000)

n = 3
n₁ = 2*n
n₂ = n
c₁ = 1.37
c₂ = 1.28
c₃ = 1.52
dx₁ = 0.25π/n₂
dx₂ = 4*(c₁-1)/(c₁^n₁-1)
dx₃ = 4*(c₁-1)/(c₁^n₁-1)*c₁^(n₁-1)
dx₄ = 5*(c₂-1)/(c₂^n₂-1)
dx₅ = 4*2^0.5*(c₃-1)/(c₃^n₁-1)
err1 = 1 - dx₂/dx₁
err2 = 1 - dx₄/dx₃
err3 = 1 - dx₅/dx₁
if abs(err1) ≤ 1e-1 && abs(err2) ≤ 1e-1 && abs(err3) ≤ 1e-1
    BenchmarkExample.PlateWithHole.generateMsh("./msh/plate_with_hole_tri3_"*string(n)*".msh", transfinite = (n₁+1,n₂+1), coef = (c₁,c₂,c₃))
    # BenchmarkExample.PlateWithHole.generateMsh("./msh/plate_with_hole_tri3_"*string(n₂)*"_"*string(n₁)*".msh", transfinite = (n₁+1,n₂+1), coef = (c₁,c₂,c₃))
    println("error_1 = $err1, error_2 = $err2, error_3 = $err3")
else
    error("coefficient = $c₁, $c₂, $c₃ is not proper!, error_1 = $err1, error_2 = $err2, error_3 = $err3")
end

# 1  -> c₁ = 1.370, c₂ =      , c₃ = 5.800
# 2  -> c₁ = 1.700, c₂ = 1.500, c₃ = 2.000
# 3  -> c₁ = 1.370, c₂ = 1.280, c₃ = 1.520
# 3-7  -> c = 0.80, β = 0.8
# 3-8  -> c = 0.85, β = 0.65
# 4-6  -> c = 0.68, β = 1.5
# 4-7  -> c = 0.75, β = 1.2
# 4  -> c₁ = 1.2500, c₂ = 1.2000, c₃ = 1.3500
# 5  -> c₁ = 1.2000, c₂ = 1.1000, c₃ = 1.2500
# 6  -> c₁ = 1.0000, c₂ = 1.0750, c₃ = 1.0000
# 7  -> c₁ = 1.0000, c₂ = 1.0750, c₃ = 1.0000
# 8  -> c₁ = 1.0000, c₂ = 1.0750, c₃ = 1.0000
# 9  -> c₁ = 1.0000, c₂ = 1.0500, c₃ = 1.0000
# 10 -> c₁ = 1.0000, c₂ = 1.0500, c₃ = 1.0000
# 11 -> c₁ = 1.0000, c₂ = 1.0500, c₃ = 1.0000
# 12 -> c₁ = 1.0000, c₂ = 1.0500, c₃ = 1.0000
# 13 -> c₁ = 1.0000, c₂ = 1.0500, c₃ = 1.0000
# 14 -> c₁ = 1.0000, c₂ = 1.0250, c₃ = 1.0000
# 15 -> c₁ = 1.0000, c₂ = 1.0250, c₃ = 1.0000
# 16 -> c₁ = 1.0000, c₂ = 1.0250, c₃ = 1.0000
# 17 -> c₁ = 1.0000, c₂ = 1.0250, c₃ = 1.0000
# 18 -> c₁ = 1.0000, c₂ = 1.0250, c₃ = 1.0000
# 19 -> c₁ = 1.0000, c₂ = 1.0250, c₃ = 1.0000
# 20 -> c₁ = 1.0000, c₂ = 1.0250, c₃ = 1.0000
# 21 -> c₁ = 1.0000, c₂ = 1.0250, c₃ = 1.0000
# 22 -> c₁ = 1.0000, c₂ = 1.0250, c₃ = 1.0000
# 23 -> c₁ = 1.0000, c₂ = 1.0200, c₃ = 1.0000
# 24 -> c₁ = 1.0000, c₂ = 1.0200, c₃ = 1.0000
# 25 -> c₁ = 1.0000, c₂ = 1.0200, c₃ = 1.0000
# 26 -> c₁ = 1.0000, c₂ = 1.0200, c₃ = 1.0000
# 27 -> c₁ = 1.0000, c₂ = 1.0200, c₃ = 1.0000
# 28 -> c₁ = 1.0000, c₂ = 1.0200, c₃ = 1.0000
# 29 -> c₁ = 1.0000, c₂ = 1.0200, c₃ = 1.0000
# 30 -> c₁ = 1.0000, c₂ = 1.0150, c₃ = 1.0000
# 31 -> c₁ = 1.0000, c₂ = 1.0150, c₃ = 1.0000
# 32 -> c₁ = 1.0000, c₂ = 1.0150, c₃ = 1.0000
# 33 -> c₁ = 1.0000, c₂ = 1.0150, c₃ = 1.0000
# 34 -> c₁ = 1.0000, c₂ = 1.0150, c₃ = 1.0000
# 35 -> c₁ = 1.0000, c₂ = 1.0110, c₃ = 1.0000
# 36 -> c₁ = 1.0000, c₂ = 1.0110, c₃ = 1.0000
# 37 -> c₁ = 1.0000, c₂ = 1.0110, c₃ = 1.0000
# 38 -> c₁ = 1.0000, c₂ = 1.0110, c₃ = 1.0000
# 39 -> c₁ = 1.0000, c₂ = 1.0110, c₃ = 1.0000
# 40 -> c₁ = 1.0000, c₂ = 1.0110, c₃ = 1.0000
# 41 -> c₁ = 1.0000, c₂ = 1.0110, c₃ = 1.0000
# 42 -> c₁ = 1.0000, c₂ = 1.0110, c₃ = 1.0000