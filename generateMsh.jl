
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

n = 42
n₁ = 14
n₂ = 8
c₁ = 1.13
c₂ = 1.04
c₃ = 1.17
dx₁ = 0.25π/n₂
dx₂ = 4*(c₁-1)/(c₁^n₁-1)
dx₃ = 4*(c₁-1)/(c₁^n₁-1)*c₁^(n₁-1)
dx₄ = 5*(c₂-1)/(c₂^n₂-1)
dx₅ = 4*2^0.5*(c₃-1)/(c₃^n₁-1)
err1 = 1 - dx₂/dx₁
err2 = 1 - dx₄/dx₃
err3 = 1 - dx₅/dx₁
if abs(err1) ≤ 1e-1 && abs(err2) ≤ 1e-1 && abs(err3) ≤ 1e-1
    # BenchmarkExample.PlateWithHole.generateMsh("./msh/plate_with_hole_tri3_"*string(n)*".msh", transfinite = (n₁+1,n₂+1), coef = (c₁,c₂,c₃))
    BenchmarkExample.PlateWithHole.generateMsh("./msh/plate_with_hole_tri3_"*string(n₂)*"_"*string(n₁)*".msh", transfinite = (n₁+1,n₂+1), coef = (c₁,c₂,c₃))
    println("error_1 = $err1, error_2 = $err2, error_3 = $err3")
else
    error("coefficient = $c₁, $c₂, $c₃ is not proper!, error_1 = $err1, error_2 = $err2, error_3 = $err3")
end

# 1    -> c₁ = 1.3700, c₂ =       , c₃ = 5.8000
# 2    -> c₁ = 1.7000, c₂ = 1.5000, c₃ = 2.0000
# 3    -> c₁ = 1.3700, c₂ = 1.2800, c₃ = 1.5200
# 3-7  -> c₁ = 1.2500, c₂ = 1.5000, c₃ = 1.3500
# 3-8  -> c₁ = 1.1800, c₂ = 1.8000, c₃ = 1.2700
# 4-6  -> c₁ = 1.5000, c₂ = 0.9000, c₃ = 1.6400
# 4-7  -> c₁ = 1.3500, c₂ = 1.0400, c₃ = 1.4600
# 4    -> c₁ = 1.2500, c₂ = 1.2000, c₃ = 1.3500
# 4-9  -> c₁ = 1.2000, c₂ = 1.3000, c₃ = 1.2800
# 4-10 -> c₁ = 1.1500, c₂ = 1.4100, c₃ = 1.2200
# 5-9  -> c₁ = 1.2500, c₂ = 1.0400, c₃ = 1.3250
# 5    -> c₁ = 1.2000, c₂ = 1.1000, c₃ = 1.2500
# 6    -> c₁ = 1.1500, c₂ = 1.0750, c₃ = 1.2000
# 7    -> c₁ = 1.1400, c₂ = 1.0750, c₃ = 1.1800
# 8-15 -> c₁ = 1.1300, c₂ = 1.0400, c₃ = 1.1700
# 8    -> c₁ = 1.1200, c₂ = 1.0750, c₃ = 1.1600
# 9    -> c₁ = 1.1000, c₂ = 1.0500, c₃ = 1.1400
# 10   -> c₁ = 1.0900, c₂ = 1.0500, c₃ = 1.1200
# 11   -> c₁ = 1.0800, c₂ = 1.0500, c₃ = 1.1100
# 12   -> c₁ = 1.0700, c₂ = 1.0500, c₃ = 1.1000
# 13   -> c₁ = 1.0650, c₂ = 1.0500, c₃ = 1.0900
# 14   -> c₁ = 1.0650, c₂ = 1.0250, c₃ = 1.0800
# 15   -> c₁ = 1.0600, c₂ = 1.0250, c₃ = 1.0800
# 16   -> c₁ = 1.0500, c₂ = 1.0250, c₃ = 1.0700
# 17   -> c₁ = 1.0500, c₂ = 1.0250, c₃ = 1.0700
# 18   -> c₁ = 1.0500, c₂ = 1.0250, c₃ = 1.0600
# 19   -> c₁ = 1.0450, c₂ = 1.0250, c₃ = 1.0600
# 20   -> c₁ = 1.0450, c₂ = 1.0250, c₃ = 1.0600
# 21   -> c₁ = 1.0400, c₂ = 1.0250, c₃ = 1.0550
# 22   -> c₁ = 1.0400, c₂ = 1.0250, c₃ = 1.0550
# 23   -> c₁ = 1.0400, c₂ = 1.0200, c₃ = 1.0500
# 24   -> c₁ = 1.0350, c₂ = 1.0200, c₃ = 1.0500
# 25   -> c₁ = 1.0350, c₂ = 1.0200, c₃ = 1.0450
# 26   -> c₁ = 1.0350, c₂ = 1.0200, c₃ = 1.0450
# 27   -> c₁ = 1.0300, c₂ = 1.0200, c₃ = 1.0400
# 28   -> c₁ = 1.0300, c₂ = 1.0200, c₃ = 1.0400
# 29   -> c₁ = 1.0300, c₂ = 1.0200, c₃ = 1.0400
# 30   -> c₁ = 1.0300, c₂ = 1.0150, c₃ = 1.0400
# 31   -> c₁ = 1.0300, c₂ = 1.0150, c₃ = 1.0350
# 32   -> c₁ = 1.0250, c₂ = 1.0150, c₃ = 1.0350
# 33   -> c₁ = 1.0250, c₂ = 1.0150, c₃ = 1.0350
# 34   -> c₁ = 1.0250, c₂ = 1.0150, c₃ = 1.0350
# 35   -> c₁ = 1.0250, c₂ = 1.0110, c₃ = 1.0300
# 36   -> c₁ = 1.0250, c₂ = 1.0110, c₃ = 1.0300
# 37   -> c₁ = 1.0250, c₂ = 1.0110, c₃ = 1.0300
# 38   -> c₁ = 1.0240, c₂ = 1.0110, c₃ = 1.0300
# 39   -> c₁ = 1.0230, c₂ = 1.0110, c₃ = 1.0300
# 40   -> c₁ = 1.0230, c₂ = 1.0110, c₃ = 1.0290
# 41   -> c₁ = 1.0220, c₂ = 1.0110, c₃ = 1.0290
# 42   -> c₁ = 1.0220, c₂ = 1.0110, c₃ = 1.0280