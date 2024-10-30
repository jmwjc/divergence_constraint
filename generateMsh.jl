
using BenchmarkExample
import Gmsh: gmsh
using BubbleMsh

n = 2

# nx = 2
# filename = "patchtest_"
# filename = "patchtest_tri6_"
# filename = "patchtest_quad_"
# for n in 1:20
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1, order = 1, quad=true)
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(nx)*"_"*string(n)*".msh", transfinite = (nx+1,n+1), order = 1)
# end

# filename = "patchtest_quad_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1, quad=true)
# filename = "patchtest_quad_"
# BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1, quad=true)

# filename = "cook_tri3_"
# BenchmarkExample.CookMembrane.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = (2*n+1,n+1), order = 1, quad=false)

# filename = "cantilever_tri6_"
# BenchmarkExample.CantileverBeam.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1, quad=false, order=2)

# filename = "cantilever_quad_"
# BenchmarkExample.CantileverBeam.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1, quad=true, order=1)

# for n in 49:56
#     filename = "cantilever_"
#     BenchmarkExample.CantileverBeam.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1, quad=false, order=1)
# end

filename = "plate_with_hole_convergence_tri6_"
BenchmarkExample.PlateWithHole.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = (n+1,2*n+1), order = 2, mode = 2)

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

# n = 31
# n₁ = 68
# n₂ = 32
# n₁ = 2*n
# n₂ = n
# c₁ = 1.023355
# c₂ = 1.01884
# c₃ = 1.03155
# dx₁ = 0.25π/n₂
# dx₂ = 4*(c₁-1)/(c₁^n₁-1)
# dx₃ = 4*(c₁-1)/(c₁^n₁-1)*c₁^(n₁-1)
# dx₄ = 5*(c₂-1)/(c₂^n₂-1)
# dx₅ = 4*2^0.5*(c₃-1)/(c₃^n₁-1)
# err1 = 1 - dx₂/dx₁
# err2 = 1 - dx₄/dx₃
# err3 = 1 - dx₅/dx₁
# if abs(err1) ≤ 1e-3 && abs(err2) ≤ 1e-3 && abs(err3) ≤ 1e-3
    # BenchmarkExample.PlateWithHole.generateMsh("./msh/plate_with_hole_tri3_"*string(n)*".msh", transfinite = (n₁+1,n₂+1), coef = (c₁,c₂,c₃))
    # BenchmarkExample.PlateWithHole.generateMsh("./msh/plate_with_hole_tri3_"*string(n₂)*"_"*string(n₁)*".msh", transfinite = (n₁+1,n₂+1), coef = (c₁,c₂,c₃))

    # BenchmarkExample.PlateWithHole.generateMsh("./msh/plate_with_hole_tri6_"*string(n)*".msh", transfinite = (n₁+1,n₂+1), coef = (c₁,c₂,c₃),order=2)
    # println("error_1 = $err1, error_2 = $err2, error_3 = $err3")
# else
    # error("coefficient = $c₁, $c₂, $c₃ is not proper!, error_1 = $err1, error_2 = $err2, error_3 = $err3")
# end

# 1    -> c₁ = 1.3700, c₂ =       , c₃ = 5.8000
# 2    -> c₁ = 1.7000, c₂ = 1.5000, c₃ = 2.0000
# 3    -> c₁ = 1.3700, c₂ = 1.2800, c₃ = 1.5200
# 3-7  -> c₁ = 1.2500, c₂ = 1.5000, c₃ = 1.3500
# 3-8  -> c₁ = 1.1800, c₂ = 1.8000, c₃ = 1.2700
# 4-6  -> c₁ = 1.5000, c₂ = 0.9000, c₃ = 1.6400
# 4-7  -> c₁ = 1.3481, c₂ = 1.0400, c₃ = 1.4630
# 4    -> c₁ = 1.2570, c₂ = 1.1690, c₃ = 1.3514
# 4-9  -> c₁ = 1.2000, c₂ = 1.3000, c₃ = 1.2800
# 4-10 -> c₁ = 1.1500, c₂ = 1.4100, c₃ = 1.2200
# 5-9  -> c₁ = 1.2500, c₂ = 1.0400, c₃ = 1.3250
# 5    -> c₁ = 1.2000, c₂ = 1.1000, c₃ = 1.2500
# 6    -> c₁ = 1.1500, c₂ = 1.0750, c₃ = 1.2000
# 7    -> c₁ = 1.13292, c₂ = 1.0754, c₃ = 1.1789
# 7-15 -> c₁ = 1.1144, c₂ = 1.1098, c₃ = 1.157
# 7-16 -> c₁ = 1.1100, c₂ = 1.0950, c₃ = 1.1500
# 7-17 -> c₁ = 1.0900, c₂ = 1.1500, c₃ = 1.1300
# 8-13 -> c₁ = 1.1748, c₂ = 0.9763, c₃ = 1.2247
# 8-14 -> c₁ = 1.1400, c₂ = 1.0400, c₃ = 1.1900
# 8-15 -> c₁ = 1.1300, c₂ = 1.0400, c₃ = 1.1700
# 8    -> c₁ = 1.1145, c₂ = 1.0634, c₃ = 1.1538
# 8-16 -> c₁ = 1.1100, c₂ = 1.0750, c₃ = 1.1400
# 8-17 -> c₁ = 1.0950, c₂ = 1.0100, c₃ = 1.1300
# 8-18 -> c₁ = 1.0900, c₂ = 1.0100, c₃ = 1.1300
# 9-17 -> c₁ = 1.1100, c₂ = 1.0500, c₃ = 1.1400
# 9    -> c₁ = 1.1000, c₂ = 1.0500, c₃ = 1.1400
# 10   -> c₁ = 1.0900, c₂ = 1.0500, c₃ = 1.1200
# 11   -> c₁ = 1.0800, c₂ = 1.0500, c₃ = 1.1100
# 12   -> c₁ = 1.0700, c₂ = 1.0500, c₃ = 1.1000
# 13   -> c₁ = 1.0650, c₂ = 1.0500, c₃ = 1.0900
# 14   -> c₁ = 1.0650, c₂ = 1.0250, c₃ = 1.0800
# 15-28-> c₁ = 1.06666, c₂ = 1.0153, c₃ = 1.08738
# 15-29-> c₁ = 1.06217, c₂ = 1.02275, c₃ = 1.08215
# 15   -> c₁ = 1.0580, c₂ = 1.0300, c₃ = 1.0774
# 15-31-> c₁ = 1.0543, c₂ = 1.0368, c₃ = 1.0730
# 15-32-> c₁ = 1.05078, c₂ = 1.04359, c₃ = 1.06897
# 15-33-> c₁ = 1.0476, c₂ = 1.0500, c₃ = 1.0652
# 16-27-> c₁ = 1.0757, c₂ = 0.9939, c₃ = 1.0970
# 16-28-> c₁ = 1.0706, c₂ = 1.0012, c₃ = 1.0912
# 16-29-> c₁ = 1.06595, c₂ = 1.0083, c₃ = 1.0858
# 16-30-> c₁ = 1.0617, c₂ = 1.0150, c₃ = 1.0809
# 16-31-> c₁ = 1.0578, c₂ = 1.0215, c₃ = 1.0764
# 16   -> c₁ = 1.0542, c₂ = 1.0279, c₃ = 1.0723
# 16-33-> c₁ = 1.0509, c₂ = 1.0339, c₃ = 1.0684
# 16-34-> c₁ = 1.0478, c₂ = 1.0398, c₃ = 1.0649
# 16-35-> c₁ = 1.0450, c₂ = 1.0455, c₃ = 1.0616
# 16-36-> c₁ = 1.04235, c₂ = 1.051, c₃ = 1.0585
# 16-37-> c₁ = 1.0399, c₂ = 1.0564, c₃ = 1.0556
# 17-32-> c₁ = 1.0550, c₂ = 1.0250, c₃ = 1.0750
# 17-33-> c₁ = 1.0500, c₂ = 1.0250, c₃ = 1.0700
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
# 29   -> c₁ = 1.0292, c₂ = 1.0145, c₃ = 1.03885
# 30   -> c₁ = 1.0282, c₂ = 1.0140, c₃ = 1.0375
# 31_60-> c₁ = 1.0273, c₂ = 1.0135, c₃ = 1.03625
# 31_61-> c₁ = 1.0273, c₂ = 1.0135, c₃ = 1.03625
# 31   -> c₁ = 1.0273, c₂ = 1.0135, c₃ = 1.03625
# 31-63-> c₁ = 1.02645, c₂ = 1.0151, c₃ = 1.03525
# 31-64-> c₁ = 1.0256, c₂ = 1.0167, c₃ = 1.0343
# 31-65-> c₁ = 1.0248, c₂ = 1.01812, c₃ = 1.03338
# 31-66-> c₁ = 1.02403, c₂ = 1.01967, c₃ = 1.0325
# 31-67-> c₁ = 1.02329, c₂ = 1.0211, c₃ = 1.03165
# 32-59-> c₁ = 1.031, c₂ = 1.00518, c₃ = 1.0403
# 32-60-> c₁ = 1.03, c₂ = 1.00683, c₃ = 1.03915
# 32-61-> c₁ = 1.02904, c₂ = 1.00845, c₃ = 1.03808
# 32-62-> c₁ = 1.02813, c₂ = 1.01002, c₃ = 1.03704
# 32-63-> c₁ = 1.02726, c₂ = 1.0115, c₃ = 1.03605
# 32   -> c₁ = 1.0264, c₂ = 1.01305, c₃ = 1.0351
# 32-65-> c₁ = 1.02561, c₂ = 1.01455, c₃ = 1.03415
# 32-66-> c₁ = 1.02483, c₂ = 1.016, c₃ = 1.03326
# 32-67-> c₁ = 1.02408, c₂ = 1.01743, c₃ = 1.0324
# 32-68-> c₁ = 1.023355, c₂ = 1.01884, c₃ = 1.03155
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