
using BenchmarkExample
import Gmsh: gmsh

n = 8
filename = "patchtest_"
BenchmarkExample.PatchTest.generateMsh("./msh/"*filename*string(n)*".msh", transfinite = n+1)