
a = 1.0;
b = 5.0;
c = 2.0;
n = 2;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {  a, 0.0, 0.0};
Point(3) = {  b, 0.0, 0.0};
Point(4) = {  b,   b, 0.0};
Point(5) = {0.0,   b, 0.0};
Point(6) = {0.0,   a, 0.0};
Point(7) = {  c, 0.0, 0.0};
Point(8) = {0.0,   c, 0.0};
Point(9) = {  b, 2^0.5/2*c, 0.0};
Point(10) = {2^0.5/2*c, b, 0.0};
Point(11) = {2^0.5/2*a, 2^0.5/2*a, 0.0};
Point(12) = {2^0.5/2*c, 2^0.5/2*c, 0.0};

Line(1) = {2,7};
Line(2) = {7,3};
Line(3) = {3,9};
Line(4) = {9,4};
Line(5) = {4,10};
Line(6) = {10,5};
Line(7) = {5,8};
Line(8) = {8,6};
Circle(9) = {6,1,11};
Circle(10) = {11,1,2};
Circle(11) = {7,1,12};
Circle(12) = {12,1,8};
Line(13) = {11,12};
Line(14) = {12,9};
Line(15) = {10,12};

Curve Loop(1) = {9,13,12,8};
Curve Loop(2) = {10,1,11,-13};
Curve Loop(3) = {-12,-15,6,7};
Curve Loop(4) = {11,14,-3,-2};
Curve Loop(5) = {14,4,5,15};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Physical Surface("Ω") = {1,2,3,4,5};
Physical Curve("Γᵗ") = {3,4,5,6,9,10};
Physical Curve("Γᵍ") = {1,2,7,8};
Transfinite Curve{1,3,6,8,9,10,11,12,13} = n+1;
Transfinite Curve{2,4,5,7,14,15} = 2*n+1;
Transfinite Surface{1,3,5};
Transfinite Surface{2,4} Right;

Mesh.Algorithm = 1;
// Mesh.MshFileVersion = 2;
// Mesh.Renumber = 0;
Mesh 2;
// RefineMesh;
// RecombineMesh;
// SetOrder 2;
// Mesh.SecondOrderIncomplete = 1;
