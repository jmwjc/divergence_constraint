
a = 1.0;
b = 5.0;
// n = 13;
// c = 0.099970;
c = 0.112;

Point(1) = {0.0, 0.0, 0.0, c};
Point(2) = {  a, 0.0, 0.0, c};
Point(3) = {  b, 0.0, 0.0, c};
Point(4) = {  b,   b, 0.0, c};
Point(5) = {0.0,   b, 0.0, c};
Point(6) = {0.0,   a, 0.0, c};

Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,5};
Line(4) = {5,6};
Circle(5) = {6,1,2};

Curve Loop(6) = {1,2,3,4,5};

Plane Surface(1) = {6};
Physical Surface("Ω") = {1};
Physical Curve("Γ") = {1,2,3,4,5};
// Transfinite Curve{5} = n+2;
// Transfinite Curve{1,4} = 4*n+1;
// Transfinite Curve{2,3} = 5*n+1;

Mesh.Algorithm = 2;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;
