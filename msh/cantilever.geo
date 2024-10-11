
a = 48.0;
b = 12.0;
n = 16;

Point(1) = {0.0, -b/2, 0.0, c};
Point(2) = {  a, -b/2, 0.0, c};
Point(3) = {  a,  b/2, 0.0, c};
Point(4) = {0.0,  b/2, 0.0, c};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Transfinite Curve{1,3} = 4*n+1;
Transfinite Curve{2,4} = n+1;
// Transfinite Curve{1,3} = 107;
// Transfinite Curve{2,4} = 28;
Transfinite Surface{1};

Physical Curve("Γᵗ") = {2};
Physical Curve("Γᵍ") = {4};
// Physical Curve("Γ") = {1,2,3,4};
Physical Surface("Ω") = {1};

Mesh.Renumber = 0;
Mesh.Algorithm = 1;
// Mesh.MshFileVersion = 2;
Mesh.SecondOrderIncomplete = 1;
Mesh 2;
RecombineMesh;
SetOrder 2;
