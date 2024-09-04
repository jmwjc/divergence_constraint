
a = 48.0;
b = 12.0;
// n = 4;
// c = 12.0/8;
// 4165:0.43
// c = 0.4285331665;
c = 12;

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

// Transfinite Curve{1,3} = 4*n+1;
// Transfinite Curve{2,4} = n+1;

Physical Curve("Γᵗ") = {2};
Physical Curve("Γᵍ") = {4};
Physical Curve("Γ₁") = {1};
Physical Curve("Γ₃") = {3};
Physical Curve("Γ") = {1,2,3,4};
Physical Surface("Ω") = {1};

// Transfinite Surface{1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;
Mesh.SecondOrderIncomplete = 1;
SetOrder 1;
//RecombineMesh;
