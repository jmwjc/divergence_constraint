a = 1.0;
b = 1.0;
c = 0.075;
// n = 17;

Point(1) = {0.0, 0.0, 0.0, c};
Point(2) = {  a, 0.0, 0.0, c};
Point(3) = {  a,   b, 0.0, c};
Point(4) = {0.0,   b, 0.0, c};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};


// Transfinite Curve{1,2,3,4} = n;
// Transfinite Surface{1};
Physical Curve("Γ¹") = {1};
Physical Curve("Γ²") = {2};
Physical Curve("Γ³") = {3};
Physical Curve("Γ⁴") = {4};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;
//Mesh.SecondOrderIncomplete = 1;
//RecombineMesh;
//SetOrder 2;