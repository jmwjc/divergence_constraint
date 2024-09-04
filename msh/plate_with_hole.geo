
a = 1.0;
b = 5.0;
// n = 8;
// c = 0.099970;
c = 0.45;

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

Plane Surface(1) = {16};
Plane Surface(2) = {17};
Plane Surface(3) = {18};
Plane Surface(4) = {19};
Plane Surface(5) = {20};


// Transfinite Curve{1,2,4,7,9,10,11,12,15} = n+1;
// Transfinite Curve{3,5,6,8,13,14} = 2*n+1;


Physical Curve("Γᵍ₃") = {4,5};
Physical Curve("Γᵗ₁") = {6,7};
Physical Curve("Γᵗ₂") = {1,10};
Physical Curve("Γᵍ₁") = {2,3};
Physical Curve("Γᵍ₂") = {8,9};
//Physical Curve("Γ") = {11,12,13,14,15};
Physical Curve("Γ") = {1,2,3,4,5};
Physical Surface("Ω") = {1,2,3,4,5};

// Transfinite Surface{1,2,3,4,5};

Mesh.Algorithm = 8;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;
