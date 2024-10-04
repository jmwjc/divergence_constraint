
a = 1.0;
b = 5.0;
n = 32;
n1 = 2*n;
n2 = n;
c1 = 1.025;
c2 = 1.015;
c3 = 1.035;
c = 0.112;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {  a, 0.0, 0.0};
Point(3) = {  b, 0.0, 0.0};
Point(4) = {  b,   b, 0.0};
Point(5) = {0.0,   b, 0.0};
Point(6) = {0.0,   a, 0.0};
Point(7) = {2^0.5/2*a, 2^0.5/2*a, 0.0};

Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,5};
Line(4) = {5,6};
Circle(5) = {6,1,7};
Circle(6) = {7,1,2};
Line(7) = {7,4};

Curve Loop(1) = {5,7,3,4};
Curve Loop(2) = {6,1,2,-7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Physical Surface("Ω") = {1,2};
Physical Curve("Γᵗ") = {2,3,5,6};
Physical Curve("Γᵍ") = {1,4};
Transfinite Curve{1} = n1+1 Using Progression c1;
Transfinite Curve{2} = n2+1 Using Progression c2;
Transfinite Curve{-3} = n2+1 Using Progression c2;
Transfinite Curve{-4} = n1+1 Using Progression c1;
Transfinite Curve{5,6} = n2+1;
Transfinite Curve{7} = n1+1 Using Progression c3;
Transfinite Surface{1};
Transfinite Surface{2} Right;

Mesh.Algorithm = 1;
// Mesh.MshFileVersion = 2;
Mesh.Renumber = 0;
Mesh 2;
RecombineMesh;
SetOrder 2;
Mesh.SecondOrderIncomplete = 1;
