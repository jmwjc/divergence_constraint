
// Geometry.MatchMeshTolerance = 1e-2;
Mesh.RecombineAll = 0;

L = 1;
n = 16;

Point(1) = {0.0,0.0,0.0};
Point(2) = {  L,0.0,0.0};
Point(3) = {  L,  L,0.0};
Point(4) = {0.0,  L,0.0};
Point(5) = {0.0,0.0,  L};
Point(6) = {  L,0.0,  L};
Point(7) = {  L,  L,  L};
Point(8) = {0.0,  L,  L};
Point(9)  = {0.5*L,  0.0,0.0};
Point(10) = {    L,0.5*L,0.0};
Point(11) = {0.5*L,    L,0.0};
Point(12) = {  0.0,0.5*L,0.0};
Point(13) = {0.5*L,0.5*L,0.0};
Point(14) = {0.5*L,  0.0,  L};
Point(15) = {    L,0.5*L,  L};
Point(16) = {0.5*L,    L,  L};
Point(17) = {  0.0,0.5*L,  L};
Point(18) = {0.5*L,0.5*L,  L};

Line(1) = {1,9};
Line(2) = {9,2};
Line(3) = {2,10};
Line(4) = {10,3};
Line(5) = {3,11};
Line(6) = {11,4};
Line(7) = {4,12};
Line(8) = {12,1};
Line(9) = {9,13};
Line(10) = {13,12};
Line(11) = {11,13};
Line(12) = {13,10};
Line(13) = {5,14};
Line(14) = {14,6};
Line(15) = {6,15};
Line(16) = {15,7};
Line(17) = {7,16};
Line(18) = {16,8};
Line(19) = {8,17};
Line(20) = {17,5};
Line(21) = {14,18};
Line(22) = {18,17};
Line(23) = {16,18};
Line(24) = {18,15};
Line(25) = {1,5};
Line(26) = {2,6};
Line(27) = {3,7};
Line(28) = {4,8};
Line(29) = {9,14};
Line(30) = {10,15};
Line(31) = {11,16};
Line(32) = {12,17};
Line(33) = {13,18};

Curve Loop(1) = {1,9,10,8};
Curve Loop(2) = {2,3,-12,-9};
Curve Loop(3) = {4,5,11,12};
Curve Loop(4) = {6,7,-10,-11};
Curve Loop(5) = {13,21,22,20};
Curve Loop(6) = {14,15,-24,-21};
Curve Loop(7) = {16,17,23,24};
Curve Loop(8) = {18,19,-22,-23};
Curve Loop(9) = {1,29,-13,-25};
Curve Loop(10) = {2,26,-14,-29};
Curve Loop(11) = {3,30,-15,-26};
Curve Loop(12) = {4,27,-16,-30};
Curve Loop(13) = {5,31,-17,-27};
Curve Loop(14) = {6,28,-18,-31};
Curve Loop(15) = {7,32,-19,-28};
Curve Loop(16) = {8,25,-20,-32};
Curve Loop(17) = {9,33,-21,-29};
Curve Loop(18) = {10,32,-22,-33};
Curve Loop(19) = {11,33,-23,-31};
Curve Loop(20) = {12,30,-24,-33};

Plane Surface(1) = {-1};
Plane Surface(2) = {-2};
Plane Surface(3) = {-3};
Plane Surface(4) = {-4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};
Plane Surface(14) = {14};
Plane Surface(15) = {15};
Plane Surface(16) = {16};
Plane Surface(17) = {17};
Plane Surface(18) = {18};
Plane Surface(19) = {19};
Plane Surface(20) = {20};

Surface Loop(1) = {1,5,9,16,17,18};
Surface Loop(2) = {2,6,10,11,-20,-17};
Surface Loop(3) = {3,7,12,13,19,20};
Surface Loop(4) = {4,8,14,15,-18,-19};

Volume(1) = {1};
Volume(2) = {2};
Volume(3) = {3};
Volume(4) = {4};

Physical Surface("Γᵗ") = {5};
Physical Surface("Γᵍ") = {1,2,3,4,6,7,8,9,10,15,16};
Physical Surface("Γʳ") = {11,12,13,14};
Physical Volume("Ω") = {1,2,3,4};

Transfinite Curve{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24} = n+1;
Transfinite Curve{25,26,27,28,29,30,31,32,33} = 2*n+1;
Transfinite Surface{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
Transfinite Volume{1,2,3,4};

Mesh 3;

// RecombineMesh;