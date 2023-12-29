// Gmsh project created on Wed Dec 27 14:13:33 2023
SetFactory("OpenCASCADE");
Mesh.MeshSizeMax = 0.036;
Mesh.MeshSizeMin = 0.036;
Mesh.ElementOrder = 1;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Circle(5) = {0.1, 0.1, 0, 0.8, 0, 2*Pi};
//+
Curve Loop(2) = {5};
//+
Plane Surface(2) = {2};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+
Physical Curve("Diri", 1) = {6, 5, 4, 3};
//+
Physical Curve("Inter", 2) = {2};
//+
Physical Surface("right", 3) = {1};
