// Gmsh project created on Wed Dec 27 14:13:33 2023
SetFactory("OpenCASCADE");
Mesh.MeshSizeMax = 0.045;
Mesh.MeshSizeMin = 0.045;
Mesh.ElementOrder = 2;
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
Physical Curve("Diri", 1) = {1, 2, 3, 4};
//+
Physical Surface("whole", 2) = {1};
