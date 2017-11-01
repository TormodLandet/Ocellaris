// Gmsh project created on Wed Nov  1 12:57:53 2017
SetFactory("OpenCASCADE");

// Create cylinder with sphere cutout
Cylinder(1) = {0, 0, 0, 0, 0, 10, 1, 2*Pi};
Sphere(2) = {0, 0, 3, 0.2, -Pi/2, Pi/2, 2*Pi};
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

// Mesh density fields
Field[1] = MathEval;
Field[1].F = "0.1";

Field[2] = MathEval;
Field[2].F = "0.05";

Field[3] = Restrict;
Field[3].IField = 2;
Field[3].FacesList = {4};

Field[4] = Min;
Field[4].FieldsList = {1, 3};
Background Field = 4;

// Physical domains for export
Physical Surface(10) = {7}; // Inlet
Physical Surface(20) = {6}; // Outlet
Physical Surface(30) = {4}; // Sphere
Physical Surface(40) = {5}; // Walls
Physical Volume(100) = {3}; // The fluid domain

