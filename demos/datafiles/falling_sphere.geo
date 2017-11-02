// Gmsh project created on Wed Nov  1 12:57:53 2017
SetFactory("OpenCASCADE");

// Create cylinder with sphere cutout
Cylinder(1) = {0, 0, 0, 0, 0, 10, 1, 2*Pi};
Sphere(2) = {0, 0, 3, 0.4, -Pi/2, Pi/2, 2*Pi};
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

// Mesh density fields
Field[1] = Ball;
Field[1].Radius = 0.6;
Field[1].ZCenter = 3.0;
Field[1].VIn = 0.1;
Field[1].VOut = 0.3;
Background Field = 1;

// Physical domains for export
Physical Surface(10) = {7}; // Inlet
Physical Surface(20) = {6}; // Outlet
Physical Surface(30) = {4}; // Sphere
Physical Surface(40) = {5}; // Walls
Physical Volume(100) = {3}; // The fluid domain

