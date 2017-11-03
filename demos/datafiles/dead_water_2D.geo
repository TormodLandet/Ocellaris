// Gmsh project created on Wed Nov  1 18:39 2017
SetFactory("OpenCASCADE");

// Create square with a partial circle cutout
// With a circle chord of c = 34 m and a sagitta of s = 5 m
// the circle radius is r = 30.9 m  since r = c^2/(8s) + s/2.
Rectangle(1) = {0, -15, 0, 150, 15};
Disk(2) = {50, 25.9, 0, 30.9};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };

// Make mesh conform to initial free surface
Point(100) = {  1, -6, 0};
Point(101) = {149, -6, 0};
Line(100) = {100, 101};
Line {100} In Surface {3}; // Conform to this line

// Mesh density fields
Field[1] = Ball;
Field[1].Radius = 33;
Field[1].XCenter = 50;
Field[1].YCenter = 25.9;
Field[1].VIn = 1.0;
Field[1].VOut = 3.0;

Field[2] = Box;
Field[2].XMin = 0;
Field[2].XMax = 150;
Field[2].YMin = -8;
Field[2].YMax = -2;
Field[2].ZMin = -100;
Field[2].ZMax = 1000;
Field[2].VIn = 1.0;
Field[2].VOut = 5.0;

Field[3] = Box;
Field[3].XMin = 50;
Field[3].XMax = 100;
Field[3].YMin = -8;
Field[3].YMax =  0;
Field[3].VIn = 1.0;
Field[3].VOut = 5.0;

Field[100] = Min;
Field[100].FieldsList = {1, 2, 3};
Background Field = 100;

// Physical domains for export
Physical Surface(10) = {1}; // Inlet
Physical Surface(20) = {5}; // Outlet
Physical Surface(30) = {3}; // Fram
Physical Surface(40) = {2, 4, 6}; // Top/bottom
Physical Volume(100) = {3}; // The fluid domain
