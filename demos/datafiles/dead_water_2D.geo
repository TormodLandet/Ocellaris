// Gmsh project created on Wed Nov  1 18:39 2017
SetFactory("OpenCASCADE");

////////////////////////////////////////////
// Parameters

// Parameters of the geometry
DefineConstant[ H = {  30.0, Name "Parameters/Total depth" } ];
DefineConstant[ h = {   6.0, Name "Parameters/Upper layer depth" } ];
DefineConstant[ L = { 200.0, Name "Parameters/Domain length" } ];
DefineConstant[ C = {  50.0, Name "Parameters/Hull inlet offset amidships" } ];
DefineConstant[ l = {  34.0, Name "Parameters/Hull length" } ];
DefineConstant[ d = {   5.0, Name "Parameters/Draught amidships" } ];

// Mesh params
DefineConstant[ lc_fine   = { 1.0, Name "Parameters/LC fine" } ];
DefineConstant[ lc_course = { 3.0, Name "Parameters/LC course" } ];

////////////////////////////////////////////
// Geometry

// Create square with a partial circle cutout
// With a circle chord of c = 34 m and a sagitta of s = 5 m
// the circle radius is r = 30.9 m  since r = c^2/(8s) + s/2.
r = l^2 / (8 * d) + d / 2;
Rectangle(1) = {0, -H, 0, L, H};
Disk(2) = {C, r - d, 0, r};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };

// Make mesh conform to initial free surface. The free surface mesh conforming
// line is inset 1 mesh cell from the boundary to avoid mesh degeneration there
Point(100) = {0 + lc_fine, -h, 0};
Point(101) = {L - lc_fine, -h, 0};
Line(100) = {100, 101};
Line {100} In Surface {3}; // Conform to this line

////////////////////////////////////////////
// Mesh cell size fields:

// Fine mesh near the hull
Field[1] = Ball;
Field[1].Radius = r + 3;
Field[1].XCenter = C;
Field[1].YCenter = r - d;
Field[1].VIn = lc_fine;
Field[1].VOut = lc_course;

// Fine mesh near the pycnocline
Field[2] = Box;
Field[2].XMin = 0;
Field[2].XMax = L;
Field[2].YMin = -h - 4;
Field[2].YMax = -h + 4;
Field[2].VIn = lc_fine;
Field[2].VOut = lc_course;

// Fine mesh in the immediate wake of the hull
Field[3] = Box;
Field[3].XMin = C;
Field[3].XMax = C + l * 2;
Field[3].YMin = -h - 4;
Field[3].YMax =  0;
Field[3].VIn = lc_fine;
Field[3].VOut = lc_course;

// Fine mesh near the inlet
Field[4] = Box;
Field[4].XMin = 0;
Field[4].XMax = C / 5;
Field[4].YMin = -H;
Field[4].YMax =  0;
Field[4].VIn = lc_fine;
Field[4].VOut = lc_course;

// Fine mesh near the surface at the outlet
Field[5] = Box;
Field[5].XMin = L - C / 10;
Field[5].XMax = L;
Field[5].YMin = -h - 2;
Field[5].YMax =  0;
Field[5].VIn = lc_fine;
Field[5].VOut = lc_course;

// The resulting mesh size field is the minimum of the above fields
Field[100] = Min;
Field[100].FieldsList = {1, 2, 3, 4, 5};
Background Field = 100;

////////////////////////////////////////////
// Physical domains

Physical Surface(10) = {1}; // Inlet
Physical Surface(20) = {5}; // Outlet
Physical Surface(30) = {3}; // Fram
Physical Surface(40) = {2, 4, 6}; // Top/bottom
Physical Volume(100) = {3}; // The fluid domain
