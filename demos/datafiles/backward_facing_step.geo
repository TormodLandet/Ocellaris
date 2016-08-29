///////////////////////////////////////////////////////////////////////////////
// Backward facing step GMSH geometry input file
//
// Mesh with
//     gmsh -2 backward_facing_step.geo -o backward_facing_step.msh
// then convert to xml with
//     dolfin-convert backward_facing_step.msh backward_facing_step.xml

///////////////////////////////////////////////////////////////////////////////
// Basic variables
// Defined only if not given on the command line as 
//     gmsh -setnumber L1 2 -setnumber L2 10 ...

If (!Exists(L1)) L1 =    2; EndIf
If (!Exists(L2)) L2 =    8; EndIf
If (!Exists(H1)) H1 =    1; EndIf
If (!Exists(H2)) H2 =    1; EndIf
If (!Exists(lc)) lc = 0.25; EndIf

///////////////////////////////////////////////////////////////////////////////
// Points
Point(1) = {-L1,  H1, 0, lc};
Point(2) = {-L1,   0, 0, lc};
Point(3) = {  0,   0, 0, lc};
Point(4) = {  0, -H2, 0, lc};
Point(5) = { L2, -H2, 0, lc};
Point(6) = { L2,  H1, 0, lc};

///////////////////////////////////////////////////////////////////////////////
// Lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

///////////////////////////////////////////////////////////////////////////////
// Surfaces
Line Loop(10) = {1, 2, 3, 4, 5, 6};
Plane Surface(11) = {10};

///////////////////////////////////////////////////////////////////////////////
// Physical entities

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Surface(1) = {11};
