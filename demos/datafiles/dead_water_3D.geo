////////////////////////////////////////////
// Gmsh project created on Wed Nov  1 18:39 2017
SetFactory("OpenCASCADE");


////////////////////////////////////////////
// Parameters

// Parameters of the geometry
DefineConstant[ H = {  35.0, Name "Parameters/Total depth" } ];
DefineConstant[ h = {   7.0, Name "Parameters/Upper layer depth" } ];
DefineConstant[ L = { 200.0, Name "Parameters/Domain length" } ];
DefineConstant[ B = {  50.0, Name "Parameters/Domain breadth" } ];
DefineConstant[ C = {  65.0, Name "Parameters/Hull inlet offset amidships" } ];
// Hull:
DefineConstant[ l = {  30.0, Name "Parameters/Hull length" } ];
DefineConstant[ w = {  11.0, Name "Parameters/Hull width" } ];
DefineConstant[ d = {   5.0, Name "Parameters/Draught amidships" } ];

// Mesh params
DefineConstant[ lc_fine   = { 2.0, Name "Parameters/LC fine" } ];
DefineConstant[ lc_med    = { 5.0, Name "Parameters/LC medium" } ];
DefineConstant[ lc_course = {15.0, Name "Parameters/LC course" } ];
DefineConstant[ fs_eps    = { 5.0, Name "Parameters/Fine mesh free surface thickness" } ];


////////////////////////////////////////////
// Geometry

// Mesh control parameter for geometry construction
Q = lc_fine * 0.8;
Qm = lc_med * 0.8;

// Create hull out of ellipses
N = Ceil(d / Q);
dz = d / (N - 0.95);
For i In {0:N+1}
    // Make ellipses
    z = dz - i * dz;
    If (i > 1)
        z = z + dz * (i / (N + 1))^3;
    EndIf
    radius_x = l / 2 * Sqrt(1 + z / d);
    radius_y = w / 2 * Sqrt(1 + z / d);
    e1 = newl;
    Ellipse(e1) = {C, 0, z, radius_x, radius_y};
    e2 = newl;
    Line Loop(e2) = {e1};
    ellipses[i] = e2;
EndFor
For i In {0:N}
    // Make hull pieces
    surfs = ThruSections{ellipses[i], ellipses[i + 1]};
    top = news; Plane Surface(top) = {ellipses[i]};
    bot = news; Plane Surface(bot) = {ellipses[i + 1]};
    pshell = news;
    Surface Loop(pshell) = {top, surfs[0], bot};
    hull_pieces[i] = newv;
    Volume(hull_pieces[i]) = {pshell};
EndFor

// Create water domain
ocean = newv; Box(ocean) = {0, -B/2, -H,   L, B, H};

For i In {0:N}
	// Delete hull shape from water domain
	ocean_new = newv;
	BooleanDifference(ocean_new) = { Volume{ocean}; Delete; }{ Volume{hull_pieces[i]}; Delete; };
	ocean = ocean_new;
EndFor

// Remove half the domain (negative y-axis part)
half = newv; Box(half) = {0, -B/2, -H,   L, B/2, H};
ocean_new = newv;
BooleanDifference(ocean_new) = { Volume{ocean}; Delete; }{ Volume{half}; Delete; };
ocean = ocean_new;

// Make mesh conform to initial free surface. The free surface mesh conforming
// line is inset 1 mesh cell from the boundary to avoid mesh degeneration there
s1 = news;
Rectangle(s1) = {Qm, Qm, -h,     L - 2 * Qm,  B / 2 - 2 * Qm};
If (h < d)
    // Remove ellipse around the hull to avoid surfce/hull intersection
    radius_x = l / 2 * Sqrt(1 - h / (d + 2 * Q));
    radius_y = w / 2 * Sqrt(1 - h / (d + 2 * Q));
    e1 = newl;
    Ellipse(e1) = {C, 0, -h, radius_x, radius_y};
    e2 = newl;
    Line Loop(e2) = {e1};
    s2 = news;
    Plane Surface(s2) = {e2};
    s3 = news;
    BooleanDifference(s3) = { Surface{s1}; Delete; }{ Surface{s2}; Delete; };
Else
    radius_x = Q / 2;
    s3 = s1;
EndIf
Surface {s3} In Volume {ocean}; // Conform to this plane surface

// Mesh conforming line in the inlet plane
p1 = newp; Point(p1) = {0,   0 + Qm, -h};
p2 = newp; Point(p2) = {0, B/2 - Qm, -h};
l1 = newl; Line(l1) = {p1, p2};
Line{l1} In Surface{1}; // FIXME: surface number extracted from GUI

// Mesh conforming line in the outlet plane
p1 = newp; Point(p1) = {L,   0 + Qm, -h};
p2 = newp; Point(p2) = {L, B/2 - Qm, -h};
l1 = newl; Line(l1) = {p1, p2};
Line{l1} In Surface{6}; // FIXME: surface number extracted from GUI

// Mesh conforming lines in the center plane
p1 = newp; Point(p1) = {0 + Qm, 0, -h};
p2 = newp; Point(p2) = {C - radius_x, 0, -h};
p3 = newp; Point(p3) = {C + radius_x, 0, -h};
p4 = newp; Point(p4) = {L - Qm, 0, -h};
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p3, p4};
Line{l1} In Surface{2}; // FIXME: surface number extracted from GUI
Line{l2} In Surface{2}; // FIXME: surface number extracted from GUI

// Mesh conforming line in the off-center plane
p1 = newp; Point(p1) = {0 + Qm, B/2, -h};
p2 = newp; Point(p2) = {L - Qm, B/2, -h};
l1 = newl; Line(l1) = {p1, p2};
Line{l1} In Surface{4}; // FIXME: surface number extracted from GUI


////////////////////////////////////////////
// Mesh cell size fields:

// Max depth of hull and depth of pycnocline
If (d < h)
  dh = h;
Else
  dh = d;
EndIf

// Approximate radius of hull
r = l^2 / (8 * d) + d / 2;

// Fine mesh near the hull
Field[10] = Box;
Field[10].XMin = C - l * 0.8;
Field[10].XMax = C + l * 0.8;
Field[10].YMin = 0 - w * 1.0;
Field[10].YMax = 0 + w * 1.0;
Field[10].ZMin = - d * 1.8;
Field[10].ZMax = 0;
Field[10].VIn = lc_fine;
Field[10].VOut = lc_course;

// Fine mesh near the pycnocline
Field[20] = Box;
Field[20].XMin = C - l * 1.5;
Field[20].XMax = C + l * 1.5;
Field[20].YMin = -B;
Field[20].YMax =  B;
Field[20].ZMin = -h - fs_eps;
Field[20].ZMax = -h + fs_eps;
Field[20].VIn = lc_fine;
Field[20].VOut = lc_course;

// Medium mesh near the pycnocline
Field[21] = Box;
Field[21].XMin = 0;
Field[21].XMax = L;
Field[21].YMin = -B;
Field[21].YMax =  B;
Field[21].ZMin = -h - 2*fs_eps;
Field[21].ZMax = -h + 2*fs_eps;
Field[21].VIn = lc_med;
Field[21].VOut = lc_course;

// Fine mesh in the immediate wake of the hull
Field[30] = Box;
Field[30].XMin = C - l * 0.25;
Field[30].XMax = C + l * 1.75;
Field[30].YMin = -1.5 * w;
Field[30].YMax =  1.5 * w;
Field[30].ZMin = -1.5 * h;
Field[30].ZMax =  0;
Field[30].VIn = lc_fine;
Field[30].VOut = lc_course;

// Overall mesh grading (lc_course on sea floor and lc_med around the free surface)
zm = -h - fs_eps * 3;
Field[40] = MathEval;
Field[40].F = Sprintf("((%g)*(%g) + (%g)*(%g) - z*((%g) - (%g)))/((%g) + (%g))", H, lc_med, lc_course, zm, lc_course, lc_med, H, zm);
// Not finer than this
Field[48] = MathEval;
Field[48].F = Sprintf("%g", lc_med);
// Take the max
Field[49] = Max;
Field[49].FieldsList = {40, 48};

// Medium/Fine mesh around the fine zones
Field[50] = Box;
Field[50].XMin = C - l * 2.0;
Field[50].XMax = C + l * 2.0;
Field[50].YMin = -B;
Field[50].YMax =  B;
Field[50].ZMin = -2.0 * d;
Field[50].ZMax =  0.0;
Field[50].VIn = (lc_fine + lc_med) / 2.0;
Field[50].VOut = lc_course;

// Overall mesh grading
Field[90] = MathEval;
Field[90].F = Sprintf("(((x - %g)/%g)^2 + (z/%g)^2)^0.5", C, L, H);
Field[99] = Threshold;
Field[99].DistMax = 2.00;
Field[99].DistMin = 0.25;
Field[99].IField = 90;
Field[99].LcMax = lc_course;
Field[99].LcMin = lc_med;

// The resulting mesh size field is the minimum of the above fields
Field[100] = Min;
Field[100].FieldsList = {10, 20, 21, 30, 50, 99};
Background Field = 100;


////////////////////////////////////////////
// Physical domains

Physical Volume(100) = { ocean }; // The fluid domain


// The end, code below is added by gmsh GUI
////////////////////////////////////////////
