////////////////////////////////////////////
// Created on 2018-04-20 by Tormod Landet
// Surface piering bottom mounted cylinder
// in a wave flume
SetFactory("OpenCASCADE");


////////////////////////////////////////////
// Parameters

// Geometrical parameters
DefineConstant[ H = { 1.00, Name "Parameters/Height" } ];
DefineConstant[ R = { 0.06, Name "Parameters/Cylinder radius" } ];
DefineConstant[ L = { 2.00, Name "Parameters/Domain length" } ];
DefineConstant[ B = { 0.50, Name "Parameters/Domain breadth" } ];
DefineConstant[ C = { 0.75, Name "Parameters/Cylinder dist from inlet" } ];

// Mesh parameters
Q = 5; // Coursener, Q=1 gives relatively fine mesh
DefineConstant[ sw_pos    = {0.60, Name "Parameters/Still water height" } ];
DefineConstant[ lc_fine   = {0.01 * Q, Name "Parameters/LC fine" } ];
DefineConstant[ lc_med    = {0.03 * Q, Name "Parameters/LC medium" } ];
DefineConstant[ lc_course = {0.07 * Q, Name "Parameters/LC course" } ];


////////////////////////////////////////////
// Geometry

// The domain
flume = newv; Box(flume) = {-C, -B/2, 0,   L, B, H};

//R = 0;
If (R > 0)
    // The cylinder
    cylinder = newv; Cylinder(cylinder) = {0, 0, 0, 0, 0, H, R};
    
    // Delete cylinder from domain
    flume_new = newv;
    BooleanDifference(flume_new) = { Volume{flume}; Delete; }{ Volume{cylinder}; Delete; };
    flume = flume_new;
EndIf

// Remove half the domain (negative y-axis part)
half = newv; Box(half) = {-C, -B/2, 0,   L, B/2, H};
flume_new = newv;
BooleanDifference(flume_new) = { Volume{flume}; Delete; }{ Volume{half}; Delete; };
flume = flume_new;

////////////////////////////////////////////
// Mesh cell size fields:

// Fine mesh near the cylinder
Field[10] = Box;
Field[10].XMin = 0 - R * 2.0;
Field[10].XMax = 0 + R * 3.0;
Field[10].YMin = 0 - R * 2.0;
Field[10].YMax = 0 + R * 2.0;
Field[10].ZMin = 0;
Field[10].ZMax = H;
Field[10].VIn = lc_fine;
Field[10].VOut = lc_course;

// Medium mesh near the cylinder
Field[20] = Box;
Field[20].XMin = 0 - R * 4.0;
Field[20].XMax = 0 + R * 6.0;
Field[20].YMin = 0 - R * 4.0;
Field[20].YMax = 0 + R * 4.0;
Field[20].ZMin = 0;
Field[20].ZMax = H;
Field[20].VIn = lc_med;
Field[20].VOut = lc_course;

// Overall mesh grading
Field[90] = MathEval;
Field[90].F = Sprintf("(((x - %g)/%g)^2 + ((z - %g)/%g)^2)^0.5", 0, 5*L, sw_pos, H);
Field[99] = Threshold;
Field[99].DistMax = 10*R;
Field[99].DistMin =  1*R;
Field[99].IField = 90;
Field[99].LcMax = lc_course;
Field[99].LcMin = lc_med;

// The resulting mesh size field is the minimum of the above fields
Field[100] = Min;
Field[100].FieldsList = {10, 20, 21, 30, 50, 99};
Background Field = 100;


////////////////////////////////////////////
// Physical domains

Physical Volume(100) = { flume }; // The fluid domain


// The end, code below is added by gmsh GUI
////////////////////////////////////////////
