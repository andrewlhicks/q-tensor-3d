W = DefineNumber[ 1, Name "Parameters/width" ];
H = DefineNumber[ 0.5, Name "Parameters/height" ];
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {W, 0, 0, 1.0};
//+
Point(3) = {W, W, 0, 1.0};
//+
Point(4) = {0, W, 0, 1.0};
//+
Point(5) = {0, 0, H, 1.0};
//+
Point(6) = {W, 0, H, 1.0};
//+
Point(7) = {W, W, H, 1.0};
//+
Point(8) = {0, W, H, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Line(9) = {1, 5};
//+
Line(10) = {2, 6};
//+
Line(11) = {3, 7};
//+
Line(12) = {4, 8};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, 10, -5, -9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {2, 11, -6, -10};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {3, 12, -7, -11};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {4, 9, -8, -12};
//+
Plane Surface(6) = {6};
//+
Physical Surface("Top/bottom", 1) = {1, 2};
//+
Physical Surface("Sides", 2) = {3, 4, 5, 6};
//+
Surface Loop(1) = {3, 1, 4, 5, 6, 2};
//+
Volume(1) = {1};
//+
Physical Volume("Volume", 3) = {1};
