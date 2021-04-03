r = DefineNumber[ 1, Name "Parameters/inner radius" ];
R = DefineNumber[ 1.5, Name "Parameters/outer radius" ];

Point(1) = {0, 0, 0, 1.0};
Point(2) = {r, 0, 0, 1.0};
Point(3) = {0, r, 0, 1.0};
Point(4) = {0, 0, r, 1.0};
Point(5) = {-r, 0, 0, 1.0};
Point(6) = {0, -r, 0, 1.0};
Point(7) = {0, 0, -r, 1.0};
Point(8) = {R, R, R, 1.0};
Point(9) = {R, R, -R, 1.0};
Point(10) = {R, -R, R, 1.0};
Point(11) = {R, -R, -R, 1.0};
Point(12) = {-R, R, R, 1.0};
Point(13) = {-R, R, -R, 1.0};
Point(14) = {-R, -R, R, 1.0};
Point(15) = {-R, -R, -R, 1.0};

Circle(1) = {5, 1, 7};
Circle(2) = {7, 1, 2};
Circle(3) = {2, 1, 4};
Circle(4) = {4, 1, 5};
Circle(5) = {5, 1, 3};
Circle(6) = {3, 1, 2};
Circle(7) = {2, 1, 6};
Circle(8) = {6, 1, 5};
Circle(9) = {4, 1, 6};
Circle(10) = {6, 1, 7};
Circle(11) = {7, 1, 3};
Circle(12) = {3, 1, 4};

Line(13) = {15, 14};
Line(14) = {14, 12};
Line(15) = {12, 13};
Line(16) = {13, 15};
Line(17) = {10, 8};
Line(18) = {8, 9};
Line(19) = {9, 11};
Line(20) = {11, 10};
Line(21) = {14, 10};
Line(22) = {12, 8};
Line(23) = {13, 9};
Line(24) = {15, 11};

Curve Loop(1) = {1, 11, -5};
Surface(1) = {1};
Curve Loop(2) = {4, 5, 12};
Surface(2) = {2};
Curve Loop(3) = {3, -12, 6};
Surface(3) = {3};
Curve Loop(4) = {2, -6, -11};
Surface(4) = {4};
Curve Loop(5) = {1, -10, 8};
Surface(5) = {5};
Curve Loop(6) = {8, -4, 9};
Surface(6) = {6};
Curve Loop(7) = {3, 9, -7};
Surface(7) = {7};
Curve Loop(8) = {2, 7, 10};
Surface(8) = {8};

Curve Loop(9) = {19, 20, 17, 18};
Surface(9) = {9};
Curve Loop(10) = {19, -24, -16, 23};
Surface(10) = {10};
Curve Loop(11) = {16, 13, 14, 15};
Surface(11) = {11};
Curve Loop(12) = {14, 22, -17, -21};
Surface(12) = {12};
Curve Loop(13) = {13, 21, -20, -24};
Surface(13) = {13};
Curve Loop(14) = {22, 18, -23, -15};
Surface(14) = {14};

Physical Surface("Inner", 1) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Surface("Outer", 2) = {9, 10, 11, 12, 13, 14};

Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Surface Loop(2) = {9, 10, 11, 12, 13, 14};

Volume(1) = {1, 2};
Physical Volume(1) = {1};
