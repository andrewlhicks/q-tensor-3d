//outer radius
r1 = DefineNumber[0.2619, Name "Parameters/r1" ];

//outer radius
r2 = DefineNumber[1, Name "Parameters/r2" ];

//x-coord of center
x = DefineNumber[ 0, Name "Parameters/x" ];

//y-coord of center
y = DefineNumber[ 0, Name "Parameters/y" ];

//z-coord of center
z = DefineNumber[ 0, Name "Parameters/z" ];

//node density
ndens = DefineNumber[ 1, Name "Parameters/ndens" ];

Point(1) = {x, y, z, ndens};
Point(2) = {x+r2, y, z, ndens};
Point(3) = {x-r2, y, z, ndens};
Point(4) = {x, y+r2, z, ndens};
Point(5) = {x, y-r2, z, ndens};
Point(6) = {x, y, z+r2, ndens};
Point(7) = {x, y, z-r2, ndens};

Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
Circle(3) = {3, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {2, 1, 6};
Circle(6) = {6, 1, 3};
Circle(7) = {3, 1, 7};
Circle(8) = {7, 1, 2};
Circle(9) = {4, 1, 6};
Circle(10) = {6, 1, 5};
Circle(11) = {5, 1, 7};
Circle(12) = {7, 1, 4};

Curve Loop(1) = {2, 7, 12};
Surface(1) = {1};

Curve Loop(2) = {2, -6, -9};
Surface(2) = {2};

Curve Loop(3) = {3, -10, 6};
Surface(3) = {3};

Curve Loop(4) = {3, 11, -7};
Surface(4) = {4};

Curve Loop(5) = {4, -8, -11};
Surface(5) = {5};

Curve Loop(6) = {4, 5, 10};
Surface(6) = {6};

Curve Loop(7) = {1, 9, -5};
Surface(7) = {7};

Curve Loop(8) = {1, -12, 8};
Surface(8) = {8};

Point(8) = {x+r1, y, z, ndens};
Point(9) = {x-r1, y, z, ndens};
Point(10) = {x, y+r1, z, ndens};
Point(11) = {x, y-r1, z, ndens};
Point(12) = {x, y, z+r1, ndens};
Point(13) = {x, y, z-r1, ndens};

Circle(13) = {8, 1, 10};
Circle(14) = {10, 1, 9};
Circle(15) = {9, 1, 11};
Circle(16) = {11, 1, 8};
Circle(17) = {8, 1, 12};
Circle(18) = {12, 1, 9};
Circle(19) = {9, 1, 13};
Circle(20) = {13, 1, 8};
Circle(21) = {10, 1, 12};
Circle(22) = {12, 1, 11};
Circle(23) = {11, 1, 13};
Circle(24) = {13, 1, 10};

Curve Loop(9) = {14, 19, 24};
Surface(9) = {9};

Curve Loop(10) = {14, -18, -21};
Surface(10) = {10};

Curve Loop(11) = {15, -22, 18};
Surface(11) = {11};

Curve Loop(12) = {15, 23, -19};
Surface(12) = {12};

Curve Loop(13) = {16, -20, -23};
Surface(13) = {13};

Curve Loop(14) = {16, 17, 22};
Surface(14) = {14};

Curve Loop(15) = {13, 21, -17};
Surface(15) = {15};

Curve Loop(16) = {13, -24, 20};
Surface(16) = {16};

Physical Surface(101) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Surface(102) = {9, 10, 11, 12, 13, 14, 15, 16};

Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Surface Loop(2) = {9, 10, 11, 12, 13, 14, 15, 16};

Volume(1) = {1, 2};
Physical Volume(103) = {1};
