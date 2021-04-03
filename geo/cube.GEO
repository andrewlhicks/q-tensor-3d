//side length
len = DefineNumber[1, Name "Parameters/len" ];

//node density
ndens = DefineNumber[ 1, Name "Parameters/ndens" ];

Point(1) = {0, 0, 0, ndens};
Point(2) = {len, 0, 0, ndens};
Point(3) = {0, len, 0, ndens};
Point(4) = {0, 0, len, ndens};
Point(5) = {len, len, 0, ndens};
Point(6) = {0, len, len, ndens};
Point(7) = {len, 0, len, ndens};
Point(8) = {len, len, len, ndens};

Line(1) = {1,2};
Line(2) = {2,7};
Line(3) = {7,4};
Line(4) = {4,1};

Line(5) = {5,3};
Line(6) = {3,6};
Line(7) = {6,8};
Line(8) = {8,5};

Line(9) = {2,5};
Line(10) = {7,8};
Line(11) = {4,6};
Line(12) = {1,3};

Curve Loop(1) = {1,2,3,4};
Surface(1) = {1};

Curve Loop(2) = {5,6,7,8};
Surface(2) = {2};

Curve Loop(3) = {-9,2,10,8};
Surface(3) = {3};

Curve Loop(4) = {-12,-4,11,-6};
Surface(4) = {4};

Curve Loop(5) = {-3,10,-7,-11};
Surface(5) = {5};

Curve Loop(6) = {-1,12,-5,-9};
Surface(6) = {6};

Physical Surface(101) = {1, 2, 3, 4, 5, 6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};

Volume(1) = {1};
Physical Volume(102) = {1};
