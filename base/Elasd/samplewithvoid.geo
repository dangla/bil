// Dimension of the sample
L = 1 ;
H = 1 ;
L05 = 0.5*L;
H05 = 0.5*H;

// Major and minor diameters of the elliptical void
a = 0.1 ;
b = 0.01 ;
// orientation of the axis
theta = 45.;
cost = Cos(theta*Pi/180);
sint = Sin(theta*Pi/180);
// Coordinates of the points in the major and minor axis
ax =  a*cost;
ay =  a*sint;
bx = -b*sint;
by =  b*cost;

// Characteristic lengths
lc1 = 0.05 ;
lc2 = 0.02 ;
lc3 = 0.005 ;

// Geometry Boundary
Point(1) = {0,0,0};
Point(2) = {L,0,0};
Point(3) = {L,H,0};
Point(4) = {0,H,0};

Point(5) = {L05   ,H05   ,0};
Point(6) = {L05-ax,H05-ay,0};
Point(7) = {L05+bx,H05+by,0};
Point(8) = {L05+ax,H05+ay,0};
Point(9) = {L05-bx,H05-by,0};


Line(10) = {1,2};
Line(20) = {2,3};
Line(30) = {3,4};
Line(40) = {4,1};


Ellipse(50) = {6,5,8,7};
Ellipse(51) = {7,5,8,8};
Ellipse(52) = {8,5,6,9};
Ellipse(53) = {9,5,8,6};


Line Loop(70) = {10,20,30,40};
Line Loop(71) = {50,51,52,53};


Characteristic Length{1,2,3,4} = lc1 ;
Characteristic Length{7,9} = lc2 ;
Characteristic Length{6,8} = lc3 ;


Plane Surface(100) = {70,71};


Physical Surface(1) = {100};
Physical Line(1) = {10,20,30,40,50,51,52,53};
Physical Point(1)= {1,4} ;
