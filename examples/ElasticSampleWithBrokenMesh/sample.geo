L = 0.1 ;
H = 0.2 ;

lc1 = 0.01 ;
lc2 = 0.001 ;


// Geometry Boundary
Point(1) = {0,0,0};
Point(2) = {L,0,0};
Point(3) = {L,H,0};
Point(4) = {0,H,0};
Point(5) = {0,0.51*H,0};
Point(6) = {0,0.49*H,0};
Point(7) = {0.2*L,0.5*H,0};


Characteristic Length{1:6} = lc1;
Characteristic Length{7} = lc2;


Line(10) ={1,2};
Line(20) ={2,3};
Line(30) ={3,4};
Line(40) ={4,5};
Line(50) ={5,7};
Line(60) ={7,6};
Line(70) ={6,1};

Line Loop(91) = {10,20,30,40,50,60,70};
Plane Surface(100) = {91};

Physical Surface(1) = {100};
Physical Line(1) = {10,30};
Physical Point(1)= {1,2} ;
