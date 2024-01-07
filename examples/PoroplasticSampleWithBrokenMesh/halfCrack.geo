L = 1. ;
H = 1. ;

lc1 = 0.1 ;
lc2 = 0.001 ;


// Geometry Boundary
Point(1) = {   0,   0,0};
Point(2) = {   L,   0,0};
Point(3) = {   L,.5*H,0};
Point(4) = {   L,   H,0};
Point(5) = {   0,   H,0};
Point(6) = {   0,.5*H,0};
Point(7) = {.5*L,.5*H,0};


Characteristic Length{1:6} = lc1;
Characteristic Length{7} = lc2;


Line(10) ={1,2};
Line(20) ={2,3};
Line(30) ={3,4};
Line(40) ={4,5};
Line(50) ={5,6};
Line(60) ={6,1};

Line(70) ={6,7}; // <-- crack top lip
Line(80) ={3,7}; // <-- dummy divider



//Make line loops and surfaces
Line Loop(91) = {10,20,80,-70,60}; // Bottom half
Line Loop(92) = {30,40,50,70,-80}; // Top half
Plane Surface(1000) = {91};
Plane Surface(2000) = {92};

//Assign physical IDs (numbers picked just so they will be obvious in .msh file)
Physical Surface(1) = {1000}; // bottom will have ID 1
Physical Surface(2) = {2000}; // top will have ID 2
Physical Line(3) = {70}; // crack will have ID 3

// more physical ids on top and bottom for boundary conditions in a solver
Physical Line(1) = {20,30,50,60};
Physical Point(1)= {1,2,6,7} ; // save node 7 at crack tip
