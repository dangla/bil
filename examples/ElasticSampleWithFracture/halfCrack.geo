// usage: gmsh -2 halfCrack.geo -format msh2
// then   gmshcrack -f halfCrack.msh -c 70 -t 7 -s 200

L = 1. ;
H = 1. ;

lc1 = 0.1 ;
lc2 = 0.001 ;


// Geometry Boundary
Point(1) = {   0,   0,0,lc1};
Point(2) = {   L,   0,0,lc1};
Point(3) = {   L,.5*H,0,lc1};
Point(4) = {   L,   H,0,lc1};
Point(5) = {   0,   H,0,lc1};
Point(6) = {   0,.5*H,0,lc1};
Point(7) = {.5*L,.5*H,0,lc2};
//Point(8) = {   0,.5*H,0,lc1};
Line(10) ={1,2};
Line(20) ={2,3};
Line(30) ={3,4};
Line(40) ={4,5};
Line(50) ={5,6};
Line(60) ={6,1};

Line(70) ={6,7}; // <-- crack top lip
Line(80) ={3,7}; // <-- dummy divider
//Line(90) ={8,7}; // <-- crack bottom lip



//Make line loops and surfaces
Line Loop(91) = {10,20,80,-70,60}; // Bottom half
Line Loop(92) = {30,40,50,70,-80}; // Top half
Plane Surface(100) = {91};
Plane Surface(200) = {92};

//Assign physical IDs (numbers picked just so they will be obvious in .msh file)
Physical Surface(1) = {100}; // bottom will have ID 1
Physical Surface(2) = {200}; // top will have ID 2
Physical Line(3) = {70}; // crack will have ID 3

// more physical ids on top and bottom for boundary conditions in a solver
Physical Line(1) = {10,40,50,60};
Physical Point(1)= {1,7} ; // save node 7 at crack tip
