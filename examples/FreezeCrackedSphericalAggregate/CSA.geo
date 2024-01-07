D = 0.008;  // Diameter of the aggregate
e = 0.0003; // Thickness of the coated paste
t = 1.396e-5; // Thickness of the crack

R  = D/2;
Re = R + e;
theta = t/R ;
costheta = Cos(theta);
sintheta = Sin(theta);

Point(1) = {0,-R,0} ;
Point(3) = {0,R,0} ;
Point(4) = {0,Re,0} ;
Point(2) = {0,-Re,0} ;

Point(5) = {R*costheta,R*sintheta,0} ;
Point(6) = {Re*costheta,Re*sintheta,0} ;
Point(7) = {R*costheta,R*(-sintheta),0} ;
Point(8) = {Re*costheta,Re*(-sintheta),0} ;

Point(10) = {0,0,0} ;


// Size of the meshing
lc1 = R/10 ;
lc2 = lc1/5 ;
Characteristic Length{1:4} = lc1;
Characteristic Length{5:8} = lc2;


Line(33) = {5,6} ;
Line(44) = {4,3} ;
Line(55) = {8,7} ;
Line(66) = {1,2} ;
Line(77) = {3,1};

Circle(111) = {7,10,5} ;
Circle(161) = {1,10,7} ;
Circle(171) = {2,10,8} ;
Circle(121) = {5,10,3} ;
Circle(131) = {6,10,4} ;


Line Loop(1000) = {77,121,111,161};

Line Loop(2000) = {-121,44,131,33};
Line Loop(3000) = {-161,55,171,66};


Surface(100) = {1000};
Surface(200) = {2000};
Surface(300) = {3000};

Physical Point(1) = {1} ;

Physical Line(1)={111,77};
Physical Line(2)={44,33,131,66,171,55};

Physical Surface(1) = {100} ;
Physical Surface(2) = {200,300} ;

Mesh.MshFileVersion =2.2;
