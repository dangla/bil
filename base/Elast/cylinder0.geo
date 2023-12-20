InnerRadius = 0.1 ;
OuterRadius = 0.2 ;

WidthLayer1 = 0.05 ;
WidthLayer2 = 0.05 ;
WidthLayer3 = 0.05 ;

r = InnerRadius ;
R = OuterRadius ;

Z0 = 0 ;
Z1 = Z0 + WidthLayer1 ;
Z2 = Z1 + WidthLayer2 ;
Z3 = Z2 + WidthLayer3 ;

lc1 = 1. ;
lc2 = 1. ;
lc3 = 2. ;

Point(1) = {r,Z0,0,lc2} ;
Point(2) = {r,Z1,0,lc1} ;
Point(3) = {r,Z2,0,lc1} ;
Point(4) = {r,Z3,0,lc2} ;
Point(5) = {R,Z3,0,lc3} ;
Point(6) = {R,Z2,0,lc3} ;
Point(7) = {R,Z1,0,lc3} ;
Point(8) = {R,Z0,0,lc3} ;

Line(10) = {1,2} ;
Line(20) = {2,3} ;
Line(30) = {3,4} ;
Line(40) = {4,5} ;
Line(50) = {5,6} ;
Line(60) = {6,7} ;
Line(70) = {7,8} ;
Line(80) = {8,1} ;
Line(90) = {2,7} ;
Line(100) = {6,3} ;

Line Loop(200) = {-10,-90,-70,-80} ;
Line Loop(300) = {-20,100,-60,90} ;
Line Loop(400) = {-30,-40,-50,-100} ;

Plane Surface(500) = {200} ;
Plane Surface(600) = {300} ;
Plane Surface(700) = {400} ;

Physical Surface(1) = {500} ;
Physical Surface(2) = {600} ;
Physical Surface(3) = {700} ;

Physical Line(1) = {-10,-70,-80} ;
Physical Line(2) = {-20,-60} ;
Physical Line(3) = {-30,-40,-50} ;

