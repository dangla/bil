R_1 = 0.1 ;
R_2 = 1. ;
lc1 = R_1/20. ;
lc2 = R_2/10 ;
lc4 = lc2*2. ;

Point(1) = {R_1,0,0,lc1} ;
Point(2) = {R_2,0,0,lc2} ;
Point(3) = {R_2,R_2,0,lc4} ;
Point(4) = {0,R_2,0,lc2} ;
Point(5) = {0,R_1,0,lc1} ;
Point(6) = {0,0,0,lc1} ;

Line(10) = {1,2} ;
Line(20) = {2,3} ;
Line(30) = {3,4} ;
Line(40) = {4,5} ;
Circle(50) = {5,6,1} ;

Line Loop(90) = {10,20,30,40,50} ;

Plane Surface(100) = {90} ;

Physical Line(1) = {10,20,30,40,50} ;
Physical Surface(1) = {100} ;