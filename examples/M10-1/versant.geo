lc1 = 0.1 ;

Point(1) = {0,0,0,lc1} ;
Point(2) = {0,0.5,0,lc1} ;
Point(3) = {0,1,0,lc1} ;
Point(4) = {-10,3,0,lc1} ;
Point(5) = {-10,2,0,lc1} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,5} ;
Line(5) = {5,1} ;

Line Loop(5) = {1,2,3,4,5} ;

Plane Surface(10) = {5} ;

Physical Surface(1) = {10} ;
Physical Line(1) = {1} ;
