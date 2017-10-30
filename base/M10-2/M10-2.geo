L  = 0.01 ;
H  = 0.2 ;
lc = 0.002 ;

Point(1) = {0,0,0,lc} ;
Point(2) = {L,0,0,lc} ;
Point(3) = {L,H,0,lc} ;
Point(4) = {0,H,0,lc} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Line Loop(5) = {1,2,3,4} ;

Plane Surface(10) = {5} ;

Physical Surface(1) = {10} ;
Physical Line(1) = {1} ;
