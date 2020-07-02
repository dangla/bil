lc1 = 0.05 ;
lc2 = 0.05 ;

L = 1 ;
R = 0.5 ;

Point(1) = {0,0,0,lc1} ;
Point(2) = {L,0,0,lc1} ;
Point(3) = {L,L,0,lc1} ;
Point(4) = {0,L,0,lc1} ;

Point(5) = {R,0,0,lc2} ;
Point(6) = {0,R,0,lc2} ;

Circle(10) = {5,1,6} ;
Line(11) = {1,5} ;
Line(12) = {5,2} ;
Line(13) = {2,3} ;
Line(14) = {3,4} ;
Line(15) = {4,6} ;
Line(16) = {6,1} ;

Line Loop(20) = {12,13,14,15,-10} ;
Line Loop(21) = {11,10,16} ;

Plane Surface(100) = {20} ;
Plane Surface(101) = {21} ;

