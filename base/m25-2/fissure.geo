E_b = 0.005 ;
L_b = 0.005 ;
e_f = 0.00025 ;
h_f = 0.0005 ;

lc1 = 0.0005 ;
lc2 = 0.0001 ;

Point(1) = {0,0,0,lc2} ;
Point(2) = {h_f,0,0,lc2} ;
Point(3) = {h_f,e_f,0,lc2} ;
Point(4) = {E_b,e_f,0,lc1} ;
Point(5) = {E_b,L_b,0,lc1} ;
Point(6) = {0,L_b,0,lc2} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,5} ;
Line(5) = {5,6} ;
Line(6) = {6,1} ;

Line Loop(10) = {1,2,3,4,5,6} ;
Plane Surface(100) = {10} ;

Physical Line(1) = {1,2,3,4,5} ;
Physical Line(2) = {6} ;
Physical Surface(1) = {100} ;
