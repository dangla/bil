// Valeurs initiales
R_a = 0.02 ;
H_a = 0.02 ;

R_t = R_a + 0.0175 ;
H_t = H_a + 0.014 ;

r   = 0.006 ;
l   = 0.03 ;

L   = 0.1 ;
D   = 0.08 ;

lc1 = 0.001 ;
lc2 = 0.003 ;

Point(1) = {0,0,0,lc1} ;
Point(2) = {0,H_a,0,lc1} ;
Point(3) = {r,H_a,0,lc1} ;
Point(4) = {R_a,H_a,0,lc1} ;
Point(5) = {R_a,0,0,lc1} ;
Point(6) = {0,H_t,0,lc1} ;
Point(7) = {r,H_t,0,lc1} ;
Point(8) = {R_t,H_t,0,lc1} ;
Point(9) = {R_t,0,0,lc1} ;
Point(10) = {0,H_t+l,0,lc1} ;
Point(11) = {r,H_t+l,0,lc1} ;

Point(12) = {0,L,0,lc2} ;
Point(13) = {D,L,0,lc2} ;
Point(14) = {D,0,0,lc2} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,5} ;
Line(5) = {5,1} ;

Line Loop(100) = {1,2,3,4,5} ;

Line(6) = {2,6} ;
Line(7) = {6,7} ;
Line(8) = {7,8} ;
Line(9) = {8,9} ;
Line(10) = {9,5} ;

Line Loop(101) = {6,7,8,9,10,-4,-3,-2} ;

Line(11) = {6,10} ;
Line(12) = {10,11} ;
Line(13) = {11,7} ;

Line Loop(102) = {11,12,13,-7} ;

Line(14) = {10,12} ;
Line(15) = {12,13} ;
Line(16) = {13,14} ;
Line(17) = {14,9} ;

Line Loop(103) = {14,15,16,17,-9,-8,-13,-12} ;

Plane Surface(1) = {100} ;
Plane Surface(2) = {101} ;
Plane Surface(3) = {102} ;
Plane Surface(4) = {103} ;

Physical Surface(1) = {1} ;
Physical Surface(2) = {2,3} ;
Physical Surface(3) = {4} ;
Physical Line(3) = {5,10,17,15} ;