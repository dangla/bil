lc1 = 0.2 ;
lc2 = 0.1 ;

L = 1.5 ;
R = 0.5 ;

Point(1) = {-L,-L,0} ;
Point(2) = {L,-L,0} ;
Point(3) = {L,L,0} ;
Point(4) = {-L,L,0} ;

Point(5) = {-R,-R,0} ;
Point(6) = {R,-R,0} ;
Point(7) = {R,R,0} ;
Point(8) = {-R,R,0} ;


Characteristic Length{1:4} = lc1;
Characteristic Length{5:8} = lc2;

Line(11) = {1,2} ;
Line(12) = {2,3} ;
Line(13) = {3,4} ;
Line(14) = {4,1} ;
Line(15) = {5,6} ;
Line(16) = {6,7} ;
Line(17) = {7,8} ;
Line(18) = {8,5} ;

Line Loop(20) = {11,12,13,14} ;
Line Loop(21) = {15,16,17,18} ;

Plane Surface(100) = {20,-21} ;
Plane Surface(101) = {21} ;

Physical Surface(1) = {100} ;
Physical Surface(2) = {101} ;

Physical Line(1) = {11,12,13,14} ;

Physical Point(1) = {1};

Mesh.MshFileVersion = 2.2;
