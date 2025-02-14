L  = 0.002 ;
H  = 0.2 ;
lc = 0.002 ;

Point(1) = {0,0,0,lc} ;
Point(2) = {L,0,0,lc} ;
Point(3) = {L,L,0,lc} ;
Point(4) = {0,L,0,lc} ;
Point(5) = {0,0,H,lc} ;
Point(6) = {L,0,H,lc} ;
Point(7) = {L,L,H,lc} ;
Point(8) = {0,L,H,lc} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Line Loop(5) = {1,2,3,4} ;
Plane Surface(100) = {5} ;

Line(11) = {5,6} ;
Line(12) = {6,7} ;
Line(13) = {7,8} ;
Line(14) = {8,5} ;
Line Loop(15) = {11,12,13,14} ;
Plane Surface(110) = {15} ;

Line(21) = {1,5} ;
Line(22) = {2,6} ;
Line(23) = {3,7} ;
Line(24) = {4,8} ;

Line Loop(25) = {1,22,-11,-21} ;
Line Loop(26) = {2,23,-12,-22} ;
Line Loop(27) = {3,24,-13,-23} ;
Line Loop(28) = {4,21,-14,-24} ;
Plane Surface(120) = {25} ;
Plane Surface(130) = {26} ;
Plane Surface(140) = {27} ;
Plane Surface(150) = {28} ;

Surface Loop(160) = {-100,110,-120,-130,-140,-150} ;

Volume(1000) = {160} ;

Physical Surface(1) = {100} ;
Physical Volume(1) = {1000} ;

Mesh.MshFileVersion = 2.2;

