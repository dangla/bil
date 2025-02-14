lc1 = 3; //0.05 ;

L = 1.5 ;


Point(1) = {-L,-L,0} ;
Point(2) = {L,-L,0} ;
Point(3) = {L,L,0} ;
Point(4) = {-L,L,0} ;

Characteristic Length{1:4} = lc1;

Line(11) = {1,2} ;
Line(12) = {2,3} ;
Line(13) = {3,4} ;
Line(14) = {4,1} ;

Line Loop(20) = {11,12,13,14} ;

Plane Surface(100) = {20} ;

Physical Surface(1) = {100} ;

Physical Line(1) = {11,12,13,14} ;

Mesh.MshFileVersion = 2.2;
