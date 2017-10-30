R_1 = 3. ;
R_2 = 60. ;
R_3 = 4. ;
lc1 = 0.15 ;
lc2 = 5. ;
lc3 = lc1 ;
lc4 = lc2*1.5 ;

Point(1) = {R_1,0,0,lc1} ;
Point(2) = {R_2,0,0,lc2} ;
Point(3) = {R_2,R_2,0,lc4} ;
Point(4) = {0,R_2,0,lc2} ;
Point(5) = {0,R_1,0,lc1} ;
Point(6) = {0,0,0,lc1} ;
Point(7) = {R_3,0,0,lc3} ;
Point(8) = {0,R_3,0,lc3} ;
Point(9) = {R_1/Sqrt(2.),R_1/Sqrt(2.),0,lc1} ;
Point(10) = {R_3/Sqrt(2.),R_3/Sqrt(2.),0,lc3} ;

Line(10) = {7,2} ;
Line(11) = {1,7} ;
Line(20) = {2,3} ;
Line(30) = {3,4} ;
Line(40) = {4,8} ;
Line(41) = {8,5} ;
Circle(50) = {5,6,9} ;
Circle(51) = {9,6,1} ;
Circle(80) = {7,6,10} ;
Circle(81) = {10,6,8} ;
Line(85) = {9,10} ;
Line(86) = {10,3} ;

Line Loop(90) = {10,20,-86,-80} ;
Line Loop(91) = {11,80,-85,51} ;
Line Loop(92) = {30,40,-81,86} ;
Line Loop(93) = {41,50,85,81} ;

Plane Surface(100) = {90} ;
Plane Surface(200) = {91} ;
Plane Surface(300) = {92} ;
Plane Surface(400) = {93} ;

Physical Line(1) = {10,11,20,30,40,41,50,51} ;
Physical Surface(1) = {100,200,300,400} ;