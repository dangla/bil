R_1 = 0.425 ;
R_2 = 1.225 ;
R_3 = 10.   ;

lc1 = 0.032 ;
lc2 = 0.032 ;
lc3 = 0.151212 ;

Point(1) = {R_1,0,0,lc1} ;
Point(2) = {R_2,0,0,lc2} ;
Point(3) = {R_3,0,0,lc3} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;

Physical Line(1) = {1,2} ;
Physical Point(1) = {3} ;