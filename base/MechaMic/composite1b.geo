Include "inclusion.geo" ;

Symmetry{1,0,0,0} {Duplicata{Surface{100,101} ;} }

Symmetry{0,1,0,0} {Duplicata{Surface{100,101,102,108} ;} }

/*
* Compound Line(200) = {13,112} ;
* Compound Line(201) = {104,122} ;
* Compound Line(300) = {14,105} ;
* Compound Line(301) = {123,113} ;
*/
/*
Compound Line(500) = {104,122} ;
Compound Line(501) = {13,112} ;
Compound Line(600) = {123,113} ;
Compound Line(601) = {105,14} ;
*/

Compound Surface(1000) = {100,102,120,110} ;
Compound Surface(1001) = {101,108,126,116} ;

Physical Point(2) = {1} ;
//Physical Line(1) = {13,112,104,122,14,105,123,113} ;
Physical Line(1) = {500,501,600,601} ;
//Physical Surface(1) = {100,102,120,110} ;
//Physical Surface(2) = {101,108,126,116} ;
Physical Surface(1) = {1000} ;
Physical Surface(2) = {1001} ;

//Translate{-L,0,0} {Symmetry{1,0,0,0} {Duplicata{Surface{100,101} ;} } }

