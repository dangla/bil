// usage: gmsh -2 composite1.geo -format msh2

Include "inclusion.geo" ;

Rotate{{0,0,1},{0,0,0},0.5*Pi}{Duplicata{Surface{100,101} ;} }
Rotate{{0,0,1},{0,0,0},Pi}{Duplicata{Surface{100,101} ;} }
Rotate{{0,0,1},{0,0,0},1.5*Pi}{Duplicata{Surface{100,101} ;} }


//Symmetry{1,0,0,0} {Duplicata{Surface{100,101} ;} }

//Symmetry{0,1,0,0} {Duplicata{Surface{100,101,102,108} ;} }

/*
* Compound Line(200) = {13,112} ;
* Compound Line(201) = {104,122} ;
* Compound Line(300) = {14,105} ;
* Compound Line(301) = {123,113} ;
*/

Physical Point(2) = {1} ;
Physical Point(1) = {2} ;
Physical Line(1) = {13,14,104,105,114,115,124,125} ;
//Physical Line(1) = {13,112,104,122,14,105,123,113} ;


Physical Surface(1) = {100,102,112,122} ;
Physical Surface(2) = {101,108,118,128} ;
//Physical Surface(1) = {100,102,120,110} ;
//Physical Surface(2) = {101,108,126,116} ;

//Translate{-L,0,0} {Symmetry{1,0,0,0} {Duplicata{Surface{100,101} ;} } }

