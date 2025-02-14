W = 0.02;
H = 0.2;
w = W*0.5;
h = H*0.5;

Point(1) = {0,0,0} ;
Point(2) = {W,0,0} ;
Point(3) = {W,H,0} ;
Point(4) = {0,H,0} ;

Point(5) = {0.5*(W-w),0.5*(H-h),0} ;
Point(6) = {0.5*(W+w),0.5*(H-h),0} ;
Point(7) = {0.5*(W+w),0.5*(H+h),0} ;
Point(8) = {0.5*(W-w),0.5*(H+h),0} ;

lc1 = w*0.2 ;
lc2 = w*0.2 ;
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
