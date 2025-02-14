// taille caractéristique des éléments
cl1 = 1;

// points
Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {1, 1, 0, cl1};
Point(4) = {0, 1, 0, cl1};

// lignes
Line(11) = {1, 2};
Line(12) = {2, 3};
Line(13) = {3, 4};
Line(14) = {4, 1};

// surfaces
Line Loop(99) = {11, 12, 13, 14};
Plane Surface(100) = {99};

// entités physiques
Physical Line(1) = {11, 12, 13, 14};
Physical Surface(1) = {100};

Recombine Surface {100};

Mesh.MshFileVersion = 2.2;
