

AttachSpec("EarlierCode/magma.spec");
Attach("Modular.m"); 

G:=GL2Cartan(-8,11);
G0:=GL2CartanNormalizer(-8,11);
G subset G0;

GL2:=GL2Ambient(8);

G:=sub<GL2|[[ 7, 7, 0, 1 ],[ 1, 2, 0, 1 ],[ 5, 0, 0, 5 ],[ 7, 6, 0, 7 ],[ 1, 4, 0, 1 ],[ 5, 0, 0, 1 ],[ 3, 2, 0, 3 ]]>;

G0:=sub<GL2|[[ 3, 2, 2, 3 ],[ 3, 0, 0, 3 ],[ 1, 0, 2, 1 ],[ 1, 4, 0, 1 ],[ 3, 3, 0, 1 ],[ 1, 0, 4, 1 ],[ 1, 2, 0, 1 ],[ 5, 4, 4, 5 ],[ 5, 0, 0, 1 ]]>;


/*
GL2:=GL2Ambient(48);
gens := [[1, 14, 46, 43], [11, 34, 2, 17], [13, 2, 14, 23], [17, 12, 12, 13], [25, 32, 8, 17], [41, 20, 44, 37]];
G:=sub<GL2|gens>;

G0:=Normalizer(GL2,G);
//gens := [[5, 36, 12, 25], [27, 16, 16, 27], [31, 30, 24, 35], [35, 32, 26, 13], [37, 14, 40, 25], [39, 2, 40, 3]];
//G0:=sub<GL2|gens>;

G0:=Normalizer(GL2,G);
for i in [1..100] do
    G1:=sub<GL2|Generators(G) join {Random(G0)}>;
    g:=GL2Genus(G1);
    X1:=CreateModularCurveRec(G1);
    X1`CPname;
    if X1`CPname eq "48B11" then break i; end if;
end for;
G0:=G1;
*/
G:=GL2Borel(23);
G0:=GL2Ambient(1);

//Z:=sub<GL2Ambient(17)| Generators(SL2Ambient(17)) join {GL2Ambient(17)![3,0,0,3]}>;
//G:=G meet Z;
//G0:=G0 meet Z;

G:=GL2Cartan(-4,13);
G0:=GL2CartanNormalizer(-4,13);
G0:=GL2Ambient(1);

X:=CreateModularCurveRec(G);
X0:=CreateModularCurveRec(G0);

X:=FindModelOfXG(X);
X0:=FindModelOfXG(X0);

// Attach("Modular.m");  FindMorphism(X,X0);
//FindRationalFunction(X,X0);
FindMorphism(X,X0);