
AttachSpec("../Modular.spec");

load "groups.m";


for r in adelicimagedata do
    GL2:=GL2Ambient(r`level);
    G:=sub<GL2|r`generators>;
    if {GL2Genus(M): M in IntermediateSubgroups(GL2,G) *} then
        r`label;
    end if;
end for;


// 24.72.2.ii.1

/*

[ 235, 5, 2, 36 ]
[ 239, 7, 10, 45 ]
[ 239, 13, 10, 45 ]
[ 282, 5, 2, 18 ]
[ 301, 3, 20, 55 ]
[ 301, 7, 20, 55 ]


for r in groups do
    GL2:=GL2Ambient(r`N);
    G:=sub<GL2|r`gens>;
    Gt:=sub<GL2|{Transpose(g): g in Generators(G)}>;
    assert IsConjugate(GL2,G,Gt);
end for;
*/

for r in groups do
    GL2:=GL2Ambient(r`N);
    G:=sub<GL2|r`gens>;
    H:=GL2CommutatorSubgroup(G);
    SL2Index(H);
end for;


for i in [1..#groups] do
    r:=groups[i];
    G:=sub<GL2Ambient(r`N)|r`gens>;
    for p in [p:p in PrimesUpTo(30) | r`N mod p ne 0] do
        if GL2PointCount(G,p) eq 0 then 
            g:=GL2Genus(G);
            [i,p,g,r`N];
        end if;
    end for;
end for;