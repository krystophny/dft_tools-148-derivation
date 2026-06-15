SeedRandom[7];
rots = {IdentityMatrix[3],
        RotationMatrix[Pi/2,{0,0,1}], RotationMatrix[Pi,{0,0,1}], RotationMatrix[3 Pi/2,{0,0,1}],
        RotationMatrix[Pi,{1,0,0}], RotationMatrix[Pi,{0,1,0}],
        RotationMatrix[Pi,{1,1,0}], RotationMatrix[Pi,{1,-1,0}]};
ng=Length[rots];
tinv[R_]:=(R[[1,1]]R[[2,2]]-R[[1,2]]R[[2,1]])<0;
Dlf[l_,R_]:=Module[{e=EulerAngles[R],a,b,g},{a,b,g}=e;
   Table[WignerD[{l,mp,m},a,b,g],{mp,l,-l,-1},{m,l,-l,-1}]];

(* DIAGNOSTIC 1: pure orbital symmetrizer (no spin, no time reversal),
   does standard WignerD form a faithful rep / idempotent projector? *)
orbS[l_,M_]:=(1/ng) Sum[Dlf[l,rots[[i]]].M.ConjugateTranspose[Dlf[l,rots[[i]]]],{i,ng}];
orbTest[l_]:=Module[{d=2l+1,A,M,Ms,Mss},A=RandomComplex[{-1-I,1+I},{d,d}];M=(A+ConjugateTranspose[A])/2;
  Ms=orbS[l,M];Mss=orbS[l,Ms];
  Print["orbital-only l=",l,"  |S^2-S|=",ScientificForm[N@Max@Abs@Flatten[Mss-Ms],3]]];

(* DIAGNOSTIC 2: does the group even close as a rep under WignerD?  D(Ri)D(Rj)=D(Ri.Rj)? *)
closure[l_]:=Module[{worst=0},Do[worst=Max[worst,
   Max@Abs@Flatten@N[Dlf[l,rots[[i]]].Dlf[l,rots[[j]]]-Dlf[l,rots[[i]].rots[[j]]]]],{i,ng},{j,ng}];
   Print["rep closure l=",l,"  max|D(Ri)D(Rj)-D(RiRj)|=",ScientificForm[worst,3]]];

orbTest[1];orbTest[2];
closure[1];closure[2];
