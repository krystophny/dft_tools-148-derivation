(* MRE for dft_tools #148 - decisive idempotency test. No Wien2k.
   The spinor matrices written by dmftproj must form a valid (anti)unitary
   rep of the magnetic point group, so the dft_tools symmetrizer
   S[M] = (1/N) sum_g  mat_g . (M or conj M) . mat_g^dag      (symmetry.py:150-162)
   must be an idempotent PROJECTOR: S[S[M]]=S[M]. A wrong spin phase breaks
   idempotency -> symmetrizing an already-symmetric M shifts its eigenvalues,
   exactly the reported failure. Whichever of {current, flipped} is idempotent
   is correct. *)
SeedRandom[20260615];
rots = {IdentityMatrix[3],
        RotationMatrix[Pi/2,{0,0,1}], RotationMatrix[Pi,{0,0,1}], RotationMatrix[3 Pi/2,{0,0,1}],
        RotationMatrix[Pi,{1,0,0}], RotationMatrix[Pi,{0,1,0}],
        RotationMatrix[Pi,{1,1,0}], RotationMatrix[Pi,{1,-1,0}]};
ng = Length[rots];
tinv[R_] := (R[[1,1]] R[[2,2]] - R[[1,2]] R[[2,1]]) < 0;       (* cos(beta)<0  setsym.f:148 *)
phaseOf[R_] := Module[{e=EulerAngles[R],a,b,g}, {a,b,g}=e; If[tinv[R], g-a, a+g]];

(* orbital block for given l: Wigner D^l(R^-1) to match dmat's inverse convention *)
Dl[l_,R_] := Module[{e=EulerAngles[R],a,b,g}, {a,b,g}=e;
   If[l==0, {{1}}, Table[Conjugate[WignerD[{l,mp,m},a,b,g]],{mp,l,-l,-1},{m,l,-l,-1}]]];

(* non-mixing output mat: spin-up block ephase*Dl, spin-dn block conj(ephase)*Dl *)
mkMat[l_,R_,sgn_] := Module[{ph=phaseOf[R], e=0, d=Dl[l,R], z},
   e=Exp[sgn I ph/2]; z=0 d; ArrayFlatten[{{e d, z},{z, Conjugate[e] d}}]];

savg[l_,sgn_,M_] := (1/ng) Sum[
   With[{g=mkMat[l,rots[[i]],sgn], ti=tinv[rots[[i]]]},
     g . If[ti, Conjugate[M], M] . ConjugateTranspose[g]], {i,ng}];

runTest[l_] := Module[{dim=2(2l+1), A, M0, res},
  A = RandomComplex[{-1-I,1+I},{dim,dim}]; M0=(A+ConjugateTranspose[A])/2;
  Print["---- l=",l,"  (spin x orbital dim = ",dim,") ----"];
  Do[Module[{Ms, Mss, idem, drift},
     Ms  = savg[l,sgn,M0];
     Mss = savg[l,sgn,Ms];
     idem  = Max@Abs@Flatten[N[Mss-Ms]];
     drift = Max@Abs[ Sort[Re@Eigenvalues[N@Ms]] - Sort[Re@Eigenvalues[N@Mss]] ];
     Print["  ", If[sgn==1,"CURRENT exp(+i ph/2)","FLIP    exp(-i ph/2)"],
           ":  idempotency |S^2-S| = ", ScientificForm[idem,3],
           "   eig-shift = ", ScientificForm[drift,3]];
     ], {sgn,{1,-1}}];
];
Print["group ops=",ng,"  time_inv=",Boole[tinv[#]&/@rots],"  phases/pi=",N[phaseOf[#]/Pi]&/@rots];
runTest[0];
runTest[1];
runTest[2];
