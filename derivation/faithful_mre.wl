(* Faithful replication of dmftproj's real construction:
   dmat + d_matrix (setsym.f:694-776), timeinv_op (timeinv.f, isrt=0),
   and outputqmc's non-mixing exported mat. Tests whether the symmetrizer
   S[M]=(1/N) sum_g mat_g.(M or conjM).mat_g^dag (symmetry.py) is idempotent
   for current ephase=exp(+i ph/2) vs proposed flip exp(-i ph/2). *)
SeedRandom[2026];
ifac[n_]:=n!;
(* d_matrix exactly as setsym.f:735-776 *)
smalld[l_,m_,n_,b_]:=Module[{f1,s=0,f2,f3,f4,t},
  f1=(ifac[l+m] ifac[l-m])/(ifac[l+n] ifac[l-n]);
  Do[If[(l-m-t)>=0&&(l-n-t)>=0&&(t+n+m)>=0,
     f2=(ifac[l+n] ifac[l-n])/(ifac[l-m-t] ifac[m+n+t] ifac[l-n-t] ifac[t]);
     f3=If[(2l-m-n-2t)==0,1,Sin[b/2]^(2l-m-n-2t)];
     f4=If[(2t+n+m)==0,1,Cos[b/2]^(2t+n+m)];
     s=s+(-1)^(l-m-t) f2 f3 f4], {t,0,2l}];
  Sqrt[f1] s];
(* dmat exactly as setsym.f:694-732 : row index m, col index n, both -l..l *)
dmatMM[l_,a_,b_,c_,det_]:=Table[
  (If[det<-0.5,(-1)^l,1]) Exp[I n a] Exp[I m c] smalld[l,m,n,b],
  {m,-l,l},{n,-l,l}];
(* orbital time reversal, isrt=0: T[m,-m]=(-1)^m ; op: mat -> T.Conjugate[mat] *)
Tflip[l_]:=Table[If[mp==-m,(-1)^m,0],{mp,-l,l},{m,-l,l}];
timeinvOrb[l_,mat_]:=Tflip[l].Conjugate[mat];

(* D4 magnetic group *)
rots={IdentityMatrix[3],RotationMatrix[Pi/2,{0,0,1}],RotationMatrix[Pi,{0,0,1}],
      RotationMatrix[3Pi/2,{0,0,1}],RotationMatrix[Pi,{1,0,0}],RotationMatrix[Pi,{0,1,0}],
      RotationMatrix[Pi,{1,1,0}],RotationMatrix[Pi,{1,-1,0}]};
ng=Length[rots];
tinv[R_]:=(R[[1,1]]R[[2,2]]-R[[1,2]]R[[2,1]])<0;
ea[R_]:=EulerAngles[R];
det3[R_]:=Det[R];

buildMat[l_,R_,sgn_]:=Module[{e=ea[R],a,b,g,ph,rotl,rr,eph,z},
  {a,b,g}=e; ph=If[tinv[R],g-a,a+g];
  rotl=dmatMM[l,a,b,g,det3[R]];
  rr=If[tinv[R],timeinvOrb[l,rotl],rotl];
  eph=Exp[sgn I ph/2]; z=0 rr;
  ArrayFlatten[{{eph rr, z},{z, Conjugate[eph] rr}}]];

symmAvg[l_,sgn_,M_]:=Module[{mats=buildMat[l,#,sgn]&/@rots},
  (1/ng) Sum[mats[[i]].If[tinv[rots[[i]]],Conjugate[M],M].ConjugateTranspose[mats[[i]]],{i,ng}]];

runTest[l_]:=Module[{d=2(2l+1),A,M0,res},
  A=RandomComplex[{-1-I,1+I},{d,d}];M0=(A+ConjugateTranspose[A])/2;
  Print["---- l=",l," (dim ",d,") ----"];
  Do[Module[{Ms,Mss,idem,drift},
    Ms=N@symmAvg[l,sgn,M0];Mss=N@symmAvg[l,sgn,Ms];
    idem=Max@Abs@Flatten[Mss-Ms];
    drift=Max@Abs[Sort[Re@Eigenvalues[Ms]]-Sort[Re@Eigenvalues[Mss]]];
    Print["  ",If[sgn==1,"CURRENT exp(+i ph/2)","FLIP    exp(-i ph/2)"],
      ": idempotency=",ScientificForm[idem,3]," eig-shift=",ScientificForm[drift,3]];
   ],{sgn,{1,-1}}]];
runTest[1]; runTest[2];
