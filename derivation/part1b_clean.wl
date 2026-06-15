$Assumptions = (a|b|g) \[Element] Reals;
sy = PauliMatrix[2]; sz = PauliMatrix[3];
Dhalf[a_,b_,g_] := MatrixExp[-I a sz/2].MatrixExp[-I b sy/2].MatrixExp[-I g sz/2];
spmt[a_,b_,g_] := {{ Exp[ I (a+g)/2] Cos[b/2],  Exp[-I (a-g)/2] Sin[b/2] },
                   {-Exp[ I (a-g)/2] Sin[b/2],  Exp[-I (a+g)/2] Cos[b/2] }};
z[m_]:=FullSimplify[m]; isZero[m_]:=FullSimplify[m]==ConstantArray[0,Dimensions[m]];

(* setsym/spinrotmat special-case spin blocks *)
srmB0[a_,g_]  := Module[{e=Exp[I (a+g)/2]},  {{e,0},{0,Conjugate[e]}}];
srmBpi[a_,g_] := Module[{e=Exp[I (g-a)/2]}, {{0,e},{-Conjugate[e],0}}];
(* outputqmc non-mixing diagonal construction (current code) *)
outB0[a_,g_]  := Module[{e=Exp[I (a+g)/2]}, {{e,0},{0,Conjugate[e]}}];
outBpi[a_,g_] := Module[{e=Exp[I (g-a)/2]}, {{e,0},{0,Conjugate[e]}}];
(* issue's proposed sign flip: exp(-i phase/2) *)
flipB0[a_,g_] := Module[{e=Exp[-I (a+g)/2]}, {{e,0},{0,Conjugate[e]}}];

Print["spmt unitary & det1: ", FullSimplify[spmt[a,b,g].ConjugateTranspose[spmt[a,b,g]]]//MatrixForm,
      "  det=",FullSimplify[Det[spmt[a,b,g]]]];
Print["B: spinrotmat==spmt  beta=0?  ", isZero[srmB0[a,g]-spmt[a,0,g]],
      "   beta=pi? ", isZero[srmBpi[a,g]-spmt[a,Pi,g]]];
Print["C: outputqmc==spmt   beta=0?  ", isZero[outB0[a,g]-spmt[a,0,g]],
      "   beta=pi? ", isZero[outBpi[a,g]-spmt[a,Pi,g]]];
Print["   outputqmc beta=pi (diag) = ", z[outBpi[a,g]]//MatrixForm,
      "   spmt beta=pi (antidiag) = ", z[spmt[a,Pi,g]]//MatrixForm];
Print["D: issue-flip==spmt  beta=0?  ", isZero[flipB0[a,g]-spmt[a,0,g]],
      "   flip==conj(spmt) beta=0? ", isZero[flipB0[a,g]-Conjugate[spmt[a,0,g]]]];
