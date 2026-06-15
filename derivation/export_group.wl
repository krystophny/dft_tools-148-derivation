rots={IdentityMatrix[3],RotationMatrix[Pi/2,{0,0,1}],RotationMatrix[Pi,{0,0,1}],
      RotationMatrix[3Pi/2,{0,0,1}],RotationMatrix[Pi,{1,0,0}],RotationMatrix[Pi,{0,1,0}],
      RotationMatrix[Pi,{1,1,0}],RotationMatrix[Pi,{1,-1,0}]};
tinv[R_]:=(R[[1,1]]R[[2,2]]-R[[1,2]]R[[2,1]])<0;
Do[Module[{R=rots[[i]],e,a,b,g,ph,ti,dt},e=EulerAngles[R];{a,b,g}=e;
   ti=tinv[R];ph=If[ti,g-a,a+g];dt=Det[R];
   Print[NumberForm[N[a,12],12]," ",NumberForm[N[b,12],12]," ",NumberForm[N[g,12],12],
         " ",NumberForm[N[dt,3],3]," ",If[ti,1,0]," ",NumberForm[N[ph,12],12]]],{i,Length[rots]}];
