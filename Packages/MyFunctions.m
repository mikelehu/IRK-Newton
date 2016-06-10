(* ::Package:: *)

(* ::Section:: *)
(*MyFunctions*)


BeginPackage["MyFunctons`"];

(*******MyFunctions-1**********)

ErrorPos2::usage="";
Chdata::usage=".......";

(*******MyFunctions-2**********)

Pert::usage="";
FunMeanHerr::usage="";
FunAllHerr::usage="";
FunHerr::usage="";
FunMeanEst::usage="";

(*******MyFunctions-3**********)

FunMeanH::usage="";


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*MyFunctions-1*)


ErrorPos2[v_,vx_ ,doi_:100]:=
Module[{q,qx,neq,d},
neq=(Length[v]-1)/2;
d=neq/2;
q=SetPrecision[v[[2;;d+1]],doi]+SetPrecision[v[[2d+2;;3d+1]]
,doi];
qx=SetPrecision[vx[[2;;d+1]],doi]+SetPrecision[vx[[2d+2;;3d+1]]
,doi];
Norm[q-qx]];




Chdata[u_,parameters_List]:=
Module[{dim,neq,d,nbody,Gm,MM,mq,mv,i,i1,i2,id,qi,vi,qnew,vnew},

dim=3;
neq=Length[u];
d=neq/2;
nbody=Length[u]/(2*dim);
Gm=parameters;
MM=Sum[Gm[[i]],{i,nbody}];
mq=Array[0&,dim];
mv=Array[0&,dim];

 Do[ 
  i1=(i-1)*dim;
  i2=d+i1;
  qi=Table[u[[i1+id]],{id,1,dim}];
  vi=Table[u[[i2+id]],{id,1,dim}];
  mq=mq+Gm[[i]]*qi;
  mv=mv+Gm[[i]]*vi,
{i,nbody}];

 mq=mq/MM;
 mv=mv/MM;

qnew={};
vnew={};

Do[ 
  i1=(i-1)*dim;
  i2=d+i1;
  AppendTo[qnew,Table[u[[i1+id]]-mq[[id]],{id,1,dim}]];
  AppendTo[vnew,Table[u[[i2+id]]-mv[[id]],{id,1,dim}]],
{i,nbody}];

Flatten[Join[qnew,vnew]]
];



(* ::Subsection::Closed:: *)
(*MyFunctions-2*)


Pert[e0_,k_]:={e0*(k*RandomReal[{-1,1}])};


FunMeanHerr[outAll_, outBAll_, nstat_,nout_,HAM_,parameters_
             ,neq_,prec_,step_Integer]:=
Module[{i,Ham0, Hamerr,Error,
             SHamerr=Array[0 &,Floor[nout/step]],
             DHamerr=Array[0 &,Floor[nout/step]],
             SError=Array[0 &,Floor[nout/step]],
             DError=Array[0 &,Floor[nout/step]],
             EndHam=Array[0&,nstat],outA,outB,MeanE,DesvE},
For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  outB=outBAll[[i]];
  Ham0 = HAM[neq,First[outA],parameters,prec];
  Hamerr = Map[HAM[neq,#,parameters,prec]/Ham0-1&,outA];
  Hamerr=  Map[First,Partition[Hamerr,step]];
  EndHam[[i]]=Last[Hamerr];
(*  SHamerr=SHamerr+Abs[Hamerr];*)
  SHamerr=SHamerr+Hamerr;
  DHamerr=DHamerr+Hamerr^2;
  Error =Chop[ MapThread[ErrorPos2[#1,#2,prec]&,{outB, outA}],10^-100];
  Error= Map[First,Partition[Error,step]];
  SError=SError+Abs[Error];
];
MeanE=SHamerr/nstat;
DesvE=Sqrt[DHamerr/nstat-MeanE^2];
{MeanE,DesvE,SError/nstat ,EndHam}
]


FunAllHerr[outAll_, nstat_,nout_,HAM_,parameters_
             ,neq_,prec_,step_Integer]:=
Module[{i,Ham0, Hamerr,      
             AllHamerr={}
},
For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  Ham0 = HAM[neq,First[outA],parameters,prec];
  Hamerr = Map[HAM[neq,#,parameters,prec]/Ham0-1&,outA];
  Hamerr=  Map[First,Partition[Hamerr,step]];
  AppendTo[AllHamerr,Hamerr];
];
AllHamerr
]


FunHerr[outAll_, nstat_,nout_,HAM_,parameters_,
        neq_,prec_,step_Integer]:=
Module[{i,Ham0, Hamerr,Hamerrdif,Error,outA},
Hamerrdif={};
For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  Ham0 = HAM[neq,First[outA],parameters,prec];
  Hamerr=Map[HAM[neq,#,parameters,prec]/Ham0-1&,outA];
  Hamerr=Map[First,Partition[Hamerr,step]];
  Hamerrdif = Join[Hamerrdif,Drop[Hamerr,1]-Drop[Hamerr,-1]];
];
Hamerrdif
]



FunMeanEst[outAll_,outBAll_,nstat_,nout_,neq_,prec_,step_Integer]:=
Module[{i,outA,outB,
             est,estQ,meanEst,desvEst,Error,
             Qty,meanQty,desvQty,
             Sout=Array[0 &,Floor[nout/step]],
             S2out=Array[0 &,Floor[nout/step]],
             SEstQty=Array[0 &,Floor[nout/step]],
             S2EstQty=Array[0 &,Floor[nout/step]]},

For[i=1, i<= nstat,i++,
outA=outAll[[i]];
outB=outBAll[[i]];
est=Take[outA, All,{2neq+2,3neq+1}];
estQ=Map[Norm,Chop[Take[est,All,{1,neq/2}],10^-100]];
estQ=Map[First,Partition[estQ,step]];
Sout=Sout+estQ;
S2out=S2out+estQ^2;
Error =Chop[ MapThread[ErrorPos2[#1,#2,prec]&,{outB, outA}],10^-100];
Error= Map[First,Partition[Error,step]];
Qty=Log[10,Quiet[estQ/Error]/.ComplexInfinity->1/.Indeterminate->1//.0->1];
SEstQty=SEstQty+Qty;
S2EstQty=S2EstQty+Qty^2;
];

meanEst=Sout/nstat;
desvEst=Sqrt[S2out/nstat-meanEst^2];
meanQty=SEstQty/nstat;
desvQty=Sqrt[S2EstQty/nstat-meanQty^2];
{meanEst,desvEst,meanQty,desvQty}
];



(* ::Subsection:: *)
(*MyFunctions-3*)


FunMeanH[outAll_, nstat_,nout_,HAM_,parameters_
             ,neq_,prec_,step_Integer]:=
Module[{i,Ham0, Hamerr,Error,
             SHamerr=Array[0 &,Floor[nout/step]],
             DHamerr=Array[0 &,Floor[nout/step]],
             EndHam=Array[0&,nstat],outA,outB,MeanE,DesvE},
For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  Ham0 = HAM[neq,First[outA],parameters,prec];
  Hamerr = Map[HAM[neq,#,parameters,prec]/Ham0-1&,outA];
  Hamerr=  Map[First,Partition[Hamerr,step]];
  EndHam[[i]]=Last[Hamerr];
(*  SHamerr=SHamerr+Abs[Hamerr];*)
  SHamerr=SHamerr+Hamerr;
  DHamerr=DHamerr+Hamerr^2;
 
];
MeanE=SHamerr/nstat;
DesvE=Sqrt[DHamerr/nstat-MeanE^2];
{MeanE,DesvE,EndHam}
]


(* ::Subsection:: *)
(*End*)


End[ ];
EndPackage[];



