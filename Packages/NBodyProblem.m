(* ::Package:: *)

(* ::Section:: *)
(*N-Body Problem *)


BeginPackage["NBodyProblem`"];


Chdata::usage=".......";
NBodyODE::usage=" 

    \!\(\*SubscriptBox[\(q\), \(i\)]\)= \!\(\*SubscriptBox[\(v\), \(i\)]\); 
    \!\(\*SubscriptBox[\(v\), \(i\)]\)= \!\(\*UnderoverscriptBox[\(\[Sum]\), \(j = 1  \\\\n  j \[NotEqual] i\), \(N\)]\)\!\(\*FractionBox[\(-Gmj\), SuperscriptBox[\(\[LeftDoubleBracketingBar]\*SubscriptBox[\(q\), \(i\)] - \*SubscriptBox[\(q\), \(j\)]\[RightDoubleBracketingBar]\), \(3\)]]\)(\!\(\*SubscriptBox[\(q\), \(i\)]\)-\!\(\*SubscriptBox[\(q\), \(j\)]\));

    parameters[[1]]=1 ----> Return (fq)
    parameters[[1]]=2 ----> Return (fp)
    otherwise        ----> Return (f=fq,fp)";
NBodyODEc::usage=" ";
NBodyODEd::usage=" ";
NBodyODEfunc::usage=" ";
NBodyHam::usage=" 
    H(q,v)=H-P= \!\(\*FractionBox[\(1\), \(2\)]\) \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(N\)]\)\!\(\*SubscriptBox[\(Gm\), \(i\)]\)\[LeftDoubleBracketingBar]\!\(\*SubscriptBox[\(v\), \(i\)]\)\!\(\*SuperscriptBox[\(\[RightDoubleBracketingBar]\), \(2\)]\) -\!\(\*UnderoverscriptBox[\(\[Sum]\), \(1 \[LessEqual] i < j \[LessEqual] N\), \(N\)]\)\!\(\*FractionBox[\(\*SubscriptBox[\(Gm\), \(i\)] \*SubscriptBox[\(Gm\), \(j\)]\), \(\(\\\\\)\(\[LeftDoubleBracketingBar]\*SubscriptBox[\(q\), \(i\)] - \*SubscriptBox[\(q\), \(j\)] \[LeftDoubleBracketingBar]\)\)]\) ";
NBodyHam2::usage="";

Begin["`Private`"];


(* ::Subsection:: *)
(*Functions*)


(* ::Subsubsection::Closed:: *)
(*Chdata*)


Chdata[u_,parameters_List]:=
Module[{dim,neq,d,nbody,Gm,MM,mq,mv,i,i1,i2,id,qi,vi,qnew,vnew},

dim=3;
neq=Length[u];
d=neq/2;
nbody=Length[u]/(2*dim);
Gm=parameters;
MM=Sum[Gm[[i]],{i,nbody}];
 mq=Array[0&,dim];
 mv=Array[0 &,dim];

 Do[ i1=(i-1)*dim;
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

Do[ i1=(i-1)*dim;
  i2=d+i1;
  AppendTo[qnew,Table[u[[i1+id]]-mq[[id]],{id,1,dim}]];
  AppendTo[vnew,Table[u[[i2+id]]-mv[[id]],{id,1,dim}]],
{i,nbody}];

Flatten[Join[qnew,vnew]]
];



(* ::Input:: *)
(**)


(* ::Subsubsection::Closed:: *)
(*Odefun*)


NBodyODE[t_,u_,parameters_List]:=
Module[{evalODE,dim,neq,d,nbody,y,r3,Gm,f,f2,i,i1,i2,j,j1,j2,id,
             qi,qj,qij,vi},
dim=3;
neq=Length[u];
d=neq/2;
nbody=Length[u]/(2*dim);
f=Array[0 &,neq];
f2=Array[0 &,d];
Gm =Drop[parameters,-1];

If[parameters[[-1]]==1, 
(*fq*)
Do[ i1=(i-1)*dim;
      i2=d+i1;
      vi=Table[u[[i2+id]],{id,1,dim}];
      Do[f2[[i1+id]]=vi[[id]],{id,1,dim}],     
{i,nbody}];
f2,

If[parameters[[-1]]==2,
(*fp*)
Do[ i1=(i-1)*dim;
      qi=Table[u[[i1+id]],{id,1,dim}];
      Do [    j1=(j-1)*dim;
                 qj=Table[u[[j1+id]],{id,1,dim}];
                 qij=qi-qj;
                 r3=Norm[qij]^3;
                 Do[ f2[[i1+id]]=f2[[i1+id]]-Gm[[j]]*qij[[id]]/r3,{id,1,dim}];
                 Do[ f2[[j1+id]]=f2[[j1+id]]+Gm[[i]]*qij[[id]]/r3,{id,1,dim}],
 {j,i+1,nbody}],
{i,nbody}];
f2,

(*else fq,fp*)
Do[ i1=(i-1)*dim;
      i2=d+i1;
      qi=Table[u[[i1+id]],{id,1,dim}];
      vi=Table[u[[i2+id]],{id,1,dim}];
      Table[f[[i1+id]]=vi[[id]],{id,1,dim}];
      Do [    j1=(j-1)*dim;
                j2=d+j1;
                 qj=Table[u[[j1+id]],{id,1,dim}];
                 qij=qi-qj;
                 r3=Norm[qij]^3;
                 Do[ f[[i2+id]]=f[[i2+id]]-Gm[[j]]*qij[[id]]/r3,{id,1,dim}];
                 Do[ f[[j2+id]]=f[[j2+id]]+Gm[[i]]*qij[[id]]/r3,{id,1,dim}],
 {j,i+1,nbody}],
{i,nbody}];
f
]]

];

NBodyODEc[t_,u_,parameters_List]:=Apply[NBodyODEfunc,Flatten[{t,u,parameters}]];

NBodyODEd[t_,u_,parameters_List]:=SetPrecision[Apply[NBodyODEfunc,Flatten[N[{t,u,parameters}]]],40];

NBodyODEfunc=Compile[{t,
q11,q12,q13,q21,q22,q23,q31,q32,q33,q41,q42,q43,q51,q52,q53,q61,q62,q63,q71,q72,q73,q81,q82,q83,q91,q92,q93,q101,q102,q103,
p11,p12,p13,p21,p22,p23,p31,p32,p33,p41,p42,p43,p51,p52,p53,p61,p62,p63,p71,p72,p73,p81,p82,p83,p91,p92,p93,p101,p102,p103,
Gma1,Gma2,Gma3,Gma4,Gma5,Gma6,Gma7,Gma8,Gma9,Gma10,par1},

Module[{dim,neq,d,nbody,y,r3,Gm,f,f2,i,i1,i2,j,j1,j2,id,
             qi,qj,qij,vi,
             u},

u={q11,q12,q13,q21,q22,q23,q31,q32,q33,q41,q42,q43,q51,q52,q53,q61,q62,q63,q71,q72,q73,q81,q82,q83,q91,q92,q93,q101,q102,q103,
p11,p12,p13,p21,p22,p23,p31,p32,p33,p41,p42,p43,p51,p52,p53,p61,p62,p63,p71,p72,p73,p81,p82,p83,p91,p92,p93,p101,p102,p103};

dim=3;
neq=60;
d=30;(*neq/2       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*)
nbody=neq/(2*dim);
f=Array[N[0] &,neq];
f2=Array[N[0] &,d];
Gm ={Gma1,Gma2,Gma3,Gma4,Gma5,Gma6,Gma7,Gma8,Gma9,Gma10};

If[par1==1, 
(*fq*)
Do[ i1=(i-1)*dim;
      i2=d+i1;
      vi=Table[u[[i2+id]],{id,1,dim}];
      Do[f2[[i1+id]]=vi[[id]],{id,1,dim}],     
{i,nbody}];
f2,

If[par1==2,
(*fp*)
Do[ i1=(i-1)*dim;
      qi=Table[u[[i1+id]],{id,1,dim}];
      Do [    j1=(j-1)*dim;
                 qj=Table[u[[j1+id]],{id,1,dim}];
                 qij=qi-qj;
                 r3=Norm[qij]^3;
                 Do[ f2[[i1+id]]=f2[[i1+id]]-Gm[[j]]*qij[[id]]/r3,{id,1,dim}];
                 Do[ f2[[j1+id]]=f2[[j1+id]]+Gm[[i]]*qij[[id]]/r3,{id,1,dim}],
 {j,i+1,nbody}],
{i,nbody}];
f2,

(*else fq,fp*)
Do[ i1=(i-1)*dim;
      i2=d+i1;
      qi=Table[u[[i1+id]],{id,1,dim}];
      vi=Table[u[[i2+id]],{id,1,dim}];
      Do[f[[i1+id]]=vi[[id]],{id,1,dim}];
      Do [    j1=(j-1)*dim;
                j2=d+j1;
                 qj=Table[u[[j1+id]],{id,1,dim}];
                 qij=qi-qj;
                 r3=Norm[qij]^3;
                 Do[ f[[i2+id]]=f[[i2+id]]-Gm[[j]]*qij[[id]]/r3,{id,1,dim}];
                 Do[ f[[j2+id]]=f[[j2+id]]+Gm[[i]]*qij[[id]]/r3,{id,1,dim}],
 {j,i+1,nbody}],
{i,nbody}];
f
]]
]
];



(* ::Subsubsection::Closed:: *)
(*Hamiltonian*)


NBodyHam[{t_,u_,e_},parameters_List, doi_:100]:=
Module[{dim,neq,d,nbody,Gm,uu,H,P,r,i,i1,i2,j,j1,j2,id,
              qi,qj,qij,vi},
dim=3;
neq=Length[u];
d=neq/2;
nbody=Length[u]/(2*dim);
Gm =SetPrecision[parameters,doi];
uu=SetPrecision[u,doi]+SetPrecision[e,doi];
H=0;
P=0;
Do[ i1=(i-1)*dim;
      i2=d+i1;
      qi=Table[uu[[i1+id]],{id,1,dim}];
      vi=Table[uu[[i2+id]],{id,1,dim}];
      H=H+Gm[[i]]*(vi.vi);
      Do [    j1=(j-1)*dim;
                j2=d+j1;
                 qj=Table[uu[[j1+id]],{id,1,dim}];
                 qij=qi-qj;
                 r=Norm[qij];
                 P=P+Gm[[i]]*Gm[[j]]/r,
 {j,i+1,nbody}],
{i,nbody}];
H/2-P
];



(* ::Subsubsection::Closed:: *)
(*Hamiltonian C*)


NBodyHam2[out_,parameters_List, doi_:100]:=
Module[{dim,neq,d,nbody,Gm,uu,H,P,r,i,i1,i2,j,j1,j2,id,
              qi,qj,qij,vi},
dim=3;
uu=SetPrecision[out[[2;;61]],doi]+SetPrecision[out[[62;;121]],doi];
neq=Length[uu];
d=neq/2;
nbody=Length[uu]/(2*dim);
Gm =SetPrecision[parameters,doi];
H=0;
P=0;
Do[ i1=(i-1)*dim;
      i2=d+i1;
      qi=Table[uu[[i1+id]],{id,1,dim}];
      vi=Table[uu[[i2+id]],{id,1,dim}];
      H=H+Gm[[i]]*(vi.vi);
      Do [    j1=(j-1)*dim;
                j2=d+j1;
                 qj=Table[uu[[j1+id]],{id,1,dim}];
                 qij=qi-qj;
                 r=Norm[qij];
                 P=P+Gm[[i]]*Gm[[j]]/r,
 {j,i+1,nbody}],
{i,nbody}];
H/2-P
];


(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];



