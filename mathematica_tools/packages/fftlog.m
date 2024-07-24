(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



(* ::Input::Initialization:: *)

BeginPackage["fftlog`"]


q0kmqmuqdecomp::usage="decompose expression into sum of powers of q0, kmq and muq ";
n1::usage="fftlog index from first power spectrum";
n2::usage="fftlog index from second power spectrum";

Begin["`Private`"]
Get[NotebookDirectory[]<>"globalvars.m"] ;
Get[NotebookDirectory[]<>"util.m"];
Get[NotebookDirectory[]<>"losintegrals.m"];
Get[NotebookDirectory[]<>"losintcalcs.m"];

(*decompose an expression into powers of k0,q0,kmq,muk,muq 
then get result of integration using master and line of sight integrals using the losint function from losintcalcs.m *)
q0kmqmuqdecomp[v0_]:=Module[{v=v0},
v=Expand[Simplfunc[v]];
(*v=If[Head[v]===Plus,List@@v,{v}];*)
Tv=Table[{ q0 D[v[[l]],q0]/v[[l]],kmq D[v[[l]],kmq]/v[[l]],muq D[v[[l]],muq]/v[[l]],v[[l]]/. {q0->1,kmq->1,muq->1}},{l,1,Length[v]}];

grouped=GroupBy[Tv,#[[{1,2,3}]]&->(#[[4]]&)];
groupedSums=KeyValueMap[#1->Total[#2]&,grouped];
Tuniq=Flatten/@List@@@Normal@groupedSums;

res=Sum[ Tuniq[[m,4]] losint[-(Tuniq[[m,1]]-Tuniq[[m,3]])/2+n1, -Tuniq[[m,2]]/2+n2,Tuniq[[m,3]]],{m,1,Length[Tuniq]}];

Collect[res,muk]
]

End[] 

EndPackage[]



