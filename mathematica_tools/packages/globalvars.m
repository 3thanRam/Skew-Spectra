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

BeginPackage["fftlogdiv`"]

k::usage="vector k";
k0::usage="Length of k";
muk::usage="cosine of k.z";

q::usage="vector q";
q0::usage="Length of q";
muq::usage="cosine of q.z";
phi::usage="last angle";

kmq::usage="Norm of k-q";

\[Epsilon]::usage="small parameter";
\[Nu]::usage="fftlog bias";
$Assumptions::usage="Assumptions";

Begin["`Private`"]
$Assumptions= k0>=0&&q0>=0&&kmq>=0&& 2 Pi>=phi>=0 && 1>=muk>=-1 && 1>=muq>=-1 && \[Epsilon] > 0&&{ k0,q0,muk, muq,phi, \[Epsilon],\[Nu]}\[Element]Reals;
k={0,k0 Sqrt[1-muk^2],k0 muk};
q={q0 Sqrt[1-muq^2] Cos[phi],q0 Sqrt[1-muq^2] Sin[phi],q0 muq};
End[] (*End Private Context*)

EndPackage[]



