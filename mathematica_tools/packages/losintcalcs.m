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
BeginPackage["losintcalcs`"]
losint::usage="integral over q of (q.z)^n3 q^-n1 kmq^-n2";

Begin["`Private`"]
Get[NotebookDirectory[]<>"globalvars.m"] ;
Get[NotebookDirectory[]<>"losintegrals.m"] ;

delta[x_,y_]/;(Order[x,y]<0):=delta[y,x];
delta[x_,x_]:=1;

(*get line of sight function by letter index and order*)
getFunction[letter_,n_]:=Module[{functionName},functionName=ToExpression[ToString[letter]<>ToString[n]]]
numberToLetter[n_]:=FromCharacterCode[64+n]

(*Function to sum an expression over all permutations of given variables*)
sumOverPermutations[expr_,vars_List]:=Module[{perms,substitutedExprs},perms=Permutations[vars];
sortedPerms=Sort/@perms;substitutedExprs=(expr/. Thread[vars->#]&)/@perms;Total[substitutedExprs]]

divideByNumberOfTerms[expr_]:=Module[{terms,numTerms,dividedExpr},
terms=If[Head[expr]===Plus,List@@expr,{expr}];
numTerms=Length[terms];
dividedExpr=expr/numTerms]


(*returns an expression with a given number of delta functions and products of k with different indices*)
Deltak[n0_,nD0_,nk0_]:=Module[{n=n0,nD=nD0,nk=nk0},
expr=1;
expr*=Product[delta[numberToLetter[2 id+1],numberToLetter[2 id+2]],{id,0,nD-1}];
expr*=Product[K[numberToLetter[2 nD+1+ik]],{ik,0,nk-1}];
	
If[n>2 && nD>0,perms=Simplify[sumOverPermutations[expr,Table[numberToLetter[i],{i,n}]]];expr=divideByNumberOfTerms[perms/perms[[1]]],expr=expr];
expr
]

(*get the expression of the line of sight integral of a given order*)
losint[n10_,n20_,n30_]:=Module[{n1=n10,n2=n20,n3=n30},
Vars=Table[getFunction[numberToLetter[i0+1],n3] ,{i0,0,Quotient[n3,2]}];
res=k0^(3-2(n1+n2));
S=Sum[k0^(-2(i0-Quotient[n3,2])) Deltak[n,-(i0-Quotient[n3,2]),n3+2(i0-Quotient[n3,2])]getFunction[numberToLetter[i0+1],n3] ,{i0,0,Quotient[n3,2]}];
S=S/.{K[x_]:>k[[3]],delta[x_,y_]:>1,getFunction[numberToLetter[1],0]->I2[n1,n2]};
For[v=1,v<Length[Vars]+1,v++,
S=S/.{Vars[[v]]->Vars[[v]][n1,n2]};
];
res*=S //Simplify
]

End[] (*End Private Context*)

EndPackage[]



