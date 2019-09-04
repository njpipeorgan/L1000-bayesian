(*
 * L1000 peak deconvolution based on Bayesian analysis
 * 
 * Copyright 2019 Yue Qiu & Tianhuan Lu
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *)

<<Developer`;<<functions.wl;
If[Head[wlldpeak]=!=LibraryFunction,Print["Failed to load library function from \"dpeak.so\""];Quit[];];
If[Length@$ScriptCommandLine==2,$folder=$ScriptCommandLine[[2]],
    Print["usage: wolframscript -file pipeline.wl <folder name>"];Quit[];];
$lxbs=FileNames[$folder<>"/*.lxb"];
Print[ToString@Length@$lxbs<>" LXB files found."];If[Length@$lxbs==0,Quit[];];
$plateid=StringTake[Last@FileNameSplit[$lxbs[[1]]],{1,-9}];
Quiet@CreateDirectory["output"];

(* read lxb files *)
lxbs=AssociationThread[(StringTake[Last@FileNameSplit[#],{1,-5}]&/@$lxbs)->(lxbtonarray@Import[#,{"FCS","Data"}]&/@$lxbs)];

(* LISS, calculate scaling parameters *)
realexps=Block[{i=invariantset[#]},bootstrap[Median,#,100]&/@i]&/@lxbs;
mappedexps=Table[(r-#1)/#2&@@LinearModelFit[Transpose@{lissexp,r[[;;,1]]},{1,x},x,Weights->(r[[;;,2]]^-2)]["BestFitParameters"],{r,realexps}];
caliblissexp=(Around[#1.#2/Total[#2],Sqrt[(#1^2).#2/Total[#2]-(#1.#2/Total[#2])^2+1/Mean[#2]]]&[#[[;;,1]],#[[;;,2]]^-2]&/@Transpose[mappedexps]);
calibparams=AssociationThread[Keys[realexps]->Table[
    Append[#1,Mean[(#2/caliblissexp[[;;,2]])^2]]&@@Quiet@LinearModelFit[
        Transpose[{r[[;;,1]],caliblissexp[[;;,1]]}],{1,x},x,Weights->(caliblissexp[[;;,2]]^-2)][{"BestFitParameters","FitResiduals"}]
    ,{r,realexps}]];

(* calculate marginal distributions and z-scores *)
LaunchKernels[4];
$alphac=0.01;
$dp52ratio=2./3;
marginals=ParallelTable[marginalprior@dpeakstrict[lxb[[2]],calibparams@lxb[[1]],$dp52ratio,$alphac],{lxb,Normal[lxbs]},Method->"FinestGrained"];
marginalsqn=marginalqn@marginals;
marginalsquantile=Transpose@ToPackedArray@Table[marginalqnquantile[Normal/@marginalsqn[[;;,g]]],{g,978}];
mqnzscore=AssociationThread[Keys[lxbs]->(N@Sqrt[2]$intinverf[Clip[2.marginalsquantile-1.,{-0.9999,0.9999}]])];

(* output results *)
Export["output/liss_parameters-"<>$plateid<>".h5",{"/colid"->Keys[lxbs],"/field"->{"intercept","slope","chi_squared"},"/data"->calibparams}];
Export["output/marginals-"<>$plateid<>".h5",{"/colid"->Keys[lxbs],"/rowid"->lmgenes,"/peakloc"->NumericArray[$peakgrid,"Real32"],"/data"->NumericArray[marginals,"Real32"]}];
Export["output/zscore_pc-"<>$plateid<>".h5",{"/colid"->Keys[lxbs],"/rowid"->lmgenes,"/data"->NumericArray[Values@mqnzscore,"Real32"]}];
