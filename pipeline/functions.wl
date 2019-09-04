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

(* import bead set and landmark genes *)
beadset=AssociationThread[{"dp52","dp53"}->Import["resources/beadset.txt","Table"]];
lmgenes=Import["resources/lmgenes.txt","Table"][[1]];
lxbtonarray[lxb_]:=NumericArray[Sort@Pick[lxb[[;;,{1,5}]],Unitize[lxb[[;;,1]]lxb[[;;,5]]],1],"UnsignedInteger32"];

(* functions for LISS *)
invariantset[lxb_List]:=Table[Log2@N@ToPackedArray[Pick[#2,#1,i]/.{0->Nothing}],{i,10}]&[lxb[[;;,1]],lxb[[;;,2]]];
invariantset[lxb_NumericArray]:=invariantset[Normal[lxb]];
bootstrap[func_,data_,n_]:=Around[func[data],Max[StandardDeviation[func/@RandomChoice[data,{n,Length@data}]],0.01]];
lissexp=Log2@ToPackedArray@N@{50,100,200,400,600,1000,2000,3000,5000,7000};

(* convert between bead set order and gene order *)
gexidxtolxbidx=Flatten[Transpose@Table[Position[lmgenes,#]/.{}->{{979}}&/@beadset[dp],{dp,{"dp52","dp53"}}]];
lxbidxtogexidx=ToPackedArray[Position[gexidxtolxbidx,#][[1,1]]&/@Range[978]];
gextolog2lxb[gex_]:=Partition[ArrayPad[Log2@gex,{0,1},-666.][[gexidxtolxbidx]],2];
gextolxb[gex_]:=Partition[ArrayPad[gex,{0,1},-666.][[gexidxtolxbidx]],2];

(* define peak location grid and background extraction functions *)
$peakgrid=ToPackedArray@N@Range[-4,20,$resolution=1/16];
$backgroundbins=Reverse@NestWhileList[#-0.7(0.15+5.3*2^(-3#/4))&,15.09,#>0&];
intbg[lxb_]:=Interpolation[Transpose@{
    ArrayPad[$backgroundbins,{0,1},100.],
    ArrayPad[0.9(#/Total[#]&@GaussianFilter[BinCounts[Log2@N@lxb[[;;,2]],{$backgroundbins}],2])/(
        Partition[$backgroundbins,2,1].{-1,1}),{1,1}]+(0.1/15.)
    },InterpolationOrder->0];

(* load cuda dpeak library *)
wlldpeak=LibraryFunctionLoad["./dpeak.so","wll_dpeak",
    {{Real,2,"Constant"},{Real,1,"Constant"},{Real,1,"Constant"},{Real,2,"Constant"}, {Real,1,"Constant"}},{Real,3}];

(* functions for likelihood calculation *)
dpeakstrict[lxb_NumericArray,other___]:=dpeakstrict[Normal[lxb],other];
dpeakstrict[lxb_List/;ArrayDepth[lxb]==2,liss:{b_Real,a_Real,_},dp52ratio_:2./3,bgratio_:0.01]:=
    Block[{g=GroupBy[lxb,First->Last],w,val,padval,nval,bgprob,l,ltotal,bgn,bgr,bgw,bg=intbg[lxb]},
        val=Log2@N@ToPackedArray[g[#]/.{0->Nothing}/.{Missing[___]->{}}]&/@Range[11,500];
        w=Max[nval=ToPackedArray[Length/@val]];
        padval=ToPackedArray[PadRight[#,w,-1.]&/@val];
        bgprob=ToPackedArray@*bg/@val;
        bgprob=ToPackedArray[PadRight[#,w,-1.]&/@bgprob];
        ltotal=ConstantArray[0.,{Length@val,Length@$peakgrid,Length@$peakgrid}];
        bgn=ToPackedArray[Max[Round[bgratio Length@#],1]+{-1,0,1}&/@val];
        bgr=Clip[ToPackedArray@N[bgn/Clip[Length/@val,{1,\[Infinity]}]],{0.,1.}];
        bgw=ToPackedArray@N@MapThread[PDF[BinomialDistribution[#2,bgratio],#1]&,{bgn,Length/@val},1];
        Do[
            l=wlldpeak[padval,ReplacePart[ConstantArray[dp52ratio,490],{1->1.,489->0.}],bgr[[;;,ib]],bgprob,InverseFunction[b+a #&][$peakgrid]];
            Do[ltotal[[il]]+=bgw[[il,ib]]Exp[Clip[l[[il]]+200.,{-650.,\[Infinity]}]];,{il,Length@val}];
        ,{ib,Length[bgn[[1]]]}];
    ltotal=Log[ltotal];ltotal[[489]]+=ltotal[[1]];Clip[ltotal[[2;;]],{-650.,\[Infinity]}]
];

(* functions for marginal distribution & z-score calculation *)
$marginalprioridx=ToPackedArray[(lxbidxtogexidx/.{21->997})-22];
marginalprior[array_]:=ArrayReshape[#,{2Length[#],Length@$peakgrid}][[$marginalprioridx]]&@
ToPackedArray[(Log@Join[Total[#,{2}],Total[#]]&@Exp@Clip[#-Max[#],{-300.,0.}])&/@array];
marginalquantile[marginals_]:=Block[{probs=ToPackedArray[(#/Total[#]&@Exp[#-Max@#])&/@marginals]},
    probs.MovingAverage[#/Last@#&@Prepend[Accumulate@Total@probs,0.0],2]];
softmax[p_List]:=#/Total[#]&@Exp[p-Max[p]];
$resamplenumber=4Length[$peakgrid]-3;
resamplebase[p_]:=ArrayResample[p,$resamplenumber,Resampling->"Quadratic"]
$resamplematrix=SparseArray[ToPackedArray@Table[N@Chop[resamplebase[SparseArray[{i->1.},Length@$peakgrid]],0.0001],{i,Length@$peakgrid}]];
marginalqn[marginals_]:=
    Table[Block[{rm=ToPackedArray[softmax/@(Normal@marginals[[i]].$resamplematrix)],intf},
        intf=Interpolation[
            Transpose@{Prepend[#/Last[#]&@Accumulate@Clip[GaussianFilter[Mean@rm,30],{1.*^-10,\[Infinity]}],0.],Range[0.,N@$resamplenumber]},
            InterpolationOrder->1];
        rm.(SparseArray[Catenate@MapIndexed[
                {#2,#1}->#3&@@@ArrayPad[#1,{{0,0},{1,0}},#2[[1]]]&,
                BlockMap[(If[#1==#3,{{#1,#4-#2}},Join[{{#1,1.-#2}},{#,1}&/@Range[#1+1,#3-1],{{#3,#4}}]]&@@Flatten@QuotientRemainder[#,1])&,
                    Clip[intf[Range[0,1,0.001]],{0.,$resamplenumber-1.*^-10}]+1.,2,1]
            ],{$resamplenumber,1000},0.])]
    ,{i,Length@marginals}];
marginalqnquantile[marginals_]:=Block[{probs=ToPackedArray[marginals]},probs.MovingAverage[#/Last@#&@Prepend[Accumulate@Total@probs,0.0],2]];
$intinverf=ListInterpolation[InverseErf[Range[-0.9999,0.9999,0.0001]],{-0.9999,0.9999}];
