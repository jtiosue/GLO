(* ::Package:: *)

(* ::Title:: *)
(*Gaussian Linear Optics (GLO) package*)


BeginPackage["GLO`"];
(*Load with Get["path/GLO.wl"] to force reload every time*)
(*or Needs["path/GLO.wl"] to load only once.*)


(* ::Section::Closed:: *)
(*Overview*)


(* ::Text:: *)
(*In this package, we provide various functions relevant to entanglement bosonic Gaussian states. We provide functions for simulating linear optical circuits, as well as functions for computing quantities after applying Haar random linear optical unitaries. Etc.*)


(* ::Section:: *)
(*Public functions*)


(* ::Subsection::Closed:: *)
(*Matrix functions*)


matrixDirectSum::usage="matrixDirectSum[args] returns the direct sum " <>
"of the matrices in the array 'args'.";


symplecticForm::usage="symplecticForm[n] gives the 2n\[Cross]2n matrix representing " <>
"the symplectic form on n modes.";


symplecticEigenvalues::usage="symplecticEigenvalues[v] for a 2n\[Cross]2n covariance " <>
"matrix v gives the n positive real symplectic eigenvalues of v.";


initialSqueezedState::usage="initialSqueezedState[squeeze] returns the " <>
"covariance matrix associated with Length[squeeze] modes and squeezing " <>
"parameter r=squeeze[[i]] on mode i for a product state of squeezed modes. " <>
"This will be such that the standard deviation of x on that mode is " <>
"exp(-r) and that of p exp(r).";


validCovarianceMatrixQ::usage="validCovarianceMatrixQ[\[Sigma]] returns a Boolean; " <>
"whether \[Sigma] is a valid covariance matrix.";


reducedCovarianceMatrix::usage="reducedCovarianceMatrix[\[Sigma],subsystem] " <>
"gives the covariance matrix of the reduced state consisting of the modes " <>
"in `subsystem`.";
reducedCovarianceMatrixOfLeftPartition::usage="reducedCovarianceMatrixOfLeftPartition" <>
"[\[Sigma],partition] gives the covariance matrix of the reduced state consisting of the " <>
"mode 1 through `partition`. I.e. to the left of or at `partition`.";
reducedCovarianceMatrixOfRightPartition::usage="reducedCovarianceMatrixOfRightPartition" <>
"[\[Sigma],partition] gives the covariance matrix of the reduced state consisting of the " <>
"modes `partition` through th end. I.e. to the right of or at `partition`.";


squeezerMatrix::usage="squeezerMatrix[r] is the 2\[Cross]2 symplectic matrix " <>
"corresponding to one mode squeezing.";


beamsplitterMatrix::usage="beamsplitterMatrix[\[Theta]] is the 4\[Cross]4 symplectic " <>
"matrix corresponding to the beamsplitter.";


phaseshiftMatrix::usage="phaseshiftMatrix[\[Theta]] is the 2\[Cross]2 symplectic " <>
"matrix corresponding to the phase shifter.";


twoModeMatrix::usage="twoModeMatrix[\[Theta],\[Phi],\[Beta]] is an arbitrary two mode " <>
"passive transformation. It is the most general 4\[Cross]4 orthogonal symplectic " <>
"matrix. It corresponds to a \[Theta],\[Phi] phase shift on mode 1,2, and a \[Beta] " <>
"beamsplitter applied between them.";


oneDimRandomBrickworkLayer::usage="oneDimRandomBrickworkLayer[n] returns " <>
"a 2n\[Cross]2n matrix corresponding to one layer of a one dimensional brickwork " <>
"on n modes, where all the angles are chosen uniformly randomly. " <>
"Note that n MUST BE EVEN.";


oneDimRandomActiveBrickworkLayer::usage="oneDimRandomActiveBrickworkLayer[n,rmin,rmax] " <>
"returns a 2n\[Cross]2n matrix corresponding to one layer of a one dimensional brickwork " <>
"on n modes, where all the angles are chosen uniformly randomly and the squeezer " <>
"is chosen uniformly in [rmin,rmax]. By default, rmin is chosen to be -1 and " <>
"rmax is chosen to be 1. Note that n MUST BE EVEN.";


oneDimFixedBrickworkLayer::usage="oneDimFixedBrickworkLayer[n,\[Theta],\[Phi]] returns " <>
"a 2n\[Cross]2n matrix corresponding to one layer of a one dimensional brickwork " <>
"on n modes, where all the beamsplitter angles are set to \[Theta] and all the " <>
"phase shift angles are set to \[Phi]. Note that n MUST BE EVEN.";


oneDimFixedBeamsplitterBrickworkLayer::usage=
"oneDimFixedBeamsplitterBrickworkLayer[n,\[Theta]] returns " <>
"a 2n\[Cross]2n matrix corresponding to one layer of a one dimensional brickwork " <>
"on n modes, where all the beamsplitter angles are set to \[Theta] and all the " <>
"phase shift angles are chosen uniformly randomly. Note that n MUST BE EVEN.";


applyByCongruence::usage="applyByCongruence[v,s] applys the orthogonal symplectic " <>
"transformation on the covariance matrix v by congruence, s.v.s\[Transpose]."


haarModeMatrix::usage="haarModeMatrix[n] gives a Haar random 2n\[Cross]2n orthogonal " <>
"symplectic matrix (via its isomorphism with the unitary group). This " <>
"corresponds to a random passive Gaussian unitary.";


haarTwoModeMatrix::usage="haarTwoModeMatrix[n,i,j] gives a 2n\[Cross]2n orthogonal " <>
"symplectic matrix. It is identity on all modes other than i and j, and applies " <>
"a Haar random 2 mode unitary on modes i and j.";


measureVacuumMode::usage="measureVacuumMode[\[Sigma],mode] measures mode `mode` and " <>
"projects the mode into the vacuum state |0>. See Serafini page 129 with " <>
"\[Sigma]m = IdentityMatrix[2].";


measureVacuumModes::usage="measureVacuumsMode[\[Sigma],modes] measures each mode in the " <>
"table `modes` and projects the mode into the vacuum state |0>. See " <>
"Serafini page 129 with \[Sigma]m = IdentityMatrix[2].";


(* ::Subsection::Closed:: *)
(*Quantities from the covariance matrix*)


renyi2EntropyFromPurity::usage="renyi2EntropyFromPurity[purity] returns the " <>
"Renyi-2 entropy from the purity.";


purity::usage="purity[v] computes the purity of the state described " <>
"by the covariance matrix v.";
renyi2Entropy::usage="renyi2Entropy[v] computes the Renyi-2 entropy of the" <>
"state described by the covariance matrix v.";
entanglementEntropy::usage="entanglementEntropy[v] computes the entanglement " <>
"entropy of the state described by the covariance matrix v.";
energy::usage="energy[v] compute the energy of the state described by " <>
"the covariance matrix v. In particular, energy[v] = Tr[v].";
averagePhotonNumber::usage="averagePhotonNumber[v] computes the expected value " <>
"of the photon number operator in the Gaussian state described by the covariance " <>
"matrix v.";


purityOfSubsystem::usage="purityOfSubsystem[v,subsystem] computes the purity " <>
"of the reduced state of the state described by the covariance matrix v, " <>
"where the reduced state is of the modes in subsystem";
renyi2EntropyOfSubsystem::usage="renyi2EntropyOfSubsystem[v,subsystem] " <>
"computes the Renyi-2 entropy of the reduced state of the state " <>
"described by the covariance matrix v, where the reduced state is of " <>
"the modes in subsystem";
entanglementEntropyOfSubsystem::usage="entanglementEntropyOfSubsystem[v,subsystem] " <>
"computes the entanglement entropy of the reduced state of the state " <>
"described by the covariance matrix v, where the reduced state is of " <>
"the modes in subsystem";
energyOfSubsystem::usage="energy[v,subsystem] computes the energy of the reduced " <>
"state described by the covariance matrix v, where the reduced state is of " <>
"the modes in subsystem";
averagePhotonNumberInSubsystem::usage="averagePhotonNumberInSubsystem[v,subsystem] " <>
"computes the expected value of the photon number operator in the reduced " <>
"state described by the covariance matrix v, where the reduced state is of " <>
"the modes in subsystem";
logNegativityOfSubsystem::usage="logNegativityOfSubsystem[v,subsystem] " <>
"computes the logorithmic negativity of the state described by the covariance " <>
"matrix v when considering subsystem A as being the modes in subsystem " <>
"and subsystem B as being the rest of " <>
"the modes. If the logorithmic negativity is larger than 0, then the state " <>
"is not separable. However, if it is 0, it may or may not be separable.";


purityOfLeftPartition::usage="purityOfLeftPartition[v,partition] " <>
"computes the purity of the reduced state of the state described by the " <>
"covariance matrix v, where the reduced state is of the modes 1 through " <>
" `partition`. In other words, the reduced state is that of all modes " <>
"to the left of or at the partition.";
renyi2EntropyOfLeftPartition::usage="renyi2EntropyOfLeftPartition" <>
"[v,partition] computes the Renyi-2 entropy of the reduced state of the " <>
"state described by the covariance matrix v, where the reduced state is " <>
"of the modes 1 through `partition`. In other words, the reduced state " <>
"is that of all modes to the left of or at the partition.";
entanglementEntropyOfLeftPartition::usage="entanglementEntropyOfLeftPartition" <>
"[v,partition] computes the entanglement entropy of the reduced state of the " <>
"state described by the covariance matrix v, where the reduced state is " <>
"of the modes 1 through `partition`. In other words, the reduced state " <>
"is that of all modes to the left of or at the partition.";
energyOfLeftPartition::usage="energyOfLeftPartition[v,partition] computes the energy of the reduced " <>
"state described by the covariance matrix v, where the reduced state is of " <>
"modes 1 through `partition`. In other words, the reduced state " <>
"is that of all modes to the left of or at the partition.";
averagePhotonNumberInLeftPartition::usage="averagePhotonNumberInLeftPartition[v,partition] " <>
"computes the expected value of the photon number operator in the reduced " <>
"state described by the covariance matrix v, where the reduced state is of " <>
"modes 1 through `partition`. In other words, the reduced state " <>
"is that of all modes to the left of or at the partition.";
logNegativityOfLeftPartition::usage="logNegativityOfLeftPartition[v,partition] " <>
"computes the logorithmic negativity of the state described by the covariance " <>
"matrix v when considering subsystem A as being the modes between " <>
"1 and partition (to the left of or at `partition`) and subsystem B as being " <>
"the rest of the modes. If the logorithmic negativity is larger than 0, then " <>
"the state is not separable. However, if it is 0, it may or may not be separable.";


purityOfRightPartition::usage="purityOfRightPartition[v,partition] " <>
"computes the purity of the reduced state of the state described by the " <>
"covariance matrix v, where the reduced state is of the modes `partition` through " <>
" the end. In other words, the reduced state is that of all modes " <>
"to the right of or at the partition.";
renyi2EntropyOfRightPartition::usage="renyi2EntropyOfRightPartition" <>
"[v,partition] computes the Renyi-2 entropy of the reduced state of the " <>
"state described by the covariance matrix v, where the reduced state is " <>
"of the modes `partition` through the end. In other words, the reduced state " <>
"is that of all modes to the right of or at the partition.";
entanglementEntropyOfRightPartition::usage="entanglementEntropyOfRightPartition" <>
"[v,partition] computes the entanglement entropy of the reduced state of the " <>
"state described by the covariance matrix v, where the reduced state is " <>
"of the modes `partition` through the end. In other words, the reduced state " <>
"is that of all modes to the right of or at the partition.";
energyOfRightPartition::usage="energyOfRightPartition[v,partition] computes " <>
"the energy of the reduced state described by the covariance matrix v, where " <>
"the reduced state is of modes `partition` though the end. In other words, the " <>
"reduced state is that of all modes to the right of or at the partition.";
averagePhotonNumberInRightPartition::usage="averagePhotonNumberInRightPartition[v,partition] " <>
"computes the expected value of the photon number operator in the reduced " <>
"state described by the covariance matrix v, where the reduced state is of " <>
"modes `partition` through the end. In other words, the reduced state " <>
"is that of all modes to the right of or at the partition.";
logNegativityOfRightPartition::usage="logNegativityOfRightPartition[v,partition] " <>
"computes the logorithmic negativity of the state described by the covariance " <>
"matrix v when considering subsystem A as being the modes between " <>
"`partition` and the end (to the right of or at `partition`) and subsystem B as being " <>
"the rest of the modes. If the logorithmic negativity is larger than 0, then " <>
"the state is not separable. However, if it is 0, it may or may not be separable.";


(* ::Subsection::Closed:: *)
(*Simulation functions*)


simulate::usage="PARAMETERS
    squeeze: array of n squeezing parameters, where n is the number of modes.
    maxdepth: simulate this many layers.
    iters: the number of times to perform the simulation in order to estimate the
        average.
    quantities: these are the quantities to simulate. quantities should be an
        array of functions of the covariance matrix. After applying each layer
        to the covariance matrix V, quantities[[i]][V] will be called and stored
        for each i. This function returns an array of all such calculated values.
    layer: layer is a function so that layer[] outputs an orthogonal
        symplectic matrix that operates on the covariance matrix at a time step.
    RETURNS
    an array with dimensions Length[quantities]\[Cross](maxdepth+1).
";


simulateWithMeasurements::usage="PARAMETERS
    rate: number between 0 and 1. After application of a layer, each mode will
        be measured and projected onto the vacuum state |0> with probability `rate`.
    squeeze: array of n squeezing parameters, where n is the number of modes.
    maxdepth: simulate this many layers.
    iters: the number of times to perform the simulation in order to estimate the
        average.
    quantities: these are the quantities to simulate. quantities should be an
        array of functions of the covariance matrix. After applying each layer
        to the covariance matrix V, quantities[[i]][V] will be called and stored
        for each i. This function returns an array of all such calculated values.
    layer: layer is a function so that layer[] outputs an orthogonal
        symplectic matrix that operates on the covariance matrix at a time step.
    RETURNS
    an array with dimensions Length[quantities]\[Cross](maxdepth+1).
";


(* ::Subsection::Closed:: *)
(*Plotting functions*)


plotLinearLinearQuantities::usage="PARAMETERS
quantities: should have dimension Length[legends]\[Cross]Length[quantityLabels]\[Cross](maxdepth+1).
quantityLabels: label for eqch quantity.
legendLabels: legends label.
depPerLayer: defaul 2; what the depth of one layer is.
RETURNS
A table (of length Length[quantityLabels]) of plots for each quantity.";


plotLinearLogQuantities::usage="PARAMETERS
quantities: should have dimension Length[legends]\[Cross]Length[quantityLabels]\[Cross](maxdepth+1).
quantityLabels: label for eqch quantity.
legendLabels: legends label.
depPerLayer: defaul 2; what the depth of one layer is.
RETURNS
A table (of length Length[quantityLabels]) of plots for each quantity.";


plotLogLinearQuantities::usage="PARAMETERS
quantities: should have dimension Length[legends]\[Cross]Length[quantityLabels]\[Cross](maxdepth+1).
quantityLabels: label for eqch quantity.
legendLabels: legends label.
depPerLayer: defaul 2; what the depth of one layer is.
RETURNS
A table (of length Length[quantityLabels]) of plots for each quantity.";


plotLogLogQuantities::usage="PARAMETERS
quantities: should have dimension Length[legends]\[Cross]Length[quantityLabels]\[Cross]maxdepth.
quantityLabels: label for eqch quantity.
legendLabels: legends label.
depPerLayer: defaul 2; what the depth of one layer is.
RETURNS
A table (of length Length[quantityLabels]) of plots for each quantity.";


(* ::Section:: *)
(*Private context*)


Begin["`Private`"]; (* Begin Private Context *) 


(* ::Subsection::Closed:: *)
(*Matrix functions*)


matrixDirectSum[args_]:=Module[{f,pargs},
    (*get rid of empty matrices *)
	pargs=Select[args, UnsameQ[#, {}] &];
	If[Length[pargs]==0, Return@{}];
    f=pargs[[1]];
	Do[
		f=ArrayFlatten[{{f,0},{0,pargs[[i]]}}];,
		{i,2,Length[pargs]}
	];
	f
];


symplecticForm[n_]:=matrixDirectSum@Table[{{0,1},{-1,0}},n];


symplecticEigenvalues[v_]:=Select[
    Eigenvalues[I symplecticForm[Length[v]/2].v]//N//Chop,
    Positive
];


initialSqueezedState[squeeze_]:=DiagonalMatrix[
    Join@@Table[{Exp[-2s],Exp[2s]}, {s,squeeze}]
];


validCovarianceMatrixQ[\[Sigma]_]:=SymmetricMatrixQ[\[Sigma]]&&
    And@@NonNegative@Eigenvalues[\[Sigma]+I symplecticForm[Length[\[Sigma]]/2]]&&
    PositiveDefiniteMatrixQ@\[Sigma];


reducedCovarianceMatrix[\[Sigma]_,subsystem_]:=With[
    {s=Flatten@Table[{2i-1,2i},{i,subsystem}]},
    \[Sigma][[s,s]]
];
reducedCovarianceMatrixOfLeftPartition[\[Sigma]_,partition_]:=reducedCovarianceMatrix[
    \[Sigma],Range[partition]
];
reducedCovarianceMatrixOfRightPartition[\[Sigma]_,partition_]:=reducedCovarianceMatrix[
    \[Sigma],Range[partition,Length[\[Sigma]]/2]
];


squeezerMatrix[r_]:=DiagonalMatrix@{Exp[r],Exp[-r]};


beamsplitterMatrix[\[Theta]_]:=({
 {Cos[\[Theta]], 0, -Sin[\[Theta]], 0},
 {0, Cos[\[Theta]], 0, -Sin[\[Theta]]},
 {Sin[\[Theta]], 0, Cos[\[Theta]], 0},
 {0, Sin[\[Theta]], 0, Cos[\[Theta]]}
});


phaseshiftMatrix[\[Theta]_]:=({
 {Cos[\[Theta]], Sin[\[Theta]]},
 {-Sin[\[Theta]], Cos[\[Theta]]}
});


twoModeMatrix[\[Theta]_,\[Phi]_,\[Beta]_]:=Dot[
    matrixDirectSum[{phaseshiftMatrix[\[Theta]],phaseshiftMatrix[\[Phi]]}],
    beamsplitterMatrix[\[Beta]]
];


(*Direct sum of m random symplectic matrices *)
ssrand[m_]:=matrixDirectSum[Table[N@haarModeMatrix@2,m]];
(*First column of brickwork, assuming n is even *)
oddSrand[n_]:=ssrand[n/2];
(*Second column of brickwork, assuming n is even and no periodic boundary conditions*)
evenSrand[n_]:=matrixDirectSum[{IdentityMatrix[2],ssrand[n/2-1],IdentityMatrix[2]}];
(*layer, assuming n is even*)
oneDimRandomBrickworkLayer[n_]:=(
    If[OddQ[n]||n\[NotElement]PositiveIntegers,Throw["n must be an even positive integer"]];
    Return[oddSrand[n].evenSrand[n]];
);


oneDimRandomActiveBrickworkLayer[n_,rmin_:-1,rmax_:1]:=Dot[
    (* Euler decompostion *)
    oneDimRandomBrickworkLayer[n],
    matrixDirectSum@Table[N@squeezerMatrix@RandomReal@{rmin,rmax},n],
    oneDimRandomBrickworkLayer[n]
];


(*Direct sum of m fixed symplectic matrices *)
ssfixed[m_,\[Theta]_,\[Phi]_]:=matrixDirectSum[Table[N@twoModeMatrix[\[Phi],\[Phi],\[Theta]],m]];
(*First column of brickwork, assuming n is even *)
oddSfixed[n_,\[Theta]_,\[Phi]_]:=ssfixed[n/2,\[Theta],\[Phi]];
(*Second column of brickwork, assuming n is even and no periodic boundary conditions*)
evenSfixed[n_,\[Theta]_,\[Phi]_]:=matrixDirectSum[{
    IdentityMatrix[2],ssfixed[n/2-1,\[Theta],\[Phi]],IdentityMatrix[2]
}];
(*layer, assuming n is even*)
oneDimFixedBrickworkLayer[n_,\[Theta]_,\[Phi]_]:=(
    If[OddQ[n]||n\[NotElement]PositiveIntegers,Throw["n must be an even positive integer"]];
    Return[oddSfixed[n,\[Theta],\[Phi]].evenSfixed[n,\[Theta],\[Phi]]];
);


(*Direct sum of m fixed symplectic matrices *)
ssfixedbeam[m_,\[Theta]_]:=matrixDirectSum[N@Table[
    twoModeMatrix[#1,#2,\[Theta]]&@@RandomReal[2Pi,2],
    m
]];
(*First column of brickwork, assuming n is even *)
oddSfixedbeam[n_,\[Theta]_]:=ssfixedbeam[n/2,\[Theta]];
(*Second column of brickwork, assuming n is even and no periodic boundary conditions*)
evenSfixedbeam[n_,\[Theta]_]:=matrixDirectSum[{
    IdentityMatrix[2],ssfixedbeam[n/2-1,\[Theta]],IdentityMatrix[2]
}];
(*layer, assuming n is even*)
oneDimFixedBeamsplitterBrickworkLayer[n_,\[Theta]_]:=(
    If[OddQ[n]||n\[NotElement]PositiveIntegers,Throw["n must be an even positive integer"]];
    Return[oddSfixedbeam[n,\[Theta]].evenSfixedbeam[n,\[Theta]]];
);


applyByCongruence[v_,s_]:=s.v.s\[Transpose];


haarModeMatrix[n_]:=Module[
    {
        u=RandomVariate@CircularUnitaryMatrixDistribution[n],
        order=Join@@Table[{i,i+n},{i,n}]
    },
    u=ArrayFlatten@{{Re@u,-Im@u},{Im@u,Re@u}};
    (*now put it in direct product form*)
    u[[order,order]]
];


haarTwoModeMatrix[n_,i_,j_]:=Module[
    {
        u=IdentityMatrix[2n],
        indices={2i-1,2i,2j-1,2j}
    },
    u[[indices,indices]]=haarModeMatrix[2];
    u
];


measureVacuumLastMode[\[Sigma]_]:=Module[(*Serafini page 129, measuring |0><0|*)
  {n=Length[\[Sigma]]/2,\[Sigma]A,\[Sigma]B,\[Sigma]AB,\[Sigma]m=IdentityMatrix@2,\[Sigma]Ap},
  \[Sigma]A=\[Sigma][[Range[2n-2],Range[2n-2]]];
  \[Sigma]B=\[Sigma][[{2n-1,2n},{2n-1,2n}]];
  \[Sigma]AB=\[Sigma][[Range[2n-2],{2n-1,2n}]];
  \[Sigma]Ap=\[Sigma]A-\[Sigma]AB.Inverse[\[Sigma]B+\[Sigma]m].\[Sigma]AB\[Transpose];
  matrixDirectSum@{\[Sigma]Ap,\[Sigma]m}
];

measureVacuumMode[\[Sigma]_,mode_]:=With[
  {
    orderswap=If[
     mode==Length[\[Sigma]]/2,
     Range[2mode],
     Join[
       Range[2mode-2],
       {Length@\[Sigma]-1,Length@\[Sigma]},
       Range[2mode+1,Length[\[Sigma]]-2],
       {2mode-1,2mode}
     ]
    ]
  },
  measureVacuumLastMode[\[Sigma][[orderswap,orderswap]]][[orderswap,orderswap]]
];

measureVacuumModes[\[Sigma]_,modes_]:=Module[{v=\[Sigma]},
  Do[v=measureVacuumMode[v,m],{m,modes}];
  v
];


(* ::Subsection::Closed:: *)
(*Quantities from the covariance matrix*)


renyi2EntropyFromPurity[purity_]:=-Log[purity]//N;


purity[v_]:=If[Length@v==0,1,1/Sqrt@Det@v//N];
renyi2Entropy[v_]:=renyi2EntropyFromPurity@purity@v;
entanglementEntropy[v_]:=Sum[
    N[(\[Nu]+1)/2 Log[(\[Nu]+1)/2]-If[\[Nu]==1,0,(\[Nu]-1)/2 Log[(\[Nu]-1)/2]]],
    {\[Nu],symplecticEigenvalues[v]}
];
energy[v_]:=Tr[v]//N;
averagePhotonNumber[v_]:=1/4 (energy[v]-Length[v]);


subv[v_,s_]:=reducedCovarianceMatrix[v,s];
purityOfSubsystem[v_,subsystem_]:=purity@subv[v,subsystem];
renyi2EntropyOfSubsystem[v_,subsystem_]:=renyi2Entropy@subv[v,subsystem];
entanglementEntropyOfSubsystem[v_,subsystem_]:=entanglementEntropy@subv[v,subsystem];
energyOfSubsystem[v_,subsystem_]:=energy@subv[v,subsystem];
averagePhotonNumberInSubsystem[v_,subsystem_]:=averagePhotonNumber@subv[v,subsystem];
logNegativityOfSubsystem[v_,subsystem_]:=With[
    {
        t=matrixDirectSum@Table[
            If[MemberQ[subsystem,i],IdentityMatrix[2],PauliMatrix[3]],
            {i,Length[v]/2}
        ]
    },
    Sum[Max[0,-Log2@x],{x,symplecticEigenvalues[t.v.t]}]
];


left[p_]:=Range[p];
purityOfLeftPartition[v_,partition_]:=purityOfSubsystem[v,left@partition];
renyi2EntropyOfLeftPartition[v_,partition_]:=renyi2EntropyOfSubsystem[v,left@partition];
entanglementEntropyOfLeftPartition[v_,partition_]:=entanglementEntropyOfSubsystem[
    v,left@partition
];
energyOfLeftPartition[v_,partition_]:=energyOfSubsystem[v,left@partition];
averagePhotonNumberInLeftPartition[v_,partition_]:=averagePhotonNumberInSubsystem[v,left@partition];
logNegativityOfLeftPartition[v_,partition_]:=logNegativityOfSubsystem[v,left@partition];


right[p_,v_]:=Range[p,Length[v]/2];
purityOfRightPartition[v_,partition_]:=purityOfSubsystem[v,partition~right~v];
renyi2EntropyOfRightPartition[v_,partition_]:=renyi2EntropyOfSubsystem[v,partition~right~v];
entanglementEntropyOfRightPartition[v_,partition_]:=entanglementEntropyOfSubsystem[
    v,partition~right~v
];
energyOfRightPartition[v_,partition_]:=energyOfSubsystem[v,partition~right~v];
averagePhotonNumberInRightPartition[v_,partition_]:=averagePhotonNumberInSubsystem[v,partition~right~v];
logNegativityOfRightPartition[v_,partition_]:=logNegativityOfSubsystem[v,partition~right~v];


(* ::Subsection::Closed:: *)
(*Simulation functions*)


simulate[squeeze_,maxdepth_,quantities_,layer_]:=(*
    PARAMETERS
    squeeze: array of n squeezing parameters, where n is the number of modes.
    maxdepth: simulate this many layers.
    quantities: these are the quantities to simulate. quantities should be an
        array of functions of the covariance matrix. After applying each layer
        to the covariance matrix V, quantities[[i]][V] will be called and stored
        for each i. This function returns an array of all such calculated values.
    layer: layer is a function so that layer[] outputs an orthogonal
        symplectic matrix that operates on the covariance matrix at a time step.
    RETURNS
    an array with dimensions Length[quantities]\[Cross](maxdepth+1).
    *)
    Module[{state,calc},
        state=initialSqueezedState[squeeze];
        calc[st_]:=Table[q[st],{q,quantities}];
        Join[
            {calc[state]},
            Table[
                state=state~applyByCongruence~layer[];
                calc[state],
                {depth,maxdepth}
            ]
        ]//Transpose
    ];


simulate[squeeze_,maxdepth_,iters_,quantities_,layer_]:=Mean[
    Table[simulate[squeeze,maxdepth,quantities,layer],iters]
];


simulateWithMeasurements[rate_,squeeze_,maxdepth_,quantities_,layer_]:=(*
    PARAMETERS
    rate: number between 0 and 1. After application of a layer, each mode will
        be measured and projected onto the vacuum state |0> with probability `rate`.
    squeeze: array of n squeezing parameters, where n is the number of modes.
    maxdepth: simulate this many layers.
    quantities: these are the quantities to simulate. quantities should be an
        array of functions of the covariance matrix. After applying each layer
        to the covariance matrix V, quantities[[i]][V] will be called and stored
        for each i. This function returns an array of all such calculated values.
    layer: layer is a function so that layer[] outputs an orthogonal
        symplectic matrix that operates on the covariance matrix at a time step.
    RETURNS
    an array with dimensions Length[quantities]\[Cross](maxdepth+1).
    *)
    Module[{state,calc,modesToMeasure},
        state=initialSqueezedState[squeeze];
        calc[st_]:=Table[q[st],{q,quantities}];
        Join[
            {calc[state]},
            Table[
                state=state~applyByCongruence~layer[];
                modesToMeasure=Select[Range[Length@squeeze],RandomReal[]<rate&];
                state=measureVacuumModes[state,modesToMeasure];
                calc[state],
                {depth,maxdepth}
            ]
        ]//Transpose
    ];


simulateWithMeasurements[rate_,squeeze_,maxdepth_,iters_,quantities_,layer_]:=Mean[
    Table[simulateWithMeasurements[rate,squeeze,maxdepth,quantities,layer],iters]
];


(* ::Subsection::Closed:: *)
(*Plotting functions*)


plotQuantities[plotFun_,quantities_,quantityLabels_,legendLabels_,depPerLayer_]:=Table[
    plotFun[
        Table[
          Table[{depPerLayer(dep-1),quantities[[j,i]][[dep]]},{dep,Dimensions[quantities][[3]]}],
          {j,Length@quantities}
        ],
        PlotLegends->legendLabels,PlotRange->All,ImageSize->Large,
        AxesLabel->{Style["depth",FontSize->16],Style[quantityLabels[[i]],FontSize->16]},
        Joined->True
    ],
    {i,1,First[quantities,0]//Length}
];


plotLinearLinearQuantities[quantities_,quantityLabels_,legendLabels_,depPerLayer_:2]:=plotQuantities[
    ListLinePlot,quantities,quantityLabels,legendLabels,depPerLayer
];


plotLinearLogQuantities[quantities_,quantityLabels_,legendLabels_,depPerLayer_:2]:=plotQuantities[
    ListLogPlot,quantities,quantityLabels,legendLabels,depPerLayer
];


plotLogLinearQuantities[quantities_,quantityLabels_,legendLabels_,depPerLayer_:2]:=plotQuantities[
    ListLogLinearPlot,quantities,quantityLabels,legendLabels,depPerLayer
];


plotLogLogQuantities[quantities_,quantityLabels_,legendLabels_,depPerLayer_:2]:=plotQuantities[
    ListLogLogPlot,quantities,quantityLabels,legendLabels,depPerLayer
];


(* ::Section:: *)
(*End*)


End[]; (* End Private Context *)


EndPackage[]
