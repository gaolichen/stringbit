ElementA[k_, l_, M_, n_, m_] := Exp[TPI*k*((m + n)/M)]*
      Sin[(m - n)*(Pi/M/2)] + Exp[TPI*(k + l)*((m + n)/M/2)]*
      Sin[(l - k - 1/2)*(m - n)*(Pi/M)]
 
TPI = (2*I)*Pi
 
MatrixA2[k_, l_, M_] := Table[ElementA[k, l, M, n, m], {n, 0, M - 1}, 
     {m, 0, M - 1}]
 
MatrixAV2[M_, L_, opts:OptionsPattern[]] := If[OptionValue[SymmetricA], 
     Table[(Sin[(m - n)*(Pi/(2*M))]/2)*(1 + Exp[TPI*L*((m + n)/M)] + 
        Exp[Pi*I*((2*L*n + m + n)/M)] + Exp[Pi*I*((2*L*m + m + n)/M)]), 
      {n, 0, M - 1}, {m, 0, M - 1}], Table[Sin[(m - n)*(Pi/(2*M))] + 
       Exp[TPI*(L + 1)*((m + n)/(2*M))]*Sin[(L + 1/2)*(m - n)*(Pi/M)], 
      {n, 0, M - 1}, {m, 0, M - 1}]]
 
Options[MatrixAV2] = {SymmetricA -> False}
 
MatrixAW2[M_, L_] := Table[(1 + Exp[I*Pi*((n + m)/M)])*
      (1 + Exp[2*I*Pi*L*((m + n)/M)])*Sin[(m - n)*(Pi/(2*M))], {n, 0, M - 1}, 
     {m, 0, M - 1}]
 
OmegaV2[M_, L_, opts:OptionsPattern[]] := 
    InverseCT2[M, L, FilterRules[{opts}, Options[InverseCT2]]] . 
     ConjugateTranspose[MatrixAV2[M, L, FilterRules[{opts}, 
        Options[MatrixAV2]]]] . Transpose[InverseCT2[M, L, 
       FilterRules[{opts}, Options[InverseCT2]]]]
 
Options[OmegaV2] = {SymmetricA -> False, Numeric -> False}
 
InverseCT2[M_, L_, opts:OptionsPattern[]] := 
    Module[{MC}, MC = MatrixCT2[M, L, FilterRules[{opts}, 
         Options[MatrixCT2]]]; Return[Inverse[MC]]; ]
 
Options[InverseCT2] = {Numeric -> False}
 
MatrixCT2[M_, L_, opts:OptionsPattern[]] := If[OptionValue[Numeric], 
     Return[N[Table[ElementCT[M, L, M - L, n, m], {n, 0, M - 1}, 
        {m, 0, M - 1}]]], Return[Table[ElementCT[M, L, M - L, n, m], 
       {n, 0, M - 1}, {m, 0, M - 1}]]]
 
Options[MatrixCT2] = {Numeric -> False}
 
ElementCT[M_, L_, K_, m_, n_] := If[m == 0, Return[Cm0T2[M, L, K, n]], 
     If[n == 0, Return[0], If[n < L, Return[Cmn1T[M, L, K, m, n]], 
       If[n < M - 1, Return[Cmn2T[M, L, K, m, n - L + 1]], 
        If[n == M - 1, Return[CmMT[M, L, K, m]], 
         Print["Error: indeices out of range. m=", m, ", n=", n]]]]]]
 
Cm0T2[M_, L_, K_, m_] := If[m == 0, Return[1], Return[0]]
 
Cmn1T[M_, L_, K_, m_, n_] := DotV1T[M, L, K, m, n]*
     Cos[Pi*(n/(2*L)) - Pi*(m/(2*M))]
 
DotV1T[M_, L_, K_, m_, n_] := If[IntegerQ[n/L - m/M], Sqrt[L/M], 
     (-Sqrt[M*L]^(-1))*((1 - Exp[(-TPI)*m*(L/M)])/
       (1 - Exp[(-TPI)*(n/L - m/M)]))]
 
Cmn2T[M_, L_, K_, m_, n_] := DotV2T[M, L, K, m, n]*
     Cos[Pi*(n/(2*K)) - Pi*(m/(2*M))]
 
DotV2T[M_, L_, K_, m_, n_] := If[IntegerQ[n/K - m/M], 
     Sqrt[K/M]*Exp[(-TPI)*n*(L/K)], (1/Sqrt[M*K])*((1 - Exp[(-TPI)*m*(L/M)])/
       (1 - Exp[(-TPI)*(n/K - m/M)]))]
 
CmMT[M_, L_, K_, m_] := DotWT[M, L, K, m]*Cos[m*(Pi/(2*M)) - Pi/4]
 
DotWT[M_, L_, K_, m_] := If[m == 0, 0, (-Sqrt[L*K]^(-1))*
      ((1 - Exp[(-TPI)*m*(L/M)])/(1 - Exp[TPI*(m/M)]))]
 
OmegaW2[M_, L_, opts:OptionsPattern[]] := 
    InverseCT2[M, L, FilterRules[{opts}, Options[InverseCT2]]] . 
     ConjugateTranspose[MatrixAW2[M, L]] . 
     Transpose[InverseCT2[M, L, FilterRules[{opts}, Options[InverseCT2]]]]
 
GammaPV2[M_, L_, opts:OptionsPattern[]] := If[OptionValue[SymmetricA], 
     -Cot[Pi/(2*M)] + (Cot[(2*L - 1)*(Pi/2/M)] - Cot[(2*L + 1)*(Pi/2/M)])/2 - 
      GammaV2[M, L, FilterRules[{opts}, Options[GammaV2]]], 
     -Cot[Pi/(2*M)] - Cot[(2*L + 1)*(Pi/2/M)] - GammaV2[M, L, 
       FilterRules[{opts}, Options[GammaV2]]]]
 
Options[GammaPV2] = {Numeric -> False, SymmetricA -> False}
 
GammaV2[M_, L_, opts:OptionsPattern[]] := 
    Tr[Conjugate[MatrixST2[M, L, FilterRules[{opts}, Options[MatrixST2]]]] . 
      InverseCT2[M, L, FilterRules[{opts}, Options[MatrixST2]]] . 
      ConjugateTranspose[MatrixAV2[M, L, FilterRules[{opts}, 
         Options[MatrixAV2]]]]]
 
Options[GammaV2] = {Numeric -> False, SymmetricA -> False}
 
MatrixST2[M_, L_, opts:OptionsPattern[]] := If[OptionValue[Numeric], 
     Return[N[Table[ElementST[M, L, M - L, n, m], {n, 0, M - 1}, 
        {m, 0, M - 1}]]], Return[Table[ElementST[M, L, M - L, n, m], 
       {n, 0, M - 1}, {m, 0, M - 1}]]]
 
Options[MatrixST2] = {Numeric -> False}
 
ElementST[M_, L_, K_, m_, n_] := If[m == 0, Return[Sm0T2[M, L, K, n]], 
     If[n == 0, Return[0], If[n < L, Return[Smn1T[M, L, K, m, n]], 
       If[n < M - 1, Return[Smn2T[M, L, K, m, n - L + 1]], 
        If[n == M - 1, Return[SmMT[M, L, K, m]], 
         Print["Error: indeices out of range. m=", m, ", n=", n]]]]]]
 
Sm0T2[M_, L_, K_, m_] := 0
 
Smn1T[M_, L_, K_, m_, n_] := DotV1T[M, L, K, m, L - n]*
     Cos[Pi*(n/(2*L)) + Pi*(m/(2*M))]
 
Smn2T[M_, L_, K_, m_, n_] := DotV2T[M, L, K, m, K - n]*
     Cos[Pi*(n/(2*K)) + Pi*(m/(2*M))]
 
SmMT[M_, L_, K_, m_] := DotWT[M, L, K, m]*Cos[m*(Pi/(2*M)) + Pi/4]
 
GammaPW2[M_, L_, opts:OptionsPattern[]] := -4*Cot[Pi/(2*M)] - 
     GammaW2[M, L, FilterRules[{opts}, Options[GammaW2]]]
 
GammaW2[M_, L_, opts:OptionsPattern[]] := 
    Tr[Conjugate[MatrixST2[M, L, FilterRules[{opts}, Options[MatrixST2]]]] . 
      InverseCT2[M, L, FilterRules[{opts}, Options[MatrixST2]]] . 
      ConjugateTranspose[MatrixAW2[M, L]]]
 
Options[GammaW2] = {Numeric -> False}
 
VerifyMatrices2[M_] := Module[{L, mat0, v0, res}, 
     mat0 = ConstantArray[0, {M, M}]; v0 = ConstantArray[0, M]; 
      For[L = 1, L < M, L++, 
       res = Chop[N[MatrixA2[M, L + 1, M] - MatrixAV2[M, L]]]; 
        Assert[res == mat0]; res = Chop[N[MatrixA2[L, L + 1, M] + 
            MatrixA2[M, 1, M] - MatrixAW2[M, L]]]; Assert[res == mat0]; ]; 
      Return[True]; ]
 
Attributes[Assert] = {HoldAllComplete}
