SignMat[x_, y_] := 0
 
BuildSignMat[maxM_] := Module[{i, j, bits, cnt}, 
     For[i = 0, i < 2^maxM, i++, bits = Reverse[IntegerDigits[i, 2]]; 
        cnt = 0; For[j = maxM, j >= 1, j--, 
         If[j <= Length[bits] && bits[[j]] == 1, cnt = cnt + 1; Continue[]]; 
          SignMat[i, j - 1] = (-1)^cnt; ]; ]; ]
 
BuildBitCount[maxM_] := Module[{i}, For[i = 0, i < 2^maxM, i++, 
       BitCount[i] = Count[IntegerDigits[i, 2], 1]]; ]
 
BitCount[n_] := 0
 
MergeSign[a_, b_] := Module[{i, bits, ret = 1, aa}, 
     If[BitAnd[a, b] != 0, Return[0]]; If[b == 0 || a == 0, Return[1]]; 
      bits = IntegerDigits[b, 2]; aa = a; For[i = Length[bits], i >= 1, i--, 
       If[bits[[i]] == 0, Continue[]]; 
        ret = ret*SignMat[aa, Length[bits] - i]; 
        aa = BitSet[aa, Length[bits] - i]; ]; Return[ret]; ]
 
OutputEigenFunction[fun_, opts:OptionsPattern[]] := 
    Module[{i, j, ret, term}, ret = 0; For[i = 0, i < Length[fun], i++, 
       If[fun[[i + 1]] === 0, Continue[]]; If[OptionValue[FullSimp], 
         term = ToRadicals[FullSimplify[fun[[i + 1]]]], 
         term = Simplify[fun[[i + 1]]]]; For[j = 0, 2^j <= i, j++, 
         If[BitGet[i, j] == 1, term = term*Subscript[\[Theta], j + 1]]; ]; 
        ret = ret + term; ]; Return[ret]; ]
 
Options[OutputEigenFunction] = {FullSimp -> False}
 
Attributes[Subscript] = {NHoldRest}
 
OutputFunctionProd[m_, prod_, opts:OptionsPattern[]] := 
    Module[{i, j, res, term1, term2}, res = 0; For[i = 0, i < Length[prod], 
       i++, If[prod[[i + 1]] === 0, Continue[]]; term1 = 1; term2 = 1; 
        For[j = 0, j < m, j++, If[BitGet[i, j] == 0, 
           term1 = term1*Subscript[\[Theta], j + 1], 
           term2 = term2*Subscript[\[Phi], j + 1]]; ]; 
        If[OptionValue[FullSimp], res += ToRadicals[FullSimplify[
             prod[[i + 1]]]]*term1*term2, res += Simplify[prod[[i + 1]]]*
           term1*term2]; ]; Return[res]]
 
Options[OutputFunctionProd] = {FullSimp -> False}
 
GrassmannDiff[poly_, vars_] := Module[{i, j, sign, list, res, bits = 0}, 
     If[Length[vars] == 0, Return[poly]]; sign = Signature[vars]; 
      If[sign == 0, Return[{0}]]; For[i = 1, i <= Length[vars], i++, 
       bits = BitSet[bits, vars[[i]] - 1]]; list = Reverse[Sort[vars]]; 
      sign = sign*Signature[list]; res = ConstantArray[0, Length[poly]]; 
      For[i = 0, i < Length[poly], i++, If[BitAnd[i, bits] != bits, 
         Continue[]]; res[[i - bits + 1]] += poly[[i + 1]]*sign*
          MergeSign[bits, i - bits]; ]; Return[res]]
 
GrassmannDot[m_, poly1_, poly2_] := Module[{i, j, res}, 
     res = ConstantArray[0, 2^m]; For[i = 0, i < 2^m && i < Length[poly1], 
       i++, If[poly1[[i + 1]] === 0, Continue[]]; 
        For[j = 0, j < 2^m && j < Length[poly2], j++, 
         If[BitAnd[i, j] != 0 || poly2[[j + 1]] === 0, Continue[]]; 
          res[[i + j + 1]] += poly1[[i + 1]]*poly2[[j + 1]]*
            MergeSign[i, j]]; ]; Return[res]; ]
 
LeftShiftBits[bits_, n_] := Module[{i, j, res}, 
     res = ConstantArray[0, BitShiftLeft[Length[bits], n]]; 
      For[i = 1, i <= Length[bits], i++, res[[BitShiftLeft[i - 1, n] + 1]] = 
         bits[[i]]; ]; Return[res]]
