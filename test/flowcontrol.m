Print["flowcontrol"];

Print["... if"]

ClearAll[a, b, c, d];

Assert[Head[If[a, b]] === If];
Assert[If[True, b] === b];
Assert[If[False, b] === Null];
Assert[Head[If[a, b, c]] === If];
Assert[If[True, b, c] === b];
Assert[If[False, b, c] === c];
Assert[If[a, b, c, d] === d];
Assert[If[True, b, c, d] === b];
Asssert[If[False, b, c, d] === c];

Print["... return"]

ClearAll[f]; f[x_] := Return[x, f]; Assert[f[5] === 5]

ClearAll[f]; f[x_] := {Return[x, f]}; Assert[f[5] === 5]

ClearAll[f]; f[x_] := Block[{t = x}, Return[t, Block]]; Assert[f[5] === 5]

ClearAll[f]; f[x_] := Module[{t = x}, Return[t, Module]]; Assert[f[5] === 5]


Print["... goto"];

t = Reap[
    i = Pi;
    Block[{i = 1},
        Label["loop"];
        Sow[i];
        MatchQ[i, _ ? (If[# > 5, Sow[i*i]; Goto["out"]]&)];
        i = i + 1;
        Goto["loop"];
    ];
    Assert[False];
    Label["out"];
    i
];
Assert[t === {Pi, {{1, 2, 3, 4, 5, 6, 36}}}];

ClearAll[f];
f[n_] :=
    Block[{i = 1, p = 1},
        Label["loop"];
        p = p*i;
        i = i + 1;
        If[i <= n, Goto["loop"]];
        p
    ];
Do[Assert[f[n] === n!], {n, 1, 10}];
