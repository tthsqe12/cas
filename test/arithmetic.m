Print["arithmetic"];

Print["... zeta"];

Assert[Zeta/@Range[-12,0] === {0, 691/32760, 0, -1/132, 0, 1/240,
                               0, -1/252, 0, 1/120, 0, -1/12, -1/2}];
Assert[Zeta/@Range[2,12] === {Pi^2/6, Zeta[3], Pi^4/90, Zeta[5],
                              Pi^6/945, Zeta[7], Pi^8/9450, Zeta[9],
                              Pi^10/93555, Zeta[11], 691*Pi^12/638512875}];

Print["... gamma"];

Assert[Gamma[1] == 1];
Assert[Gamma[2] == 1];
Assert[Gamma[3] == 2];
Assert[Gamma[4] == 6];
Assert[Gamma[-1/2] == -2*Sqrt[Pi]];
Assert[Gamma[1/2] == Sqrt[Pi]];
Assert[Gamma[3/2] == 1/2*Sqrt[Pi]];

Assert[Factorial[0] == 1];
Assert[Factorial[1] == 1];
Assert[Factorial[2] == 2];
Assert[Factorial[3] == 6];
Assert[Factorial[-3/2] == -2*Sqrt[Pi]];
Assert[Factorial[-1/2] == Sqrt[Pi]];
Assert[Factorial[1/2] == 1/2*Sqrt[Pi]];

Print["... compare"];

Assert[N[1] == N[1+2^-46]];
Assert[N[1+10I] == N[1+2^-46+10I]];
Assert[Not[N[1] == N[1+2^-45]]];
Assert[Not[N[1+10I] == N[1+2^-45+10I]]];
Assert[Not[2.00006 - 2.00005 == 0.00001]];


Print["... integer"];

Assert[1 + 2 == 3];

Do[
    Assert[Sum[n^0, {n, 1, m}] == m];
    Assert[Sum[n^1, {n, 1, m}] == m*(m + 1)/2];
    Assert[Sum[n^2, {n, 1, m}] == m*(m + 1)*(2*m + 1)/6];
    Assert[Sum[n^3, {n, 1, m}] == m^2*(m + 1)^2/4];
, {m, 1, 50}];

Print["... rational"];

Do[
    t = Sum[1/n, {n, 1, m}] + Sum[1/n, {n, m + 1, 2*m}];
    s = Sum[1/n, {n, 1, 2*m}];
    Assert[t == s]
, {m, 1, 50}];

Do[
    Assert[m/(m + 1) == Sum[1/(n*(n+1)), {n, 1, m}]];
    Assert[m + 1 == Product[1 + 1/n, {n, 1, m}]];
    Assert[Fibonacci[m]/Fibonacci[m + 1] == Nest[1/(1 + #)&, 0, m]]
, {m, 1, 50}];

Print["... elementary functions"]

Do[
    Assert[Exp[Log[x]] == x];
    Assert[Sin[ArcSin[x]] == x];
    Assert[Cos[ArcCos[x]] == x];
    Assert[Tan[ArcTan[x]] == x];
, {x, {1.2 + I*3.4, 1.2`100 + I*3.4`100}}]