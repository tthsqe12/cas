Print["dirichlet"];

Do[
    ms = Select[Range[q], CoprimeQ[q, #]&];
    Do[
        chi1 = DirichletCharacter[q, m1];
        chi2 = DirichletCharacter[q, m2];
        chi3 = DirichletCharacter[q, m1*m2];
        Do[
            Assert[chi1[n]*chi2[n] == chi3[n]];
        , {n, -q, q}];
    , {m1, ms}, {m2, ms}];
, {q, 2, 15}]