Print["map"];

Print["... Map"];

Assert[Map[f, {{a, b}, {c, d, e}}] === {f[{a, b}], f[{c, d, e}]}];

Assert[Map[f, {{a, b}, {c, d, e}}, {2}] === {{f[a], f[b]}, {f[c], f[d], f[e]}}];

Assert[Map[f, {{a, b}, {c, d, e}}, 2] === {f[{f[a], f[b]}], f[{f[c], f[d], f[e]}]}];

Assert[Map[f, h1[h2[h3[x]]], -1] === h1[f[h2[f[h3[f[x]]]]]]];

Assert[Map[f, h1[h2[h3[x]]], {0, -1}] === f[h1[f[h2[f[h3[f[x]]]]]]]];

Print["... Apply"];

Assert[Apply[f, {a, b, c, d}] === f[a, b, c, d]];

Assert[Apply[f, {{a, b, c}, {d, e}}] === f[{a, b, c}, {d, e}]];

Assert[Apply[f, {{a, b, c}, {d, e}}, {1}] === {f[a, b, c], f[d, e]}];

Assert[Apply[f, {{a, b, c}, {d, e}}, {0, 1}] === f[f[a, b, c], f[d, e]]];

Assert[Apply[f, {{{{{a}}}}}, 2] === {f[f[{{a}}]]}];

Assert[Apply[f, {{{{{a}}}}}, {0, 2}] === f[f[f[{{a}}]]]];

Assert[Apply[f, {{{{{a}}}}}, Infinity] === {f[f[f[f[a]]]]}];

Assert[Apply[f, {{{{{a}}}}}, {0, Infinity}] === f[f[f[f[f[a]]]]]];

Assert[Apply[f, {{{{{a}}}}}, -1] === {f[f[f[f[a]]]]}];

Assert[Apply[f, {{{{{a}}}}}, -2] === {f[f[f[f[a]]]]}];

Assert[Apply[f, {{{{{a}}}}}, -3] === {f[f[f[{a}]]]}];

Assert[Apply[f, {{{{{a}}}}}, {2, -3}] === {{f[f[{a}]]}}];

Assert[Apply[f, h0[h1[h2[h3[h4[a]]]]], {2, -3}] === h0[h1[f[f[h4[a]]]]]];

Assert[Apply[f, a] === a];

Assert[Apply[f, {a, "string", 3}, {-1}] === {a, "string", 3}];
