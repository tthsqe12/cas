$AssertFunction = Exit[17]&;
Print["passed in ", First@Timing[
(*  Get["test/hypergeometric.m"];*)
  Get["test/flowcontrol.m"];
  Get["test/arithmetic.m"];
  Get["test/ContinuedFraction.m"];
  Get["test/dirichlet.m"];
  Get["test/sort.m"];
  Get["test/part.m"];
  Get["test/map.m"];
  Get["test/function.m"];
  Get["test/algebra.m"];
  Get["test/derivative.m"];
  Get["test/pattern.m"];
  Get["test/parse.m"];
  Get["test/modular_subgroups.m"];
]];
Exit[0];
