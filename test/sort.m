Print["string"]

Print["... Sort"];

Assert[Sort[{"ab", "abc", "bc", "bca"}] === {"ab", "abc", "bc", "bca"}];
Assert[Sort[{"abc", "ab", "bc", "bca"}] === {"ab", "abc", "bc", "bca"}];
Assert[Sort[{"abc", "bc", "ab", "bca"}] === {"ab", "abc", "bc", "bca"}];
Assert[Sort[{"abc", "bc", "bca", "ab"}] === {"ab", "abc", "bc", "bca"}];
Assert[Sort[{"bc", "abc", "bca", "ab"}] === {"ab", "abc", "bc", "bca"}];

Print["... StringJoin"];

Assert[Sort[Map[StringJoin, Tuples[{"a", "b", "A", "B"}, 2]]] === {"AA", "AB",
    "Aa", "Ab", "BA", "BB", "Ba", "Bb", "aA", "aB", "aa", "ab", "bA", "bB", "ba", "bb"}];
