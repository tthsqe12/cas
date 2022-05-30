Print["hypergeometric"];

Print["  terminating"];

tab = {
    {{1, -3}, {2}, 10} -> -164,
    {{1/2, 1/3, -20}, {25/2}, -1} -> 6.07920215418685724687919*^8
};
Do[Assert[N[HypergeometricPFQ@@a[[1]]] == a[[2]]], {a, tab}];

Print["  singular points"];

tab = {
    {{1/2, 1/3}, {9/4}, 1/10} -> 1.007646373450970945394293353,
    {{1/2, 1/3}, {9/4}, 9/10} -> 1.102850891223585999379,
    {{1/2, 1/3}, {9/4}, 1} -> 1.129450433651130012503,
    {{1, 2}, {4}, 9/10} -> 2.1789423102929665152,
    {{1, 2}, {4}, 1} -> 3,
    {{1, 2, 3}, {4, 5}, 0} -> 1,
    {{1, 2, 3}, {4, 5}, 1} -> N[120 - 12*Pi^2]
};
Do[Assert[N[HypergeometricPFQ@@a[[1]]] == a[[2]]], {a, tab}];

Print["  algebraic case"];

Do[
    t = N[(Sqrt[1 - z] + Sqrt[-z])^(1/3)/2 + (Sqrt[1 - z] - Sqrt[-z])^(1/3)/2];
    Assert[N[HypergeometricPFQ[{-1/6, 1/6}, {1/2}, z]] == t];
, {z, {-1/2, 1/2}}];

Do[
    f1 = HypergeometricPFQ[{-1/12, 1/4}, {2/3}, 1728/j];
    f2 = HypergeometricPFQ[{1/4, 7/12}, {4/3}, 1728/j];
    t = N[j/3^3*f1^3/f2^3];
    Assert[j == 27*t*(t + 8)^3/(t - 1)^3];

    f1 = HypergeometricPFQ[{-1/24, 7/24}, {3/4}, 1728/j];
    f2 = HypergeometricPFQ[{5/24, 13/24}, {5/4}, 1728/j];
    t = N[j/2^4*f1^4/f2^4];
    Assert[j == 16*(t^2 + 14*t + 1)^3/(t*(t - 1)^4)];

    f1 = HypergeometricPFQ[{-1/60, 19/60}, {4/5}, 1728/j];
    f2 = HypergeometricPFQ[{11/60, 31/60}, {6/5}, 1728/j];
    t = N[j*f1^5/f2^5];
    Assert[j == (t^4 + 228*t^3 + 494*t^2 - 228*t + 1)^3/(t*(t^2 - 11*t - 1)^5)];

    f1 = HypergeometricPFQ[{-1/42, 13/42, 9/14}, {4/7, 6/7}, 1728/j];
    f2 = HypergeometricPFQ[{5/42, 19/42, 11/14}, {5/7, 8/7}, 1728/j];
    f3 = HypergeometricPFQ[{17/42, 31/42, 15/14}, {9/7, 10/7}, 1728/j];
    t = N[j*f1*f2^2/f3^3];
    Assert[j == (t^2 - t +1)^3*(t^6 + 229*t^5 + 270*t^4 - 1695*t^3 + 1430*t^2 - 235*t + 1)^3/((t^2 - t)*(t^3 - 8*t^2 + 5*t + 1)^7)];

    f1 = HypergeometricPFQ[{-1/10, 1/10, 1/2}, {1/5, 4/5}, -64/j];
    f2 = HypergeometricPFQ[{1/10, 3/10, 7/10}, {2/5, 6/5}, -64/j];
    f3 = HypergeometricPFQ[{7/10, 9/10, 13/10}, {8/5, 9/5}, -64/j];
    t = N[j*f2^2/f3*(j*f1^2 + f2*f3)/(j*f2^3 + f1*f3^2)];
    Assert[j == (t^2 - 1)*(t^2 - 4*t - 1)^5/(t*(t^2 + t - 1)^5)];
, {j, {-20000, 20000}}];

Print["  divergent case"];

tab = {
    {{1, 1}, {}, -1} -> 0.5963473623231940743410784993,
    {{1, 1}, {}, 1} -> 0.69717488323506606876547868 - 1.1557273497909217179100931833*I,
    {{1, 1, 1}, {}, -1} -> 0.6680913263777777654330149875,
    {{1, 1, 1}, {}, 1} -> 0.8513164493201240765836883157 - 0.7156163078376499733488632214*I,
    {{1, 1, 1, 1}, {},-1} -> 0.7237499441815530522203241574,
    {{1, 1, 1, 1}, {}, 1} -> 0.90488007818086251696673303 - 0.515351906643764121580160246*I
};
Do[Assert[N[HypergeometricPFQ@@a[[1]]] == a[[2]]], {a, tab}];

Print["  entire case"];

Do[N[HypergeometricPFQ[{1},{2}, z]] == N[(Exp[z]-1)/z], {z, -9, 9, 2}];

Print["  benchmark in ",First@Timing[
tab = {
  {{1/2,1/3,1/5,1/4},{8/7,1/6,1/9},       1/10,   20, 1.0416551010780920877},
  {{1/2,1/3,1/5,1/4},{8/7,1/6,1/9},   -1+I/300,  200, 0.736146410756869965365418193202119572901652552025166999461296551468894556787397814580414025634098153403217178511117034131271820480437262088020734160211512374994493992189454246946088680872 + I*0.000609114243450659361950770848608680476215098076015633291634716235158853489834358733499703366390313679697687566349542370727179353462509062589882297673226996739782668404393491789738044058793522},
  {{1/2,1/3,1/5,1/4},{8/7,1/6,1/9},    1/300+I,  200, 0.853498211060579954490363364784481268197877903709746131856530610099120791429240016786586006359060701719614321941609976993732551828969874878632105881859561180925930687147324142274043333887472 + I*0.30223571830517708188043830846390180594396734422086540281089516218535003720494168386082577514241218819831483834674824765948731963332253786613021854910934261385527579972720139449948486556402},
  {{1/2,1/3,1/5,1/4},{8/7,1/6,1/9},      11/10,  100, 2.0947383700681176102853448748295509569772539992373820829947310008697960661994302886089470591379582 - I*1.3206342448888487914568920364082049536583031035313367092197531894227907419533993812995218656908316},
  {{1/2,1/3,1/5,1/4},{8/7,1/6,409/420},  11/10,  100, 1.10012992690145570229907705218466524693847390580609878681247630305435529779555060922581931979395195 - I*0.019650631936840220248345681395791981619880186753148199866219643526635725512689324474868785073820769},
  {{1/2,1/3,1/5,1/4},{8/7,1/6,1/9},    21/10+I,  100, 0.7147492149025488915684098194900081106205006718857531504256433965697117195412174890267458603786 + I*0.978124944753727608286669747283128480165412084816829506605452158840356988903421048385043},
  {{1/2,1/3,1/5,1/4},{8/7,1/6,10/9},     21/10,  100, 1.09611960714180653165793929034538750471702782740589204327790995048489804128469673421588286818944088 - I*0.11033813903294061048496297271513399387602730602212159996056871315399179910457530832536941039103546},
  {{1/2,1/3,3/2,3/2},{8/7,1/6,409/420}, -21/10,  100, -0.08111719109576625004949013645065663426},
  {{1/2,1/3,3/2,3/2},{8/7,1/6,409/420}, 21I/10,  100, -0.383693615635771510953247000877434203008702287833914190219349246387863784807276497366211439242383 + I*0.351056068887940720554612223365826032178226364387669038498056038396666040714820367798417390995224},
  {{1/2,1/3,3/2,3/2},{8/7,1/6,7/2+1/4}, 21I/10,  100, 0.49828627736249915445467228767693397760685372220619431985808462963351614294504694016696025797 + I*0.5388483851104630403622090287525048639997504587807914458817597010957463072020568473516818},
  {{1/2,1/3,3/2,3/2},{8/7,1/6,409/420},  21/10,  100, -0.3778888881782393774553127309816076760288807629794011361623259687233181125952890271460921393790923 + I*1.032405719773183086686241272266191251045059248944307109581453130336189134197752486814029528854289},
  {{1/2,1/3,3/2,3/2},{8/7,1/6,7/2+1/4},  21/10,  100, 0.5587476752887922757498928650974538228321146670224016380361419696854589138256157940556936141300072 - I*2.64978703162468283617876316171841523330924854565456964046456128264004644150597368088432880196093},
  {{1/2,1/3,3/2,5/2},{8/7,1/6,7/2+1/4},  21/10,  100, -1.799478585201052861552266948753576293798081804909492183540383822249557260558139389282027739306081 - I*2.055021309529954251848254071722409177646100548975120412982678947474050608126033726808850703354638},
  {{1,1,2,2,3},{3/2,3/2,3/2,3/2},        21/10,  100, -0.515856214128177889536364351387403447925456482215887373490756429690251426189984609386387046002781 - I*0.054708919318247836303446919070054608880399867448242344451898096159820146328422450193688791295906},
  {{-2,1,2,2,3},{4/3,4/3,4/3,4/3},       21/10,  100, 25.72463010204081632653061224489795918367346938775510204081632653061224489795918367346938775510204},
  {{1/2,1/3,1/3,3/2,3/2},{8/7,1/6,7/2+1/5,1/4},         21/10, 100, -0.116857050064488324880005256242593338635544257506370503116752352065017363357264920954229325317 - I*3.5154068199906069488964490621819761650472811211473757611117105966445694384432626973651155973030},
  {{1/2,1/3,1/3,3/2,3/2,4/3},{8/7,1/6,7/2,1/4,1/5},     21/10, 100, -19.94044851545386332011632936129290609826931703330950304355642151716710042503908899870498108255 + I*17.002267047972761110922409197388381021873913435295734935475234691998423671940643618802726394435},
  {{1/2,1/3,1/3,3/2,3/2,4/3},{8/7,1/6,7/2,1/4,5+1/5}, 101/100, 100, 1.2614237787431717747560268320458380373654348418286269662499700441213544811397018495438863637237191417 - I*2.5917416945228937312557772395492111006139098181686942376029167147937383597296691524489068756*^-9},
  Null};
Do[
(*Print[{a[[1]], a[[2]], a[[3]], a[[4]]}];*)
    r = N[HypergeometricPFQ[a[[1]], a[[2]], a[[3]]], a[[4]]];
    Assert[r == a[[5]]]
, {a, Most[tab]}];
]];
