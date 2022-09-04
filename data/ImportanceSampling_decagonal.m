(* ImportanceSampling_decagonal.m
   V.O. UdeM january 2004
*)

(*----------------------------- params -----------------------------*)
lbl = "decagonal";

firstLUTdimDepth = 6;  (* importance index *)
secondLUTdimDepth = 6; (* structural index *)

firstLUTdimSize = Fibonacci[firstLUTdimDepth+2];
secondLUTdimSize = Fibonacci[secondLUTdimDepth+2];
(*******************************************************************
firstLUTdimDepth secondLUTdimDepth   lut
4                       6            8x21
6                       6           21x21
*******************************************************************)

nSubdivisionTiles = 4; (* tile types in [1..nSubdivisionTiles] *)
nSamplingTiles = 2; (* tile types in [nSubdivisionTiles+1..nSubdivisionTiles+nSamplingTiles] *)

(****************** end of params *******************)

(****************** constants *******************)
nFold = 20; (* basic characteristic number *)
phi = tau = 1.618033988749895; (* phi: Knuth, tau: Grunbaum & Shephard; phi == (1+Sqrt[5])/2 *)
tau2 = tau^2;

phi1 = 1/phi;
phi2 = phi^2;
zphi = (1 + Sqrt[5])/4//N; (* intersection[vv[[20]],vv[[16]],zero,vv[[18]]]//euclidlen//FullSimplify *)

hv = Table[rotatedaround[{1,0},zero, (i) 2 Pi/nFold], {i,nFold}]//N;
vv = Table[rotatedaround[{0,1},zero, (i) 2 Pi/nFold], {i,nFold}]//N;


pentagonSizeFactor = 1/phi^3;
pentaCoefs = Table[pentagonSizeFactor,{i,20}];
Do[pentaCoefs[[i]] = zphi pentagonSizeFactor, {i,3,20,4}];
starCoefs = Table[pentagonSizeFactor,{i,20}];
Do[starCoefs[[i]] = zphi pentagonSizeFactor, {i,3,20,4}]; (* pentagons *)

scaletab = Table[phi^(-i+1),{i,100}]//N;

phitab = Table[phi^(-i+3),{i,100}];
fibotab=Table[Fibonacci[100-i],{i,100}];

lut = Table[{0,0},{firstLUTdimSize},{secondLUTdimSize}]; (* 0-init *)
If[firstLUTdimDepth == 4,
  lut={{{116,907},{-94,-144},{-34,-3176},{-144,-2187},{-97,1069},{2846,3981},{-1113,4651},{-3006,3872},{367,8875},{0,0},{0,0},
    {0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},{{-684,-127},{3063,-3111},{-452,-6140},{-251,-3387},
    {-2334,-1309},{4820,625},{-6566,-1778},{-4766,190},{576,9038},{4360,3418},{-206,-50},{937,3954},{-306,1135},{0,0},{0,0},
    {0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},{{-830,-318},{4105,-3750},{-1280,-7124},{-203,-1985},{-2929,-34},
    {4822,-350},{-5960,-3103},{-5213,-2903},{731,8002},{4335,3371},{-457,29},{1492,7984},{-988,3047},{367,138},{224,253},{0,0},
    {0,0},{0,0},{0,0},{0,0},{0,0}},{{-820,-14},{4531,-3856},{-1814,-7053},{533,1239},{-2338,1878},{3743,-924},{-5529,-3770},
    {-5700,-4782},{1361,5736},{3928,2680},{-168,625},{99,9036},{-2282,6750},{1110,606},{676,756},{192,352},{-51,-162},{0,0},
    {0,0},{0,0},{0,0}},{{-1001,-75},{4198,-3668},{-2734,-5900},{2023,3151},{-423,3727},{2142,-169},{-4913,-3479},{-6415,-5450},
    {2029,3703},{2867,1477},{-343,807},{-911,8510},{-3182,7640},{2486,1295},{1434,1616},{436,968},{-180,-166},{-232,205},
    {104,127},{0,0},{0,0}},{{-1297,-238},{2765,-3184},{-2831,-3885},{2663,3705},{1938,3466},{-120,1884},{-5058,-2537},{-6504,-4580},
    {3318,2423},{2379,659},{-238,714},{-1382,6326},{-3161,6368},{2926,1822},{1547,1919},{886,1999},{-392,-24},{-614,711},
    {241,370},{-149,201},{-152,462}},{{-1341,-47},{1847,-2594},{-2427,-1753},{2294,3522},{2393,3008},{-1282,4173},{-4520,-999},
    {-6045,-2691},{3634,2254},{1972,194},{-113,55},{-1066,3112},{-2860,3808},{2735,1290},{1344,1545},{633,1917},{-478,1120},
    {-1210,1726},{113,667},{-455,460},{-305,1038}},{{-914,636},{1123,-2108},{-1501,-401},{1103,2531},{1990,1909},
    {-1266,6167},{-3000,912},{-4962,-916},{3315,2772},{1364,113},{-199,-35},{-536,1271},{-2563,2245},{2242,906},{1179,1351},
    {413,1327},{-377,1978},{-1142,2466},{-440,638},{-1011,1021},{-678,2306}}}/10000.
];
If[firstLUTdimDepth == 6,
  lut={{{155,418},{-252,139},{-8,-3117},{-15,-1263},{-430,122},{922,4061},{-586,4382},{-1391,3924},{210,10795},{0,0},{0,0},
    {0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},{{107,168},{143,-779},{-20,-4261},{-127,-2825},{-657,91},
    {3183,3020},{-2968,3845},{-3483,3091},{203,10073},{499,505},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
    {0,0},{0,0}},{{59,-50},{462,-1352},{-140,-4881},{-214,-3471},{-787,-650},{4231,2348},{-4373,2718},{-4275,2489},{292,9763},
    {1350,1388},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},{{-38,-383},{1130,-2182},
    {-470,-5747},{-431,-4202},{-1093,-1541},{5233,1516},{-5854,587},{-5080,1345},{455,9418},{3145,3248},{-63,8},{0,0},{0,0},{0,0},
    {0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},{{-262,-763},{2160,-2948},{-1122,-6591},{-874,-4499},{-1578,-2258},
    {5666,878},{-6628,-2055},{-5416,-314},{694,9104},{4001,4229},{-225,69},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
    {0,0},{0,0}},{{-463,-985},{3510,-3381},{-2097,-7420},{-1542,-4401},{-2387,-2424},{5792,352},{-6763,-3608},
    {-5081,-2556},{919,8654},{4296,4676},{-567,275},{200,837},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
    {{-664,-977},{4544,-3582},{-3249,-8031},{-2381,-3605},{-3081,-2133},{5745,-158},{-6468,-4341},{-4359,-4634},{1064,7977},
    {3988,4567},{-853,699},{569,2646},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},{{-632,-625},{4980,-3694},
    {-4285,-8300},{-3031,-1992},{-3674,-1249},{5507,-666},{-6042,-4603},{-3627,-6035},{1139,7042},{3624,4302},{-979,1302},
    {1273,6249},{-250,849},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},{{-591,-175},{5024,-3749},{-4865,-8226},
    {-3395,388},{-3893,279},{5011,-937},{-5561,-4741},{-3391,-6684},{1272,5916},{3307,3796},{-854,1863},{1380,8390},{-760,2616},
    {338,-25},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},{{-566,366},{4939,-3652},{-4857,-8026},{-3119,2558},{-3651,2281},
    {4028,-1082},{-5117,-4783},{-3781,-6967},{1573,4769},{3113,3160},{-598,2302},{1095,9323},{-1765,6188},{1049,-29},{0,0},{0,0},
    {0,0},{0,0},{0,0},{0,0},{0,0}},{{-681,677},{4890,-3402},{-4432,-7805},{-2007,3767},{-2979,4140},{2681,-1017},
    {-4853,-4710},{-4773,-7238},{2063,3703},{3015,2515},{-367,2590},{348,9164},{-2228,8206},{2457,-1},{170,52},{0,0},{0,0},{0,0},
    {0,0},{0,0},{0,0}},{{-775,670},{4763,-3165},{-4029,-7402},{-346,4158},{-1975,4773},{1417,-600},{-4901,-4524},{-6123,-7589},
    {2657,2798},{2967,1920},{-232,2870},{-431,8848},{-2322,9198},{3216,245},{523,241},{214,127},{0,0},{0,0},{0,0},{0,0},
    {0,0}},{{-904,376},{4283,-3092},{-3806,-6743},{1241,4059},{-557,4380},{184,130},{-5321,-4178},{-7165,-7732},{3322,2122},
    {2899,1257},{-149,3183},{-1099,8523},{-2136,8868},{3477,584},{1189,647},{627,373},{0,0},{0,0},{0,0},{0,0},{0,0}},
    {{-965,-135},{3418,-3148},{-3564,-5872},{2270,3908},{1076,3337},{-1145,1095},{-5954,-3645},{-7592,-7290},{4021,1575},{2962,660},
    {-1,3359},{-1568,8092},{-2068,7992},{3248,958},{1462,1148},{1414,874},{-104,138},{0,0},{0,0},{0,0},{0,0}},
    {{-1088,-729},{2518,-2933},{-3381,-4910},{2891,3817},{2148,2580},{-2598,2159},{-6373,-3031},{-7561,-6215},{4667,1094},{3107,219},
    {246,3206},{-1902,7202},{-2174,6607},{2879,1189},{1398,1460},{1659,1123},{-302,461},{-95,233},{0,0},{0,0},{0,0}},
    {{-1104,-1118},{1775,-2528},{-3218,-3980},{3224,3897},{2683,2270},{-3720,3015},{-6578,-2425},{-7460,-5018},{5067,676},{3209,-61},
    {396,2524},{-2107,5667},{-2071,5042},{2448,1169},{1107,1532},{1573,1272},{-690,1095},{-304,696},{95,120},{0,0},{0,0}},
    {{-1148,-919},{1221,-1997},{-3025,-2980},{3258,3935},{2818,2473},{-4012,3680},{-6555,-1877},{-7465,-4140},{5141,414},{3193,-227},
    {449,1665},{-2205,3579},{-1838,3421},{2021,933},{838,1424},{1242,1300},{-817,1563},{-730,1639},{212,363},{-186,278},{0,0}},
    {{-923,-138},{869,-1536},{-2624,-2019},{2879,3663},{2661,2741},{-3535,4336},{-6118,-1204},{-7454,-3451},{4883,402},{2987,-293},
    {294,634},{-1948,1930},{-1633,2093},{1698,721},{680,1296},{984,1288},{-799,1830},{-1014,2164},{372,833},{-519,817},{0,0}},
    {{-624,550},{648,-1249},{-1878,-1077},{2054,2901},{2170,2514},{-2421,5208},{-5217,-332},{-7114,-2667},{4378,650},{2609,-290},
    {107,-48},{-1442,889},{-1688,1389},{1500,683},{659,1247},{822,1274},{-617,2029},{-1143,2536},{63,1039},{-1158,1862},
    {-112,327}},{{-385,910},{497,-1281},{-1069,-493},{1189,2124},{1289,1868},{-1501,5888},{-3524,1138},{-5931,-1164},{3463,1666},
    {2003,-269},{-1,-387},{-873,454},{-1646,954},{1477,772},{801,1205},{775,1211},{-414,2182},{-1035,2656},{-528,1023},
    {-1252,2162},{-313,757}},{{-158,976},{367,-1440},{-520,-270},{585,1582},{714,1489},{-1007,6257},{-2159,2331},{-4290,506},
    {2567,2790},{1474,-262},{44,-471},{-480,313},{-1591,705},{1553,905},{974,1143},{700,1089},{-256,2268},{-811,2672},
    {-1113,876},{-1186,2139},{-695,1682}}}/10000.
];

tileShapes = 
   {{zero,vv[[15]], phi vv[[19]]},        (* 1 *)
    {zero,vv[[15]], phi vv[[19]]},        (* 2 *)
    {zero,phi^2 vv[[15]], phi vv[[17]]},  (* 3 *)
    {zero,phi^2 vv[[15]], phi vv[[17]]},  (* 4 *)
    {zero},                               (* 5 *)
    {zero}                                (* 6 *)
}

tileFillCol={
  RGBColor[.8, 1., 1.],  (* Light Bluish *)
  RGBColor[.0, 1., 1.],  (* Cyan *)
  RGBColor[1., 1, .5],  (* Yellowish *)
  RGBColor[1., .7, 0.],  (* Yellow *)
  RGBColor[1., .6, .6],  (* Reddish *)
  RGBColor[.5, 1, .5]  (* Greenish *)
};

(****************** end of constants *******************)

(****************** decagonal-specific Procedures *******************)

fcodeValue[fcode_]:=Plus@@(fcode Take[fibotab,{-Length[fcode]-2,-3}]) (* fibotab must be known *)

getNofBits[decimal_]:=Max[1,Floor[Log[phi,(decimal  sqrt5)+1]]-1]
getNofSubdivisions[decimal_]:=Ceiling[(getNofBits[decimal])/2.]

getLUTIndex[decimal_]:= 1+Floor[firstLUTdimSize*FractionalPart[Log[phi2,(decimal sqrt5)+1]]] 
  (* firstLUTdimSize comprises number of subdivisions of one octave *)

getTileShape[fig_]:=Module[{tileType,tileRefPt,tileVLst,tileDir,LOS,fcode,res,sz,v1,v2,v3},
  {tileType,tileRefPt,tileVLst,tileDir,LOS,fcode} = fig;
  If[tileType <= nSubdivisionTiles, {v1,v2,v3} = tileVLst];
  sz = scaletab[[LOS]];
  If[shapeTrianlesOnly,
    If[tileType <= nSubdivisionTiles,
      Return[{v1,v2,v3,v1}]
    ,(*ELSE*)
      Return[{}]
    ]
  ,(*ELSE*)
    res = Switch[tileType
    ,1,
       Join[
         Plus[v3,#]& /@ (sz Table[rotatedaround[{pentaCoefs[[i]],0},{0,0}, (i+tileDir+1) Pi/10], {i,13,15,2}]),
         Plus[v2,#]& /@ (sz Table[rotatedaround[{pentaCoefs[[i]],0},{0,0}, (i+tileDir+3) Pi/10], {i,3,7,2}]),
         Plus[v1,#]& /@ (sz Table[rotatedaround[{starCoefs[[i]],0},{0,0}, (i+tileDir+13) Pi/10], {i,7,11,2}])]
    ,2,
       Join[
         Plus[v3,#]& /@ (sz Table[rotatedaround[{pentaCoefs[[i]],0},{0,0}, (i+tileDir+3) Pi/10], {i,11,13,2}]),
         Plus[v2,#]& /@ (sz Table[rotatedaround[{starCoefs[[i]],0},{0,0}, (i+tileDir-1) Pi/10], {i,7,11,2}]),
         Plus[v1,#]& /@ (sz Table[rotatedaround[{pentaCoefs[[i]],0},{0,0}, (i+tileDir-3) Pi/10], {i,3,7,2}])]
    ,3,
       Join[
         Plus[v3,#]& /@ (sz Table[rotatedaround[{pentaCoefs[[i]],0},{0,0}, (i+tileDir+1) Pi/10], {i,11,17,2}]),
         Plus[v2,#]& /@ (sz Table[rotatedaround[{starCoefs[[i]],0},{0,0}, (i+tileDir+1) Pi/10], {i,7,9,2}]),
         Plus[v1,#]& /@ (sz Table[rotatedaround[{pentaCoefs[[i]],0},{0,0}, (i+tileDir-1) Pi/10], {i,1,3,2}])]
    ,4,
       Join[
         Plus[v3,#]& /@ (sz Table[rotatedaround[{pentaCoefs[[i]],0},{0,0}, (i+tileDir-1) Pi/10], {i,13,20,2}]),
         Plus[v2,#]& /@ (sz Table[rotatedaround[{pentaCoefs[[i]],0},{0,0}, (i+tileDir+1) Pi/10], {i,7,9,2}]),
         Plus[v1,#]& /@ (sz Table[rotatedaround[{starCoefs[[i]],0},{0,0}, (i+tileDir-1) Pi/10], {i,1,3,2}])]
    ,5, Plus[tileRefPt,#]& /@ (sz Table[rotatedaround[{pentaCoefs[[i]],0},{0,0}, (i+tileDir) Pi/10], {i,1,20,2}])
    ,6, Plus[tileRefPt,#]& /@ (sz Table[rotatedaround[{starCoefs[[i]],0},{0,0}, (i+tileDir) Pi/10], {i,1,20,2}])
    ];
    Return[Append[res,First[res]]]
  ]
] (* getTileShape *)

recursiveSubdiv[fig_]:=Module[{tileType,tileRefPt,tileVLst,tileDir,LOS,i,v1,v2,v3,
fcode,restypeseq,fig1,res={},scheme,types,midpt1,midpt2,r,curLOS},
  {tileType,tileRefPt,tileVLst,tileDir,LOS,fcode} = fig;
  If[tileType <= nSubdivisionTiles, {v1,v2,v3} = tileVLst];
  If[cutOutOfRangeFlag && !figInRange[fig,scaletab[[LOS]]], Return[{}]];
  If[adaptiveFlag,
    If[tileType <= nSubdivisionTiles, (* generator *)
      If[LOS > getMaxLOSinTile[{v1,v2,v3}], Return[{fig}] ] (* getMaxLOSinTile or getMaxLOSinTrianle *)
    ,(*ELSE: sampler *)
      If[LOS > getLOS[tileRefPt], Return[{fig}] ]
    ];
  ,(*ELSE: non-adaptive*)
    If[LOS > maxLOS, Return[{fig}]]
  ];
  Which[
  tileType==1,
    midpt1 = (v2 + phi1 v3) / phi;
    res = Join[res, recursiveSubdiv[{4,v3,{v3,v1,midpt1},modN[tileDir+14],LOS+1,Join[{0,0},fcode]}] ];
    res = Join[res, recursiveSubdiv[{1,v2,{v2,midpt1,v1},modN[tileDir+6],LOS+1,Join[{1,0},fcode]}] ];
    res = Join[res, recursiveSubdiv[{5,midpt1,{midpt1},modN[tileDir+1],LOS+1,Join[{1,0},fcode]}] ];
  ,tileType==2,
    midpt1 = (v1 + phi1 v3) / phi;
    res = Join[res, recursiveSubdiv[{3,v2,{v2,v3,midpt1},modN[tileDir+6],LOS+1,Join[{0,0},fcode]}] ];
    res = Join[res, recursiveSubdiv[{2,midpt1,{midpt1,v1,v2},modN[tileDir+14],LOS+1,Join[{1,0},fcode]}] ];
  ,tileType==3,
    midpt1 = (v2 + phi1 v1) / phi;
    midpt2 = (phi1 v1 + v3) / phi;
    res = Join[res, recursiveSubdiv[{4,v1,{v1,midpt1,midpt2},modN[tileDir+0],LOS+1,Join[{0,0},fcode]}] ];
    res = Join[res, recursiveSubdiv[{1,v3,{v3,midpt2,midpt1},modN[tileDir+12],LOS+1,Join[{1,0},fcode]}] ];
    res = Join[res, recursiveSubdiv[{3,v2,{v2,v3,midpt1},modN[tileDir+8],LOS+1,Join[{0,1},fcode]}] ];
    res = Join[res, recursiveSubdiv[{5,midpt2,{midpt2},modN[tileDir+7],LOS+1,Join[{1,0},fcode]}] ];
  ,tileType==4,
    midpt1 = (v1 + phi1 v2) / phi;
    midpt2 = (phi1 v2 + v3) / phi;
    res = Join[res, recursiveSubdiv[{3,midpt1,{midpt1,v2,midpt2},modN[tileDir+0],LOS+1,Join[{0,0},fcode]}] ];
    res = Join[res, recursiveSubdiv[{2,midpt2,{midpt2,v3,midpt1},modN[tileDir+8],LOS+1,Join[{1,0},fcode]}] ];
    res = Join[res, recursiveSubdiv[{4,v3,{v3,v1,midpt1},modN[tileDir+12],LOS+1,Join[{0,1},fcode]}] ];
    res = Join[res, recursiveSubdiv[{5,midpt1,{midpt1},modN[tileDir+15],LOS+1,Join[{0,1},fcode]}] ];
  ,tileType==5,
    res = Join[res, recursiveSubdiv[{6,tileRefPt,{tileRefPt},modN[tileDir],LOS+1,Join[{0,0},fcode]}] ];
  ,tileType==6,
    res = Join[res, recursiveSubdiv[{5,tileRefPt,{tileRefPt},modN[tileDir+10],LOS+1,Join[{0,0},fcode]}] ];
  ];
  Return[res]
] (* recursiveSubdiv *)

getFLstSeed[dz_]:=
  {getFigure[3,{-1.30902,0}+dz,20,1,{}],
   getFigure[4,{ 1.30902,0}+dz,10,1,{}]}
(****************** end of decagonal-specific Procedures *******************)

(****************** init *******************)
f0:=getFLstSeed[{0,0}];
If[unknown[outputAspectRatio], outputAspectRatio = 1];
If[outputAspectRatio == 1,
  {ysize,xsize} = {256,256};
  area34 = 1.24495;  (*triangles 3,4 Det[{tileShapes[[3,2]]-tileShapes[[3,1]],tileShapes[[3,3]]-tileShapes[[3,1]]}]/2. *)
  workingArea  = area34/phi^2;
  workingSqSize = Sqrt[workingArea];
  {{xmin, xmax},{ymin,ymax}}={{-workingSqSize/2,workingSqSize/2},{.005,.005+workingSqSize}}; 
];

tileShapesTh = Thickness[.0002]; tileShapesCol = Red; 
samplingPtSize = PointSize[.0075]; samplingPtCol = Black;

dir = nFold; iter = 0; threshold = 1;
arrows=txtgl = {}; basicgl={};
sgl = {samplingPtSize,samplingPtCol}; (* sampling points *)
shapesgl =  {tileShapesTh,tileShapesCol}; (* tile shapes *)
frame={Gray,Thickness[.005],Line[{{xmin, ymin},{xmin, ymax},{xmax, ymax},{xmax, ymin},{xmin, ymin}}]};
rng = {{xmin, xmax},{ymin, ymax}};
font = "Courier"; textsz = 10;
adaptiveFlag=False;
maxLOS = 5;
dbg = False;
SetOptions[Graphics, ImageSize -> {800, Automatic}];
  showTileType = False;
  showFcode = True;
  showShortFcode = False;
  fillFigure = True;
  showSamplingPt = False;
  showdir = True;
  showTileShapes = True;
  showTileShapesWithDashing = False;
  shapeTrianlesOnly = False;
  showiSimplifiedTileShapes = False;

cutOutOfRangeFlag = True;

 (****************** end of init *******************)

 (****************** init late version *******************)

Print["ImportanceSampling_decagonal.m loaded."];
