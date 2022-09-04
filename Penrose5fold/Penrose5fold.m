(* Penrose5fold.m
   V.O. based on version 2002/08/28
   2022/09/04
   
   showGrowthAmmannBars[]
   dbgInflationsAmmannBars[]
   
   visuPenroseNara[]
*)
 
(****************** parameters *******************)

SetDirectory[ToFileName[$HomeDirectory,"Penrose5fold/"]];
SetOptions[Graphics, ImageSize -> {Automatic, 3/2 1024}];

(*Get["data/ImportanceSampling_common.m"];
Get["data/ImportanceSampling_decagonal.m"];
*)
showdir = False;
showFcode = False;
showShortFcode = False;
cutOutOfRangeFlag = False;
rng = All;
font = "Courier-Bold"; textsz = 18;

shapeTrianlesOnly = False;
shapeTrianlesOnly = True;

showAmmannBars = True;
{ammanBarCol,ammanBarTh} = {Black,AbsoluteThickness[1]};

tileShapesTh = Thickness[0];
tileShapesCol = Yellow;

PI = Pi//N;
(****************** end of parameters *******************)


showGrowthAmmannBars[inniters_:8, dbg_:False] :=
    Module[ {},
    	showTileType = False;
    	showFcode = False;
    	showSamplingPt = True;

        niters = inniters;
        prevsgl = {};
        Do[
              maxLOS = iter;
              flst = Flatten[recursiveSubdiv /@ f0, 1];
              txtgl = shapesgl = sgl = {};
              basicgl = {getGLst[flst]};
              gl1 = {shapesgl,basicgl,txtgl};
              g = Graphics[{gl1,AbsolutePointSize[10],sgl,Red,prevsgl}, ImageSize -> {Automatic, 3/2 1024}];
              g//Print;
              
              prevsgl = sgl;
              If[dbg, Export["g_AmmannBars_iter"<>ToString[iter]<>".pdf",g]];
         ,{iter,0,niters}];
    ]


dbgInflationsAmmannBars[inniters_:4, dbg_:False] :=
    Module[ {},
        niters = inniters;
        {px,py} = {2.8,-1.2};
        n = 2;
        Do[
          maxLOS = iter;
          plbl = ToString[lbl]<>" iter="<>ToString[iter]<>" "<>date;
          shapesgl = gl = {};
          Do[
            type = i;
            {x,y} = {px,py} {(i-1) - n Floor[(i-1)/n], Floor[(i-1)/n]};
            fig = getFigure[type,{x,y},dir,1,{}];
            flst = recursiveSubdiv[fig];
            AppendTo[gl, getGLst[flst]];
          ,{i, Length[tileShapes] }];
          g = Graphics[{PointSize[.003],gl,shapesgl}];
          g//Print;
          If[dbg, Export["g_dbgInflations_AmmannBars_iter"<>ToString[iter]<>".pdf",g] ];
        ,{iter,0,niters}];
    ] (* dbgInflations *)
	
(*
(****************** below: visuPenroseNara *******************)
(****************** constants *******************)
tau = GoldenRatio//N;
tau2 = tau^2;

r23Dummy = 0;
r14Dummy = -100;
r23U = 1;
r23D = 2;
r23L = 3;
r14U = -3;
r14D = -4;

tileset = {
  r23U,
  r23D,
  r23L,
  r14U,
  r14D
}


zero = {0,0};
uvect = Table[rotatedaround[{0,1},zero,(i-1) PI/5]//N//Chop, {i,10}];
uvectrot = RotateRight[uvect,2];
hueTable={3/6., 5/6., 4/3., 0/3., 2/3., 1/6.}; (* C M G R B Y *)

col0 = Yellow;
col1 = Red;
col2 = Black;
col3 = Cyan;

col02 = Red;
col13 = Blue;

colA = Yellow;
colB = Cyan;

  {fibonplus1,fibon} = {1,1}; (* fibo[n+1],fibo[n] *)
  fibolbl = ToString[fibonplus1]<>"_"<>ToString[fibon];
(****************** end of constants *******************)

(**************** System-dependent setup ******************)
If[$System == "Microsoft Windows",
  SetDirectory["C:\\0_tilings\\optim"];
  (*epsOut:=Identity[#1]&;*)
,(*ELSE: Linux, Unix *)
  AppendTo[$Path,$HomeDirectory<>"/math"];
];
pdfOut:=Export[#1,#2,"PDF"]&; (* pdfOut[fname,graph] *)
pictOut:=Export[#1,#2,"PICT"]&; (* pictOut[fname,graph] *)
tiffOut:=Export[#1,#2,"TIFF",ImageSize->{512,512}]&; (* tiffOut[fname,graph] *)
date = ToString[Date[][[1]]]<>"/"<>ToString[Date[][[2]]]<>"/"<>
        ToString[Date[][[3]]]<>"_"<>ToString[Date[][[4]]]<>":"<>
        ToString[Date[][[5]]];
(**************** end of System-dependent setup ******************)


(****************** procedure *******************)
mf := MatrixForm

mod10[x_]:=Mod[x-1,10]+1

euclidlen[z_] :=
    Sqrt[z[[1]]^2 + z[[2]]^2];

rotatedaround[pt1_, pt2_, angle_] :=
 (* gives point which is pt1 rotated around pt2 by angle *)
    Block[ {res, a},
        If[ pt1 == pt2,
            Return[pt1]
        ];
        a = ArcTan[(pt1 - pt2)[[1]], (pt1 - pt2)[[2]]];
        res = pt2 + euclidlen[pt1 - pt2] {Cos[a + angle], Sin[a + angle]};
        Return[res]
    ]; (* rotatedaround *)

rotatedaroundandscaled[pt1_, pt2_, angle_, k_] :=
 (* pt1 rotated around pt2 by angle then vector {pt2,pt1} is scaled k \
times *)
    Block[ {res, a},
        If[ pt1 == pt2,
            Return[pt1]
        ];
        a = ArcTan[(pt1 - pt2)[[1]], (pt1 - pt2)[[2]]];
        res = pt2 + k euclidlen[pt1 - pt2] {Cos[a + angle], Sin[a + angle]};
        Return[res]
    ]; (* rotatedaround *)

arrow[from_, to_] :=
    Module[ {v1, v2, len, ka1 = .07, ka2 = .025, ka3 = .06, ka},
        len = euclidlen[to - from] // N;
        v1 = (to - from)/len;
        v2 = {-v1[[2]], v1[[1]]};
        ka = 1/2;
        {Line[{to - .05*ka*v1, from}], 
        Polygon[{to, to - ka1 ka v1 - ka2 ka v2, to - ka3 ka v1, 
        to - ka1 ka v1 + ka2 ka v2}]}
    ]

arr[{from_,to_}]:= (*** arrow, compatible with Line[{pt1,pt2}] ***)
Block[{v1,v2,len,ka1=.07,ka2=.025,ka3=.06,ka},
    len = euclidlen[to - from]//N;
    v1 = (to - from) / len;
    v2 = {-v1[[2]], v1[[1]]};
    ka = 5 curmag;
    {Line[{to-.05*ka*v1,from}],
       Polygon[{to,to-ka1 ka v1 - ka2 ka v2, to-ka3 ka v1,
                to -ka1 ka v1 + ka2 ka v2}]}
]; (* arr *)
invArr[{to_,from_}]:=arr[{from,to}]

getFigure[{type_,z0_,dir_,mag_,extra_}]:=Block[{z2,z3},
  Return[{type, z0, dir, mag, extra}]
] (* getFigure *)


getFigureVertices[{type_,z0_,dir_,mag_,extra_}]:=Block[{z1,z2,z3},
  If[type >= 0, (* tThickXXX *)
    z1 = z0 + mag uvect[[mod10[dir-1]]];
    z3 = z0 + mag uvect[[mod10[dir+1]]];
  ,(*ELSE*)
    z1 = z0 + mag uvect[[mod10[dir-2]]];
    z3 = z0 + mag uvect[[mod10[dir+2]]];
  ];
  z2 = z1 + (z3-z0);
  Return[{z0,z1,z2,z3}]
] (* getFigureVertices *)


petal1[z0_:{0,0},scale_:1,dir_:1,rand_:True] :=
    Module[ {controlpts,npts = 40,f,w = {2.5,0},h = {0,3},zero = {0,0}, randfactor =.5 },
        controlpts = Plus[#,-(h/2)]& /@ { zero,.5w,.25w+.2h,.25w+.8h,.5w+h,-.5w+h,-.25w+.8h,-.25w+.2h,-.5w,zero};
        controlpts = rotatedaround[#,zero,(dir-1)2 PI/40]&  /@ controlpts;
        If[ rand, Do[ controlpts[[i]] += {RandomReal[]-.5,RandomReal[]-.5} randfactor ,{i,2,Length[controlpts]-1}] ];
        f = BSplineFunction[controlpts];
        Table[z0 + scale f[i/npts],{i,0,npts}]
    ]
 petal2[z0_:{0,0},scale_:1,dir_:1,rand_:True] :=
    Module[ {controlpts,npts = 40,f,w = {1,0},h = {0,.8},zero = {0,0}, randfactor =.1 },
        controlpts = Plus[#,-(h/2)]& /@ { zero,.25w,.7w-.4h,.8w-.3h,w-.2h,.25w+.3h,.3w+.8h,.2w+.9h,h,-.2w+.9h,-.3w+.8h,-.25w+.3h, -w-.2h,-.8w-.3h,-.7w-.4h,-.25w, zero};
        controlpts = rotatedaround[#,zero,(dir-1)2 PI/40]&  /@ controlpts;
        If[ rand, Do[ controlpts[[i]] += {RandomReal[]-.5,RandomReal[]-.5} randfactor ,{i,2,Length[controlpts]-1}] ];
        f = BSplineFunction[controlpts];
        Table[z0 + scale f[i/npts],{i,0,npts}]
    ]
 

petal3[z0_:{0,0},scale_:1,dir_:1,rand_:True] :=
    Module[ {controlpts,npts = 50,f,w = {1,0},h = {0,2},zero = {0,0}, randfactor =.2 },
        controlpts = Plus[#,-(h/2)]& /@ { zero,.3w+.001h,.8w+.25h,w+.5h,.7w+.9h,h,-.2w+h,-.4w+.95h,.4w+.65h,-.2w+.5h,.3w+.35h,.4w+.2h,-.35w+.15h,-.3w+.05h,-.2w,-.1w,zero};
        controlpts = rotatedaround[#,zero,(dir-1)2 PI/40]&  /@ controlpts;
        If[ rand, Do[ controlpts[[i]] += {RandomReal[]-.5,RandomReal[]-.5} randfactor ,{i,2,Length[controlpts]-1}] ];
        f = BSplineFunction[controlpts];
        Table[z0 + scale f[i/npts],{i,0,npts}]
    ]

getGLst[flst_]:=Map[getGFigure, flst, {1}]

getGFigure[fig_]:=Block[{gl={},type, z0,z1,z2,z3, dir, mag, extra,subdivLevel,symb,z01,z12,z23,z30,z02,thval,z,x1,x3,x0,r},
  If[dbgLst, Print["getGFigure ",++count,"/",lstlen," iter=",iter]];
  {type, z0, dir, mag, extra} = fig;
  {subdivLevel,symb,thval} = extra;
  {z0,z1,z2,z3} = getFigureVertices[fig];
  If[showOrigin,
      AppendTo[gl,Point[(10z0+z2)/11] ]
  ]; (* showOrigin *)
  If[showThValues,
    If[symbolicForm,
      AppendTo[marks,Text[ToString[thval//InputForm],getPivot[fig],{-1,0}] ];
    ,(*ELSE*)
      AppendTo[marks,Text[ToString[N[thval,4],TotalWidth->11],getPivot[fig],{-1,0}] ];
    ];
  ]; (* showThValues *)
  If[showThDisks,
    AppendTo[marks,{GrayLevel[thval/.{t:>tau}],Point[getPivot[fig]]} ];
  ]; (* showThDisks *)
  If[showP1marks,
    If[type == r23U || type == r23D || type == r23L,
      z01=(tau z0 + z1)/(tau+1);
      z12=( z2 + z1)/(1+1);
      z23=( z3 + z2)/(1+1);
      z30=(tau z0 + z3)/(tau+1);
      z02= (tau z2 + z0)/(tau+1);
      x1 = z1 - 1/tau (z2-z1);
      x3 = z3 - 1/tau (z2-z3);
      r = euclidlen[z1-z3];
      AppendTo[marks,{p1marksth,LightGray,Polygon[(*{z0,z01,z02,z30,z0}*)Join[{z0},Plus[#,x3]& /@ (r RotateLeft[ngon,(dir-1)*4][[1;;1+4]]),Plus[#,x1]& /@ (r RotateLeft[ngon,(dir-1)*4][[1+20-4;;1+20]])]]
        ,Red,(*Line[{z12,z02,z23}],*)Line[Plus[#,x3]& /@ (r RotateLeft[ngon,(dir-1)*4][[1;;1+6]])],Line[Plus[#,x1]& /@ (r RotateLeft[ngon,(dir-1)*4][[1+20-6;;1+20]])]
        } ];
    ];
    If[type == r14U || type == r14D,
      z01=(tau z0 + z1)/(tau+1);
      z12=( z2 + z1)/(1+1);
      z23=( z3 + z2)/(1+1);
      z30=(tau z0 + z3)/(tau+1);
      z02= (tau z2 + z0)/(tau+1);
      x0 = z0 - tau (z2-z0);
      x1 = z2 + 1/tau (z2-z1);
      x3 = z2 + 1/tau (z2-z3);
      r = euclidlen[x0-z01];
      AppendTo[marks,{p1marksth,LightGray,Polygon[Join[{z0},Plus[#,x0]& /@ (r RotateLeft[ngon,(dir-1)*4][[1+8;;1+12]]) ]]
		,Red(*,Line[{z01,z12}],Line[{z23,z30}]*)
			,Line[Join[Plus[#,x0]& /@ (r RotateLeft[ngon,(dir-1)*4][[1+8;;1+12]]) ]]
			,Line[Join[Plus[#,x1]& /@ (r RotateLeft[ngon,(dir-1)*4][[37;;39]]) ]]
			,Line[Join[Plus[#,x3]& /@ (r RotateLeft[ngon,(dir-1)*4][[23;;25]]) ]]
		} ];
    ];
  ]; (* If[showP1marks, *)
  If[showNaramarks,
    If[type == r23U || type == r23D || type == r23L,
      z01=(tau z0 + z1)/(tau+1);
      z12=( z2 + z1)/(1+1);
      z23=( z3 + z2)/(1+1);
      z30=(tau z0 + z3)/(tau+1);
      z02= (tau z2 + z0)/(tau+1);
      x1 = z1 - 1/tau (z2-z1);
      x3 = z3 - 1/tau (z2-z3);
      r = euclidlen[z0-z1]/2;
      rvert = r/10;
      AppendTo[marks,{Green,Thickness[.005]
         ,Line[Plus[#,z0]& /@ (r RotateLeft[ngon,(dir-1)*4][[7;;15]])]
         ,Line[Plus[#,z1]& /@ (r RotateLeft[ngon,(dir-1)*4][[15;;27]])]
         ,Line[Plus[#,z2]& /@ (r RotateLeft[ngon,(dir-1)*4][[20+7;;20+15]])]
         ,Line[Plus[#,z3]& /@ (r RotateLeft[ngon,(dir-1)*4+20][[15;;27]])]
         ,Black,Polygon[Join[{z0},Plus[#,z0]& /@ (rvert RotateLeft[ngon,(dir-1)*4][[7;;15]]) ]]
         , Table[refpt =  (z0 + .9 r RotateLeft[ngonPlusHalf,(dir-1)*4][[i]] ); Polygon[petal1[refpt,rvert/2.5,(dir-1)*4 + i-9.5] ] , {i,7,14}] 
         , refpt =  (z0 + .2 (z2-z0) ); Polygon[petal2[refpt,rvert,(dir-1)*4 ] ]
        
         
        } ];
    ];
    If[type == r14U || type == r14D,
      z01=(tau z0 + z1)/(tau+1);
      z12=( z2 + z1)/(1+1);
      z23=( z3 + z2)/(1+1);
      z30=(tau z0 + z3)/(tau+1);
      z02= (tau z2 + z0)/(tau+1);
      x0 = z0 - tau (z2-z0);
      x1 = z2 + 1/tau (z2-z1);
      x3 = z2 + 1/tau (z2-z3);
      r = euclidlen[x0-z01];
      AppendTo[marks,{p1marksth,LightGray,Polygon[Join[{z0},Plus[#,x0]& /@ (r RotateLeft[ngon,(dir-1)*4][[1+8;;1+12]]) ]]
		,Red(*,Line[{z01,z12}],Line[{z23,z30}]*)
			,Line[Join[Plus[#,x0]& /@ (r RotateLeft[ngon,(dir-1)*4][[1+8;;1+12]]) ]]
			,Line[Join[Plus[#,x1]& /@ (r RotateLeft[ngon,(dir-1)*4][[37;;39]]) ]]
			,Line[Join[Plus[#,x3]& /@ (r RotateLeft[ngon,(dir-1)*4][[23;;25]]) ]]
		} ];
    ];
  ]; (* If[showNaramarks, *)
  
  
  If[showSepArrows,
    Switch[type
    ,r23U,
      AppendTo[arrgl,{If[EvenQ[dir],col02,col13],arr[{z0,z0 + tau mag uvect[[dir]] }]} ];
    ,r23D,
      AppendTo[arrgl,{If[EvenQ[dir],col02,col13],arr[{z0 + tau mag uvect[[dir]],z0 }]} ];
    ,r23L,
      AppendTo[arrgl,{If[OddQ[dir],col02,col13],
        Line[{z0+mag uvect[[mod10[dir-1]]],z0+mag/tau uvect[[mod10[dir+0]]]}],
        arr[{z0+mag/tau uvect[[mod10[dir-0]]],z0+mag uvect[[mod10[dir+1]]]}]}];
    ,r14U,
      AppendTo[arrgl,{If[OddQ[dir],col02,col13],arr[{z0,z0 + 1/tau mag uvect[[dir]] }]} ];
    ,r14D,
      AppendTo[arrgl,{If[OddQ[dir],col02,col13],arr[{z0 + 1/tau mag uvect[[dir]],z0 }]} ];
    ]; (*switch*)
  ]; (* showSepArrows *)
  If[showSFC,
    Switch[type
    ,r23U,
      AppendTo[sfcgl,{ Line[{(z0+z1)/2.,(z1+z2)/2.}],Line[{(z2+z3)/2.,(z3+z0)/2.}] } ];
    ,r23D,
      AppendTo[sfcgl,{ Line[{(z0+z1)/2.,(z1+z2)/2.}],Line[{(z2+z3)/2.,(z3+z0)/2.}] } ];
    ,r23L,
      AppendTo[sfcgl,{ Line[{(z0+z1)/2.,(z0+z3)/2.}],Line[{(z2+z3)/2.,(z1+z2)/2.}] } ];
    ,r14U,
       AppendTo[sfcgl,{ Line[{(z0+z1)/2.,(z1+z2)/2.}],Line[{(z2+z3)/2.,(z3+z0)/2.}] } ];
    ,r14D,
       AppendTo[sfcgl,{ Line[{(z0+z1)/2.,(z1+z2)/2.}],Line[{(z2+z3)/2.,(z3+z0)/2.}] } ];
    ]; (*switch*)
  ]; (* showSepArrows *)
  If[showEdgeArrows,
    ar:=Line;
    Switch[type
    ,r23U,
      AppendTo[arrgl,{colA,ar[{z0,z1}],ar[{z0,z3}],colB,ar[{z1,z2}],ar[{z3,z2}]}];
    ,r23D,
      AppendTo[arrgl,{colA,ar[{z0,z1}],ar[{z0,z3}],colB,ar[{z1,z2}],ar[{z3,z2}]}];
    ,r14U,
      AppendTo[arrgl,{colA,ar[{z0,z1}],ar[{z0,z3}],colB,ar[{z2,z1}],ar[{z2,z3}]}];
    ,r14D,
      AppendTo[arrgl,{colA,ar[{z0,z1}],ar[{z0,z3}],colB,ar[{z2,z1}],ar[{z2,z3}]}];
    ,r23L,
      AppendTo[arrgl,{colA,ar[{z0,z1}],ar[{z0,z3}],colB,ar[{z1,z2}],ar[{z3,z2}]}];
    ]; (*switch*)
  ]; (* showEdgeArrows *)
 
  If[borderFlag,
    AppendTo[bordergl,Line[{z0,z1,z2,z3,z0}] ];
  ]; (* If[borderFlag *)
  If[showVoronoi,
    If[type==r23U || type==r23D || type==r23L,
      z002 = (tau z0 + z2)/tau2;
      z022 = (z0 + tau z2)/tau2;
      AppendTo[voronoigl,Line[{(z0+z1)/2,z002,z022,z002,(z0+z3)/2}] ];
      AppendTo[voronoigl,Line[{(z2+z1)/2,z022,(z2+z3)/2}] ];
    ];
    If[type==r14U || type==r14D,
      z113 = (tau2 z1 + z3)/(1+tau2);
      z133 = (z1 + tau2 z3)/(1+tau2);
      AppendTo[voronoigl,Line[{(z0+z1)/2,z113,z133,z113,(z2+z1)/2}] ];
      AppendTo[voronoigl,Line[{(z0+z3)/2,z133,(z2+z3)/2}] ];
    ];
  ]; (* If[borderFlag *)
  If[showDeBruijnIndices,
    If[EvenQ[dir],
      Switch[type
      ,r23U,AppendTo[txtgl,{Text["2",z2]}]
      ,r23D,AppendTo[txtgl,{Text["0",z0]}]
      ,r23L,AppendTo[txtgl,{Text["1",z3]}]
      ,r14U,AppendTo[txtgl,{Text["1",z2]}]
      ,r14D,AppendTo[txtgl,{Text["3",z0]}]
     ]; (*switch*)
    ,(*ELSE*)
      Switch[type
      ,r23U,AppendTo[txtgl,{Text["1",z2]}]
      ,r23D,AppendTo[txtgl,{Text["3",z0]}]
      ,r23L,AppendTo[txtgl,{Text["2",z3]}]
      ,r14U,AppendTo[txtgl,{Text["2",z2]}]
      ,r14D,AppendTo[txtgl,{Text["0",z0]}]
      ]; (*switch*)
    ]; (*If[EvenQ[dir]*)
    txt = ToString[Simplify[symb]];
  ]; (* If[showDeBruijnIndices, *)
  If[showXDeBruijnIndices,
    r = euclidlenN[z1-z0]/6;
    r = euclidlenN[z1-z0]^(1/2)/50.;
    r = euclidlenN[z1-z0]^(2/3)/30.;
    If[EvenQ[dir],
      Switch[type
      ,r23U,AppendTo[marks,{{col2,Disk[z2,r]}}]
      ,r23D,AppendTo[marks,{{col0,Disk[z0,r]}}]
      ,r23L,AppendTo[marks,{{col1,Disk[z3,r]}}]
      ,r14U,AppendTo[marks,{{col1,Disk[z2,r]}}]
      ,r14D,AppendTo[marks,{{col3,Disk[z0,r]}}]
     ]; (*switch*)
    ,(*ELSE*)
      Switch[type
      ,r23U,AppendTo[marks,{{col1,Disk[z2,r]}}]
      ,r23D,AppendTo[marks,{{col3,Disk[z0,r]}}]
      ,r23L,AppendTo[marks,{{col2,Disk[z3,r]}}]
      ,r14U,AppendTo[marks,{{col2,Disk[z2,r]}}]
      ,r14D,AppendTo[marks,{{col0,Disk[z0,r]}}]
      ]; (*switch*)
    ]; (*If[EvenQ[dir]*)
  ]; (* If[showXDeBruijnIndices, *)
  Return[gl]
] (* getGFigure *)

ptinrange[{x_,y_}]:=(xMinCutRng-margin < x < xMaxCutRng+margin && yMinCutRng-margin < y < yMaxCutRng+margin)

inrange[lst_]:= ptinrange[lst];

decomposeFLst[flst_]:=Block[{res={},i,len=Length[flst]},
  Do[
    If[cutOutOfRangeFlag,
      If[inrange[flst[[i,2]]],
 	    If[twostepsFlag,
	      res = Join[res, decomposeFig2steps[flst[[i]] ] ]
	    ,(*ELSE*)
	      res = Join[res, decomposeFig[flst[[i]] ] ]
	    ];
	  ,(*ELSE*)
        If[dbgList,Print[(Plus@@flst[[i,2]])/Length[flst[[i,2]]],"---Eliminating ",flst[[i,2]] ]];
      ];
    ,(*ELSE*)
      If[twostepsFlag,
		res = Join[res, decomposeFig2steps[flst[[i]] ] ]
      ,(*ELSE*)
		res = Join[res, decomposeFig[flst[[i]] ] ]
      ];
    ];
  ,{i,len}];
  Return[res];
] (* decomposeFLst *)

decomposeFig2steps[fig_]:=
Block[{res={},type,z0,z1,z2,z3,subdivLevel,symb,s0,thval,z02,z12,z23},
  {type, z0, dir, mag, extra} = fig;
  {subdivLevel,symb,thval} = extra;
  If[OddQ[dir], sign=-1, sign=1];
  If[dbgLst, Print["decomposeFig2steps ",++count,"/",lstlen,{type, extra}]];
  {z0,z1,z2,z3} = getFigureVertices[fig];
  Switch[type(*wwwwwwwwwwwwwwwwwwwwwwwww*)
  ,r23U,
    If[OddQ[dir], s0=0, s0=3];
    z02 = (tau z0+z2)/tau2; z12 = (tau z1+z2)/tau2; z23 = (tau z3+z2)/tau2;
    AppendTo[res,getFigure[{r23U,z0,mod10[dir+1],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23U,z02,mod10[dir+0],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23D,z02,mod10[dir+2],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23L,z02,mod10[dir-2],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z02,mod10[dir+5],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23U,z12,mod10[dir+1],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14D,z12,mod10[dir+4],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z23,mod10[dir-4],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
  ,r23D,
    If[OddQ[dir], s0=1, s0=2];
    z02 = (tau z0+z2)/tau2; z12 = (tau z1+z2)/tau2; z23 = (tau z3+z2)/tau2;
    AppendTo[res,getFigure[{r23D,z0,mod10[dir+1],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23U,z02,mod10[dir+0],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23L,z02,mod10[dir+2],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23D,z02,mod10[dir-2],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z02,mod10[dir+5],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23D,z12,mod10[dir+1],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z12,mod10[dir+4],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z23,mod10[dir-4],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
  ,r23L,
    If[OddQ[dir], s0=3, s0=0];
    z02 = (tau z0+z2)/tau2; z12 = (tau z1+z2)/tau2; z23 = (tau z3+z2)/tau2;
    AppendTo[res,getFigure[{r23U,z0,mod10[dir+1],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23U,z02,mod10[dir+0],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23U,z02,mod10[dir+2],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23D,z02,mod10[dir-2],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z02,mod10[dir+5],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23D,z12,mod10[dir+1],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z12,mod10[dir+4],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z23,mod10[dir-4],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
  ,r14U,
    If[OddQ[dir], s0=3, s0=0];
    z12 = (z1+tau z2)/tau2; z23 = (z3+tau z2)/tau2;
    AppendTo[res,getFigure[{r23U,z0,mod10[dir+0],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23U,z0,mod10[dir+2],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z12,mod10[dir+4],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23D,z23,mod10[dir+3],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z23,mod10[dir-4],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
  ,r14D,
    If[OddQ[dir], s0=2, s0=1];
    z12 = (z1+tau z2)/tau2; z23 = (z3+tau z2)/tau2;
    AppendTo[res,getFigure[{r23D,z0,mod10[dir+0],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23U,z0,mod10[dir+2],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z12,mod10[dir+4],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r23D,z23,mod10[dir+3],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z23,mod10[dir-4],mag/tau2,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
  ]; (*switch*)
  Return[res]
] (* decomposeFig2steps *)

decomposeFig[fig_]:=
Block[{res={},type,z0,z1,z2,z3,subdivLevel,symb,s0,thval},
  {type, z0, dir, mag, extra} = fig;
  {subdivLevel,symb,thval} = extra;
  If[OddQ[dir], sign=-1, sign=1];
  If[dbgLst, Print["decomposeFig ",++count,"/",lstlen,{type, extra}]];
  {z0,z1,z2,z3} = getFigureVertices[fig];
  Switch[type(*zzzzzzzzzz*)
  ,r23U,
    If[OddQ[dir], s0=0, s0=3];
    AppendTo[res,getFigure[{r23L,z3,mod10[dir-4],mag/tau,{curSubdivLevel,s0-2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z1,mod10[dir+1],mag/tau,{curSubdivLevel,s0-sign,thval/t^2 + 1/t }}] ];
    AppendTo[res,getFigure[{r23D,z2,mod10[dir+5],mag/tau,{subdivLevel,s0,thval/t^2 + 1/t^2 }}] ];
    If[redundantFlag,
      AppendTo[res,getFigure[{r23Dummy,z1,mod10[dir+4],mag/tau,{xxx}}] ];
      AppendTo[res,getFigure[{r14Dummy,z3,mod10[dir-1],mag/tau,{xxx}}] ];
    ]; (* If[redundantFlag, *)
  ,r23D,
    If[OddQ[dir], s0=1, s0=2];
    AppendTo[res,getFigure[{r23U,z3,mod10[dir-4],mag/tau,{subdivLevel,s0,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z1,mod10[dir+1],mag/tau,{curSubdivLevel,s0,thval/t^2 + 1/t }}] ];
    AppendTo[res,getFigure[{r23U,z2,mod10[dir+5],mag/tau,{curSubdivLevel,s0-sign,thval/t^2 + 1/t^2 }}] ];
    If[redundantFlag,
      AppendTo[res,getFigure[{r23Dummy,z1,mod10[dir+4],mag/tau,{xxx}}] ];
      AppendTo[res,getFigure[{r14Dummy,z3,mod10[dir-1],mag/tau,{xxx}}] ];
    ]; (* If[redundantFlag, *)
  ,r23L,
    If[OddQ[dir], s0=3, s0=0];
    AppendTo[res,getFigure[{r23D,z3,mod10[dir-4],mag/tau,{subdivLevel,s0,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z1,mod10[dir+1],mag/tau,{curSubdivLevel,s0+2 sign,thval/t^2 + 1/t }}] ];
    AppendTo[res,getFigure[{r23U,z2,mod10[dir+5],mag/tau,{curSubdivLevel,s0+sign,thval/t^2 + 1/t^2 }}] ];
    If[redundantFlag,
      AppendTo[res,getFigure[{r23Dummy,z1,mod10[dir+4],mag/tau,{xxx}}] ];
      AppendTo[res,getFigure[{r14Dummy,z3,mod10[dir-1],mag/tau,{xxx}}] ];
    ]; (* If[redundantFlag, *)
  ,r14U,
    If[OddQ[dir], s0=3, s0=0];
    AppendTo[res,getFigure[{r23L,z3,mod10[dir-3],mag/tau,{curSubdivLevel,s0+2 sign,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14D,z2,mod10[dir+3],mag/tau,{subdivLevel,s0,thval/t^2 + 1/t }}] ];
    If[redundantFlag,
      AppendTo[res,getFigure[{r23Dummy,z1,mod10[dir+3],mag/tau,{xxx}}] ];
      AppendTo[res,getFigure[{r14Dummy,z2,mod10[dir-3],mag/tau,{xxx}}] ];
    ]; (* If[redundantFlag, *)
  ,r14D,
    If[OddQ[dir], s0=2, s0=1];
    AppendTo[res,getFigure[{r23U,z3,mod10[dir-3],mag/tau,{subdivLevel,s0,thval/t^2 }}] ];
    AppendTo[res,getFigure[{r14U,z2,mod10[dir+3],mag/tau,{curSubdivLevel,s0+sign,thval/t^2 + 1/t }}] ];
    If[redundantFlag,
      AppendTo[res,getFigure[{r23Dummy,z1,mod10[dir+3],mag/tau,{xxx}}] ];
      AppendTo[res,getFigure[{r14Dummy,z2,mod10[dir-3],mag/tau,{xxx}}] ];
    ]; (* If[redundantFlag, *)
  ]; (*switch*)
  Return[res]
] (* decomposeFigure *)

writePts[pts_,fname_]:=
Block[{i,x,y,res={},len=Length[pts],out,count=0,countmarg=0},
  out = OpenWrite[fname];
  Do[
    {x,y} = pts[[i]]//Chop;
    Print[i,"/",len," -> ",{x,y}];
    If[rngminx+outmarg <x< rngmaxx-outmarg && 
       rngminy+outmarg <y< rngmaxy-outmarg,
      count++;
      WriteString[out,ToString[CForm[x]]<>" "<>ToString[CForm[y]]<>" 0\n"];
      AppendTo[res,Point[{x,y}]];
    ,(*ELSE*)
      If[rngminx < x < rngmaxx && rngminy < y < rngmaxy,
        countmarg++;
        WriteString[out,ToString[CForm[x]]<>" "<>ToString[CForm[y]]<>" -1\n"];
        AppendTo[res,{Red,Point[{x,y}]} ];
      ];
    ]
  ,{i,len}];
   Close[out];
  Print["Written: ",count,"+",countmarg,"=",countmarg+count];
  Return[res]
]


addImmovablePts[pts_,fname_]:=
Block[{i,x,y,res={},len=Length[pts],out,count=0,countmarg=0},
  out = OpenAppend[fname];
  Do[
    {x,y} = pts[[i]]//Chop;
      Print["++ ",i,"/",len," -> ",{x,y}];
      If[rngminx < x < rngmaxx && rngminy < y < rngmaxy,
        countmarg++;
        WriteString[out,ToString[CForm[x]]<>" "<>ToString[CForm[y]]<>" -1\n"];
        AppendTo[res,{Red,Point[{x,y}]} ];
      ];
  ,{i,len}];
   Close[out];
  Print["Written: ",count,"+",countmarg,"=",countmarg+count];
  Return[res]
]

selectVisible[lst_]:=Block[{i,x,y,res={},len=Length[lst]},
  Do[
    {x,y} = lst[[i]];
    If[-dr+globaldelta[[1]] < x < dr+globaldelta[[1]] &&
       -dr+globaldelta[[2]] < y < dr+globaldelta[[2]],
      AppendTo[res,Chop[{x,y}] ];
    ];
  ,{i,len}];
  Return[res]
] (* selectVisible *)

findclosestpt[pt_,lst_]:=Block[{i,x,y,d,res,mindist=1000000,len=Length[lst]},
  Do[
    {x,y} = lst[[i]];
    d = euclidlen[{x,y}-pt];
    If[euclidlen[{x,y}-pt] < mindist,
      res = {x,y}//Chop;
      mindist = d;
    ];
  ,{i,len}];
  Return[res]
] (* findclosestpt *)

grada[flst_]:=Block[{res={},i,len=Length[flst],z0,z,inputVal,fig},
  Do[
    fig=flst[[i]];
    {type, z0, dir, mag, extra} = fig;
    {subdivLevel,symb,thval} = extra;
    z = getPivot[fig];
    inputVal = (z[[1]] - xMinCutRng)/(xMaxCutRng-xMinCutRng);
    If[inputVal > thval, AppendTo[res,Point[z] ] ];
  ,{i,len}];
  Return[res]
] (*grada*)

flat[flst_,inputVal_]:=Block[{res={},i,len=Length[flst],z0,z},
  Do[
    fig=flst[[i]];
    {type, z0, dir, mag, extra} = fig;
    {subdivLevel,symb,thval} = extra;
    If[inputVal > thval, AppendTo[res,Point[getPivot[fig]] ] ];
  ,{i,len}];
  Return[res]
] (*flat*)

showValues[flst_]:=Block[{res={},i,len=Length[flst],z0,z,inputVal},
  Do[
    fig=flst[[i]];
    {type, z0, dir, mag, extra} = fig;
    {subdivLevel,symb,thval} = extra;
    AppendTo[res,thval];
  ,{i,len}];
  Return[res]
] (*showValues*)

getFinalSampling[flst_]:=Block[{res={},i,len=Length[flst],z0,z,inputVal},
  Do[
    fig=flst[[i]];
    AppendTo[res,getPivot[fig]];
  ,{i,len}];
  Return[res]
] (*getFinalSampling*)

getPivot[fig_]:=Block[{z0,z1,z2,z3,type=fig[[1]]},
   {z0,z1,z2,z3} = getFigureVertices[fig];
    Switch[type
    ,r23U, z2
    ,r23D, z0
    ,r23L, z3
    ,r14U, z2
    ,r14D, z0
    ]
] (* getPivot *)
(****************** end of procedure *******************)

(************* prog starts here *************)

visuPenroseNara[inniters_:5, dbg_:False] :=
    Module[ {},
		SetOptions[Graphics, ImageSize -> {1024,Automatic},AspectRatio->Automatic, PlotRange->All];
		dbgFlag = dbg;
		niters = inniters;
		twostepsFlag = False;
		If[twostepsFlag, niters = niters/2];
		If[twostepsFlag, lbl = "5fold_2steps", lbl = "5fold_1step"];
		
		showSepArrows = False;
		showSFC = False;
		sfcCol = Green;
		showEdgeArrows = False; (* separate marking for 01 and 12 edges *)
		arrTh = Thickness[.005];
		showVoronoi = False;
		voronoiColor = GrayLevel[.85];
		voronoiTh = Thickness[.004];
		
		{sfcCol,sfcTh} = {Green,Thickness[.001]};
		
		nngon = 40;
		ngon = Table[rotatedaround[{1,0},{0,0},2PI i/nngon],{i,0, nngon-1}];
		ngon20 = Table[rotatedaround[{1,0},{0,0},2PI i/nngon],{i,0, 20}];
		marksz = .01;
		showThDisks = False;
		borderFlag = True;
		showDeBruijnIndices = False;
		showXDeBruijnIndices = False;(*extended DeBruijnIndices*)
		
		showP1marks = True; (* marking like in Penrose P1 set *)
		showNaramarks = False; (* marking like in Penrose P1 set *)
		p1marksth = Thickness[.001];
		
		showOrigin = False;
		showGraphics = True;
		  borderColor = Cyan;
		  borderTh = .001;
		  cutOutOfRangeFlag = False;
		  symbolicForm = True; t = .; (* t in symbolic form *)
		  symbolicForm = False; t = GoldenRatio//N;
		  showThValues = True;
		  showThValues = False;
		redundantFlag = False;		
		globaldelta = {0,0};
		dbgLst = False;
		
        flst = {};
        curSubdivLevel = 0;
        curmag = 1;
        z0 = {0,0} + globaldelta;
        zstar = Table[rotatedaround[{0,tau},z0,(i-2/2) 2Pi/5]//N, {i,5}];
        initgl = txtgl = {};
        arrgl = {arrTh};
        sfcgl = {sfcCol, sfcTh};
        p1marks = marks = {PointSize[marksz]};
        prev1bordergl = prev2bordergl = {};
        prev1marks = prev2marks = {};
        iter = 0;
        bordergl = {PointSize[.005],borderColor, Thickness[borderTh]};
        voronoigl = {voronoiColor,voronoiTh};
        rng = All;
        n = 3;
        {px,py} = {1.2,1.9};
        flst = {};
        Do[
          {x,y} = {px,py 1.1^(-Floor[(i-1)/n])} {(i-1) - n Floor[(i-1)/n], Floor[(i-1)/n]};
          type = tileset[[i]];
          AppendTo[flst,getFigure[{type,z0+{x,y},2,curmag,{curSubdivLevel,0,0}}]];
        ,{i,If[dbgFlag,Length[tileset],1]}];
        Print[flst//mf];
        
        gl0 = getGLst[flst];
        plbl = lbl<>" "<>date;
        p = Graphics[{bordergl,gl0,voronoigl,arrgl,sfcgl,marks,txtgl}];
        p//Print;
       
        gl0 = bordergl/.{Line->Polygon,Cyan->LightYellow};
        bordergl = {PointSize[.005],borderColor, Thickness[borderTh]};
        Do[
          If[ twostepsFlag,
              curmag = curmag/tau2,
              curmag = curmag/tau
          ];
          curSubdivLevel++;
          textsz--;
          Print[">>iter=",iter," curSubdivLevel=",curSubdivLevel];
          plbl = lbl<>" "<>date;
          flst = decomposeFLst[flst];
          gl = {};
          marks = {PointSize[marksz]};
          txtgl = {};
          arrgl = {arrTh};
          sfcgl = {sfcCol, sfcTh};
          voronoigl = {voronoiColor,voronoiTh};
          AppendTo[gl,getGLst[flst]];
          plbl = lbl<>" "<>date;
          p = Graphics[{CapForm["Round"],gl0,voronoigl,arrgl,sfcgl,marks,txtgl, bordergl,gl}, PlotRange->If[dbgFlag,All,All(*{{-.3,.3},{.5,1.1}}*)]];
          p//Print;
          bordergl = {borderColor};
          voronoigl = {voronoiColor,voronoiTh};
        ,{iter,niters}];
    ] (* visuPenroseNara *)
    *)