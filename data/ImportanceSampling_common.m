(* ImportanceSampling_common.m
   V.O. UdeM january 2004
*)


(**************** System-dependent setup ******************)
unknown[name_]:=(Head[name] === Symbol)

If[unknown[workingDirectory], workingDirectory = "c:\\Penrose5fold" ];

date:= ToString[Date[][[1]]]<>"/"<>ToString[Date[][[2]]]<>"/"<>
        ToString[Date[][[3]]]<>"_"<>ToString[Date[][[4]]]<>":"<>
        ToString[Date[][[5]]]<>":"<>ToString[Date[][[6]]];
Off[Solve::"ifun"]; Off[InverseFunction::"ifun"];
SetOptions["stdout", PageWidth -> 100];
pid=ToString[$ProcessID];
(**************** end of System-dependent setup ******************)

(****************** constants *******************)
colorTable={Red,Green,Blue,Cyan,Magenta,Yellow,Gray,Red,Green,Blue,Cyan,Magenta,Yellow};
eps = 10.^(-5);
epsilon = 10.^(-15);
comaprisonepsilon = 10.^(-15); (* for comparisons L2N(v1-v2) *)
zero = {0,0};
(****************** end of constants *******************)

(****************** Procedures *******************)
modN[x_]:=Mod[x,nFold,1] (* nFold must be known *)

dPrint[par___]:=If[dbg,Print[par]]
suppressedText[___]:={};
toString[n_]:=If[0<=n<=9,ToString[n], FromCharacterCode[ToCharacterCode["A"]+n-10] ];

uPerpVectorCCW[{x_,y_}]:= uVector[{-y,x}]
uPerpVectorCW[{x_,y_}]:= uVector[{y,-x}]
uVector[{x_,y_}]:= ({x,y}/euclidlen[{x,y}])

euclidlen[z_]:= Sqrt[ Plus@@(z^2) ]

getAngle[z0_, z1_, z2_] := Block[{res},
  If[z0==z1 || z2==z1, res = 0, res = ArcTan@@(z0-z1) - ArcTan@@(z2-z1) ];
  If[res < 0, res += 2Pi];
  Return[N[res]]
] (* getAngle *)

rotatedaround[pt1_, pt2_, alpha_]:=Block[{a},
(* gives a point which is pt1 rotated around pt2 by angle alpha *)
    If [pt1 == pt2, Return[pt1]];
    a = ArcTan[(pt1-pt2)[[1]],(pt1-pt2)[[2]]];
    Return[ pt2 + euclidlen[pt1-pt2] {Cos[a+alpha], Sin[a+alpha]} ]
]; (* rotatedaround *)

reorgFFT[lst2D_]:= Block[{ydim,xdim,res=lst2D}, (*  2D reorgFFT after FFT *)
  {ydim,xdim} = Dimensions[lst2D];
  Print["reorgFFT xdim=",xdim," ydim=",ydim];
  If [Dimensions[res] != {},
    res = RotateLeft[res,Floor[ydim/2]  ]//Transpose;
    res = RotateLeft[res,Floor[xdim/2]  ]//Transpose
  ];
  Return[res]
] (* reorgFFT *)

getFigure[tileType_,tileRefPt_,tileDir_,LOS_,fcode_]:=Module[
  {shiftedcont = Plus[#,tileRefPt]& /@ tileShapes[[tileType]]},
  {tileType,tileRefPt,rotatedaround[#,tileRefPt,modN[tileDir] 2PI/nFold]& /@ shiftedcont,tileDir,LOS,fcode}
]

    
getGLst[flst_]:=Map[getGFigure, flst, {1}]

getGFigure[fig_]:=Block[{gl={},tileType,tileRefPt,tileVLst,tileDir,LOS,fcode,center,index,subseq,shape},
  If[dbgLst, Print["getGFigure ",++count," iter=",iter]];
  {tileType,tileRefPt,tileVLst,tileDir,LOS,fcode} = fig;
  If[fillFigure, 
    shape = getTileShape[fig];
    If[shape =!= {}, AppendTo[gl, {tileFillCol[[tileType]], Polygon[getTileShape[fig]]} ] ]
  ];
  If[showTileShapes,
    shape = getTileShape[fig];
    If[shape =!= {}, AppendTo[shapesgl, {tileShapesTh,tileShapesCol,Line[getTileShape[fig]]} ] ];
  ];
  If[showiSimplifiedTileShapes&&tileType<=nSubdivisionTiles,AppendTo[shapesgl,Line[Append[tileVLst,First[tileVLst]]]]];
  If[showTileShapesWithDashing,
      If[tileType <= nSubdivisionTiles,
	AppendTo[gl, {GrayLevel[.7], Polygon[Append[tileVLst,First[tileVLst]] ]} ];
        AppendTo[gl, {Dashing[{.003,.003}],Line[Append[tileVLst,First[tileVLst]] ]} ];
      ];
      AppendTo[gl, {tileShapesTh,tileShapesCol, Line[getTileShape[fig]]} ];
  ];
  If[showTileType,
    center = getCenter[fig];
    AppendTo[gl, {Black,Text[ToString[tileType],center,{0,1}] } ];
  ];
  If[showFcode,
    center = getCenter[fig];
    str={};Do[str=str<>ToString[fcode[[i]]],{i, Length[fcode]  }];
    AppendTo[gl, {Red,Text[str,center,{0,-1}] } ];
    AppendTo[gl, {Blue,Text[fcodeValue[fcode],center,{0,1}] } ];
  ];
  If[showShortFcode && tileType > nSubdivisionTiles,
    center = getCenter[fig];
    subseq = Take[fcode,Min[secondLUTdimDepth,Length[fcode]]];
    index = fcodeValue[subseq];
    str={};Do[str=str<>ToString[subseq[[i]]],{i, Length[subseq]  }];
    AppendTo[gl, {Black,Text[ToString[index],center,{0,0}] } ];
    
  ];
  (*If[showSamplingPt, If[tileType > nSubdivisionTiles, AppendTo[sgl, Point[tileVLst[[1]]] ] ] ];*)
  If[showSamplingPt, AppendTo[sgl, Point[tileVLst[[1]]] ] ];
  If[showdir,
      center = getCenter[fig];
      vect2 = scaletab[[LOS+3]] vv[[tileDir]];
      vect1 = {vect2[[2]],-vect2[[1]]};
      AppendTo[gl, {Blue,Line[{center+vect2,center,center+vect1}]} ];
  ];
  If[showAmmannBars && tileType <= nSubdivisionTiles,
    {z0,z1,z2} = tileVLst;
    If[tileType == 1,
      z01 = (z0+3 z1)/4;
      z12 = (z1+z2)/2;
      z220 = z2 + phi^-2 (z0-z2) / 2;
      z122 = z2 + phi1 (z1-z2) / 4;
      AppendTo[gl, {ammanBarCol,ammanBarTh,Line[{z01,z12,z220,z122}]} ];
    ];
    If[tileType == 2,
      z01 = (3 z0+z1)/4;
      z20 = (z0+z2)/2;
      z122 = z2 + phi^-2 (z1-z2) / 2;
      z220 = z2 + phi1 (z0-z2) / 4;
      AppendTo[gl, {ammanBarCol,ammanBarTh,Line[{z01,z20,z122,z220}]} ];
    ];
    If[tileType == 3,
      z001 = (phi z0 + (z0 + phi (z0+z1)/2)/phi2)/phi2;
      z20 = (z0+z2)/2;
      z200 = z0 + phi1 (z2-z0) / 4;
      z122 = z2 + phi^-2 (z1-z2) / 2;
      z01 = z0 + (phi1-phi^-5/4) (z1-z0) ;
      AppendTo[gl, {ammanBarCol,ammanBarTh,Line[{z200,z001,z20,z122,z01}]} ];
    ];
    If[tileType == 4,
      z100 = (phi z1 + (z1 + phi (z0+z1)/2)/phi2)/phi2;
      z12 = (z2+z1)/2;
      z112 = z1 + phi1 (z2-z1) / 4;
      z220 = z2 + phi^-2 (z0-z2) / 2;
      z01 = z1 + (phi1-phi^-5/4) (z0-z1) ;
      AppendTo[gl, {ammanBarCol,ammanBarTh,Line[{z112,z100,z12,z220,z01}]} ];
    ];
(*
Print[euclidlen[z200-z0]/euclidlen[z2-z0] ];
*)
  ]; (* showAmmannBars *)
  Return[gl]
] (* getGFigure *)

 
getMaxLOSinTrianle[{v1_,v2_,v3_}]:=Max[getLOS[v1],getLOS[v2],getLOS[v3],getLOS[(v1+v2+v3)/3.],
    getLOS[(v1+v2)/2],getLOS[(v2+v3)/2],getLOS[(v3+v1)/2],
    getLOS[(v1+3v2)/4],getLOS[(v2+3v3)/4],getLOS[(v3+3v1)/4],
    getLOS[(3v1+v2)/4],getLOS[(3v2+v3)/4],getLOS[(3v3+v1)/4],
    getLOS[(v1+v2+2v3)/4],getLOS[(v2+v3+2v1)/4],getLOS[(v3+v1+2v2)/4]
] (* getMaxLOSinTrianle slightly more sophisticated version *)
 
getMaxLOSinTile[vlst_]:=Max[(getLOS /@ vlst),(Plus @@ vlst)/Length[vlst]] (* the simplest version *)

getImportance[{x0_,y0_}]:= (* x,y between 0 and 1; provided global importanceImage od size {xsize,ysize} *)
Module[{ix,iy,x,y},
  {x,y} = {(x0-xmin)/(xmax-xmin),(y0-ymin)/(ymax-ymin)};
  {ix,iy} = Mod[Round[{x,y} {xsize,ysize}], {xsize,ysize}, 1];
  importanceImage[[iy,ix]]
] (* getImportance *)

getLOS[{x_,y_}]:= (* provided global imageLOS of known size {xsize,ysize} *)
Module[{ix,iy},
  {ix,iy} = Ceiling[{(x-xmin)/(xmax-xmin),(y-ymin)/(ymax-ymin)} {xsize,ysize}];
  If[1 <= ix <= xsize && 1 <= iy <= ysize,
    imageLOS[[iy,ix]]
  ,(*ELSE*)
    imageLOS[[Max[Min[iy,ysize],1],Max[Min[ix,xsize],1] ]]
    (*maxLOS*)
  ]
] (* getLOS *)

ptinrangePlusMargin[{x_,y_},margin_]:=(xmin-margin < x < xmax+margin && ymin-margin < y < ymax+margin)
ptinrangeMinusMargin[{x_,y_},margin_]:=(xmin+margin < x < xmax-margin && ymin+margin < y < ymax-margin)
ptinrange[{x_,y_}]:= (xmin < x < xmax && ymin < y < ymax)
ptinrange[x_,y_]:= (xmin < x < xmax && ymin < y < ymax)

figInRange[{tileType_,tileRefPt_,tileVLst_,tileDir_,LOS_,fcode_},margin_]:= Or @@ Flatten[{
  {((xmin-margin < #[[1]] < xmax+margin && ymin-margin < #[[2]] < ymax+margin)& /@ tileVLst)},
  (xmin-margin < (Plus@@tileVLst/Length[tileVLst])[[1]] < xmax+margin &&  ymin-margin < (Plus@@tileVLst/Length[tileVLst])[[2]] < ymax+margin)}]

(***** 2 variants:
figInRange[{tileType_,tileRefPt_,tileVLst_,tileDir_,LOS_,fcode_},margin_]:= Or @@ Flatten[{
  {((xmin-margin < #[[1]] < xmax+margin && ymin-margin < #[[2]] < ymax+margin)& /@ tileVLst)},
  (xmin-margin < (Plus@@tileVLst/Length[tileVLst])[[1]] < xmax+margin &&  ymin-margin < (Plus@@tileVLst/Length[tileVLst])[[2]] < ymax+margin)}]
figInRange[{tileType_,tileRefPt_,tileVLst_,tileDir_,LOS_,fcode_},margin_]:= Or @@ 
  ((xmin-margin < #[[1]] < xmax+margin && ymin-margin < #[[2]] < ymax+margin)& /@ tileVLst)
***********)
  
getSamplingPtUsingLUT[fig_,lutIndex1_]:=Module[{tileType,tileRefPt,tileVLst,tileDir,LOS,fcode,k21,k31,lutindex2,significantSeq,vect1,vect2},
  {tileType,tileRefPt,tileVLst,tileDir,LOS,fcode} = fig;
  vect2 = vv[[tileDir]];
  vect1 = {vect2[[2]],-vect2[[1]]};
  significantSeq = Take[fcode,Min[secondLUTdimDepth,Length[fcode]]];
  lutindex2 = fcodeValue[significantSeq]+1;
  {k21,k31} = lut[[lutIndex1,lutindex2]];
  tileRefPt + k21 scaletab[[LOS]] vect1 + k31 scaletab[[LOS]] vect2
] (* getSamplingPtUsingLUT *)

getCenter[fig_]:= (Plus @@ fig[[3]])/Length[fig[[3]]]

dumpFlst[fname_]:=Block[{i,out,fig,countmarg=0,tileType,tileRefPt,tileVLst,tileDir,LOS,fcode},
  out = OpenWrite[fname];
  Do[
    fig = {tileType,tileRefPt,tileVLst,tileDir,LOS,fcode} = flst[[i]];
      WriteString[out,ToString[CForm[i]]<>" "<>ToString[Chop[fig]]<>"\n"];
  ,{i,Length[flst]}];
   Close[out];
] (* dumpFlst *)

showCoefs[lut_,lbl_]:=Module[{val,indexmin,indexmax,g1,g2,j,tab},
  indexmin = 1;
  indexmax = Length[lut];
  g1 = Table[
    tab = Table[{i,lut[[i,j,1]]}, {i,indexmin,indexmax}];
    ListPlot[tab,PlotStyle->getColor[j],DisplayFunction:>Identity]
  ,{j,secondLUTdimSize}];
  g2 = Table[
    tab = Table[{i,lut[[i,j,2]]}, {i,indexmin,indexmax}];
    ListPlot[tab,PlotStyle->getColor[j],DisplayFunction:>Identity]
  ,{j,secondLUTdimSize}];
  Show[g1,g2,DisplayFunction->$DisplayFunction,PlotRange->All,PlotLabel->lbl];
]; (* showCoefs *)

smoothLUT[lut_]:=Module[{newlut=lut,level,firstLevel=1,lastLevel=firstLUTdimSize+1},
   newlut[[firstLevel+1]] = (2*lut[[firstLevel]]+4*lut[[firstLevel+1]]+2*lut[[firstLevel+2]]+lut[[firstLevel+3]])/9.;
   Do[
     newlut[[level]] = (lut[[level-2]]+2*lut[[level-1]]+4*lut[[level]]+2*lut[[level+1]]+lut[[level+2]])/10.;
  ,{level,firstLevel+2,lastLevel-2}];
   newlut[[lastLevel-1]] = (2*lut[[lastLevel]]+4*lut[[lastLevel-1]]+2*lut[[lastLevel-2]]+lut[[lastLevel-3]])/9.;
  newlut
] (* smoothLUT *)

(******
*********)

writeSamplingPtsOverThreshold[flst_,level_,fname_,margin_]:=
Block[{i,x,y,gl={},innerTab={},outerTab={},len=Length[flst],out,count=0,fig,countmarg=0,tileType,tileRefPt,tileVLst,tileDir,LOS,fcode,threshold},
  out = OpenWrite[fname];
  (* first pass: inner points only *)
  Do[
    fig = {tileType,tileRefPt,tileVLst,tileDir,LOS,fcode} = flst[[i]];
    If[tileType <= nSubdivisionTiles, Continue[] ];
    {x,y} = getSamplingPtUsingLUT[fig,getLUTIndex[level]] + rtab[[Mod[i,1024,1] ]];
    If[ptinrangeMinusMargin[{x,y},margin] && level >=  fcodeValue[fcode],
      count++;
      WriteString[out,ToString[CForm[x]]<>" "<>ToString[CForm[y]]<>" 0\n"];
      AppendTo[gl,{Red,Point[{x,y}]}];
      AppendTo[innerTab,i];
    ];
  ,{i,len}];
  (* second pass: outer points only *)
  Do[
    fig = {tileType,tileRefPt,tileVLst,tileDir,LOS,fcode} = flst[[i]];
    If[tileType <= nSubdivisionTiles, Continue[] ];
    {x,y} = getSamplingPtUsingLUT[fig,getLUTIndex[level]] + rtab[[Mod[i,1024,1] ]];
    If[!ptinrangeMinusMargin[{x,y},margin] && ptinrange[{x,y}] && level >= fcodeValue[fcode],
      countmarg++;
      WriteString[out,ToString[CForm[x]]<>" "<>ToString[CForm[y]]<>" -1\n"];
      AppendTo[gl,{Black,Point[{x,y}]} ];
      AppendTo[outerTab,i];
    ];
  ,{i,len}];
   Close[out];
  dPrint["Written into ",fname," : ",count,"+",countmarg,"=",count+countmarg];
  Return[{innerTab,outerTab,gl}] (* indTab is the table of indeces: sampling points -> flst *)
] (* writeSamplingPtsOverThreshold *)



showOneOctave[lbl_]:=Module[{ix,iy,importance,samplingPt,selectedPts,count=0},
  {ysize,xsize} = {256,256};
  importanceImage = Table[Round[Fibonacci[2*maxLOS]+ (ix-1)*(Fibonacci[2*maxLOS+2]-Fibonacci[2*maxLOS])/256],{iy,256},{ix,256}];
  area12 = 0.769421; (*triangles 1,2 Det[{tileShapes[[1,2]]-tileShapes[[1,1]],tileShapes[[1,3]]-tileShapes[[1,1]]}]/2. *)
  area34 = 1.24495;  (*triangles 3,4 Det[{tileShapes[[3,2]]-tileShapes[[3,1]],tileShapes[[3,3]]-tileShapes[[3,1]]}]/2. *)
  workingArea  = area34/phi^2;
  workingSqSize = Sqrt[workingArea];
  workingRectSize = Sqrt[workingArea/2];
  If[xsize == ysize,{{xmin, xmax},{ymin,ymax}}={{1.30902-workingSqSize/2,1.30902+workingSqSize/2},{.005,.005+workingSqSize}}]; 
  If[xsize == 2*ysize,{{xmin, xmax},{ymin,ymax}}={{1.30902-workingRectSize,1.30902+workingRectSize},{.005,.005+workingRectSize}}]; 
  rng = {{xmin, xmax},{ymin, ymax}};
  frame={Yellow,Thickness[.005],Line[{{xmin, ymin},{xmin, ymax},{xmax, ymax},{xmax, ymin},{xmin, ymin}}]};

  selectedPts = {Red,PointSize[.0075]}; count=0;
  Print["showOneOctave: Processing flst of length ",Length[flst]];
  Do[
    fig={tileType,tileRefPt,tileVLst,tileDir,LOS,fcode} = flst[[i]];
    If[tileType <= nSubdivisionTiles, Continue[] ];
    importance = getImportance[tileRefPt];
    If[ptinrange[tileRefPt] && importance >= fcodeValue[fcode],
        samplingPt = getSamplingPtUsingLUT[fig, getLUTIndex[importance]];
        AppendTo[selectedPts,Point[samplingPt] ];
        count++;
    ];
  ,{i,Length[flst]}];
  Show[Graphics[{selectedPts,frame}],PlotRange->rng,PlotLabel->lbl];
]; (* showOneOctave *)


findk21k31[{x_,y_},{x21_,y21_},{x31_,y31_}]:={(x31*y-x*y31)/(x31*y21-x21*y31),(x21*y-x*y21)/(-(x31*y21)+x21*y31)}
(* find {k21,k31} such that {x,y}==k21{x21,y21}+k31{x31,y31} *)


getRandomOnUnitCircle[]:=With[{r=1,theta=2 Pi Random[]}, r {Sin[theta], Cos[theta]} ]
getRandomInUnitDisc[]:=With[{r=Random[],theta=2 Pi Random[]}, r {Sin[theta], Cos[theta]} ]

rtab={100, 82, 58, -14, 116, -132, 222, 226, 100, 19, 201, -222, 100, 248, -43, 
 154, 172, 93, -7, 161, 14, 112, 126, -196, -84, 99, -89, -239, 97, -145, 
 134, 104, 177, 141, -142, 212, 227, 70, -119, -183, -241, 39, 7, -156, -10, 
 185, -250, -203, -102, -224, 214, 218, 102, -73, 131, -84, -239, 230, -15, 
 -150, -247, 19, 101, 208, 19, 215, 206, 144, 207, 120, 147, -57, -120, -214, 
 100, -133, -48, -191, 108, -159, -16, -55, 46, -174, 212, -165, 159, 74, 
 164, 238, -78, 165, 127, 166, 135, -232, -143, 67, -150, -208, -181, -37, 
 -147, 139, 176, -150, -233, 78, 99, -198, 87, 98, 36, 95, 33, 173, -43, 
 -215, -16, -181, 10, -26, 169, -167, 96, -30, -189, 7, -25, 3, -91, -231, 
 199, 251, -164, -91, 145, -221, -120, -95, -50, -108, -237, 172, -233, 125, 
 -78, -206, 168, 189, -166, -11, 52, -204, 118, 142, 210, 186, 232, -81, -74, 
 -8, 100, 193, -201, 190, -108, -209, 2, -64, -68, 139, -20, 19, -202, -173, 
 -28, -114, -206, -244, -248, 171, 161, 186, 24, 117, -112, 101, -49, -220, 
 234, -43, 34, -69, 115, 162, 74, 116, 27, 11, -97, 60, 153, -127, -157, -28, 
 189, -70, -151, 111, 53, 170, 71, 211, -177, 71, 168, 9, 225, -56, -215, 
 105, 118, 113, -96, -236, -246, 235, -193, 197, 225, 9, 5, -169, -242, -35, 
 48, 128, -26, 244, 226, -158, 50, 96, 195, -192, -69, 40, 255, -37, 149, 
 222, -104, 127, 11, 20, 56, -54, -80, -251, 245, 8, -77, -196, 171, 43, -15, 
 -133, 203, -65, 126, -64, -82, -78, 191, -53, 208, 19, -250, 194, -46, 254, 
 -187, -131, 225, 27, -13, -232, -217, -234, -171, -134, 249, -112, 38, -68, 
 184, 162, -27, -23, 27, -163, -199, -80, 232, -48, -30, -233, 158, -160, 
 -187, -60, 111, 67, 84, -207, 156, -222, 117, -48, -136, 144, 122, -246, 
 -225, 224, 163, 235, -42, -118, 253, 3, -107, 240, 199, 249, -71, 30, -132, 
 -188, -18, 60, -145, -142, -162, -81, -201, -180, -131, -172, -84, 247, -63, 
 -33, -124, -27, 159, 9, -163, 11, -55, -219, -50, -204, -70, -28, 174, 149, 
 175, -198, -203, -205, 253, 5, 178, 152, 136, -40, 225, -54, -256, -152, 
 -108, -35, 232, -119, 3, 172, -205, 4, -24, -79, 236, 217, -247, 102, 164, 
 -214, -242, -34, 7, -206, -34, 233, -121, 192, 15, 146, 247, -119, -152, 
 130, 230, 195, 67, -113, -15, -174, 107, 195, 65, 198, 166, -43, 46, 65, 80, 
 17, -152, -154, 170, -12, -93, 228, 80, -163, 122, -2, -165, -230, 121, 
 -232, -46, -79, 244, 242, -14, -20, -238, 20, -256, 217, -43, -236, 45, 44, 
 -118, -112, -198, -105, 187, -161, -31, -76, 130, -124, 213, 174, -79, -166, 
 -151, -95, 7, -75, -253, -223, -59, 133, 67, -165, 7, -223, 87, 129, -248, 
 87, 141, -169, 247, 21, -214, 222, -207, 197, -134, -241, -81, 104, 37, 
 -211, -91, 11, 196, 46, 11, 41, -218, -35, -225, -134, -124, -56, -116, 19, 
 157, 139, -9, -28, 201, 29, -212, -60, -174, -139, -105, 138, -104, 35, 234, 
 -12, 163, 90, 179, 166, 218, 148, -184, -215, -86, -190, -190, -241, 238, 
 -194, -81, -47, -79, -14, 182, 16, 222, -61, 4, 84, 0, 218, 66, -216, -206, 
 43, -48, 128, 128, -218, -26, 136, -28, -101, -32, 32, -62, 198, -40, -212, 
 148, 237, 38, -253, -202, -232, 170, 254, -154, -210, -244, 106, -46, -193, 
 -194, -116, -164, -11, 83, -84, -156, -94, 159, 0, 49, -109, -143, -73, 160, 
 -97, 162, 42, 83, -191, 231, 221, -172, -98, 192, -103, 79, 41, 77, 18, 
 -208, 246, 57, -13, -64, 182, 163, -163, -153, 210, 66, -81, 188, 157, 222, 
 83, 255, 114, -200, 25, -212, 5, -39, -243, -15, -254, -32, -109, 33, 183, 
 -128, -74, -2, 202, 51, 172, 109, 83, -234, -142, -135, -238, 35, -93, -98, 
 -89, -89, 245, 214, -169, 134, -191, 192, -253, 83, -130, -48, 164, -71, 
 182, 155, 84, -130, -69, 164, 48, -198, -169, -28, 119, -239, -52, 43, -88, 
 -105, 224, -31, 26, 115, -80, 73, -179, -105, 58, 177, -152, -197, 108, 214, 
 152, -64, 43, -196, 119, -216, -12, -233, -28, 111, 20, -202, 148, -107, 
 171, 115, 190, 140, 6, 153, -57, -6, -137, -98, -145, 130, -148, -142, 247, 
 -4, 186, -153, 136, 206, -83, 42, -89, -233, -81, -128, 220, -211, 84, 91, 
 157, -139, 132, 192, 5, 151, 14, 208, -151, 249, 233, -44, 211, -40, -141, 
 -75, -205, -147, 62, -55, 79, -188, 2, 44, 155, -169, 207, -245, -36, -128, 
 37, 186, -131, -72, 94, 176, -225, 224, -95, 103, -89, -52, -147, -169, 
 -239, -254, -151, 52, -77, -162, 240, -240, 212, -61, -221, 219, 77, 243, 
 -25, -180, 44, -112, 76, -225, -56, 43, 247, 243, 106, 9, -123, -168, 245, 
 -189, 55, -174, 211, 238, 184, 17, 198, 158, -159, -110, -199, -248, -70, 
 128, 136, 10, 2, -66, 249, 176, -115, 231, -224, -112, 153, 194, -201, -168, 
 -83, 6, 51, 196, 149, -157, -77, 158, 38, 224, 228, 161, 60, -166, -44, 136, 
 -200, 203, -5, 151, -217, -25, -97, 159, 30, 94, 108, -174, -57, 199, 69, 
 -58, 170, 32, 119, -86, 124, -88, -159, 43, 251, -119, -95, 158, 203, 189, 
 135, -123, -6, -171, -19, 57, 177, 216, -203, -82, 235, -94, 93, -68, 53, 
 -129, 154, -18, 204, 240, 118, 197, -60, -62, -145, -15, -157, -82, -254, 
 -67, -44, -74, -63, -98, -45, -198, -90, 90, -7, -195, -19, 148, -217, 86, 
 -91, 239, 20, -76, -99, -13, -158, -227, 95, -225, -183, 173, -141, 236, 
 -177, -94, -253, -58, -95, 228, 201, -98, 243, -70, -199, -50, 14, 23, -138, 
 -52, 101, 156, 87, -151, -121, -157, -116, -22, -92, 59, 14, 236, -58, 60, 
 -250, -43, -246, 189, 127, -33, -49, -199, 115, 207, 172, -1, -129, 193, 
 231, -199, -42, 165, -96, 226, -164, 224, 216, 173, -24, 15, -77, 98, -169, 
 -94, -127, -67, -125, -149, -178, 77, 39}/10^10//N;
(****************** end of Procedures *******************)
Print["ImportanceSampling_common.m loaded."];

