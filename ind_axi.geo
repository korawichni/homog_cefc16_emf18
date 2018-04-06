Include "ind_axi_dat.pro";

p=.07e-3;
pc = Rc/8;

pr = 0.3e-3*1.5;
pr2 = 1e-3;
pr3 = 4e-3;
pax = 0.5e-3;

prW = pr*0.7*MD ;

dPo[]+=newp ; Point(dPo[0])  = {X1, Y2, 0 , pax};
dPo[]+=newp ; Point(dPo[1])  = {X2, Y2, 0 , pr2};
dPo[]+=newp ; Point(dPo[2])  = {X3, Y2, 0 , pr2};
dPo[]+=newp ; Point(dPo[3])  = {X4, Y2, 0 , pr2};

dPo[]+=newp ; Point(dPo[4])  = {X4, Y3, 0 , pr2};
dPo[]+=newp ; Point(dPo[5])  = {X4, Y4, 0 , pr};
dPo[]+=newp ; Point(dPo[6])  = {X4,  0, 0 , pr};

dPo[]+=newp ; Point(dPo[7])  = {X3,  0, 0 , pr};
dPo[]+=newp ; Point(dPo[8])  = {X3, Y4, 0 , pr};
dPo[]+=newp ; Point(dPo[9])  = {X3, Y3, 0 , pr};

dPo[]+=newp ; Point(dPo[10]) = {X2, Y3, 0 , pr};
dPo[]+=newp ; Point(dPo[11]) = {X2, Y4, 0 , pr*0.7};
dPo[]+=newp ; Point(dPo[12]) = {X1, Y4, 0 , pax};
dPo[]+=newp ; Point(dPo[13]) = {X1, Y3, 0 , pax};

dPg[]+=newp ; Point(dPg[0]) = {X2, 0, 0 , pr*0.7};
dPg[]+=newp ; Point(dPg[1]) = {X1, 0, 0 , pax};


For k1 In {0:#dPo[]-1}
  k2 = (k1 != #dPo[]-1) ? k1+1 : 0 ;
  dLo[]+=newl ; Line(dLo[k1]) = {dPo[k1], dPo[k2]};
EndFor
lliron=newll ; Line Loop(lliron) = {dLo[]};
surfiron[] += news ; Plane Surface(surfiron[0]) = {lliron};

dLg[]+=newl ; Line(dLg[0]) = {dPo[11],dPg[0]};
dLg[]+=newl ; Line(dLg[1]) = {dPg[0],dPg[1]};
dLg[]+=newl ; Line(dLg[2]) = {dPg[1],dPo[12]};

llgap1 = newll ;
Line Loop(llgap1) = {dLg[], -dLo[11]};
surfgap1[] += news ; Plane Surface(surfgap1[0]) = {llgap1};

dPw[]+=newp ; Point(dPw[0])  = {Xw1, 0*Yw1, 0., prW};
dPw[]+=newp ; Point(dPw[1])  = {Xw1, Yw2, 0., prW};
dPw[]+=newp ; Point(dPw[2])  = {Xw2, Yw2, 0., prW};
dPw[]+=newp ; Point(dPw[3])  = {Xw2, 0*Yw1, 0., prW};

dLw[]+=newl ; Line(dLw[0]) = {dPg[0], dPw[0]};
For k1 In {0:#dPw[]-2}
  k2 = (k1 != #dPw[]-1) ? k1+1 : 0 ;
  dLw[]+=newl ; Line(dLw[1+k1]) = {dPw[k1], dPw[k2]};
EndFor
dLw[]+=newl ; Line(dLw[4]) = {dPw[3],dPo[7]};

llgap2 = newll ; Line Loop(llgap2) = {dLg[0],dLw[],dLo[{7:10}]};
surfgap2[] += news ; Plane Surface(surfgap2[0]) = {llgap2};

surfhomo[]={}; sw[]={}; is[]={};

If (Fine)
  // Auxiliar elementary cell
  cen0 = newp ; Point(cen0)  = { 0, 0, 0 , pc};
  pw[]+=newp ; Point(pw[0])  = { Rc,  0, 0 , pc};
  pw[]+=newp ; Point(pw[1])  = {  0, Rc, 0 , pc};
  pw[]+=newp ; Point(pw[2])  = {-Rc,  0, 0 , pc};
  pw[]+=newp ; Point(pw[3])  = {  0,-Rc, 0 , pc};

  cw[]+=newl ; Circle(cw[0])  = {pw[0], cen0, pw[1]};
  cw[]+=newl ; Circle(cw[1])  = {pw[1], cen0, pw[2]};
  cw[]+=newl ; Circle(cw[2])  = {pw[2], cen0, pw[3]};
  cw[]+=newl ; Circle(cw[3])  = {pw[3], cen0, pw[0]};

  llcw0 = newll ; Line Loop(llcw0) = {cw[]};
  sw0 = news ; Plane Surface(sw0) = {llcw0};

  pb[]+=newp ; Point(pb[0])  = { Dex/2,  Dey/2, 0 , p};
  pb[]+=newp ; Point(pb[1])  = {-Dex/2,  Dey/2, 0 , p};
  pb[]+=newp ; Point(pb[2])  = {-Dex/2, -Dey/2, 0 , p};
  pb[]+=newp ; Point(pb[3])  = { Dex/2, -Dey/2, 0 , p};

  lb[]+=newl ; Line(lb[0])  = {pb[0], pb[1]};
  lb[]+=newl ; Line(lb[1])  = {pb[1], pb[2]};
  lb[]+=newl ; Line(lb[2])  = {pb[2], pb[3]};
  lb[]+=newl ; Line(lb[3])  = {pb[3], pb[0]};

  llb0 = newll ; Line Loop(llb0) = {lb[]};
  is0 = news ; Plane Surface(is0) = {llb0,llcw0};

  gx = (Xw2-Xw1)-NbrLayersX * Dex ;
  gy = 0.5*((Yw2-Yw1)-NbrLayersY * Dey) ;
  NbrLayersY2 = 0.5*NbrLayersY;

  For ix In {0:NbrLayersX-1}
    For iy In {0:NbrLayersY2-1}
      xaux = Xw1+gx/2+ix*Dex+Dex/2 ;
      yaux = 0*(Yw1+gy/2)+iy*Dey+Dey/2 ;
      surf[]=Translate {xaux, yaux, 0} {Duplicata{ Surface{sw0};}}; sw[] += surf[0] ;
      surf[]=Translate {xaux, yaux, 0} {Duplicata{ Surface{is0};}}; is[] += surf[0] ;
      If(ix==0)
        bnd[] = Boundary{Surface{is[ix*NbrLayersY2+iy]};}; bndi[] += bnd[{1}];
      EndIf
      If(ix==NbrLayersX-1)
        bnd[] = Boundary{Surface{is[ix*NbrLayersY2+iy]};}; bndi[] += bnd[{3}];
      EndIf
      If(iy==0)
        bnd[] = Boundary{Surface{is[ix*NbrLayersY2+iy]};}; axx[] += bnd[{2}];
      EndIf
      If(iy==NbrLayersY2-1)
        bnd[] = Boundary{Surface{is[ix*NbrLayersY2+iy]};}; bndi[] += bnd[{0}];
      EndIf
    EndFor
  EndFor
  //Removing auxiliar elementary cell
  Delete{Surface{sw0,is0};Line{cw[],lb[]};Point{pw[],pb[]};}
  bnd[] = Boundary{Line{axx[0]};}; dPw[]+= bnd[0];
  bnd[] = Boundary{Line{axx[{#axx[]-1}]};}; dPw[]+= bnd[1];

  dLw2[]+=newl ; Line(dLw2[0]) = {dPw[4], dPw[0]};
  dLw2[]+=newl ; Line(dLw2[1]) = {dPw[3], dPw[5]};

  llisol = newll ; Line Loop(llisol) = {bndi[],dLw2[0], dLw[{1:3}],dLw2[1]};
  surfgap3[] += news ; Plane Surface(surfgap3[0]) = {llisol};
EndIf


If(!Fine)
  /*
  gx = (Xw2-Xw1)-NbrLayersX * Dex ;
  gy = (Yw2-Yw1)-NbrLayersY * Dey ;
  Xw1_ = Xw1+gx/2 ;
  Xw2_ = Xw2-gx/2 ;
  Yw1_ = 0 ;
  Yw2_ = Yw2-gy/2 ;
  */
  dPis[]+=newp ; Point(dPis[0])  = {Xw1_, Yw1_, 0 , prW};
  dPis[]+=newp ; Point(dPis[1])  = {Xw1_, Yw2_, 0 , prW};
  dPis[]+=newp ; Point(dPis[2])  = {Xw2_, Yw2_, 0 , prW};
  dPis[]+=newp ; Point(dPis[3])  = {Xw2_, Yw1_, 0 , prW};

  For k1 In {0:#dPw[]-1}
    k2 = (k1 != #dPw[]-1) ? k1+1 : 0 ;
    bndi[]+=newl ; Line(bndi[k1]) = {dPis[k1], dPis[k2]};
  EndFor
  llhomo = newll ; Line Loop(llhomo) = {bndi[]};
  surfhomo[] += news ; Plane Surface(surfhomo[0]) = {llhomo};

  dLw2[]+=newl ; Line(dLw2[0]) = {dPis[0], dPw[0]};
  dLw2[]+=newl ; Line(dLw2[1]) = {dPw[3], dPis[3]};
  llisol = newll ; Line Loop(llisol) = {bndi[{0:2}], -dLw2[1], -dLw[{3:1}],-dLw2[0]};
  surfgap3[] += news ; Plane Surface(surfgap3[0]) = {llisol};
EndIf


line_out[] = {dLo[{0:5,12:13}],dLg[2]};

If(!Flag_HalfModel)
  line_out[]+=Symmetry {0,1,0,0} { Duplicata{Line{line_out[]};} }; // For convenience
  surfiron[]+=Symmetry {0,1,0,0} { Duplicata{Surface{surfiron[0]};} };
  surfgap1[]+=Symmetry {0,1,0,0} { Duplicata{Surface{surfgap1[0]};} };
  surfgap2[]+=Symmetry {0,1,0,0} { Duplicata{Surface{surfgap2[0]};} };
  surfgap3[]+=Symmetry {0,1,0,0} { Duplicata{Surface{surfgap3[0]};} };

  If(Fine)
    nn = #sw[];
    For k In {0:nn-1}
      sw[] += Symmetry {0,1,0,0} { Duplicata{Surface{sw[k]};} };
    EndFor
    is[] += Symmetry {0,1,0,0} { Duplicata{Surface{is[]};} };
  Else
    surfhomo[]+=Symmetry {0,1,0,0} { Duplicata{Surface{surfhomo[0]};} };
  EndIf
EndIf


// Physical regions
//===================================
Physical Line(OUTBND) = {line_out[]};
Physical Surface(IRON) = {surfiron[]};
Physical Surface(AIRGAP)  = {surfgap1[]};
Physical Surface(AIR)  = {surfgap2[],surfgap3[]};

Testing = 0;
If(Testing==0)
  If(Fine)
    nn = #sw[];
    For k In {0:nn-1}
      Physical Surface(iCOND+k) = {sw[k]};
    EndFor
    Physical Surface(INSULATION) = {is[]};

    Physical Surface(ALLCOND) = {sw[]};
  Else
    Physical Surface(iCOND) = {surfhomo[]};
  EndIf
EndIf

If(Testing==1)
  //  gmsh ind_axi.geo -setnumber Flag_HomogenisedModel 0  -2 -o homog_fine.msh -v 3
  If(Fine)
    Physical Surface(iCOND) = {sw[],is[]};
  Else
    Physical Surface(iCOND) = {surfhomo[]};
  EndIf
EndIf




Mesh.Light = 0;

Recursive Color SteelBlue{
  Surface{surfiron[]};
}
Recursive Color Red{
  Surface{sw[], surfhomo[]};
}
Recursive Color Gold{
  Surface{is[]};
}
Recursive Color SkyBlue{
  Surface{surfgap1[],surfgap2[],surfgap3[]};
}

//Hide{Point{Point '*'};}
