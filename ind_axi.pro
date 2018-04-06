// Time domain + fine model:
//   gmsh ind_axi.geo -setnumber Flag_HomogenisedModel 0  -2 -o fine.msh -v 3
//   getdp ind_axi -msh fine.msh -setnumber Flag_HomogenisedModel 0 -setnumber Flag_FD 0 -setnumber Flag_MH 0 -setnumber Flag_SrcType 2 -setnumber NbT 1 -sol Analysis

// saved res file for additional post-processing... ind_axi_refTD_ph88.1.res

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
// Harmonic balance + fine model
// gmsh ind_axi.geo -setnumber Flag_HomogenisedModel 0  -2 -o fine.msh -v 3
// getdp_mh ind_axi -msh fine.msh -setnumber Flag_HomogenisedModel 0 -setnumber Flag_FD 1 -setnumber Flag_MH 1 -setnumber NbHars 6 -setnumber Flag_SrcType 4 -sol Analysis_HB

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
// Harmonic balance + homog:
//   gmsh ind_axi.geo -setnumber Flag_HomogenisedModel 1 -setnumber Flag_HalfModel 1 -2 -o homog.msh
//   getdp_mh ind_axi -msh homog.msh -setnumber Flag_HomogenisedModel 1 -setnumber Flag_FD 1 -setnumber Flag_SrcType 2 -setnumber NbHars 1 -sol Analysis_HB

// Resistance (Ohms): p 'SF_s0.dat' u 1:2 w l lw 2, 'SH.dat' u 1:2 w p ps 2 pt 6 lw 2
// Inductance (mH):   p 'SF_s0.dat' u 1:($3*1e3) w l lw 2, 'SH.dat' u 1:($3*1e3) w p ps 2 pt 6 lw 2


Include "ind_axi_dat.pro"
Include "BH.pro"

//la = 0.64909 ; // Fill factor
la = Fill;

// new files, finer mesh, max X=8
file_ZSkinRe  = Sprintf("coeff/pI_RS_la%.2g.dat", la);
file_ZSkinIm  = Sprintf("coeff/qI_RS_la%.2g.dat", la);
file_NuProxRe = Sprintf("coeff/qB_RS_la%.2g.dat", la);
file_NuProxIm = Sprintf("coeff/pB_RS_la%.2g.dat", la);

file_TDcoeffs = Sprintf("coeff/CoeffsTD_RS_la%.2g.pro", la);
Include file_TDcoeffs;


DirRes = "res/";
po = "{Output/";

addtofile="_aaa";

DefineConstant[
  visu = {0, Choices{0, 1}, AutoCheck 0,
    Name StrCat[mfem,"Visu/Real-time visualization"], Highlight "LightPink"}

  Flag_IronCore    = {1, Choices{0,1}, Name StrCat[mfem, "3Core/Iron core?"], Highlight Str[col1]}
  Flag_CoreWithGap = {0, Choices{0,1}, Name StrCat[mfem, "3Core/Iron core with gap?"],
    Highlight Str[col1], Visible Flag_IronCore}

  Flag_NL  = {0, Choices{0,1}, Name StrCat[mfem, "3Core/Nonlinear bh-curve?"],
    Highlight Str[col1],Visible Flag_IronCore}
  Flag_NewtonRaphson  = {1, Choices{0,1}, Name StrCat[mfem, "3Core/Nonlinear solver/0Newton Raphson?"],
    Highlight Str[col1],Visible Flag_IronCore && Flag_NL}

  Nb_max_iter = {30, Name StrCat[mfem, "3Core/Nonlinear solver/Max. num. iterations"],
    Visible Flag_NL, Highlight Str[col3]}
  iter_max = Nb_max_iter
  relaxation_factor = {1., Name StrCat[mfem, "3Core/Nonlinear solver/Relaxation factor"],
    Visible Flag_NL, Highlight Str[col3]}
  stop_criterion = {1e-4, Name StrCat[mfem, "3Core/Nonlinear solver/Stopping criterion"],
    Visible Flag_NL, Highlight Str[col3]}

  Flag_FD = { 0, Choices{0,1}, Name StrCat[mfem,"001Frequency domain analysis?"], Highlight "AliceBlue" }
  Flag_TD = !Flag_FD
  Flag_MH = { 0, Choices{0,1}, Name StrCat[mfem,"002Multi-harmonic analysis?"], Visible Flag_FD, Highlight "AliceBlue"}
  // Flag_MH Needed for running fine model with MH

  Flag_ECeffects = { 1, Choices{0,1}, Name StrCat[mfem,"2Eddy current effects (skin and proximity)?"], Visible Flag_FD, Highlight "AliceBlue"}

  // Frequency domain ==> sinusoidal source
  // Time domain ==> any source, only fine reference model
  Flag_SrcType = { 4, Choices{0="sinusoidal", 1="triangular", 2="pulse", 3="two sin", 4="PWM"},
                   Name StrCat[mfem, "Source/003Shape"], Highlight Str[col2],
    Help Str["- Use 'sinusoidal' for frequency domain",
      "- Use any source in time domain or multi-harmonic cases"] }

  Flag_ImposedVoltage = {1, Choices{0,1}, Name StrCat[mfem,"Source/001Imposed Voltage?"]}
  Flag_Circuit = {Flag_ImposedVoltage, Choices{0,1}, Name StrCat[mfem,"Source/002Use circuit"],
    ReadOnly (Flag_ImposedVoltage==1)}

  NbHars = { 1,
    Name StrCat[mfem,"Source/4Harmonics in HBFE"],
    Help Str["- Use '1' for classical harmonic analysis with GetDP",
      "- Use '>1' for HBFE with GetDP_mh"], Visible Flag_MH }

  Flag_imposedRr   = {!(Flag_SrcType==4), Choices{0,1},
    Name StrCat[mfem, "Source/1Imposed reduced frequency X"], Highlight Str[col1]}
  Flag_imposedFreq = !Flag_imposedRr

  ORDER = {1, Choices{
      1="order 1",
      2="order 2",
      3="order 3",
      4="order 4"},
           Name StrCat[mfem,"002Homog. approx. (TD)"], Highlight Str[col2],
    Visible (Flag_TD==1 && Flag_HomogenisedModel==1),
    ReadOnly !(Flag_TD==1 && Flag_HomogenisedModel==1)}
];



If(Flag_imposedRr) // Reduced frequency
  DefineConstant[
    Rr = {4, Min 0.1, Max 5, Step 0.1, Name StrCat[mfem,"Source/1Reduced frequency"], ReadOnly 0, Highlight "Ivory"}
    delta = {Rc/Rr, Name StrCat[mfem,"Source/2Skin depth [m]"], ReadOnly 1, Highlight "LightGrey"}
    Freq  = {1/(delta*delta*mu0*SigmaCu*Pi), Name StrCat[mfem,"Source/3Frequency [Hz]"], ReadOnly 1, Highlight "LightGrey"}
  ];
Else
  DefineConstant[
    Freq = { (Flag_SrcType==4)?1000:50000, Min 0.1, Max 500e3, Name StrCat[mfem,"Source/3Frequency [Hz]"], ReadOnly 0, Highlight "Ivory"}
    delta = {1/Sqrt[mu0*SigmaCu*Freq*Pi], Name StrCat[mfem,"Source/2Skin depth [m]"], ReadOnly 1, Highlight "LightGrey"}
    Rr = {Rc/delta, Name StrCat[mfem,"Source/1Reduced frequency"], ReadOnly 1, Highlight "LightGrey"}
  ];
EndIf

//Printf("Flag_imposedRr %g Freq %g", Flag_imposedRr, Freq);

DefineConstant[
  Omega  = 2*Pi*Freq,
  Period = 1./Freq,

  //------------------------------------------------------
  // Pulse width modulation
  //------------------------------------------------------
  // mf= frequency index (ratio between two frequencies)
  // ma= amplitude index (ratio between the amplitudes, up to 1 or bigger (overmodulation))

  fm = Freq  // fundamental of the output voltage
  fc = 20000+Freq // switching or carrier frequency
  mf = fc/fm
  ma = 0.8

  Ac = 1
  Am = ma*Ac

  Pm = Period
  Pc = 1/fc
  //------------------------------------------------------

  // Time stepping variables...
  NbT = {1, Name StrCat[mfem,"Source/5Number of periods to simulate"],  Highlight "Ivory", Visible Flag_TD}//5
  NbSteps = {(Flag_SrcType==4) ? Ceil[mf*250]:120, Name StrCat[mfem,"Source/6Number of steps per period"],
    Highlight "Ivory", Visible Flag_TD, ReadOnly (Flag_SrcType==4)}
  tinit = {0.,  Name StrCat[mfem,"Source/7Initial time"], Highlight "Ivory", Visible Flag_TD}
  tend = {NbT*Period ,  Name StrCat[mfem,"Source/8Final time"], ReadOnly 1, Highlight "LightGrey", Visible Flag_TD}
  deltat = {Period/NbSteps, Name StrCat[mfem,"Source/9Step size"], ReadOnly 1, Highlight "LightGrey", Visible Flag_TD}
  thetav = 1.

  Flag_AdaptiveTime = {0, Choices {0,1}, Name StrCat[mfem,"Source/9Adaptive time?"]}
  // -------- from helix.pro (superconductors)
  // tol_abs = {1e-9, Name "Input/Solver/3Absolute tolerance on nonlinear residual"},
  // tol_rel = {1e-6, Name "Input/Solver/3Relative tolerance on nonlinear residual"},
  // --------
  rel_tol = stop_criterion
  abs_tol = 1e-3*stop_criterion

  dt_init = deltat/4
  dt_min = 1e-18
  dt_max = deltat*2
  TimeIntMethod = "Gear_5" // Euler, Trapezoidal, Gear_2, Gear_3, Gear_4, Gear_5, Gear_6

  ExtGmsh = ".pos"
  ExtGnuplot = ".dat"
];

// List of time points to be met [s]
// Breakpoints = {};
Breakpoints = {tinit:tend:2*deltat};


Group{
  Air  = Region[{AIR}];
  Insulation = Region[{INSULATION}];

  If(Flag_IronCore)
    Iron = Region[{IRON}];
    If (Flag_CoreWithGap)
       Air  += Region[{AIRGAP}];
    Else
       Iron += Region[{AIRGAP}];
    EndIf

  Else
    Iron = Region[{}];
    Air  += Region[{IRON, AIRGAP}];
  EndIf

  OuterBoundary = Region[{OUTBND}]; // including symmetry

  Winding =  Region[{}] ;

  DomainCC = Region[{Air, Insulation, Iron}] ;

  SymFactor = Flag_HalfModel ? 2.:1. ; //half inductor with axisymmetry

  nbturns = (Flag_HomogenisedModel==0) ? NbrCond/SymFactor : 1  ; // number of turns

  DefineConstant[// ****
    Flag_stranded = {0, Choices{0,1}, Name StrCat[mfem,"Winding/Each turn stranded?"], Visible (Flag_HomogenisedModel==0) }
  ];

  If (Flag_HomogenisedModel==0) // Fine case
    For iF In {1:nbturns}
      Turn~{iF} = Region[{(iCOND+iF-1)}] ;
      Winding  += Region[{(iCOND+iF-1)}] ;
    EndFor

    If(!Flag_stranded) // Normal case
      DomainC = Region[{Winding}] ;
      DomainS = Region[{}] ;
    Else  // Neglecting eddy current effects in fine model
      DomainC = Region[{}] ;
      DomainS = Region[{Winding}] ;
    EndIf
  EndIf

  If (Flag_HomogenisedModel==1) //Homogenised case
    Turn~{1} = Region[{iCOND}] ;
    Winding = Region[{iCOND}];
    DomainC = Region[{}] ;
    DomainS = Region[{Winding}] ;
  EndIf

  DomainCC += Region[{DomainS}] ;

  If(Flag_NL)
    Domain_Lin = Region[{Air, Insulation, Winding}];
    Domain_Lin_NoS = Region[{Air, Insulation}];
    Domain_NonLin = Region[{Iron}];
  Else
    Domain_Lin = Region[{Air, Insulation, Winding, Iron}];
    Domain_Lin_NoS = Region[{Air, Insulation, Iron}];
    Domain_NonLin = Region[{}];
  EndIf


  Domain = Region[{DomainC, DomainCC}] ;

  //--------------------------------------------------------
  //--------------------------------------------------------

  DummyDomain = #123456;
  // Groups related to the circuit
  Input = # 12345 ;
  iZH    = 10000 ;
  iLH    = 20000 ; // Zskin in homog. coil
  iLHp   = 30000 ;

  For k In {1:ORDER}
    Zh~{k} = Region[{(iZH+k)}];
    Lh~{k} = Region[{(iLH+k)}];
  EndFor
  For k In {2:ORDER}
    Lhp~{k-1} = Region[{(iLHp+k-1)}];
  EndFor

  Resistance_Cir  = Region[{}];
  If(Flag_FD && Flag_HomogenisedModel==1)
    Resistance_Cir += Region[{Zh~{1}}]; // Frequency domain
  EndIf
  If(Flag_TD && Flag_HomogenisedModel==1) // Time domain
    For k In {1:ORDER}
      Resistance_Cir += Region[{Zh~{k}}];
    EndFor
  EndIf

  Inductance_Cir  = Region[{}];
  If(Flag_TD && Flag_HomogenisedModel==1) // Time domain
    For k In {1:ORDER}
      Inductance_Cir += Region[{Lh~{k}}];
    EndFor
    For k In {2:ORDER}
      Inductance_Cir += Region[{Lhp~{k-1}}];
    EndFor
  EndIf

  Capacitance1_Cir = Region[ {} ] ;
  Capacitance2_Cir = Region[ {} ] ;
  Capacitance_Cir  = Region[ {Capacitance1_Cir, Capacitance2_Cir} ] ;
  Diode_Cir  = Region[ {} ] ;

  SourceV_Cir = Region[ {Input} ] ;
  SourceI_Cir = Region[ {} ] ;

  DomainZ_Cir = Region[ {Resistance_Cir, Inductance_Cir, Capacitance_Cir,
                         Diode_Cir} ] ;

  DomainSource_Cir = Region[ {SourceV_Cir, SourceI_Cir} ] ;
  DomainZt_Cir = Region[ {DomainZ_Cir, DomainSource_Cir} ] ;

}


Function {
  CoefGeo = 2*Pi*SymFactor ; // axisymmetry + symmetry factor

  sigma[#{Winding}] = SigmaCu ;
  sigma[#{Air, Insulation, Iron}] = 0.;
  rho[] = 1/sigma[];

  nu[#{Air, Insulation}] = nu0;

  If (!Flag_NL)
    nu[#{Iron}]   = nu0/1000;
  Else
    nu[ #{Iron} ] = nu_3kW[$1] ;
    h[ #{Iron} ]  = h_3kW[$1];
    dhdb_NL[ #{Iron} ]= dhdb_3kW_NL[$1] ;
    dhdb[ #{Iron} ]   = dhdb_3kW[$1] ;
  EndIf

  //----------------------------------------------
  // Parameters for MH-integration
  MH_SamplesRHS = 16 ; MH_SamplesJac = 16 ; FreqOffSet = 0 ;
  RelaxFactors1 = LinSpace[1, 0.1, 10];
  RelaxFactors2 = LogSpace[0,-4,20];
  TestAllFactors = 0;
  //----------------------------------------------

  //==================================================================

  Val_EE = 1.; // 200 (nonlinear case)

  If(Flag_FD) // Frequency domain
    FSinusoidal[] = Complex_MH[1,0]{Freq} ; //Cos F_Cos_wt_p[]{2*Pi*Freq, 0};
  Else // Time domain
    FSinusoidal[] = Complex_MH[1,0]{Freq} ; //Sin F_Sin_wt_p[]{2*Pi*Freq, 0}; //Complex_MH[0,1]{Freq}
  EndIf

  FSaw[] = // starting at zero
  ( F_Period[$Time]{Period} <    Period/4. )?     4.* F_Period[$Time]{Period}/Period :
  ( F_Period[$Time]{Period} < 3.*Period/4. )?  2.-4.* F_Period[$Time]{Period}/Period :
  ( F_Period[$Time]{Period} <    Period    )? -4.+4.* F_Period[$Time]{Period}/Period :
  4.* F_Period[$Time]{Period}/Period ;


  FSaw2[] = // starting at one
  ( F_Period[$Time]{Period} <    Period/2. )?  1-4.* F_Period[$Time]{Period}/Period :
                                              -3+4.* F_Period[$Time]{Period}/Period ;

  FPulse1[]= (F_Period[$Time]{Period} < Period/2.) ? 1. : -1.;
  // shift of a quarter of a period
  // shift = Period/4 ;

  // shift = Period/4 - 0.235741582372815e-6 ; // Trying to reduce the transient
  // shift = Period/4 + 0.227409987348189e-6 ;
  shift = Period/4 + 0.019085210124540e-5;
  // shift = Period/4 + 0.194101454849494e-6;

  //FPulse2[]= ($Time < shift) ? 1. : ((F_Period[$Time-shift]{Period} < Period/2)? -1. : 1.);
  FPulse2[]= (F_Period[$Time]{Period} < shift || F_Period[$Time]{Period} > (Period/2+shift)) ? 1. : -1.;



  //------------------------------------------------------
  // Pulse width modulation
  //------------------------------------------------------
  FSaw_fc[] = // starting at zero
  ( F_Period[$Time]{Pc} <    Pc/4. )?     4.* F_Period[$Time]{Pc}/Pc :
  ( F_Period[$Time]{Pc} < 3.*Pc/4. )?  2.-4.* F_Period[$Time]{Pc}/Pc :
  ( F_Period[$Time]{Pc} <    Pc    )? -4.+4.* F_Period[$Time]{Pc}/Pc :
  4.* F_Period[$Time]{Pc}/Pc ;

  FCarrier[] = Ac * FSaw_fc[]; // with the frequency of the carrier


  // additional phase shift for reducing the transient

  // Lf1 = 0.002727400387396;
  // tau = Lf1/Rdc; // original defined in data file

  phase_m = 88.1*Pi/180; // Zero crossing (determined via interpolation in Matlab)
  shift_t = phase_m/(2*Pi*Freq);
  Printf("Phase shift %g [rad/s]; shift_t %g [s]", phase_m, shift_t);

  FV_out[]   = Am * F_Sin_wt_p[]{2*Pi*Freq, phase_m}; // modulator
  Vdc = 100;
  If (ma<=1)
    Fpwm[] = (FV_out[] >= FCarrier[]) ? -Vdc/2 : Vdc/2;
  Else
    Fpwm[] = (FV_out[] < FCarrier[])  ? -Vdc/2 : Vdc/2;
  EndIf

  //------------------------------------------------------
  allFreqs() = {};// Frequency spectrum
  allOmegas()= {};// pulsation
  list_EE()  = {};// Amplitudes different harmonics (HBFE)

  If(Flag_SrcType==0)
    If(NbHars!=1)
      Printf("*** Warning: sinusoidal source ===> NbHars forced to 1, not %g", NbHars);
      NbHars = 1;
    EndIf
    Fct_Src[] = FSinusoidal[];
    f~{1} = Freq ;
    Om~{1} = 2*Pi*f~{1} ;

    allFreqs() += f~{1};
    allOmegas() += Om~{1};
    list_EE()  = {1.,0.};
  EndIf

  If(Flag_SrcType==1)
    Fct_Src[] = FSaw2[];

    For k In {1:NbHars}
      f~{k} = Freq*(2*k-1);
      allFreqs() += f~{k};
      Om~{k} = 2*Pi*f~{k} ;
      allOmegas() += Om~{k};
      //list_EE() += {4*(1-(-1)^(2*k-1))/((2*k-1)*Pi)^2, 0.} ;
      list_EE() += {8./((2*k-1)*Pi)^2, 0.} ;
    EndFor
  EndIf

  /*
  // coefficients of harmonic components given by python
  fft_pulse() = ListFromFile[ "coef_fft_pulse_f2000.txt" ];
  nn = #fft_pulse(); // Size of list
  freq_content() = fft_pulse({0:nn-1:2}); // First column in file
  coefs() = fft_pulse({1:nn-1:2}); // Second column in file
  */

  If(Flag_SrcType==2)
    Fct_Src[] = FPulse2[];
    For k In {1:NbHars}
      f~{k} = Freq*(2*k-1);
      allFreqs() += f~{k};

      Om~{k} = 2*Pi*f~{k} ;
      allOmegas() += Om~{k};

      list_EE() += {4*(-1)^(k-1)/(2*k-1)/Pi, 0.} ; // cos at different frequencies
      //list_EE() += {coefs(k-1), 0.}; // This works too! (from python)
    EndFor
  EndIf

  If(Flag_SrcType==3)
    If(NbHars!=2)
      Printf("*** Warning: source comprises only two harmonics ===> NbHars forced to 2 not %g", NbHars);
      NbHars = 2;
    EndIf

    a1 = 1. ;
    a2 = 0.5;
    Fct_Src[] = a1*Complex_MH[1,0]{Freq} + a2*Complex_MH[1,0]{3*Freq};

    f~{1} = Freq;
    allFreqs() += f~{1};
    Om~{1} = 2*Pi*f~{1} ;
    allOmegas() += Om~{1};

    list_EE() += {a1, 0.}; //addition of two sin signals

    If(NbHars>=2)
      If(NbHars>2)
        Printf("Warning: this source has only two harmonics: %g>2", NbHars);
      EndIf
      NbHars = 2;
      f~{2} = 3*Freq;
      allFreqs() += f~{2};
      Om~{2} = 2*Pi*f~{2} ;
      allOmegas() += Om~{2};
      list_EE()  += {a2, 0.}; //addition of two sin signals
    EndIf
  EndIf

  If(Flag_SrcType==4)
    Fct_Src[] = Fpwm[];
    // coefficients of harmonic components
    fft_pwm() = ListFromFile[ Sprintf("coef_fft_pwm_fm%g_fc%g_ma%g_with_phaseshift_sorted.txt", fm, fc, ma) ];
    nn = #fft_pwm();
    // Printf("Size of list fft_pwm() = %g", nn);

    freq_content() = fft_pwm({0:nn-1:3}); // First column in file contains frequency
    coefs() = fft_pwm({1:nn-1:3}); // Second column in file
    phas()  = fft_pwm({2:nn-1:3}); // Third column in file

    If (NbHars>nn)
      NbHars = nn ;
      Printf("===> Warning: PWM harmonic data available for NbHars=%g",nn);
    EndIf
    For k In {1:NbHars}
      f~{k} = freq_content(k-1);
      allFreqs() += f~{k};

      Om~{k} = 2*Pi*f~{k} ;
      allOmegas() += Om~{k};

      // Harmonics considered in matlab
      // g0 = 1;    g1 = 1:6 ; g2 = 1:12;
      // g3 = 1:19; g4 = 1:27; g5 = 1:36;
      // Fi_sam = a_sam.*exp(1i*phi_sam); % voltage from harmonics
      list_EE() += {coefs(k-1)*Cos[phas(k-1)], coefs(k-1)*Sin[phas(k-1)]} ; // cos, sin at different frequencies
    EndFor
  EndIf

  //==================================================================

  // Homogenization coefficients: round conductor & square packing
  // Frequency domain
  skin_rhor_list() = ListFromFile[ file_ZSkinRe ];
  skin_rhoi_list() = ListFromFile[ file_ZSkinIm ];
  prox_nur_list()  = ListFromFile[ file_NuProxRe ];
  prox_nui_list()  = ListFromFile[ file_NuProxIm ];

  skin_rhor[] = InterpolationLinear[$1]{ skin_rhor_list() };
  skin_rhoi[] = InterpolationLinear[$1]{ skin_rhoi_list() };

  prox_nur[]  = InterpolationLinear[$1]{ prox_nur_list() } ;
  prox_nui[]  = InterpolationLinear[$1]{ prox_nui_list() } ;

  // both expressions for nu are equivalent (1 harmonic), generalisation hereafter
  // nu[Winding] = Complex[ prox_nur[Rr]*nu0 , prox_nui[Rr]*Fill*SigmaCu*Rc^2/4*Omega];
  // nu[Winding] = nu0*Complex[prox_nur[Rr], prox_nui[Rr]*Fill*Rr^2/2];

  // Complex impedance -- skin effect (1 harmonic), generalisation hereafter
  // Zskin[] = 1/SymFactor*Complex[ skin_rhor[Rr]*Rdc, 2*Pi*Freq*skin_rhoi[Rr]*mu0*Len/(8*Pi*Fill)];

  For k In {1:NbHars}
    delta~{k} = 1/Sqrt[mu0*SigmaCu*f~{k}*Pi];
    Rr~{k} = Rc/delta~{k} ;
    nu_re~{k}[] = prox_nur[Rr~{k}];
    nu_im~{k}[] = prox_nui[Rr~{k}]*Fill*Rr~{k}^2/2;

    zskin_re~{k}[] = skin_rhor[Rr~{k}] * Rdc ; // Rdc and Len are computed for all turns (no symmetry)
    zskin_im~{k}[] = 2*Pi*f~{k} * skin_rhoi[Rr~{k}]*mu0*Len/(8*Pi*Fill) ;
    Printf("HB ===> freq_%g %g Hz, reduced freq Rr_%g %g; delta_%g %g mm", k, f~{k}, k, Rr~{k}, k, delta~{k}*1e3);
  EndFor

  // Multiharmonic nu, all values together
  str_nu = "";
  For k In {1:NbHars}
    If(k < NbHars)

      str_nu = StrCat[str_nu, Sprintf("$nur~{%g}, $nui~{%g}, ", k, k)];
    Else
      str_nu = StrCat[str_nu, Sprintf("$nur~{%g}, $nui~{%g} ", k, k)];
    EndIf
  EndFor
  str = StrCat["nu_mh[] = nu0 * Complex[ ", str_nu, " ];"];
  Parse[str];

  str_nu = "";
  For k In {1:NbHars}
    If(k < NbHars)
      str_nu = StrCat[str_nu, Sprintf("Om~{%g}*$nui~{%g}, -$nur~{%g}, ", k, k, k)];
    Else
      str_nu = StrCat[str_nu, Sprintf("Om~{%g}*$nui~{%g}, -$nur~{%g} ", k, k, k)];
    EndIf
  EndFor
  str = StrCat["nuOm_mh[] = nu0 * Complex[ ", str_nu, " ];"];
  Parse[str];

  str_z = "";
  For k In {1:NbHars}
    If(k < NbHars)
      str_z = StrCat[str_z, Sprintf("$zr~{%g}, $zi~{%g}, ", k, k)];
    Else
      str_z = StrCat[str_z, Sprintf("$zr~{%g}, $zi~{%g} ", k, k)];
    EndIf
  EndFor
  str = StrCat["Zskin_mh[] =  1/SymFactor * Complex[ ", str_z, " ];"];
  Parse[str];

  If(Flag_HomogenisedModel==0)
    nu[Winding] = nu0 ;
  Else
  If (Flag_ECeffects==1) // Default
      nu[Winding] = nu_mh[] ; // Valid also for NbHar=1
    Else
      nu[Winding] = nu0 ; // No EC effects
    EndIf
  EndIf

  //==================================================================
  // Initialise to zero
  // zeros() = {} ;
  // For k In {1:NbHars}
  //   zeros() += {0.,0.} ;
  //   select_har~{k}() = {};
  // EndFor
  // // Select harmonic to do harmonic to time
  // For k In {1:NbHars}
  //   select_har~{k}() += zeros();
  //   select_har~{k}(2*k-2) = 1.;
  //   select_har~{k}(2*k-1) = 0.;
  //   har~{k}[] = Complex[]{select_har~{k}()};
  // EndFor

  //==================================================================
  // Time domain
  // Prox effect
  For k In {1:ORDER}
    Ca~{k}[] = nu0 * Ca_Sq~{ORDER}~{k} ;
    Cb~{k}[] = SigmaCu * Fill * Rc^2/4 * Cb_Sq~{ORDER}~{k} ;
  EndFor
  For k In {2:ORDER}
    Cc~{k-1}[] = SigmaCu * Fill * Rc^2/4 * Cc_Sq~{ORDER}~{k-1} ;
  EndFor
  // Skin effect
  For k In {1:ORDER}
    Sa~{k}[] = Rdc/1 * Sa_Sq~{ORDER}~{k} ; // Rdc/0.9768
    Sb~{k}[] = Rdc/1 * SigmaCu * mu0* Rc^2/(8*Fill) * Sb_Sq~{ORDER}~{k} ; //8*Fill
  EndFor
  For k In {2:ORDER}
    Sc~{k-1}[] = Rdc/1 * SigmaCu * mu0 * Rc^2/(8*Fill) * Sc_Sq~{ORDER}~{k-1} ; //8*Fill
  EndFor


  //==================================================================




  If (Flag_FD && Flag_ECeffects==1) // Default
    Zskin[] = Zskin_mh[] ;    // Valid also for NbHar=1
  Else
    If(Flag_ECeffects==0) // No EC effects
      //Zskin0[] = 1/SymFactor*Complex[ Rdc, 2*Pi*Freq*mu0*Len/(8*Pi*Fill)];
      Zskin[] = 1/SymFactor*Rdc ;
      Lskin[] = 1/SymFactor*mu0*Len/(8*Pi*Fill) ;
    Else //Time domain
      For k In {1:ORDER}
        Zskin~{k}[] = 1/SymFactor * Sa~{k}[];
        Lskin~{k}[] = 1/SymFactor * Sb~{k}[];
      EndFor
      For k In {2:ORDER}
        Lskin_p~{k-1}[] = 2*1/SymFactor * Sc~{k-1}[]; //2
      EndFor

      Zskin[] = Zskin~{1}[];
      Lskin[] = Lskin~{1}[];
    EndIf
  EndIf

  If(Flag_ECeffects==0)
    addtofile=StrCat[addtofile,"_noECs"];
  EndIf

  //==================================================================
  // Auxiliary functions for post-processing (valid with only one freq)
  // nuOm[] = Complex[ Omega * Im[nu[]], -Re[nu[]] ];
  // kkk[] =  skin_rhor[Rr] / Fill /SigmaCu ;

  // Generalisation to mh case
  // nuOm[#{Air, Insulation, Iron}] = Complex[ 0., -nu0 ];

  list_1i() = {} ;
  For k In {1:NbHars}
    list_1i() += {0.,1.} ;
  EndFor

  nuOm[#{Air, Insulation}] = -nu[]*Complex[]{list_1i()};
  nuOm[#{Iron}] = -nu[$1]*Complex[]{list_1i()};

  If (Flag_ECeffects==1) // Default
    nuOm[#{Winding}] = nuOm_mh[];
    kkk[] =  SymFactor*(Re[Zskin_mh[]]/Rdc)/ Fill/ SigmaCu;
  Else
    nuOm[#{Winding}] = -nu0*Complex[]{list_1i()};
    kkk[] =  SymFactor/ Fill/ SigmaCu;
  EndIf



  DefineFunction[
    Resistance, Inductance, Capacitance
  ];

  Ns[] = (Flag_HomogenisedModel==0) ? 1 : NbrCond/SymFactor ;
  Sc[] =  SurfaceArea[]{iCOND};

  If (Flag_HomogenisedModel==1)
    // Accounting for eddy currents (homogenization)
    If(Flag_FD)
      Resistance[Zh~{1}] = Zskin[] ;
    EndIf

    If(Flag_TD)
      For k In {1:ORDER}
        Resistance[Zh~{k}] = Zskin~{k}[] ;
        Inductance[Lh~{k}] = Lskin~{k}[] ;
      EndFor
      For k In {2:ORDER}
        Inductance[Lhp~{k-1}] = Lskin_p~{k-1}[] ;
      EndFor
    EndIf
  EndIf

  // List of nodes related to circuit
  N1() = {1:nbturns};   // Node 1 for each turn
  N2() = {2:nbturns+1}; // Node 2 for each turn

  relaxation_function[] = ($Iteration<Nb_max_iter/2) ? relaxation_factor : relaxation_factor/10 ;

}


Constraint {

  { Name MVP_2D ;
    Case {
      { Region OuterBoundary ; Type Assign ;  Value 0. ; }
    }
  }

  // Massive/stranded conductor constraints
  { Name Current_2D ;
    Case {
      If(Flag_Circuit==0)
        If(Flag_HomogenisedModel==0 && Flag_MH == 0)
          { Region Winding ; Value Val_EE; TimeFunction Fct_Src[] ; }
        Else
          { Region Winding ; Value Val_EE * Complex[]{list_EE()} ; }
        EndIf
      EndIf
    }
  }

  { Name Voltage_2D ;
    Case{
    }
  }

  { Name Voltage_Cir ;
    Case {
      If(Flag_Circuit && Flag_ImposedVoltage)
      //If(Flag_HomogenisedModel==0 && Flag_MH == 0)
        If(Flag_MH == 0)
          { Region Input ; Value Val_EE; TimeFunction Fct_Src[] ; }
        Else
          { Region Input ; Value Val_EE * Complex[]{list_EE()} ; }
        EndIf
      EndIf
    }
  }
  { Name Current_Cir ;
    Case {
      If(Flag_Circuit && !Flag_ImposedVoltage)
        //If(Flag_HomogenisedModel==0 && Flag_MH == 0)
        If(Flag_MH == 0)
          { Region Input ; Value -Val_EE; TimeFunction Fct_Src[] ; }
          Else
          { Region Input ; Value -Val_EE * Complex[]{list_EE()} ; }
        EndIf
      EndIf
    }
  }

  { Name ElectricalCircuit ; Type Network ;
    Case Circuit1 { // Common to fine and homogenised models
      If(Flag_HomogenisedModel==0)
        { Region Input ;  Branch {N1(0), N2(nbturns-1)} ; }
      EndIf
      If(Flag_HomogenisedModel==1 && Flag_FD)
          { Region Input  ; Branch {777, N2(nbturns-1)} ; }
          { Region Zh~{1} ; Branch {777, N1(0)}; } // Complex impedance
      EndIf
      If(Flag_HomogenisedModel==1 && Flag_TD)
        { Region Input ; Branch {777, N2(nbturns-1)} ; }
        If(ORDER==1)
          { Region Zh~{1}; Branch {777, 800}; }
          { Region Lh~{1}; Branch {800, N1(0)}; }
        EndIf
        If(ORDER==2)
          { Region Zh~{1} ; Branch {777, 800}; }
          { Region Lh~{1} ; Branch {800, 801}; }

          { Region Lhp~{1}; Branch {801, N1(0)}; }

          { Region Zh~{2} ; Branch {801, 802}; }
          { Region Lh~{2} ; Branch {802, N1(0)}; }
        EndIf
        If(ORDER==3)
          { Region Zh~{1} ; Branch {777, 800}; }
          { Region Lh~{1} ; Branch {800, 801}; }

          { Region Lhp~{1}; Branch {801, N1(0)}; }

          { Region Zh~{2} ; Branch {801, 802}; }
          { Region Lh~{2} ; Branch {802, 803}; }

          { Region Lhp~{2}; Branch {803, N1(0)}; }

          { Region Zh~{3} ; Branch {803, 804}; }
          { Region Lh~{3} ; Branch {804, N1(0)}; }
        EndIf
        If(ORDER==4)
          { Region Zh~{1} ; Branch {777, 800}; }
          { Region Lh~{1} ; Branch {800, 801}; }

          { Region Lhp~{1}; Branch {801, N1(0)}; }

          { Region Zh~{2} ; Branch {801, 802}; }
          { Region Lh~{2} ; Branch {802, 803}; }

          { Region Lhp~{2}; Branch {803, N1(0)}; }

          { Region Zh~{3} ; Branch {803, 804}; }
          { Region Lh~{3} ; Branch {804, 805}; }

          { Region Lhp~{3}; Branch {805, N1(0)}; }

          { Region Zh~{4} ; Branch {805, 806}; }
          { Region Lh~{4} ; Branch {806, N1(0)}; }
        EndIf
      EndIf

      For k In {0:nbturns-1} // list indexes start at 0
        { Region Turn~{k+1} ; Branch {N1(k), N2(k)} ; }
      EndFor
    }
  }




}

//-----------------------------------------------------------------------------

Jacobian {
  { Name Vol ; Case { { Region All ; Jacobian VolAxiSqu ; } } }
  { Name Sur ; Case { { Region All ; Jacobian SurAxi ; } } }
}

Integration {
  { Name II ; Case {
      { Type Gauss ; Case {
          { GeoElement Triangle ;    NumberOfPoints 4 ; }
          { GeoElement Quadrangle  ; NumberOfPoints 4 ; }
        } }
    } }
}

//-----------------------------------------------------------------------------

FunctionSpace {
  /*
  { Name Hcurl_a_2D ; Type Form1P ;
    BasisFunction {
      { Name se ; NameOfCoef ae ; Function BF_PerpendicularEdge ;
        Support Domain ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef ae  ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
    }
  }
  */

  { Name Hcurl_a_2D ; Type Form1P ; // Split for TD + homog: subspace isolating unknowns
    BasisFunction {
      { Name se ; NameOfCoef ae ; Function BF_PerpendicularEdge ;
        Support Domain ; Entity NodesOf[ All, Not Winding ] ; }
      { Name seh ; NameOfCoef aeh ; Function BF_PerpendicularEdge ;
        Support Domain ; Entity NodesOf[ Winding ] ; }
    }

    SubSpace {
      { Name aH ; NameOfBasisFunction {seh}; }
    }

    Constraint {
      { NameOfCoef ae  ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
      { NameOfCoef aeh  ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
    }
  }


  { Name Hregion_u_2D ; Type Form1P ; // Gradient of Electric scalar potential (2D)
    BasisFunction {
      { Name sr ; NameOfCoef ur ; Function BF_RegionZ ;
        Support DomainC ; Entity DomainC ; }
    }
    GlobalQuantity {
      { Name U ; Type AliasOf        ; NameOfCoef ur ; }
      { Name I ; Type AssociatedWith ; NameOfCoef ur ; }
    }
    Constraint {
      { NameOfCoef U ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef I ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }


  { Name Hregion_i_2D ; Type Vector ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_RegionZ ;
        Support DomainS ; Entity DomainS ; }
    }
    GlobalQuantity {
      { Name Is ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Us ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Us ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef Is ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }


  // For circuit equations
  { Name Hregion_Z ; Type Scalar ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_Region ;
        Support DomainZt_Cir ; Entity DomainZt_Cir ; }
    }
    GlobalQuantity {
      { Name Iz ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Uz ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Uz ; EntityType Region ; NameOfConstraint Voltage_Cir ; }
      { NameOfCoef Iz ; EntityType Region ; NameOfConstraint Current_Cir ; }
    }
  }

  // Time domain basis functions with homogenisation
  For i In {2:ORDER}
    { Name Hcurl_a~{i} ; Type Form1P ;
      BasisFunction {
        { Name se ; NameOfCoef ae ; Function BF_PerpendicularEdge ;
          Support Winding ; Entity NodesOf[ All ] ; }
      }
      Constraint {
        { NameOfCoef ae ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
      }
    }
  EndFor

  For k In {2:ORDER}
    { Name Hb~{k} ; Type Vector;
      BasisFunction {
        { Name sex ; NameOfCoef aex ; Function BF_VolumeX ; Support Winding ; Entity VolumesOf[ All ] ; }
        { Name sey ; NameOfCoef aey ; Function BF_VolumeY ; Support Winding ; Entity VolumesOf[ All ] ; }
        { Name sez ; NameOfCoef aez ; Function BF_VolumeZ ; Support Winding ; Entity VolumesOf[ All ] ; }
      }
    }
  EndFor

}


Formulation {
  { Name MagDyn_a ; Type FemEquation ;
    Quantity {
       { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }

       { Name ur ; Type Local  ; NameOfSpace Hregion_u_2D  ; }
       { Name I  ; Type Global ; NameOfSpace Hregion_u_2D[I] ; }
       { Name U  ; Type Global ; NameOfSpace Hregion_u_2D[U] ; }

       { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
       { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
       { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }

       { Name Uz ; Type Global ; NameOfSpace Hregion_Z [Uz] ; }
       { Name Iz ; Type Global ; NameOfSpace Hregion_Z [Iz] ; }
    }

    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ]  ;
        In Domain_Lin ; Jacobian Vol ; Integration II ; }

      Galerkin { [ nu[{d a}] * Dof{d a} , {d a} ]  ;
        In Domain_NonLin ; Jacobian Vol ; Integration II ; }

      If(Flag_NewtonRaphson)
        Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
          In Domain_NonLin ; Jacobian Vol ; Integration II ; }
      EndIf

      Galerkin { DtDof [ sigma[] * Dof{a} , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { [ sigma[] * Dof{ur}/CoefGeo , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }

      Galerkin { DtDof [ sigma[] * Dof{a} , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { [ sigma[] * Dof{ur}/CoefGeo , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      GlobalTerm { [ Dof{I}, {U} ] ; In DomainC ; }

      Galerkin { [ -Ns[]/Sc[] * Dof{ir}, {a} ] ;
        In DomainS ; Jacobian Vol ; Integration II ; }
      Galerkin { DtDof [ CoefGeo*Ns[]/Sc[] * Dof{a}, {ir} ] ; // CoefGeo?
        In DomainS ; Jacobian Vol ; Integration II ; }

      Galerkin { [ Ns[]/Sc[] / sigma[] * Ns[]/Sc[]* Dof{ir} , {ir} ] ; // resistance term
        In DomainS ; Jacobian Vol ; Integration II ; }
      //GlobalTerm { [ Rdc * Dof{Is} , {Is} ] ; In DomainS ; } // or this resistance term
      GlobalTerm { [ Dof{Us}/CoefGeo , {Is} ] ; In DomainS ; }

      If(Flag_Circuit)
        // Circuit equations
        GlobalTerm { NeverDt[ Dof{Uz}                , {Iz} ] ; In Resistance_Cir ; }
        GlobalTerm { NeverDt[ Resistance[] * Dof{Iz} , {Iz} ] ; In Resistance_Cir ; }
        GlobalTerm { DtDof [ Inductance[]  * Dof{Iz} , {Iz} ] ; In Inductance_Cir ; }

	GlobalTerm { [ 0. * Dof{Iz} , {Iz} ] ; In DomainZt_Cir ; }

        GlobalEquation {
          Type Network ; NameOfConstraint ElectricalCircuit ;
          { Node {I};  Loop {U};  Equation {I};  In DomainC ; }
          { Node {Is}; Loop {Us}; Equation {Us}; In DomainS ; }
          { Node {Iz}; Loop {Uz}; Equation {Uz}; In DomainZt_Cir ; }
        }
      EndIf
    }
  }

  { Name MagDyn_a_Homog ; Type FemEquation ; // FD multiharmonic or not
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }

      { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
      { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
      { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }

      { Name Uz ; Type Global ; NameOfSpace Hregion_Z [Uz] ; }
      { Name Iz ; Type Global ; NameOfSpace Hregion_Z [Iz] ; }
    }

    Equation {

      Galerkin { [ nu[] * Dof{d a} , {d a} ]  ;
        In Domain_Lin ; Jacobian Vol ; Integration II ; }

      If(Flag_NL)
        If(!Flag_MH)
          Galerkin { [ nu[{d a}] * Dof{d a} , {d a} ]  ;
            In Domain_NonLin ; Jacobian Vol ; Integration II ; }
          Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
            In Domain_NonLin ; Jacobian Vol ; Integration II ; }
        EndIf
        If(Flag_MH)
          Galerkin { [ MHTransform[ h[{d a}] ]{MH_SamplesRHS} , {d a} ]  ;
            In Domain_NonLin ; Jacobian Vol ; Integration II ; }
          Galerkin { JacNL[ MHJacNL[ dhdb[{d a}] ]{MH_SamplesJac, FreqOffSet} * Dof{d a} , {d a} ] ;
            In Domain_NonLin ; Jacobian Vol ; Integration II ; }
        EndIf
      EndIf

      Galerkin { [ -1/AreaCell * Dof{ir}, {a} ] ; // +++ check sign
        In DomainS ; Jacobian Vol ; Integration II ; }
      Galerkin { DtDof [ 1/AreaCell * Dof{a}, {ir} ] ; // CoefGeo? ++++ Ruth
        In DomainS ; Jacobian Vol ; Integration II ; }

      // Galerkin { [ Zskin[]/AreaCell * Dof{ir} , {ir} ] ; // Rdc => Zskin
      //  In DomainS ; Jacobian Vol ; Integration II ; }  // Resistance term (use one of the two!)
      // GlobalTerm { [ Zskin[] * Dof{Is} ,  {Is} ] ; In DomainS ; } // Or in Resistance_Cir
      GlobalTerm { [ Dof{Us}/CoefGeo, {Is} ]     ; In DomainS ; }

      // Circuit equations
      If(Flag_Circuit)
        GlobalTerm { NeverDt[ Dof{Uz}                , {Iz} ] ; In Resistance_Cir ; }
        GlobalTerm { NeverDt[ Resistance[] * Dof{Iz} , {Iz} ] ; In Resistance_Cir ; }// Zskin considered here

        GlobalTerm { [ 0. * Dof{Iz}, {Iz} ] ; In DomainZt_Cir ; }

        GlobalEquation {
          Type Network ; NameOfConstraint ElectricalCircuit ;
          { Node {Is};  Loop {Us};  Equation {Us}; In DomainS ; }
          { Node {Iz};  Loop {Uz};  Equation {Uz}; In DomainZt_Cir ; }
        }
      EndIf
    }
  }

  { Name MagDyn_a_Homog_time ; Type FemEquation ;
    Quantity {
      { Name a   ; Type Local  ; NameOfSpace Hcurl_a_2D ; }
      { Name a_1 ; Type LocalQuantity ; NameOfSpace Hcurl_a_2D[aH] ; }

      For i In {2:ORDER}
        { Name a~{i}  ; Type Local ; NameOfSpace Hcurl_a~{i} ; }
      EndFor

      { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
      { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
      { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }

      { Name Uz ; Type Global ; NameOfSpace Hregion_Z [Uz] ; }
      { Name Iz ; Type Global ; NameOfSpace Hregion_Z [Iz] ; }
    }

    Equation {

      Galerkin { [ nu[] * Dof{d a} , {d a} ]  ;
        In Region[{Domain_Lin_NoS}] ; Jacobian Vol ; Integration II ; }

      For i In {1:ORDER}
        Galerkin { [ Ca~{i}[] * Dof{d a~{i} } , {d a~{i}} ]  ;
          In DomainS; Jacobian Vol; Integration II ; }
        Galerkin { DtDof [ Cb~{i}[] * Dof{d a~{i}} , {d a~{i} } ] ;
          In DomainS; Jacobian Vol; Integration II ; }
      EndFor
      For i In {2:ORDER}
        Galerkin { DtDof [ Cc~{i-1}[] * Dof{d a~{i-1}} , {d a~{i}} ] ;
          In DomainS ; Jacobian Vol ; Integration II ; }
        Galerkin { DtDof [ Cc~{i-1}[] * Dof{d a~{i}} , {d a~{i-1}} ] ;
          In DomainS ; Jacobian Vol ; Integration II ; }
      EndFor

      If(Flag_NL) // Never MH
          Galerkin { [ nu[{d a}] * Dof{d a} , {d a} ]  ;
            In Domain_NonLin ; Jacobian Vol ; Integration II ; }
          Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
            In Domain_NonLin ; Jacobian Vol ; Integration II ; }
      EndIf

      Galerkin { [ -1/AreaCell * Dof{ir}, {a} ] ; // +++ check sign
        In DomainS ; Jacobian Vol ; Integration II ; }
      Galerkin { DtDof [ 1/AreaCell * Dof{a}, {ir} ] ;
        In DomainS ; Jacobian Vol ; Integration II ; }

      //Zskin0[] = 1/SymFactor*Complex[ Rdc, 2*Pi*Freq*mu0*Len/(8*Pi*Fill)];

      // Galerkin { [ 1/AreaCell/sigma[]*1/AreaCell * Dof{ir} , {ir} ] ; // Rdc => Zskin
      //  In DomainS ; Jacobian Vol ; Integration II ; }  // Resistance term (use one of the two!)
      // GlobalTerm { [ Zskin[] * Dof{Is} ,  {Is} ] ; In DomainS ; } // Or in Resistance_Cir
      GlobalTerm { [ Dof{Us}/CoefGeo, {Is} ]     ; In DomainS ; }

      // Circuit equations
      If(Flag_Circuit)
        GlobalTerm { NeverDt[ Dof{Uz}                , {Iz} ] ; In Resistance_Cir ; }
        GlobalTerm { NeverDt[ Resistance[] * Dof{Iz} , {Iz} ] ; In Resistance_Cir ; } // Zskin considered here

        GlobalTerm { [ Dof{Uz}                        , {Iz} ] ; In Inductance_Cir ; }// Lskin considered here
        GlobalTerm { DtDof [ Inductance[]  * Dof{Iz}  , {Iz} ] ; In Inductance_Cir ; }
        GlobalTerm { [ 0. * Dof{Iz}, {Iz} ] ; In DomainZt_Cir ; }

        GlobalEquation {
          Type Network ; NameOfConstraint ElectricalCircuit ;
          { Node {Is};  Loop {Us};  Equation {Us}; In DomainS ; }
          { Node {Iz};  Loop {Uz};  Equation {Uz}; In DomainZt_Cir ; }
        }
      EndIf
    }
  }

  { Name MagDyn_a_Homog_time_b ; Type FemEquation ; // NOT Valid in multiharmonic
    Quantity {
      { Name a   ; Type Local  ; NameOfSpace Hcurl_a_2D ; }
      For i In {2:ORDER}
        { Name b~{i}  ; Type Local ; NameOfSpace Hb~{i} ; }
      EndFor

      { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
      { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
      { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }

      { Name Uz ; Type Global ; NameOfSpace Hregion_Z [Uz] ; }
      { Name Iz ; Type Global ; NameOfSpace Hregion_Z [Iz] ; }
    }

    Equation {

      Galerkin { [ nu[] * Dof{d a} , {d a} ]  ;
        In Region[{Domain_Lin_NoS}] ; Jacobian Vol ; Integration II ; }

      // (ORDER >= 1)
      Galerkin { [ Ca_1[] * Dof{d a}, {d a} ]  ;  In DomainS ; Jacobian Vol ; Integration II ; }
      Galerkin { DtDof [ Cb_1[] * Dof{d a} , {d a} ] ;  In DomainS ; Jacobian Vol ; Integration II ; }

      /*
      If(ORDER>=2)
        Galerkin { [ Ca~{2}[] * Dof{b~{2}} , {b~{2}} ]  ; In DomainS; Jacobian Vol; Integration II ; }
        Galerkin { DtDof [ Cb~{2}[] * Dof{b~{2}} , {b~{2} } ] ; In DomainS; Jacobian Vol; Integration II ; }
        Galerkin { DtDof [ Cc_1[] * Dof{d a}, {b_2} ] ; In DomainS ; Jacobian Vol ; Integration II ; }
        Galerkin { DtDof [ Cc_1[] * Dof{b_2}, {d a} ] ; In DomainS ; Jacobian Vol ; Integration II ; }
      EndIf
      If(ORDER>=3)
        Galerkin { [ Ca~{3}[] * Dof{b~{3}} , {b~{3}} ]  ; In DomainS; Jacobian Vol; Integration II ; }
        Galerkin { DtDof [ Cb~{3}[] * Dof{b~{3}} , {b~{3} } ] ; In DomainS; Jacobian Vol; Integration II ; }
        Galerkin { DtDof [ Cc_2[] * Dof{b_2}, {b_3} ] ; In DomainS ; Jacobian Vol ; Integration II ; }
        Galerkin { DtDof [ Cc_2[] * Dof{b_3}, {b_2} ] ; In DomainS ; Jacobian Vol ; Integration II ; }
      EndIf
      If(ORDER>=4)
        Galerkin { [ Ca~{4}[] * Dof{b~{4}} , {b~{4}} ]  ; In DomainS; Jacobian Vol; Integration II ; }
        Galerkin { DtDof [ Cb~{4}[] * Dof{b~{4}} , {b~{4} } ] ; In DomainS; Jacobian Vol; Integration II ; }
        Galerkin { DtDof [ Cc_3[] * Dof{b_3}, {b_4} ] ; In DomainS ; Jacobian Vol ; Integration II ; }
        Galerkin { DtDof [ Cc_3[] * Dof{b_4}, {b_3} ] ; In DomainS ; Jacobian Vol ; Integration II ; }
      EndIf
      */

      For i In {2:ORDER}
        Galerkin { [ Ca~{i}[] * Dof{b~{i}} , {b~{i}} ]  ; In DomainS; Jacobian Vol; Integration II ; }
        Galerkin { DtDof [ Cb~{i}[] * Dof{b~{i}} , {b~{i} } ] ; In DomainS; Jacobian Vol; Integration II ; }
      EndFor
      If (ORDER >= 2)
        Galerkin { DtDof [ Cc_1[] * Dof{d a}, {b_2} ] ;
          In DomainS ; Jacobian Vol ; Integration II ; }
        Galerkin { DtDof [ Cc_1[] * Dof{b_2}, {d a} ] ;
          In DomainS ; Jacobian Vol ; Integration II ; }
      EndIf
      For i In {3:ORDER}
        Galerkin { DtDof [ Cc~{i-1}[] * Dof{b~{i-1}} , {b~{i}} ] ;
          In DomainS ; Jacobian Vol ; Integration II ; }
        Galerkin { DtDof [ Cc~{i-1}[] * Dof{b~{i}} , {b~{i-1}} ] ;
          In DomainS ; Jacobian Vol ; Integration II ; }
      EndFor

      If(Flag_NL) // Never MH
          Galerkin { [ nu[{d a}] * Dof{d a} , {d a} ]  ;
            In Domain_NonLin ; Jacobian Vol ; Integration II ; }
          Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
            In Domain_NonLin ; Jacobian Vol ; Integration II ; }
      EndIf

      Galerkin { [ -1/AreaCell * Dof{ir}, {a} ] ; // +++ check sign
        In DomainS ; Jacobian Vol ; Integration II ; }
      Galerkin { DtDof [ 1/AreaCell * Dof{a}, {ir} ] ;
        In DomainS ; Jacobian Vol ; Integration II ; }

      // Galerkin { [ 1/AreaCell/(Fill*sigma[])*1/AreaCell * Dof{ir}, {ir} ] ; // Rdc => Zskin
      //  In DomainS ; Jacobian Vol ; Integration II ; }  // Resistance term (use one of the two!)
      // GlobalTerm { [ Zskin[] * Dof{Is} ,  {Is} ] ; In DomainS ; } // Or in Resistance_Cir
      // GlobalTerm { DtDof [ Lskin[]*Dof{Is} , {Is} ] ; In DomainS ; }    // total inductance

      GlobalTerm { [ Dof{Us}/CoefGeo, {Is} ]     ; In DomainS ; }

      // Circuit equations
      If(Flag_Circuit)
        GlobalTerm { NeverDt[ Dof{Uz}                , {Iz} ] ; In Resistance_Cir ; }
        GlobalTerm { NeverDt[ Resistance[] * Dof{Iz} , {Iz} ] ; In Resistance_Cir ; }// Zskin considered here

        GlobalTerm { [ Dof{Uz}                        , {Iz} ] ; In Inductance_Cir ; }
        GlobalTerm { DtDof [ Inductance[]  * Dof{Iz}  , {Iz} ] ; In Inductance_Cir ; }

        GlobalTerm { [ 0. * Dof{Iz}, {Iz} ] ; In DomainZt_Cir ; }

        GlobalEquation {
          Type Network ; NameOfConstraint ElectricalCircuit ;
          { Node {Is};  Loop {Us};  Equation {Us}; In DomainS ; }
          { Node {Iz};  Loop {Uz};  Equation {Uz}; In DomainZt_Cir ; }
        }
      EndIf
    }
  }



}


Resolution {
  { Name Analysis ;
    System {
      If(Flag_HomogenisedModel==0)
        If(Flag_FD) // Frequency domain
          { Name A ; NameOfFormulation MagDyn_a ; Type ComplexValue ; Frequency Freq ; }
        Else // Time domain
          { Name A ; NameOfFormulation MagDyn_a ; }
        EndIf
      EndIf
      If(Flag_HomogenisedModel==1)
        If(Flag_FD) // Frequency domain
          { Name A ; NameOfFormulation MagDyn_a_Homog ; Type ComplexValue ; Frequency Freq ; }
        Else // Time domain
        //{ Name A ; NameOfFormulation MagDyn_a_Homog_time ; } // aux BFs based on MVP
        { Name A ; NameOfFormulation MagDyn_a_Homog_time_b ; } // aux BFs directly b
        EndIf
      EndIf
    }
    Operation {
      CreateDir[DirRes];
      InitSolution[A]; SaveSolution[A];

      If(Flag_FD) // Frequency domain
// added 29 Mar 2018	    
		If(Flag_HomogenisedModel==1) // Frequency domain
          For k In {1:NbHars}
            Evaluate[
              $nur~{k} = nu_re~{k}[],
              $nui~{k} = nu_im~{k}[],
              $zr~{k}  = zskin_re~{k}[],
              $zi~{k}  = zskin_im~{k}[]
            ];
          EndFor
		EndIf
		
        If(Flag_imposedRr && Flag_SrcType==0 && NbHars==1)
          SetTime[Rr];
        EndIf
        Generate[A] ; Solve[A] ; SaveSolution[A] ;
        PostOperation [Map_local];
		If(Flag_HomogenisedModel==1) 
		  PostOperation [Get_global_Homog];
        Else  	  
          PostOperation [Get_global];
	    EndIf
      Else //Fine reference in time domain
        If(!Flag_AdaptiveTime)
          TimeLoopTheta[tinit, tend, deltat, thetav]{
            If(!Flag_NL)
              Generate[A] ; Solve[A] ;
            EndIf
            If(Flag_NL)
              //IterativeLoop[Nb_max_iter, stop_criterion, relaxation_function[]] {
              IterativeLoop[Nb_max_iter, stop_criterion, relaxation_factor] {
                GenerateJac[A]; SolveJac[A]; }
            EndIf
            SaveSolution[A] ;
            Test[ GetNumberRunTime[visu]{StrCat[mfem,"Visu/Real-time visualization"]} ]
            {
              If(Flag_HomogenisedModel==0)
                PostOperation[Map_local];
              Else
                PostOperation[Map_local_Homog_time];
              EndIf
            }
          }
        EndIf
        If(Flag_AdaptiveTime)
          TimeLoopAdaptive[ tinit, tend, dt_init, dt_min, dt_max, TimeIntMethod, List[Breakpoints],
              System { { A, rel_tol, abs_tol, MeanL1Norm} } ]
            // norm-type: L1Norm, MeanL1Norm, L2Norm, MeanL2Norm, LinfNorm
            {
              If(!Flag_NL)
                Generate[A]; Solve[A];
                Else
                IterativeLoop[Nb_max_iter, stop_criterion, relaxation_factor] {
                  GenerateJac[A] ; SolveJac[A] ; }
              EndIf
            }
            {
              Test[$Breakpoint >= 0]{
                SaveSolution[A];
                Test[ GetNumberRunTime[visu]{StrCat[mfem,"Visu/Real-time visualization"]} ]
                { PostOperation[Map_local]; }
                //PostOperation[Map_local];
              }
            }
        EndIf
      EndIf
    }
  }

  { Name Analysis_HB ; // Homogenisation with NbHar >=1
    System {
      If(Flag_HomogenisedModel==1)
        { Name A ; NameOfFormulation MagDyn_a_Homog ; Type ComplexValue; Frequency allFreqs() ; }
      Else
        { Name A ; NameOfFormulation MagDyn_a ; Type ComplexValue; Frequency allFreqs() ; }
      EndIf
    }
    Operation {
      CreateDir[DirRes];

      If(Flag_imposedRr && Flag_SrcType==0 && NbHars==1)
        SetTime[Rr];
      EndIf

      If(Flag_HomogenisedModel==1)
        For k In {1:NbHars}
          Evaluate[
            $nur~{k} = nu_re~{k}[],
            $nui~{k} = nu_im~{k}[],
            $zr~{k}  = zskin_re~{k}[],
            $zi~{k}  = zskin_im~{k}[]
          ];
        EndFor
      EndIf

      InitSolution[A];
      If(!Flag_NL)
        Generate[A]; Solve[A];
      Else
        IterativeLoop[Nb_max_iter, stop_criterion, relaxation_factor] {
          GenerateJac[A] ; SolveJac[A] ;
          //SolveJac_AdaptRelax [A, RelaxFactors1{}, TestAllFactors] ;
        }
      EndIf
      SaveSolution[A] ;
      //SaveSolutionExtendedMH [A, 1, "next"];

      If(Flag_HomogenisedModel==1)
        PostOperation [Map_local_Homog];
        PostOperation [Get_global_Homog];
      Else
      //PostOperation [Map_local_mh];
        PostOperation [Get_global_mh];
      EndIf
    }
  }

  //=====================================================
  // Experimental...
  { Name Analysis_fine_manual_NL ; // Fixed point
    System {
      If(Flag_FD) // Frequency domain
        { Name A ; NameOfFormulation MagDyn_a ; Type ComplexValue ; Frequency Freq ; }
        Else // Time domain
        { Name A ; NameOfFormulation MagDyn_a ; }
      EndIf
    }
    Operation {
      CreateDir[DirRes];
      InitSolution[A]; SaveSolution[A];

      If(Flag_FD) // Frequency domain
        If(Flag_imposedRr && Flag_SrcType==0 && NbHars==1)
          SetTime[Rr];
        EndIf
        Generate[A] ; Solve[A] ; SaveSolution[A] ;
        PostOperation [Map_local];
        PostOperation [Get_global];
      Else //Fine reference in time domain

        // set a runtime variable to count the number of linear system solves (to
        // compare the performance of adaptive vs. non-adaptive time stepping
        // scheme)
        Evaluate[ $syscount = 0 ];

        // enter implicit Euler time-stepping loop
        TimeLoopTheta[tinit, tend, deltat, thetav]{
          If(!Flag_NL)
            Generate[A] ; Solve[A] ;
          EndIf
          If(Flag_NL)
            // compute first solution guess and residual at step $TimeStep
            Generate[A]; Solve[A]; Evaluate[ $syscount = $syscount + 1 ];
            Generate[A]; GetResidual[A, $res0]; Evaluate[ $res = $res0, $iter = 0 ];
            Print[{$iter, $res, $res / $res0},
              Format "Residual %03g: abs %14.12e rel %14.12e"];

            // iterate until convergence
            While[$res > abs_tol && $res / $res0 > rel_tol &&
              $res / $res0 <= 1 && $iter < iter_max]{
              Solve[A]; Evaluate[ $syscount = $syscount + 1 ];
              Generate[A]; GetResidual[A, $res]; Evaluate[ $iter = $iter + 1 ];
              Print[{$iter, $res, $res / $res0},
                Format "Residual %g: abs %14.12e rel %14.12e"];
            }

            // save and visualize the solution if converged...
            Test[ $iter < iter_max && $res / $res0 <= 1 ]{
              SaveSolution[A];
              Test[ GetNumberRunTime[visu]{"Input/Solver/Visu"} ]{
                PostOperation[Map_local];
              }
              // increase the step if we converged sufficiently "fast"
              Test[ $iter < iter_max / 4 && $DTime < dt_max ]{
                Evaluate[ $dt_new = Min[$DTime * 1.5, dt_max] ];
                Print[{$dt_new},
                  Format "*** Fast convergence: increasing time step to %g"];
                SetDTime[$dt_new];
              }
            }
            // ...otherwise reduce the time step and try again
            {
              Evaluate[ $dt_new = $DTime / 3 ];
              Print[{$iter, $dt_new},
                Format "*** Non convergence (iter %g): recomputing with reduced step %g"];
              SetTime[$Time - $DTime];
              SetTimeStep[$TimeStep - 1];
              RemoveLastSolution[A];
              SetDTime[$dt_new];
            }
          EndIf //Flag_NL
        }// End TimeLoop
        Print[{$syscount}, Format "Total number of linear systems solved: %g"];
      EndIf //Flag_FD

    }// Operation
  } // Analysis_fine_manual_NL

}// Resolution

PostProcessing {
  //deltaI = 1.185565943533830e-04;

  { Name MagDyn_a ; NameOfFormulation MagDyn_a ; NameOfSystem A;
    PostQuantity {
      { Name a ; Value { Term { [ {a} ] ; In Domain ; Jacobian Vol  ;} } }
      { Name az ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name raz ; Value { Term { [ CompZ[{a}]*X[] ] ; In Domain ; Jacobian Vol  ;} } }

      { Name b ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name h ; Value { Term { [ nu[{d a}]*{d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name j ; Value { Term { [ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ; In DomainC ; Jacobian Vol ; } } }
      { Name jz ; Value {
          Term { [ CompZ[ -sigma[]*(Dt[{a}]+{ur}/CoefGeo) ] ] ; In DomainC ; Jacobian Vol ; }
          Term { [ CompZ[ {ir}/AreaCond ] ] ; In DomainS ; Jacobian Vol ; }
        } }

      { Name b2av ; Value { Integral { [ CoefGeo*{d a}/AreaCell ] ;
            In Domain ; Jacobian Vol  ; Integration II ; } } }


      { Name j2F ; Value { Integral {
            [ CoefGeo*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] ] ;
            In DomainC ; Jacobian Vol ; Integration II ; } } }

      { Name b2F ; Value { Integral { [ CoefGeo*nu[{d a}]*SquNorm[{d a}] ] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }

      { Name SoF ; Value {
          Integral { // generalisation MH case
            [ CoefGeo * (sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] + Complex[]{list_1i()}*nu[{d a}]*SquNorm[{d a}]) ] ;
            //[ CoefGeo * Complex[ sigma[] * SquNorm[(Dt[{a}]+{ur}/CoefGeo)], nu[]*SquNorm[{d a}] ] ] ;
            In Domain ; Jacobian Vol ; Integration II ; }
        } }//Complex power

      { Name U ; Value {
          Term { [ {U} ]   ; In DomainC ; }
          Term { [ {Us} ]   ; In DomainS ; }
          Term { [ {Uz} ]  ; In DomainZt_Cir ; }
        } }
      { Name I ; Value {
          Term { [ {I} ]   ; In DomainC ; }
          Term { [ {Is} ]   ; In DomainS ; }
          Term { [ {Iz} ]  ; In DomainZt_Cir ; }
        } }


      { Name MagEnergy ; Value {
          Integral { [ CoefGeo*nu[{d a}]*({d a}*{d a})/2 ] ;
            In Domain ; Jacobian Vol ; Integration II ; } } }
     }
  }


  { Name MagDyn_a_Homog ; NameOfFormulation MagDyn_a_Homog ;
    PostQuantity {
      { Name a   ; Value { Term { [ {a} ] ; In Domain ; Jacobian Vol  ;} } }
      { Name az  ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol  ;} } }
      { Name raz ; Value { Term { [ CompZ[{a}]*X[] ] ; In Domain ; Jacobian Vol  ;} } }

      { Name b ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }

      For k In {1:NbHars} // Selecting only one of the harmonics
        // { Name b~{k} ; Value { Term { [ har~{k}[]*{d a} ] ; In Domain ; Jacobian Vol ; } } }
        { Name nui~{k} ; Value { Term { Type Global; [ nu_im~{k}[] ] ; In DummyDomain ; } } }
        { Name nur~{k} ; Value { Term { Type Global; [ nu_re~{k}[] ] ; In DummyDomain ; } } }
        { Name nuiOm~{k} ; Value { Term { Type Global; [ Om~{k}*nu_im~{k}[] ] ; In DummyDomain ; } } }
      EndFor



      { Name j  ; Value { Term { [ -1/AreaCell*{ir} ] ; In DomainS ; Jacobian Vol ; } } }
      { Name jz ; Value { Term { [ CompZ[ -1/AreaCell*{ir} ] ] ; In DomainS ; Jacobian Vol ; } } }
      { Name int_jz ; Value { Integral { [ CompZ[ -1/AreaCell*{ir} ] ] ; In DomainS ; Jacobian Vol ; Integration II ; } } }

      { Name j2H ; Value { Integral { // Joule losses
            //[ SymFactor*CoefGeo*kkk[]*SquNorm[js[]] ] ;
            [ CoefGeo*(Re[{d a}*Conj[nuOm[]*{d a}]]+kkk[]*SquNorm[-1/AreaCell*{ir}]) ] ;
            In DomainS ; Jacobian Vol ; Integration II ; } } }

      // nuOm~{%g}[] = nu0*Complex[ Om~{%g}*$nui~{%g}, -$nur~{%g}];
      { Name SoH ; Value { Integral { // Complex power = Active power +j * Reactive power => S = P+j*Q
            //[ SymFactor*CoefGeo* ({d a}*Conj[nuOm[]*{d a}] + kkk[]*SquNorm[js[]]) ] ;
            [ CoefGeo * ({d a}*Conj[nuOm[{d a}]*{d a}] + kkk[]*SquNorm[-1/AreaCell*{ir}]) ] ;
            In Domain ; Jacobian Vol ; Integration II ; } } } //Complex power

      { Name j2H_Zskin ; Value { // Joule Losses -- total skin losses
          Term { [ SymFactor*Zskin[] * SquNorm[{Iz}] ]; In Zh~{1}; }
        }
      }
      { Name Zskin ; Value { Term { Type Global; [ SymFactor*Zskin[] ] ; In DummyDomain ; } } }

      { Name U ; Value {
          Term { [ {Us} ]  ; In DomainS ; }
          Term { [ {Uz} ]  ; In DomainZt_Cir ; }
        } }
      { Name I ; Value {
          Term { [ {Is} ]  ; In DomainS ; }
          Term { [ {Iz} ]  ; In DomainZt_Cir ; }
        } }
      // For checking multi-harmonic case...
      // Conj, Re and Im act per frequency
      { Name conjI ; Value {
          Term { [ Conj[{Is}] ]  ; In DomainS ; }
          Term { [ Conj[{Iz}] ]  ; In DomainZt_Cir ; }
        } }
      { Name squnormI ; Value {
          Term { [ SquNorm[{Is}] ]  ; In DomainS ; }
          Term { [ SquNorm[{Iz}] ]  ; In DomainZt_Cir ; }
        } }
      { Name reI ; Value {
          Term { [ Re[{Is}] ]  ; In DomainS ; }
          Term { [ Re[{Iz}] ]  ; In DomainZt_Cir ; }
        } }
      { Name imI ; Value {
          Term { [ Im[{Is}] ]  ; In DomainS ; }
          Term { [ Im[{Iz}] ]  ; In DomainZt_Cir ; }
        } }

    }
  }

  // Time domain + homog
  { Name MagDyn_a_Homog_time ; NameOfFormulation MagDyn_a_Homog_time ;
    PostQuantity {
      { Name a  ; Value { Term { [ {a} ]           ; In Domain ; Jacobian Vol; } } }
      { Name az ; Value { Term { [ CompZ[{a}] ]    ; In Domain ; Jacobian Vol; } } }
      { Name raz; Value { Term { [ CompZ[{a}]*X[]] ; In Domain ; Jacobian Vol; } } }
      For i In {1:ORDER}
        { Name a~{i}  ; Value { Term { [ {a~{i}} ]        ; In DomainS ; Jacobian Vol; } } }
        { Name az~{i} ; Value { Term { [ CompZ[{a~{i}}] ] ; In DomainS ; Jacobian Vol; } } }
        { Name raz~{i}; Value { Term { [ CompZ[{a~{i}}]*X[] ] ; In DomainS ; Jacobian Vol; } } }
      EndFor

      { Name b ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name j  ; Value { Term { [ -1/AreaCell*{ir} ] ; In DomainS ; Jacobian Vol ; } } }
      { Name jz ; Value { Term { [ CompZ[ -1/AreaCell*{ir} ] ] ; In DomainS ; Jacobian Vol ; } } }

      { Name jI ; Value { Term { [ Rdc*{Is}^2 ]  ; In DomainS ; } } }
      { Name j2HH ; Value { // Joule Losses -- total proximity losses
          For i In {1:ORDER}
            Integral { [ CoefGeo * Cb~{i}[] * Dt[{d a~{i}}] * Dt[{d a~{i}}] ];
              In Winding ; Jacobian Vol ; Integration II ; }
          EndFor
          For i In {2:ORDER}
            Integral { [ CoefGeo * 2 * Cc~{i-1}[] * Dt[{d a~{i-1}}] * Dt[{d a~{i}}] ];
              In Winding ; Jacobian Vol ; Integration II ; }
          EndFor
        }
      }

      { Name j2HH_skin ; Value { // Joule Losses -- total skin losses
          If(Flag_FD)
              Term { [ SymFactor*Zskin[] * SquNorm[{Iz}] ]; In Zh~{1}; }
          EndIf
          If(Flag_TD)
            For k In {1:ORDER}
              Term { [ SymFactor*Zskin~{k}[] * SquNorm[{Iz}] ]; In Zh~{k}; }
            EndFor
          EndIf
        }
      }

      { Name j2F ; Value { // Joule losses
          //          Integral { [ CoefGeo*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] ] ;
          //  In DomainC ; Jacobian Vol ; Integration II ; }
          Integral { [ CoefGeo *SquNorm[1/AreaCell*{ir}]/(Fill*sigma[]) ] ;
            In Winding ; Jacobian Vol ; Integration II ; }
        } }

      { Name MagEnergy ; Value {
          For i In {1:ORDER}
            Integral { [ CoefGeo * nu0/2 * {d a~{i}}*{d a~{i}} ] ;
              In Winding ; Jacobian Vol ; Integration II ; }
          EndFor
        }
      }

       { Name U ; Value {
          Term { [ {Us} ]  ; In DomainS ; }
          Term { [ {Uz} ]  ; In DomainZt_Cir ; }
        } }
      { Name I ; Value {
          Term { [ {Is} ]  ; In DomainS ; }
          Term { [ {Iz} ]  ; In DomainZt_Cir ; }
        } }


    }
  }

  { Name MagDyn_a_Homog_time_b ; NameOfFormulation MagDyn_a_Homog_time_b ;
    PostQuantity {
      { Name a  ; Value { Term { [ {a} ]           ; In Domain ; Jacobian Vol; } } }
      { Name az ; Value { Term { [ CompZ[{a}] ]    ; In Domain ; Jacobian Vol; } } }
      { Name raz; Value { Term { [ CompZ[{a}]*X[]] ; In Domain ; Jacobian Vol; } } }
      For i In {2:ORDER}
        { Name b~{i}  ; Value { Term { [ {b~{i}} ]        ; In DomainS ; Jacobian Vol; } } }
      EndFor

      { Name b ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name j  ; Value { Term { [ -1/AreaCell*{ir} ] ; In DomainS ; Jacobian Vol ; } } }
      { Name jz ; Value { Term { [ CompZ[ -1/AreaCell*{ir} ] ] ; In DomainS ; Jacobian Vol ; } } }

      { Name jI ; Value { Term { [ Rdc*{Is}^2 ]  ; In DomainS ; } } }

      { Name j2HH ; Value { // Joule Losses -- total proximity losses
          Integral { [ CoefGeo * Cb~{1}[] * Dt[{d a}] * Dt[{d a}] ];
            In Winding ; Jacobian Vol ; Integration II ; }
          For i In {2:ORDER}
            Integral { [ CoefGeo * Cb~{i}[] * Dt[{b~{i}}] * Dt[{b~{i}}] ];
              In Winding ; Jacobian Vol ; Integration II ; }
          EndFor
          If(ORDER>1)
            Integral { [ CoefGeo * 2 * Cc~{1}[] * Dt[{d a}] * Dt[{b~{2}}] ];
              In Winding ; Jacobian Vol ; Integration II ; }
          EndIf
          For i In {3:ORDER}
            Integral { [ CoefGeo * 2 * Cc~{i-1}[] * Dt[{b~{i-1}}] * Dt[{b~{i}}] ];
              In Winding ; Jacobian Vol ; Integration II ; }
          EndFor
        }
      }

      { Name j2F ; Value { // Joule losses
          //          Integral { [ CoefGeo*sigma[]*SquNorm[(Dt[{a}]+{ur}/CoefGeo)] ] ;
          //  In DomainC ; Jacobian Vol ; Integration II ; }
          Integral { [ CoefGeo *SquNorm[1/AreaCell*{ir}]/(Fill*sigma[]) ] ;
            In Winding ; Jacobian Vol ; Integration II ; }
        } }

      { Name j2HH_skin ; Value { // Joule Losses -- total skin losses
          If(Flag_FD)
              Term { [ SymFactor*Zskin[] * SquNorm[{Iz}] ]; In Zh~{1}; }
          EndIf
          If(Flag_TD)
            For k In {1:ORDER}
              Term { [ SymFactor*Zskin~{k}[] * SquNorm[{Iz}] ]; In Zh~{k}; }
            EndFor
          EndIf
        }
      }

      { Name MagEnergy ; Value {
          Integral { [ CoefGeo * nu0/2 * {d a}*{d a} ] ;
              In Winding ; Jacobian Vol ; Integration II ; }
          For i In {2:ORDER}
            Integral { [ CoefGeo * nu0/2 * {b~{i}}*{b~{i}} ] ;
              In Winding ; Jacobian Vol ; Integration II ; }
          EndFor
        }
      }

       { Name U ; Value {
          Term { [ {Us} ]  ; In DomainS ; }
          Term { [ {Uz} ]  ; In DomainZt_Cir ; }
        } }
      { Name I ; Value {
          Term { [ {Is} ]  ; In DomainS ; }
          Term { [ {Iz} ]  ; In DomainZt_Cir ; }
        } }


    }
  }

}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

PostOperation Map_local UsingPost MagDyn_a {
  Print[ jz, OnElementsOf Region[{DomainC,DomainS}], File StrCat[DirRes, "jz", ExtGmsh], LastTimeStepOnly ] ;
  Print[ b,  OnElementsOf Domain,  File StrCat[DirRes, "b", ExtGmsh],  LastTimeStepOnly ] ;
  Print[ raz,OnElementsOf Domain,  File StrCat[DirRes, "a", ExtGmsh], LastTimeStepOnly ] ;
  // Print[ raz,OnElementsOf #{Air, Insulation, Winding},  File StrCat[DirRes, "a_in", ExtGmsh], LastTimeStepOnly ] ;
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=3;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
  // RangeType = 1; // Value scale range type (1=default, 2=custom, 3=per time step)
  // IntervalsType = 2; // Type of interval display (1=iso, 2=continuous, 3=discrete, 4=numeric)

  If(!Flag_Circuit)
    Print[ I, OnRegion Winding, Format TimeTable, >File Sprintf("res/I_s%g_f%g.dat", Flag_SrcType, Freq),
      SendToServer StrCat[po,"I [A]"], Color "Pink", LastTimeStepOnly];
    Print[ U, OnRegion Winding, Format TimeTable, >File Sprintf("res/U_s%g_f%g.dat", Flag_SrcType, Freq),
      SendToServer StrCat[po,"V [V]"], Color "LightGreen", LastTimeStepOnly];
  Else
    Print[ I, OnRegion Input, Format TimeTable, File >Sprintf("res/I_s%g_f%g.dat", Flag_SrcType, Freq),
      SendToServer StrCat[po,"I [A]"], Color "Pink", LastTimeStepOnly];
    Print[ U, OnRegion Input, Format TimeTable, File >Sprintf("res/U_s%g_f%g.dat", Flag_SrcType, Freq),
      SendToServer StrCat[po,"V [V]"], Color "LightGreen", LastTimeStepOnly];
  EndIf
  // Echo[ Str["View[PostProcessing.NbViews-1].LineWidth = 2;
  //            View[PostProcessing.NbViews-1].IntervalsType=2;
  //            View[PostProcessing.NbViews-1].NbIso = 25;"],
  //          File "res/optionVI.pos"];
}

PostOperation Get_global UsingPost MagDyn_a {
  Print[ j2F[ Winding ], OnGlobal, Format TimeTable, File > Sprintf("res/j2F_s%g_iron%g.dat", Flag_SrcType, Flag_IronCore)] ;// Joule losses
  Print[ SoF[ Domain ], OnGlobal, Format TimeTable,  File > Sprintf("res/SF_s%g_iron%g.dat", Flag_SrcType, Flag_IronCore),
         SendToServer "Output/Complex power" {0,1} ] ; // Complex power
}


PostOperation Get_allTS UsingPost MagDyn_a {
  // Print[ jz, OnElementsOf DomainC, File StrCat[DirRes, "jz", ExtGmsh] ] ;
  //Print[ b,  OnElementsOf Domain,  File StrCat[DirRes, "b", ExtGmsh]  ] ;
  /*
  Print[ raz,OnElementsOf Domain,  File StrCat[DirRes, "a", ExtGmsh] ] ;
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
  */

  If(!Flag_Circuit)
    Print[ I, OnRegion Winding, Format TimeTable, File Sprintf("res/I_s%g_f%g.dat", Flag_SrcType, Freq) ];
    Print[ U, OnRegion Winding, Format TimeTable, File Sprintf("res/U_s%g_f%g.dat", Flag_SrcType, Freq) ];
  Else
    Print[ I, OnRegion Input, Format TimeTable, File Sprintf("res/I_s%g_f%g.dat", Flag_SrcType, Freq),
      SendToServer StrCat[po,"I [A]"]{0}, Color "Pink"];
    Print[ U, OnRegion Input, Format TimeTable, File Sprintf("res/U_s%g_f%g.dat", Flag_SrcType, Freq),
      SendToServer StrCat[po,"V [V]"]{0}, Color "LightGreen"];
  EndIf

  Print[ j2F[Winding], OnGlobal, Format TimeTable, File Sprintf("res/jl_s%g_f%g.dat", Flag_SrcType, Freq) ] ; // Joule losses
  // Print[ MagEnergy[#{Winding,Insulation}], OnGlobal, Format TimeTable, File StrCat[DirRes,"ME",ExtGnuplot]] ;
}


NbT0 = NbT-1;
NbT1 = NbT;
list_steps()= (NbT0*NbSteps-1):(NbT1*NbSteps-1) ;
NbHars_out = 20;
PostOperation Get_global_Fourier UsingPost MagDyn_a{
  Print[ I, OnRegion Input, FourierTransform, Format FrequencyTable, TimeStep{list_steps()}, TimeToHarmonic NbHars_out,
    File Sprintf("res/Ifft_s%g_f%g.dat", Flag_SrcType, Freq)];
  Print[ U, OnRegion Input, FourierTransform, Format FrequencyTable, TimeStep{list_steps()}, TimeToHarmonic NbHars_out,
    File Sprintf("res/Ufft_s%g_f%g.dat", Flag_SrcType, Freq)];

  // Print[ j2F[Winding], OnGlobal, FourierTransform, Format FrequencyTable, TimeStep{list_steps()}, TimeToHarmonic NbHars_out,
  //  File Sprintf("res/jl_fft_s%g_f%g.dat", Flag_SrcType, Freq) ] ;

  // Generated files contain: Freq, abs(fft(quantity)), phase(fft(quantity))
}

PostOperation Map_local_mh UsingPost MagDyn_a {
  Print[ jz, OnElementsOf Region[{DomainC,DomainS}], File Sprintf("res/jz_s%g_f%g_mh%g.pos", Flag_SrcType, Freq, NbHars), LastTimeStepOnly ] ;
  Print[ b,  OnElementsOf Domain,  File Sprintf("res/b_s%g_f%g_mh%g.pos", Flag_SrcType, Freq, NbHars),  LastTimeStepOnly ] ;
  Print[ raz,OnElementsOf Domain,  File Sprintf("res/a_s%g_f%g_mh%g.pos", Flag_SrcType, Freq, NbHars), LastTimeStepOnly ] ;
  Print[ raz,OnElementsOf #{Air, Insulation, Winding},  File Sprintf("res/a_in_s%g_f%g_mh%g.pos", Flag_SrcType, Freq, NbHars), LastTimeStepOnly ] ;
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
  // RangeType = 1; // Value scale range type (1=default, 2=custom, 3=per time step)
  // IntervalsType = 2; // Type of interval display (1=iso, 2=continuous, 3=discrete, 4=numeric)
}

NbSteps2 = 2*NbSteps;

PostOperation Get_global_mh UsingPost MagDyn_a {// Fine model with multiharmonic approach
  //Print[ j2F[ Winding ], OnGlobal, Format TimeTable, File Sprintf("res/j2F_s%g_mh%g.dat", Flag_SrcType, NbHars)] ;// Joule losses
  Print[ SoF[ Winding ], OnGlobal, Format TimeTable,
    File Sprintf(StrCat["res/SF_s%g_f%g_mh%g",addtofile,".dat"], Flag_SrcType, Freq, NbHars)] ; // Complex power

  If(!Flag_Circuit)
    Print[ I, OnRegion Winding, Format TimeTable, File Sprintf("res/If_s%g_f%g_mh%g.dat", Flag_SrcType, Freq, NbHars) ];
    Print[ U, OnRegion Winding, Format TimeTable, File Sprintf("res/Uf_s%g_f%g.dat", Flag_SrcType, Freq, NbHars) ];
  Else
    Print[ I, OnRegion Input, Format TimeTable, File Sprintf(StrCat["res/If_s%g_f%g_mh%g",addtofile,".dat"], Flag_SrcType, Freq, NbHars) ];
    Print[ U, OnRegion Input, Format TimeTable, File Sprintf(StrCat["res/Uf_s%g_f%g_mh%g",addtofile,".dat"], Flag_SrcType, Freq, NbHars) ];

    Print[ U, OnRegion Input, Format TimeTable,
      File Sprintf(StrCat["res/Uf_s%g_f%g_mh%g_t",addtofile,".dat"], Flag_SrcType, Freq, NbHars), HarmonicToTime NbSteps2 ] ;
    Print[ I, OnRegion Input, Format TimeTable,
      File Sprintf(StrCat["res/If_s%g_f%g_mh%g_t",addtofile,".dat"], Flag_SrcType, Freq, NbHars), HarmonicToTime NbSteps2 ] ;
  EndIf
}


//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

PostOperation Map_local_Homog UsingPost MagDyn_a_Homog {
  // Print[ j,   OnElementsOf DomainS, File StrCat[DirRes,"jH_mh",ExtGmsh], FrequencyLegend ] ;
  Print[ b,  OnElementsOf Winding, File Sprintf("res/bH_mh%g.pos", NbHars), FrequencyLegend ] ;
  Print[ jz, OnElementsOf Winding, File Sprintf("res/jzH_mh%g.pos", NbHars), FrequencyLegend ] ;
  Print[ int_jz[Winding], OnGlobal, Format TimeTable,
    File Sprintf("res/int_jzH_mh%g.dat", NbHars) ] ;

  For k In {1:NbHars} // Selecting only one of the harmonics
    Print[ nui~{k}, OnRegion DummyDomain, Format Table, File > Sprintf("res/nui_%g.dat", NbHars)] ;
    Print[ nur~{k}, OnRegion DummyDomain, Format Table, File > Sprintf("res/nur_%g.dat", NbHars)] ;
    Print[ nuiOm~{k}, OnRegion DummyDomain, Format Table, File > Sprintf("res/nuiOm_%g.dat", NbHars)] ;
  EndFor

  Print[ raz, OnElementsOf Domain, File StrCat[DirRes,"aH_mh",ExtGmsh], FrequencyLegend ] ;
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
}



PostOperation Get_global_Homog UsingPost MagDyn_a_Homog {
  // Complex power: S = P + i*Q, P=active power, Q=reactive power
  // Print[ SoH[Domain], OnGlobal, Format Table, Format TimeTable] ;

  Print[ SoH[Domain], OnGlobal, Format TimeTable, //SendToServer "Output/Complex power Homog." {0:2*NbHars},
    File Sprintf(StrCat["res/SH_s%g_f%g_mh%g",addtofile,".dat"], Flag_SrcType, Freq, NbHars) ] ;

  Print[ j2H[Winding], OnGlobal, Format Table,
    File Sprintf(StrCat["res/j2H_s%g_f%g_mh%g",addtofile,".dat"], Flag_SrcType, Freq, NbHars) ] ;

  If(!Flag_Circuit)
    Print[ U, OnRegion Winding, Format Table, File Sprintf("res/U_s%g_f%g_mh%g.dat", Flag_SrcType, Freq, NbHars) ] ;
    Print[ I, OnRegion Winding, Format Table, File Sprintf("res/I_s%g_f%g_mh%g.dat", Flag_SrcType, Freq, NbHars) ] ;

    Print[ U, OnRegion Winding, Format TimeTable, File Sprintf("res/U_s%g_f%g_mh%g_t.dat", Flag_SrcType, Freq, NbHars), HarmonicToTime NbSteps2 ] ;
    Print[ I, OnRegion Winding, Format TimeTable, File Sprintf("res/I_s%g_f%g_mh%g_t.dat", Flag_SrcType, Freq, NbHars), HarmonicToTime NbSteps2 ] ;
  Else
    Print[ U, OnRegion Input, Format Table, File Sprintf(StrCat["res/U_s%g_f%g_mh%g",addtofile,".dat"], Flag_SrcType, Freq, NbHars) ];
    Print[ I, OnRegion Input, Format Table, File Sprintf(StrCat["res/I_s%g_f%g_mh%g",addtofile,".dat"], Flag_SrcType, Freq, NbHars) ];

    Print[ U, OnRegion Input, Format TimeTable,
      File Sprintf(StrCat["res/U_s%g_f%g_mh%g",addtofile,"_t.dat"], Flag_SrcType, Freq, NbHars), HarmonicToTime NbSteps2 ] ;
    Print[ I, OnRegion Input, Format TimeTable,
      File Sprintf(StrCat["res/I_s%g_f%g_mh%g",addtofile,"_t.dat"], Flag_SrcType, Freq, NbHars), HarmonicToTime NbSteps2 ] ;

    //
    Print[ j2H_Zskin, OnRegion Zh~{1}, Format TimeTable,
      File Sprintf(StrCat["res/j2H_zskin_s%g_f%g_mh%g",addtofile,".dat"], Flag_SrcType, Freq, NbHars) ] ;
  EndIf

  Print[ Zskin, OnRegion DummyDomain, Format TimeTable,
    File Sprintf(StrCat["res/zskin_s%g_f%g_mh%g",addtofile,".dat"], Flag_SrcType, Freq, NbHars) ] ;


}


//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

PostOperation Map_local_Homog_time UsingPost MagDyn_a_Homog_time {
  Print[ j,   OnElementsOf DomainS, File StrCat[DirRes,"jH_time",ExtGmsh], LastTimeStepOnly ] ;
  Print[ b,   OnElementsOf Domain, File StrCat[DirRes,"bH_time",ExtGmsh], LastTimeStepOnly ] ;
  Print[ raz, OnElementsOf Domain, File StrCat[DirRes,"aH_time",ExtGmsh], LastTimeStepOnly ] ;
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
}



PostOperation Get_global_Homog_time UsingPost MagDyn_a_Homog_time {
  Print[ j2HH[Winding], OnGlobal, Format Table, LastTimeStepOnly,
    File Sprintf(StrCat["res/j2H_time_s%g_f%g",addtofile,".dat"], Flag_SrcType, Freq) ] ;
  Print[ MagEnergy[Winding], OnGlobal, Format Table, LastTimeStepOnly,
    File Sprintf(StrCat["res/ME_time_s%g_f%g",addtofile,".dat"], Flag_SrcType, Freq) ] ;

  If(!Flag_Circuit)
    Print[ U, OnRegion Winding, Format Table, File Sprintf("res/Uh_time_s%g_f%g.dat", Flag_SrcType, Freq), LastTimeStepOnly ] ;
    Print[ I, OnRegion Winding, Format Table, File Sprintf("res/Ih_time_s%g_f%g.dat", Flag_SrcType, Freq), LastTimeStepOnly ] ;
  Else
    Print[ U, OnRegion Input, Format Table, File Sprintf(StrCat["res/Uh_time_s%g_f%g",addtofile,".dat"], Flag_SrcType, Freq), LastTimeStepOnly ];
    Print[ I, OnRegion Input, Format Table, File Sprintf(StrCat["res/Ih_time_s%g_f%g",addtofile,".dat"], Flag_SrcType, Freq), LastTimeStepOnly ];
  EndIf

}


PostOperation Get_allTS_Homog_time UsingPost MagDyn_a_Homog_time {
  /*
  Print[ j,   OnElementsOf DomainS, File StrCat[DirRes,"jH_time",ExtGmsh] ] ;
  Print[ b,   OnElementsOf Domain, File StrCat[DirRes,"bH_time",ExtGmsh] ] ;
  Print[ raz, OnElementsOf Domain, File StrCat[DirRes,"aH_time",ExtGmsh] ] ;
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
  */

  If(!Flag_Circuit)
    Print[ I, OnRegion Winding, Format TimeTable, File Sprintf(StrCat["res/Ih_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER) ];
    Print[ U, OnRegion Winding, Format TimeTable, File Sprintf(StrCat["res/Uh_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER) ];
  Else
    Print[ I, OnRegion Input, Format TimeTable, File Sprintf(StrCat["res/Ih_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER),
      SendToServer StrCat[po,"I [A]"]{0}, Color "Pink"];
    Print[ U, OnRegion Input, Format TimeTable, File Sprintf(StrCat["res/Uh_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER),
      SendToServer StrCat[po,"V [V]"]{0}, Color "LightGreen"];
  EndIf

  Print[ j2HH[Winding], OnGlobal, Format TimeTable, File Sprintf(StrCat["res/jlh_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER) ] ; // Joule losses (proximity)
  Print[ j2F[Winding],  OnGlobal, Format TimeTable, File Sprintf(StrCat["res/jlsh_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER) ] ; // Joule losses (skin)
  // Print[ MagEnergy[Winding], OnGlobal, Format TimeTable, File Sprintf("res/meh_time_s%g_f%g_o%g.dat", Flag_SrcType, Freq, ORDER) ] ;

}


PostOperation Get_allTS_Homog_time_b UsingPost MagDyn_a_Homog_time_b {
  /*
  Print[ j,   OnElementsOf DomainS, File StrCat[DirRes,"jH_time",ExtGmsh] ] ;
  Print[ b,   OnElementsOf Domain, File StrCat[DirRes,"bH_time",ExtGmsh] ] ;
  Print[ raz, OnElementsOf Domain, File StrCat[DirRes,"aH_time",ExtGmsh] ] ;
  Echo[ Str["View[PostProcessing.NbViews-1].Light=0;
             View[PostProcessing.NbViews-1].LineWidth = 2;
             View[PostProcessing.NbViews-1].RangeType=3;
             View[PostProcessing.NbViews-1].IntervalsType=1;
             View[PostProcessing.NbViews-1].NbIso = 25;"],
           File "res/option.pos"];
  */

  If(!Flag_Circuit)
    Print[ I, OnRegion Winding, Format TimeTable, File Sprintf(StrCat["res/Ih_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER) ];
    Print[ U, OnRegion Winding, Format TimeTable, File Sprintf(StrCat["res/Uh_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER) ];
  Else
    Print[ I, OnRegion Input, Format TimeTable, File Sprintf(StrCat["res/Ih_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER),
      SendToServer StrCat[po,"I [A]"]{0}, Color "Pink"];
    Print[ U, OnRegion Input, Format TimeTable, File Sprintf(StrCat["res/Uh_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER),
      SendToServer StrCat[po,"V [V]"]{0}, Color "LightGreen"];
  EndIf

  Print[ j2HH[Winding], OnGlobal, Format TimeTable,
    File Sprintf(StrCat["res/jlh_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER) ] ; // Joule losses (proximity)
  Print[ j2HH_skin, OnRegion Resistance_Cir, Format TimeTable,
    File Sprintf(StrCat["res/jlsh2_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER) ] ;

  // Print[ j2F[Winding],  OnGlobal, Format TimeTable, File Sprintf(StrCat["res/jlsh_time_s%g_f%g_o%g",addtofile,".dat"], Flag_SrcType, Freq, ORDER) ] ; // Joule losses (skin)

  // Print[ MagEnergy[Winding], OnGlobal, Format TimeTable, File Sprintf("res/meh_time_s%g_f%g_o%g.dat", Flag_SrcType, Freq, ORDER) ] ;

}




// This is only for Onelab
// Careful, it will run the last GetDP charged with the corresponding formulation
DefineConstant[
  R_ = {"Analysis", Name "GetDP/1ResolutionChoices", Visible 1, Closed 1}
  C_ = {"-solve -v 4 -v2 -bin", Name "GetDP/9ComputeCommand", Visible 1}
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 1}
];
DefineConstant[
  RH_ = {"Analysis_HB", Name "GetDP_MH/1ResolutionChoices", Visible 1, Closed 1}
  CH_ = {"-solve -v 4 -v2 -bin", Name "GetDP_MH/9ComputeCommand", Visible 1}
  PH_ = {"", Name "GetDP_MH/2PostOperationChoices", Visible 1}
];


//PrintConstants;
