! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
!   CALCULATION OF ELECTRON TRANSPORT COEFFICIENTS IN MAGNETIC FIELDS  *
!*      Remarks and suggestions are welcome. Please send them to        *
!*       Alexander Potekhin <palex@astro.ioffe.ru>                      *
!*   For theoretical background and references see:                     *
!*           http://www.ioffe.ru/astro/conduct/                         *
!* ---------------------------------------------------
!*   This code stems from conduct08.f v.12.05.2011.
!* Differences: 1. ion thermal conduction (CONDI) is disabled by default,
!*   because its implemented approximation might cause problems (2013);
!* 2. frequency of e-e collisions is not reduced at T << T_plasma (2016);
!* 3. relax.time of e-e collisions (TAUEE) is switched to asymptote at
!*   \theta=\sqrt{3} T_{pe}/T = Y > 1.d8, to avoid underflow (2016).
!* 4. criterion for inclusion of e-e scattering is modified 07.12.17.
!* 5. call of COUL99I fixed 10.03.18.
!* 6. CONDI fixed 10.09.19 (typo in Eq.(26) of Chugunov & Haensel 2007)
!* 7. 17.09.19: Updated thermal averaging to conform to Fortran 2018
!*      and to take finite size of nuclei into account;
!*      replaced Coulomb logarithm subroutines COUL01, COUL99I, COULAN
!*      by their improved versions COUL19, COUL18I, COULAN3
!*      to fix numerical stability in limiting cases;
!*      improved applicability condition for e-e-scattering correction
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **

!*   ----------------------   MAIN block   -----------------------      *
!*       This is auxiliary MAIN program for input/output purposes.      *
!*It may be commented out (then uncomment it in order to use it "as is")*
!*                You can change it or write your own.                  *
!*     Calculations are performed in the subroutine CONDUCT (below),    *
!*     where most of internal quantities are in the relativistic units  *
!*         (\hbar = m_e = c = 1)                                        *
!*     Subroutine CONDCONV converts the conductivities to CGS units.    *
!* In this MAIN program we also demonstrate the convertions to SI units *
!* for electrical conductivity and to cm^2/g for thermal opacity        *
! ---------------------------------------------------------------    *
 
C%C      implicit double precision (A-H), double precision (O-Z)
C%C      character KEY
C%C      parameter (CONVSIGM=8.99d9) ! conversion S/m (SI) to 1/s (CGS)
C%C      parameter (CONVOPAC=3.02242e-4)
C%C* conv.th.conductivity [erg cm^{-1} s^{-1} K^{-1}] to opac. [cm^2/g]
C%C      write(*,'('' Charge and atomic mass of ions (Z,A): ''$)')
C%C      read*,Zion,CMI
C%C      write(*,'('' Impurity parameter (effective Z): ''$)')
C%C      read*,Zimp
C%C      print*,'Note: if B is nonzero, then ee-collisions are ignored'
C%C      write(*,'('' Magnetic field B (in Gauss): ''$)')
C%C      read*,B
C%C   50 continue
C%C      write(*,'('' Temperature T6 (in 10^6 K): ''$)')
C%C      read*,T6
C%C  100 continue
C%C      write(*,'('' Density of ions rho (in g/cm^3): ''$)')
C%C      read*,RHO
C%C      call CONDCONV(T6,RHO,B,Zion,CMI,Zimp, ! input
C%C     *  SIGMA,CKAPPA,QJ,SIGMAT,CKAPPAT,QJT,SIGMAH,CKAPPAH,QJH)
!C%C*   ----------------------   OUTPUT:
C%C      T=T6*1.D6
C%C      if (B.gt.0.) then
C%C         write(*,113)
C%C         write(*,111) RHO,T,SIGMA,CKAPPA,QJ,
C%C     *      SIGMAT,CKAPPAT,QJT,
C%C     *      SIGMAH,CKAPPAH,QJH
C%C      else
C%C         write(*,114)
C%C         write(*,111) RHO,T,SIGMA,CKAPPA,QJ
C%C      endif
C%C      OPAClg=dlog10(CONVOPAC/CKAPPA/RHO*T**3)
C%C      write(*,'('' log_{10}(opacity[cm^2/g]) ='',F8.3)') OPAClg
C%C      SIGMAlg=dlog10(SIGMA/CONVSIGM)
C%C      write(*,'('' log_{10}(el.conductivity[S/m] ='',F8.3)') SIGMAlg
C%C      write(*,'('' New density? (Y/N) ''$)')
C%C      read(*,'(A)') KEY
C%C      if (KEY.ne.'n'.and.KEY.ne.'N') goto 100
C%C      write(*,'('' New temperature? (Y/N) ''$)')
C%C      read(*,'(A)') KEY
C%C      if (KEY.ne.'n'.and.KEY.ne.'N') goto 50
C%C      stop
C%C  111 format(1PE12.4,E10.3,3(2X,3E10.2))
C%C  113 format(27X,'l o n g i t u d i n a l',
C%C     *10X,'t r a n s v e r s e',10X,'o f f - d i a g o n a l'/
C%C     *'      rho       T  ',3('        sigma     kappa      J_Q'))
C%C  114 format(20X,'electron-electron and electron-ion',
C%C     *10X,'electron-ion only'/
C%C     *'      rho       T          sigma     kappa      J_Q')
C%C      end



!*   ---------------------------------------------------------------    *
      subroutine CONDCONV(T6,RHO,B,Zion,CMI,Zimp,
     *  SIGMA,CKAPPA,QJ,SIGMAT,CKAPPAT,QJT,SIGMAH,CKAPPAH,QJH)
!* This subroutine serves as interface for conversion of the input
!* into relativistic units and the output into CGS units
!* Input: T6 - temperature in MK (10^6 K)
!*        RHO - density in g/cm^3
!*        B - magnetic field in Gauss
!*        Zion - charge number of plasma ions
!*        CMI - mass number of plasma ions
!*        Zimp - impurity parameter Z_{imp}^2 = < n_j (Z-Z_j)^2 > / n,
!*         where Z_j, n_j are partial charges and densities of impurities
!* Output: longitudinal components of the tensors --
!*         SIGMA - electrical conductivity in s^{-1}
!*         CKAPPA - thermal conductivity in erg cm^{-1} s^{-1} K^{-1}
!*         QJ=Q_J - thermopower coefficient Q in units of k/e;
!*     transverse components of the same tensors --
!*         SIGMAT,CKAPPAT,QJT,
!*     off-diagonal (Hall) components of the same tensors --
!*         SIGMAH,CKAPPAH,QJH
!*     (e.g., see Eq.(17) of Potekhin 1999, Astron. Astrophys. 351, 787)
!*                                                       Version 17.02.08


      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (AUM=1822.9,AUD=15819.4,BOHR=137.036)
!* NOTATIONS:
!*        AUM - atomic mass unit divided by the electron mass
!*        AUD - relativistic unit of density in g/cm^3
!*        BOHR - inverse fine structure constant (Bohr radius in rel.un.)
      parameter (UNISIG=7.763d20,UNIKAP=2.778d15)
!* These are rel.units of SIGMA and KAPPA expressed in CGS.
      parameter (UNITEMP=5930.,UNIB=4.414d13)
!* These are rel.units of temperature (in MK) and magnetic field (in G)
      TEMR=T6/UNITEMP ! Temperature in mc^2
      B0=B/UNIB ! B0 is the magnetic field in relativistic units
      DENRI=RHO/(AUD*AUM*CMI) ! number density of ions
!*  Call for the central subroutine which calculates the transport      *
!*  coefficients SIGMA,KAPPA,Q (in Relativistic units):                 *
      call CONDUCT(TEMR,DENRI,B0,Zion,CMI,Zimp,
     & RSIGMA,RLET,RLTT,RQ,RKAPPA,
     & RTSIGMA,RTLET,RTLTT,RTQ,RTKAPPA,
     & RHSIGMA,RHLET,RHLTT,RHQ,RHKAPPA)
      RKAPi=0.
!* If you want to allow for ion thermal conduction in the approximation
!*  of Chugunov & Haensel (2007), then uncomment the next line:
C      call CONDI(TEMR,DENRI,Zion,CMI,RKAPi)
      RKAPPA=RKAPPA+RKAPi
      RTKAPPA=RTKAPPA+RKAPi
!*   -------  CONVERSION TO ORDINARY PHYSICAL (CGSE) UNITS:   --------  
      SIGMA=RSIGMA*UNISIG ! SIGMA in s^{-1}
      CKAPPA=RKAPPA*UNIKAP ! KAPPA in erg/(K cm s)
      QJ=RQ/sqrt(BOHR) ! Q in units of k_B/e, =J_Q
      SIGMAT=RTSIGMA*UNISIG
      CKAPPAT=RTKAPPA*UNIKAP
      QJT=RTQ/sqrt(BOHR)
      SIGMAH=RHSIGMA*UNISIG
      CKAPPAH=RHKAPPA*UNIKAP
      QJH=RHQ/sqrt(BOHR)
      return
      end


!* * * * * * * * * * * * *   Block CONDUCT   * * * * * * * * * * * * *  *
!*           This is the central subroutine which calculates            *
!*              the electron transport coefficients                     *
!*                                                       Version 12.09.19
!* Stems from v.17.02.08. Diff.: more correct account of THtoEL in EECOR 


      subroutine CONDUCT(TEMR,DENRI,B0,Zion,CMI,Zimp,
     & RSIGMA,RLET,RLTT,RQ,RKAPPA,
     & RTSIGMA,RTLET,RTLTT,RTQ,RTKAPPA,
     & RHSIGMA,RHLET,RHLTT,RHQ,RHKAPPA)

!* Input: TEMR - temperature, DENRI - number density of ions, B0- magn.f.
!*        (all in the relativistic units: \hbar = m_e = c = 1), 
!*        Zion and CMI - ion charge and mass numbers,
!*        Zimp - impurity parameter: Z_{imp}^2 = < n_j (Z-Z_j)^2 > / n,
!*         where Z_j, n_j are partial charges and densities of impurities
!* Output: Rx, RTx, RHx - longitudinal, transverse and off-diagonal
!*        transport coefficients in the relativistic units, where 'x' is:
!*         SIGMA - electrical conductivity,
!*         LET - thermoelectric coefficient,
!*         LTT - thermal transport coefficient,
!*         Q - thermopower,
!*         KAPPA - thermal conductivity.

      implicit double precision (A-H), double precision (O-Z)
      save
      data PI/3.14159265d0/
      data AUM/1822.9/,BOHR/137.036/

!*   -------------------   RESTRICTIONS:   --------------------------   *
      if (TEMR.le.0..or.DENRI.le.0..or.B0.lt.0..or.Zion.le.0.
     & .or.CMI.le.0.) then
         print*,'TEMR,DENRI,B0,Zion,CMI:',TEMR,DENRI,B0,Zion,CMI
        stop'CONDUCT: Non-positive input parameter'
      endif
      if (Zion.lt..5) stop'Too small ion charge'
      if (CMI.lt.1.) stop'Too small ion mass'
      if (DENRI.gt.1.d5) stop'Too high density' ! 1.d4 -> 1.d5: 17.02.04
!*   -----------------   PLASMA PARAMETERS   ------------------------   *
      DENR=DENRI*Zion ! DENR - number density of electrons
      SPHERION=(.75/PI/DENRI)**.3333333 ! Ion sphere radius
      GAMMA=Zion**2/BOHR/TEMR/SPHERION ! Ion coupling parameter
      XSR=(3.*PI**2*DENR)**.3333333 ! special relativity parameter 
!*   XSR equals the non-magnetic Fermi momentum
      if (B0.gt.0.) then
         CNL=((sqrt(1.d0+XSR**2)+3.d0*TEMR)**2-1.d0)/2.d0/B0
      else
         CNL=10000.
      endif
!*   CNL estimates maximum contributing Landau number
      if (CNL.lt.200.) then ! 15.07.97
         MAGNET=1
      else
         MAGNET=0
      endif
!*   -------------------  Chemical potential   ----------------------   *
      CMU=CHEMPOT(B0,DENR,DDENR,TEMR)
      Q2e=4.d0*PI/BOHR*DDENR ! squared electron screening momentum
!*   -------------------  Thermal averaging    ----------------------   *
      xnuc=0. ! (1.1d-7*DENR)**.3333333d0 ! r_nuc=1.15 A^{1/3} fm (outer cr.)
      xnuct=xnuc/1.1 ! fiducial correction for thermal conductivity
      xnimp=xnuc
      call ThAv18(MAGNET,CMU,B0,TEMR,DENRI,
     *  XSR,GAMMA,Zion,CMI,Zimp,Q2e,
     *  xnuc,xnuct,xnimp,
     *  F0B,F1B,F2B,F0C,F1C,F2C,F0D,F1D,F2D)
!*   ----------------------------------------------------------------   *
!*   Longitudinal transport coefficients:
      RSIGMA=F0B/BOHR ! el. conductivity
      RLET=F1B/sqrt(BOHR) ! thermoel. coeff.
      RLTT=F2B*TEMR ! thermal transport coeff.
      RQ=RLET/RSIGMA ! thermopower
      RKAPPA=RLTT-TEMR*RLET**2/RSIGMA ! thermal conductivity
!*   Transverse transport coefficients:
!*     ... dissipative:
      RTSIGMA=F0C/BOHR
      RTLET=F1C/sqrt(BOHR)
      RTLTT=F2C*TEMR
!*     ... non-dissipative:
      RHSIGMA=F0D/BOHR
      RHLET=F1D/sqrt(BOHR)
      RHLTT=F2D*TEMR
      D=RTSIGMA**2+RHSIGMA**2 ! deteminant
!*     ... thermopower:
      RTQ=(RTSIGMA*RTLET+RHSIGMA*RHLET)/D ! diag.
      RHQ=(RTSIGMA*RHLET-RHSIGMA*RTLET)/D ! non-diag.
!*     ... thermal conductivity
      RTKAPPA=RTLTT-TEMR*(RTLET*RTQ-RHLET*RHQ) ! diag.
      RHKAPPA=RHLTT-TEMR*(RTLET*RHQ+RHLET*RTQ) ! non-diag.
!*   ------  Include nonmagnetic corrections for thermal conductivity:
      if (B0.lt.TEMR*CMU.or.B0.lt.CMU) then ! corrected 12.09.19
!* (1) correction for e-e collisions:
         CTH=TEMR/9.d0*XSR**3/sqrt(1.d0+XSR**2)
         TAUlongT=RKAPPA/CTH ! eff.thermal.e-i relax.time
         call TAUEESY(XSR,TEMR,TAUEE) ! eff.e-e relax.time
!* (2) correction for ionic quantum effects (W-F law violation at low T)
         VCL=XSR/dsqrt(1.d0+XSR**2)
         TRP=Zion/GAMMA*dsqrt(CMI*AUM*SPHERION/3.d0/BOHR) ! =T/T_p
         BORNCOR=VCL*Zion*PI/BOHR ! first non-Born correction
         T0=.19/Zion**.16667 ! critical T/T_p parameter
         G0=TRP/sqrt(TRP**2+T0**2)*(1.d0+(Zion/125.)**2)
         G2=TRP/sqrt(.0081+TRP**2)**3
         THtoEL=1.d0+G2/G0*(1.d0+BORNCOR*VCL**3)*.0105*(1.d0-1.d0/Zion)
!* result:
         EECOR=TAUEE/(TAUlongT+TAUEE*THtoEL)
         RKAPPA=RKAPPA*EECOR
         RTKAPPA=RTKAPPA*EECOR
         RHKAPPA=RHKAPPA*EECOR
      endif
      return
      end

      subroutine ThAv18(MAGNET,CMU,B,TEMR,DENRI,
     *  XSR,GAMMA,Zion,CMI,Zimp,Q2e,
     *  xnuc,xnuct,xnimp,
     *  F0,F1,F2,F0C,F1C,F2C,F0D,F1D,F2D)
!*   Thermal Averaging                                   Version 09.07.18
!* Stems from ThAv16 v.17.05.17. Difference - incl.xnimp.
!* Diff.from ThAv99 v.26.01.01/12.05.11:
!*   - incl. xnuc,xnuct;
!*   - got rid of "assign MARK"
!*    MAGNET=1 - magnetic case, 0 - non-magnetic case
!*    CMU - chemical potential, B -magn.field, TEMR -temperature (mc^2)
!*    XSR - relativity (density) param., GAMMA - ion coupling param.
!*    Zion and CMI - ion charge and mass numbers
!*    Q2e - e-screening wavenumber squared
!*   Output: Fn=\int [(E-\mu)/T]^n [N_e(E)/E] \tau_{eff}(E) df_0(E)
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(BOHR=137.036d0,PI=3.14159265d0)
      dimension AXI(4),AWI(4)
      data TAIL1/8.d0/,TAIL2/40.d0/,NM/128/
!*        Nodes and weights of the Gauss formula:
      data AXI/.32425 34234 0381d0, .61337 14327 0059d0,
     *         .83603 11073 2664d0, .96816 02395 0763d0/,
     *  AWI/.31234 70770 4   d0, .26061 06964 0294d0,
     *      .18064 81606 9486d0, .08127 43883 6157d0/,NXI/4/
      if (MAGNET.ne.0.and.MAGNET.ne.1) stop'ThAv: incorrect MAGNET'
CC      CST=XSR**3/.75d0/PI*Zion/BOHR**2 ! OLD: valid only if Zion=Zcell
      CST=4.d0*PI*DENRI*(Zion/BOHR)**2 ! =4\pi n_i(Ze^2)^2
      EMIN=dmax1(1.d0,CMU-TAIL2*TEMR)
      EMAX=CMU+TAIL2*TEMR
      E1=dmax1(1.d0,CMU-TAIL1*TEMR)
      E2=dmax1(1.d0,CMU+TAIL1*TEMR)
      if (MAGNET.eq.1) then
         NLMIN=(EMIN**2-1.d0)/2.d0/B
         NLMAX=(EMAX**2-1.d0)/2.d0/B
         NL1=(E1**2-1.d0)/2.d0/B
         NL2=(E2**2-1.d0)/2.d0/B
      else
         NLMIN=0
         NLMAX=0
         NL1=0
         NL2=0
      endif
      F=0. ! correction 12.05.11
      FC=0. ! correction 12.05.11
      FD=0. ! correction 12.05.11
      S0=0.
      S1=0.
      S2=0.
      S0C=0.
      S1C=0.
      S2C=0.
      S0D=0.
      S1D=0.
      S2D=0.
      do 50 NL=NLMIN,NLMAX
      if (MAGNET.eq.1) then
         EN=dsqrt(1.d0+2.d0*B*NL)
         ENP=dsqrt(EN**2+2.d0*B)
      else ! non-quantizing field
         EN=EMIN
         ENP=EMAX
         goto 20
      endif
      if ((NL.ge.NL1-1.and.NL.le.NL2+1).and.NL2.lt.NL1+5) goto 20
      EFIN=ENP
!*  -----------   The first type of integration: Gauss   -------------  *
      SCALE=(EFIN-EN)/TEMR
      DS0=0.
      DS1=0.
      DS2=0.
      DS0C=0.
      DS1C=0.
      DS2C=0.
      DS0D=0.
      DS1D=0.
      DS2D=0.
      do 12 IX=1,NXI ! Abramowitz & Stegun 25.4.34.
      XI=AXI(IX)
      WI=AWI(IX)
      X=XI**2
      W=2.d0*XI**2*WI
      E=EN+(EFIN-EN)*X
      call ThAvI18(E,CMU,TEMR,CST, ! input
     *  XSR,GAMMA,B,Zion,CMI,Zimp,Q2e,xnuc,xnuct,xnimp, ! input
     *  T,F,FC,FD,MARK) ! output
      if (MARK.eq.1) goto 12
   11 continue
!*        F:=F1*PHI, T:=(E-CMU)/TEMR
      F=F/XI ! improvement of 01.09.96
      DDS=W*F*SCALE
      DS0=DS0+DDS
      DS1=DS1+DDS*T
      DS2=DS2+DDS*T**2
      FC=FC/XI ! improvement of 01.09.96
      DDSC=W*FC*SCALE
      DS0C=DS0C+DDSC
      DS1C=DS1C+DDSC*T
      DS2C=DS2C+DDSC*T**2
      FD=FD/XI ! improvement of 01.09.96
      DDSD=W*FD*SCALE
      DS0D=DS0D+DDSD
      DS1D=DS1D+DDSD*T
      DS2D=DS2D+DDSD*T**2
   12 continue ! next integration step
      S0=S0+DS0
      S1=S1+DS1
      S2=S2+DS2
      S0C=S0C+DS0C
      S1C=S1C+DS1C
      S2C=S2C+DS2C
      S0D=S0D+DS0D
      S1D=S1D+DS1D
      S2D=S2D+DS2D
      goto 50
!*  --------   The second type of integration: Simpson   -------------  *
!*INTERVALS: NL=NL1=NL2 => max(EN,EMIN)..(1)..E1..(2)..E2..(3)..ENP
!*    NL2 > NL1, NL=NL1 => max(EN,EMIN)..(1)..E1..(2)..ENP
!*    NL2 > NL1, NL=NL2 =>           EN..(1)..E2..(2)..min(ENP,EMAX)
!*       NL1 < NL < NL2 =>           EN..(1)..ENP
   20 continue
      ESTART=dmax1(EN,EMIN)
      IK=0
   13 IK=IK+1
      if (IK.gt.1) ESTART=EFIN
      EFIN=dmin1(ENP,EMAX)
      if (IK.eq.1) then ! 1st interval of 2 or 3
        if (NL.eq.NL1) EFIN=E1
        if (NL.eq.NL2.and.NL.ne.NL1) EFIN=E2
      endif
      if (IK.eq.2.and.NL.eq.NL1.and.NL.eq.NL2) EFIN=E2 ! 2nd int. of 3
      DE=EFIN-ESTART
      SCALE=DE/TEMR
      H=1.d0/NM
      X=0.
      BLINK=1
      W=H/3.d0*SCALE
      E=ESTART
      call ThAvI18(E,CMU,TEMR,CST, ! input
     *  XSR,GAMMA,B,Zion,CMI,Zimp,Q2e,xnuc,xnuct,xnimp, ! input
     *  T,F,FC,FD,MARK) ! output
   14 continue
!*        F:=F1*PHI, T:=(E-CMU)/TEMR
      DDS=W*F
      DS0=DDS
      DS1=DDS*T
      DS2=DDS*T**2
      DDSC=W*FC
      DS0C=DDSC
      DS1C=DDSC*T
      DS2C=DDSC*T**2
      DDSD=W*FD
      DS0D=DDSD
      DS1D=DDSD*T
      DS2D=DDSD*T**2
      do 16 K=1,NM
      X=X+H
      E=ESTART+DE*X
      call ThAvI18(E,CMU,TEMR,CST, ! input
     *  XSR,GAMMA,B,Zion,CMI,Zimp,Q2e,xnuc,xnuct,xnimp, ! input
     *  T,F,FC,FD,MARK) ! output
      if (MARK.eq.1) goto 16
   15 continue
!*        F:=F1*PHI, T:=(E-CMU)/TEMR
      DDS=W*F*(3.d0+BLINK)
      DS0=DS0+DDS
      DS1=DS1+DDS*T
      DS2=DS2+DDS*T**2
      DDSC=W*FC*(3.d0+BLINK)
      DS0C=DS0C+DDSC
      DS1C=DS1C+DDSC*T
      DS2C=DS2C+DDSC*T**2
      DDSD=W*FD*(3.d0+BLINK)
      DS0D=DS0D+DDSD
      DS1D=DS1D+DDSD*T
      DS2D=DS2D+DDSD*T**2
      BLINK=-BLINK
   16 continue
      S0=S0+DS0-DDS/2.d0
      S1=S1+DS1-DDS*T/2.d0
      S2=S2+DS2-DDS*T**2/2.d0
      S0C=S0C+DS0C-DDSC/2.d0
      S1C=S1C+DS1C-DDSC*T/2.d0
      S2C=S2C+DS2C-DDSC*T**2/2.d0
      S0D=S0D+DS0D-DDSD/2.d0
      S1D=S1D+DS1D-DDSD*T/2.d0
      S2D=S2D+DS2D-DDSD*T**2/2.d0
      if (IK.eq.1.and.(NL.eq.NL1.or.NL.eq.NL2))goto 13 ! to 2nd interval
      if (NL1.eq.NL2.and.IK.eq.2) goto 13 ! to the 3rd interval
   50 continue
      F0=S0
      F1=S1
      F2=S2
      F0C=S0C
      F1C=S1C
      F2C=S2C
      F0D=S0D
      F1D=S1D
      F2D=S2D
      return
      end

      subroutine ThAvI18(E,CMU,TEMR,CST, ! input
     *  XSR,GAMMA,B,Zion,CMI,Zimp,Q2e,xnuc,xnuct,xnimp, ! input
     *  T,F,FC,FD,MARK) ! output
!* Calculation of the integrand for thermal averaging S/R ThAv18
!*                                                       Version 08.07.18
!* Stems from ThAvI v.10.03.18.
!* Diff.: xnimp
!* Input: E - electron longitudinal kinetic energy /mc^2
!*        CMU - electron chem.pot./mc^2
!*        TEMR - temperature [rel.un.]
!*        CST=XSR**3/.75/PI*Zion/BOHR**2*Zion/Zcell=4\pi n_i(Ze^2)^2
!*        XSR - non-magnetic relativity parameter at given density n_e
!*        GAMMA - ion Coulomb coupling parameter
!*        B - magnetic field,
!*        Zion - mean charge of the ion,
!*        CMI - mean atomic weight,
!*        Zimp - impurity parameter (eff.charge)
!*        Q2e - squared electron screening wavenumber
!*        xnuc=r_nuc/a_i=sqrt{5<r^2>/3}/a_i - nucl.size parameter
!*        xnuct - auxiliary nuc.size param. for thermal conductivity
!*        xnimp - xnuc for impurities
!* Output: T=(E-CMU)/TEMR
!*         F, FC, FD - integrand for longit., transv., off-diag.component
!*         MARK - auxiliary index, equal to 1 when the integrand =0.
!* Everything is in the relativistic units.
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(PI=3.14159265d0)
      T=(E-CMU)/TEMR
      PCL=dsqrt(E**2-1.d0) ! non-magnetic momentum
C      if (dabs(T).gt.40.d0.or.E.eq.EN) then ! EN excluded 06.11.16
      if (dabs(T).gt.40.d0) then
         F=0.
         FC=0. ! corr.06.11.16
         FD=0. ! corr.06.11.16
         MARK=1
        return
      else
         EX=dexp(T)
         F1=EX/(EX+1.d0)**2
      endif
!*        F1 = -df_0/d(E/T)
      if (E-1.d0.lt.1.d-10) then
         F=0.
         FC=0. ! corr.06.11.16
         FD=0. ! corr.06.11.16
         MARK=2
        return
      endif
      call COUL19(PCL,XSR,GAMMA,B,Zion,CMI,Q2e,xnuc,xnuct,
     *   CLeff,CLlong,CLtran,SN,THtoEL)
      GYROM=B/E ! gyrofrequency
      TAUlong=PCL**3/E/CST/CLlong/SN
      TAUt=PCL**3/E/CST/CLtran*SN
      if (Zimp.gt.0.) then ! include impurity scattering
         call COUL18I(PCL,XSR,GAMMA,B,Q2e,xnimp,
     *   CLeffI,CLlongI,CLtranI,SNDUMMY) ! corrected 10.03.18 (DUMMY)
         TAUlong=TAUlong/(1.d0+CLlongI/CLlong*(Zimp/Zion)**2)
         TAUt=TAUt*CLtran*Zion**2/(CLtran*Zion**2+CLtranI*Zimp**2)
      endif
      TAUtran=TAUt/(1.d0+(TAUt*GYROM)**2)
      TAUhall=TAUt*GYROM*TAUtran
      C=SN*PCL**3/3.d0/PI**2/E ! common factor
      F=F1*C*TAUlong
      FC=F1*C*TAUtran
      FD=F1*C*TAUhall
      MARK=2
      return
      end

      subroutine COUL19(PCL,XSR,GAMMA,B,Zion,CMI,Q2e,xnuc,xnuct,
     *   CLeff,CLlong,CLtran,SN,THtoEL)
!*                                                       Version 17.09.19
!* Origins: COUL01 v.26.01.01, COULINsc v.24.02.00, and COULIN v.23.05.07
!* Predecessor: COUL15 v.05.12.17. Difference - explicit PI,AUM,BOHR.
!*   Input: PCL - non-magnetic electron momentum \equiv \sqrt(E^2-1),
!*          Q2e - squared electron screening wavenumber (IN REL.UNITS)
!*          XSR = p_F/mc - relativity (density) parameter,
!*          GAMMA - Coulomb coupling parameter of ions,
!*          B - magnetic field,
!*          Zion - mean charge of the ion,
!*          CMI - mean atomic weight,
!*          xnuc - eff. ion radius divided by Wigner-Seitz cell radius
!*          xnuct - a second xnuc parameter for use in a quantum crystal
!*   Output: CLlong, CLtran - eff.Coulomb log.,
!*           SN = N_e(E)/N_0(E) = (3/2)(eB\hbar/c)\sum_{ns} p_n/p_0^3
!*           THtoEL - ratio of thermal to electric eff.relax.times
!* Dimensional quantities are in the relativistic units (m_e=\hbar=c=1)
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (PI=3.14159265359d0)
      parameter(AUM=1822.88848d0) ! a.m.u./m_e
      parameter (BOHR=137.035999d0) ! inverse fine-structure constant
      parameter (Uminus1=2.79855,Uminus2=12.972) ! bcc-lattice param.
      parameter(BIG=1.d6)
!*   ----------------------   Preliminaries   -----------------------   *
      DENR=XSR**3/3.d0/PI**2 ! number density of electrons
      DENRI=DENR/Zion ! number density of ions (rel.)
      SPHERION=(.75/PI/DENRI)**.3333333 ! Ion sphere radius
      Q2icl=3.d0*GAMMA/SPHERION**2 ! squared Debye screening momentum
      ECL=sqrt(1.d0+PCL**2) ! Energy
      VCL=PCL/ECL ! Velocity
      PM2=(2.d0*PCL)**2 ! squared max.momentum transfer
      TRP=Zion/GAMMA*dsqrt(CMI*AUM*SPHERION/3.d0/BOHR) ! =T/T_p
      BORNCOR=VCL*Zion*PI/BOHR ! first non-Born correction
!*   ---------------------   Non-magnetic fit   ---------------------   *
      C=(1.d0+.06*GAMMA)*dexp(-dsqrt(GAMMA))
      Q2s=(Q2icl*C+Q2e)*dexp(-BORNCOR) ! eff.scr.wavenumber in rel.un.
      XS=Q2s/PM2 ! eff.screening param.
      R2W=Uminus2/Q2icl*(1.d0+.333333*BORNCOR)
      XW=R2W*PM2 ! eff. Debye-Waller param.
!** Modification WITH FINITE SIZES OF NUCLEI; xnuc=r_{nuc}/a_i
      XW1=14.7327d0*xnuc**2 ! =4(9\pi/4)^{2/3} x_{nucl}^2 =coeff.at q^2
      XW1=XW1*(1.d0+C13*BORNCOR)*(1.d0+Zion/13.d0*dsqrt(xnuc))
      XW1=XW1*(PCL/XSR)**2 ! correction of 05.12.17
      call COULAN3(XS,XW,PCL,XW1,CL,DLeff)
      A0=1.683*sqrt(PCL/CMI/Zion) ! zero-vibr.param.(Baiko&Yakovlev95)
      VIBRCOR=exp(-A0/4.*Uminus1*exp(-9.1*TRP)) ! corr.for zero-vibr.
      T0=.19/Zion**.16667 ! critical T/T_p parameter
      G0=TRP/sqrt(TRP**2+T0**2)*(1.+(Zion/125.)**2) ! mod.10.01.99
      GW=G0*VIBRCOR
      CLeff=CL*GW ! 1st FIT (for non-magnetic electrical conductivity)
      G2=TRP/sqrt(.0081+TRP**2)**3
      THtoEL=1.+G2/G0*(1.+BORNCOR*VCL**3)*.0105*(1.-1./Zion)*
     *  (1.d0+xnuct**2*dsqrt(2.d0*Zion))
      if (PCL**2.gt.4.d2*B) then ! Non-magnetic case
         CLlong=CLeff
         CLtran=CLeff
         SN=1.d0
         goto 50
      endif
!*   -----------------------   Magnetic fit   -----------------------   *
      ENU=PCL**2/2.d0/B
      NL=ENU
!*   Cutoff according to Hernquist'84:
      TEMR=(Zion/BOHR)**2/SPHERION/GAMMA ! T[rel.un.] from a_i,Z,\Gamma
      if (GAMMA.lt.175.d0) then
         GAMMAN1=dsqrt(TEMR/CMI/AUM)*2.*PCL ! ei-inelasticity (liquid)
         TAUcl=PCL**2*VCL*BOHR**2/(4.*PI*Zion**2*DENRI)
      else
         GAMMAN1=TEMR/TRP ! =T_p ! ei(e-ph)-inelasticity (solid)
         TAUcl=VCL*BOHR/(TEMR*(2.-VCL**2)*Uminus2)
      endif
      GAMMAN2=1./TAUcl ! collisional broadening
      GAMMAN=dmax1(GAMMAN1,GAMMAN2)
      SN=0.
      do N=0,NL
         PB=dsqrt(ENU-N) ! =p_n/sqrt(2b)
         SN=SN+PB
        if (N.ne.0) SN=SN+PB
      enddo
      SN=SN*1.5d0*B*dsqrt(2.d0*B)/PCL**3
      if (ENU.le.1.d0) then ! Exact calculation
         Xis=Q2s/2.d0/B ! Screening parameter, scaled magnetically
         ZETA=R2W*2.d0*B ! magn.scaled exponential coefficient
         Xi=2.d0*PCL**2/B
         Xsum=Xi+Xis
         Q2M=(EXPINT(Xsum,1)-
     -     dexp(-ZETA*Xi)*EXPINT((1.d0+ZETA)*Xsum,1))/Xsum
         CLlong=(PCL*VCL/B)**2*Q2M/1.5d0*GW
        if (Xsum.lt.BIG) then ! condition added 26.06.17
           QTM1=(1.d0+Xsum)*EXPINT(Xsum,0)-1.d0
           QTM2=(1.d0+(1.d0+ZETA)*Xsum)*EXPINT((1.d0+ZETA)*Xsum,0)-1.d0
        else
           QTM1=1.d0/Xsum**2
           QTM2=1.d0/((1.d0+ZETA)*Xsum)**2
        endif
         QtranM=QTM1-dexp(-ZETA*Xi)*QTM2
        if (Xis.lt.BIG) then ! condition added 26.06.17
           QTP1=(1.d0+Xis)*EXPINT(Xis,0)-1.d0
           QTP2=(1.d0+(1.d0+ZETA)*Xis)*EXPINT((1.d0+ZETA)*Xis,0)-1.d0
        else
           QTP1=1.d0/Xis**2
           QTP2=1.d0/((1.d0+ZETA)*Xis)**2
        endif
         QtranP=QTP1-QTP2
         Q=(ECL**2*QtranP+QtranM)*B/PCL**2 ! Q(E,b)
         CLtran=.375d0*Q/ECL**2*GW
      else
!*   Preliminaries:
         DNU=ENU-NL
         DNU=dmax1(DNU,GAMMAN) ! cutoff
         XS1=(dsqrt(XS)+1.d0/(2.d0+XW/2.d0))**2
         PN=dsqrt(2.d0*B*DNU)
         SQB=dsqrt(B)
         X=dmax1(PN/SQB,1.d-10)
!*   Longitudinal:
        if (XW.lt..01) then
           EXW=1.d0
        elseif (XW.gt.50.) then
           EXW=1.d0/XW
        else
           EXW=(1.d0-dexp(-XW))/XW
        endif
         A1=(30.-15.*EXW-(15.-6.*EXW)*VCL**2)/
     /     (30.-10.*EXW-(20.-5.*EXW)*VCL**2)
         Q1=.25d0*VCL**2/(1.d0-.6666667d0*VCL**2)
         DLT=SQB/PCL*(A1/X-sqrt(X)*(1.5d0-.5d0*EXW+Q1)+
     +     (1.d0-EXW+.75d0*VCL**2)/(1.d0+VCL**2)*(X-sqrt(X))/NL)
         Y1=1.d0/(1.d0+DLT)
         CL0=dlog(1.d0+1.d0/XS1)
         P2=CL0*(.07+.2*EXW)
         Y2=1.5*CL0*(X**3-X/3.)/(NL+.75/(1.+2.*B)**2*X**2)+P2*X
         PY=1.+.06*CL0**2/NL**2
         DT=dsqrt(PY*Y1**2+Y2**2) ! ratio of relax.times
         CLlong=CLeff/DT
!*   Transverse:
         DB=1./(1.+.5/B)
         CL1=XS1*CL0
         P1=.8*(1.+CL1)+.2*CL0
         P2=1.42-.1*DB+sqrt(CL1)/3.
         P3=(.68-.13*DB)*CL1**.165
         P4=(.52-.1*DB)*sqrt(sqrt(CL1))
         DLT=SQB/PCL*(P1/X**2*SQB/PCL+P3*alog(NL+0.)/X-
     -     (P2+P4*alog(NL+0.))*sqrt(X))
         CLtran=CLeff*(1.+DLT)
      endif
   50 return
      end

      subroutine COUL18I(PCL,XSR,GAMMA,B,Q2e,xnimp,
     *   CLeff,CLlong,CLtran,SN) ! IMPURITY
!*                                                       Version 09.07.18
!*  This is a simplified version of COUL19 for the el.-impurity sc.
!*  Diff.from COUL99I - incl.xnimp
!*   Input: XSR = p_F/mc - relativity (density) parameter,
!*          PCL - non-magnetic electron momentum \equiv \sqrt(E^2-1),
!*          GAMMA - Coulomb coupling parameter of ions,
!*          B - magnetic field,
!*          Q2e - squared electron screening wavenumber
!*    (ALL IN THE RELATIVISTIC UNITS)
!*          xnimp - eff. ion radius divided by Wigner-Seitz cell radius
!*   Output: CLlong, CLtran - eff.Coulomb log.,
!*           SN = N_e(E)/N_0(E) = (3/2)(eB\hbar/c)\sum_{ns} p_n/p_0^3
      implicit double precision (A-H), double precision (O-Z)
      save
!* Dimensional quantities are in the relativistic units (m_e=\hbar=c=1)
!*        BOHR - radius of the first Bohr orbit in the rel.units
      data PI/3.14159265/,XW/1.d99/ ! XW=infinity
!*   ----------------------   Preliminaries   -----------------------   *
      DENR=XSR**3/3./PI**2 ! number density of electrons
      ECL=sqrt(1.+PCL**2) ! Energy
      VCL=PCL/ECL ! Velocity
      PM2=(2.*PCL)**2 ! squared max.momentum transfer
!*   ---------------------   Non-magnetic fit   ---------------------   *
      C=(1.+.06*GAMMA)*dexp(-dsqrt(GAMMA))
      Q2s=Q2e ! eff.scr.wavenumber in rel.un.
      XS=Q2s/PM2 ! eff.screening param.
      XW1=14.7327d0*xnimp**2 ! =4(9\pi/4)^{2/3} x_{nucl}^2 =coeff.at q^2
      XW1=XW1*(PCL/XSR)**2
      call COULAN3(XS,XW,PCL,XW1,CLeff,DLeff) ! nonmagn.Coulomb log.
      if (PCL**2.gt.4.d2*B) then ! Non-magnetic case
         CLlong=CLeff
         CLtran=CLeff
         SN=1.d0
         goto 50
      endif
!*   -----------------------   Magnetic fit   -----------------------   *
      ENU=PCL**2/2.d0/B
      NL=ENU
      SN=0.
      do N=0,NL
         PB=dsqrt(ENU-N) ! =p_n/sqrt(2b)
         SN=SN+PB
        if (N.ne.0) SN=SN+PB
      enddo
      SN=SN*1.5d0*B*dsqrt(2.d0*B)/PCL**3
      if (ENU.le.1.d0) then ! Exact calculation
         Xis=Q2s/2./B ! Screening parameter, scaled magnetically
         Xi=2.*PCL**2/B
         Xsum=Xi+Xis
         Q2M=EXPINT(Xsum,1)
         CLlong=(PCL*VCL/B)**2*Q2M/1.5
         QtranM=(1.+Xsum)*EXPINT(Xsum,0)-1.
         QtranP=(1.+Xis)*EXPINT(Xis,0)-1.
         Q=(ECL**2*QtranP+QtranM)*B/PCL**2 ! Q(E,b)
         CLtran=.375*Q/ECL**2
      else
!*   Preliminaries:
         DNU=ENU-NL
         XS1=(dsqrt(XS)+1./(2.+XW/2.))**2
         PN=dsqrt(2.*B*DNU)
         SQB=dsqrt(B)
         X=dmax1(PN/SQB,1.d-10)
!*   Longitudinal:
        if (XW.lt..01) then
           EXW=1.d0
        elseif (XW.gt.50.) then
           EXW=1.d0/XW
        else
           EXW=(1.d0-dexp(-XW))/XW
        endif
         A1=(30.-15.*EXW-(15.-6.*EXW)*VCL**2)/
     /     (30.-10.*EXW-(20.-5.*EXW)*VCL**2)
         Q1=.25*VCL**2/(1.-.667*VCL**2)
         DLT=SQB/PCL*(A1/X-sqrt(X)*(1.5-.5*EXW+Q1)+
     +     (1.-EXW+.75*VCL**2)/(1.+VCL**2)*(X-sqrt(X))/NL)
         Y1=1./(1.+DLT)
         CL0=dlog(1.d0+1.d0/XS1)
         P2=CL0*(.07+.2*EXW)
         Y2=1.5*CL0*(X**3-X/3.)/(NL+.75/(1.+2.*B)**2*X**2)+P2*X
         PY=1.+.06*CL0**2/NL**2
         DT=dsqrt(PY*Y1**2+Y2**2) ! ratio of relax.times
         CLlong=CLeff/DT
!*   Transverse:
         DB=1./(1.+.5/B)
         CL1=XS1*CL0
         P1=.8*(1.+CL1)+.2*CL0
         P2=1.42-.1*DB+sqrt(CL1)/3.
         P3=(.68-.13*DB)*CL1**.165
         P4=(.52-.1*DB)*sqrt(sqrt(CL1))
         DLT=SQB/PCL*(P1/X**2*SQB/PCL+P3*alog(NL+0.)/X-
     -     (P2+P4*alog(NL+0.))*sqrt(X))
         CLtran=CLeff*(1.+DLT)
      endif
   50 return
      end

      subroutine COULAN3(XS,XW0,PCL,XW1,CLeff,DLeff)
!*  ------    Analytic expression for Coulomb logarithm  Version 19.11.17
!* Stems from function COULAN2 v.23.05.00.
!*  Implements expressions from:
!*    Appendix A1.1 of Gnedin, Potekhin, Yakovlev [MNRAS 324, 725 (2001)]
!*    and in Appendix of Potekhin et al. [A&A 346, 345 (1999)
!*    with a patch and a supplement of 2017 (pp.67,70 of WRITTEN/CONDUCT)
!* The patch of 18.05.17 fixes an accuracy loss.
!* Input: XS=(q_s/2p_0)^2, where p_0 - cl.momentum, q_s - eff.scr.mom.
!*   XW0=u_{-2} (2p_0/\hbar q_D)^2, where u_{-2}=13, q_D^2=3\Gamma/a_i^2
!*   PCL=p_0/mc^2 - classical momentum in rel.un.
!*   XW1=s1*(2p_0/\hbar)^2, where s1 \approx r_{nuc}^2
!* Output: CLeff - eff. Coulomb logarithm \Lambda
!*         DLeff - its derivative over energy d\Lambda / dE [rel.un.]
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(EPS=1.d-2,EPS1=1.d-3,TINY=1.d-9,EULER=0.5772156649d0)
      parameter(BIG=1.d0/EPS)
      if(XS.lt.0..or.XW0.lt.0..or.V.lt.0..or.XW1.lt.0.)stop'COULAN3'
      E=dsqrt(1.d0+PCL**2)
      V=PCL/E ! V=p/E - cl.velocity [rel.un.]
      KEY=0
      MI=1
      if (XW1.lt.TINY) MI=0
!* ----------------------------------------------------------------------
      do I=0,MI
        if (I.eq.0) then
           XW=XW0+XW1
           B=XS*XW
        else ! to do the 2nd term
          XW=XW1
          B=XS*XW
        endif
        if (I.eq.0.or.KEY.eq.2) then
!* Check applicability of asymptotes:
          if (XW.lt.EPS) then
             KEY=1
             goto 50
          endif
          if (XW.gt.BIG.and.B.gt.BIG) then
             KEY=2
          elseif (XS.lt.EPS1.and.B.lt.EPS1/(1.d0+XW)) then
             KEY=3
          else
             KEY=4
          endif
        endif
   50   continue
         EA=dexp(-XW)
         E1=1.d0-EA
        if (XW.lt.EPS1) E1=XW*(1.d0-XW/2.d0)
        if (KEY.ne.1) then
          if (XW.gt.EPS1) then
             E2=(XW-E1)/XW ! (e^{-w}-1+w)/w
          else
             E2=.5d0*XW*(1.d0-XW/3.d0)
          endif
        endif
        XS1=XS+1.d0
        if (KEY.eq.1) then
           CL0=dlog(XS1/XS)
           CL1=.5d0*XW*(2.d0-1.d0/XS1-2.d0*XS*CL0)
           CL2=.5d0*XW*(1.5d0-3.d0*XS-1.d0/XS1+3.d0*XS**2*CL0)
        elseif (KEY.eq.2) then
           CL0=dlog(XS1/XS)
           CL1=(CL0-1.d0/XS1-1.d0/B**2)/2.d0 ! -1/B^2 added wrt COULAN2
           CL2=(2.d0*XS+1.d0)/(2.d0*XS1)-XS*CL0
        elseif (KEY.eq.3) then
           CL1=.5d0*(EA*EXPINT(XW,0)+dlog(XW)+EULER)
           CL2=.5d0*E2
        elseif (KEY.ge.4) then
           CL0=dlog(XS1/XS)
           EL=EXPINT(B,0)-EXPINT(B+XW,0)*EA
           CL1=.5d0*(CL0+XS/XS1*E1-(1.d0+B)*EL)
           CL2=.5d0*(E2-XS*XS/XS1*E1-2.d0*XS*CL0+XS*(2.d0+B)*EL)
        else
           stop'COULAN3: invalid KEY'
        endif
        if (I.eq.0) then ! 1st term calculated
           CLeff=CL1-V**2*CL2
           CL10=CL1
           CL20=CL2
        else ! 2nd term calculated
           CLeff=CLeff-(CL1-V**2*CL2)
           CL1=CL10-CL1
           CL2=CL20-CL2
        endif
      enddo
      if (XW0.gt.EPS) then
         EXW0=1.d0-dexp(-XW0)
      else
         EXW0=XW0
      endif
      DL1=EXW0/(PCL*V)/XS1**2*dexp(-XW1)
      DLeff=DL1/E**2+2.d0*V**2/E*CL2
      return
      end

!*   =========  CHEMICAL POTENTIAL and DENRITY DERIVATIVE  ==========   *
      function CHEMPOT(B0,DENR,DDENR,TEMR)
      implicit double precision (A-H), double precision (O-Z)
      save
!*  Determines chemical potential of free electron gas
!*  All quantities are by default in relativistic units
!*  DENR - electron density, TEMR - temperature (input),
!*  DDENR - partial derivative of DENR with respect to the chem.pot.
!*  B0 - magnetic field (=h\omega_B/mc^2)
!*                                                       Version 03.05.07
      data PI/3.14159265/,EPS/1.D-3/,NI/50/
      PF0=(3.*PI**2*DENR)**.3333333 ! Classical Fermi momentum
      if (PF0.gt..001) then
         THETA=TEMR/(sqrt(1.+PF0**2)-1.) ! degeneracy
      else
         THETA=2.*TEMR/PF0**2
      endif
      if (B0.gt.0.) then
         CNL=((sqrt(1.+PF0**2)+6.*TEMR)**2-1.)/2./B0 ! 25.11.96
      else
         CNL=10000.
      endif
!*   CNL estimates maximum contributing Landau number
      if (CNL.lt.200.) then ! 25.11.96
         MAGNET=1
         B=B0
         ZET=(dmin1(1.d0,PF0**2/B/1.5)*PF0)**2
      else
!*   ------------------  Fit to chemical potential  -----------------   *
         CHEMPOT=CHEMFIT(DENR,DDENR,TEMR)
         return
      endif
!*   ----------  Calculation of chemical potential CMU   ------------   *
      Kup=0.
      Klow=0.
      POWER=.666667
      if (CNL.lt.1.) POWER=2.
      DENR1=0. ! 22.01.97
      do 5 I=1,NI
      if (ZET.gt.3.*TEMR) then
!*  Degenerate plasma
         SHIFTMU=0.
      else
!*  Non-degenerate plasma
         SHIFTMU=1.5*TEMR*dlog(3.*TEMR/ZET)
         if (CNL.lt.2.) SHIFTMU=SHIFTMU/3.
      endif
      if (SHIFTMU.lt..5) then
         CMU=sqrt(1.+ZET)*(1.-SHIFTMU)
      else
         CMU=sqrt(1.+ZET)-SHIFTMU
      endif
      DENR1old=DENR1
      call DENRFIT(B,TEMR,CMU,DENR1,DDENR)
      if (DENR1.gt.DENR) then
        if (Kup.eq.0) Zup=ZET
         Kup=1
         Zup=dmin1(Zup,ZET)
      else
        if (Klow.eq.0) Zlow=ZET
         Klow=1
         Zlow=dmax1(Zlow,ZET)
      endif
!* DENR - el.density, DENRI - ion
      ZET1=ZET
      DENRMIN=dmin1(1.d-19,DENR/100.) ! 23.08.97
      DENR1=dmax1(DENRMIN,DENR1) ! 25.11.96
      ZET=ZET1*(DENR/DENR1)**POWER
      if (Kup.ne.0.and.(ZET.gt.Zup.or.(I.gt.NI/2.and.DENR1.lt.DENR)))
     &   ZET=(ZET1+Zup)/2.
      if (Klow.ne.0.and.(ZET.lt.Zlow.or.(I.gt.NI/2.and.DENR1.gt.DENR)))
     &   ZET=(ZET1+Zlow)/2.
      if (abs(DENR/DENR1-1.).lt.EPS) goto 6
    5 continue
   55 print*,'Warning CHEMPOT: INCORRECT density:',DENR1,'  i.o.',DENR
    6 continue
      CHEMPOT=CMU
      return
      end

      function CHEMFIT(DENR,DDENR,TEMR)
      implicit double precision (A-H), double precision (O-Z)
      save
!*  New fit to the chemical potential of free electron gas
!*  All quantities are by default in relativistic units
!*  DENR - electron density, TEMR - temperature (input),
!*  DDENR - partial derivative of DENR with respect to the chem.pot.
!*  B0 - magnetic field (=\hbar\omega_c/mc^2)
!*                                                       Version 10.10.96
!*                                                       updated 15.03.99
      data PI/3.14159265/
      data AICH/.25954/,BICH/.072/,B1ICH/.858/
!*  These Ichimaru(93 Rev.Mod.Phys.65 N2 255) fit parameters "A,B,b" 
!*   are used here only for the derivative dn/d\mu, not for \mu.
      XSR=(3.d0*PI**2*DENR)**.33333333d0 ! Classical Fermi momentum
      if (XSR.gt..001) then
         THETA=TEMR/(dsqrt(1.d0+XSR**2)-1.d0) ! degeneracy parameter
      else
         THETA=2.d0*TEMR/XSR**2
      endif
!*   -------------- NEW FIT to the chemical potential ---------------   *
      THETA32=THETA*dsqrt(THETA)
      P01=12.+8./THETA32
      T=dexp(dmin1(THETA,40.d0))
      Y=(1.d0/T+1.612*T)/
     /(6.192*THETA**.0944/T+5.535*T*THETA**.698)
      P02=1.3656-Y ! 1.3656=2/\pi^{1/3}
      if (THETA.gt.1.e-5) then 
         P03=1.5d0/(T-1.d0)
      else
         P03=1.5d0/THETA
      endif
      UP=1.d0+P01*TEMR*P02+P03*sqrt(TEMR)
      DN=(1.d0+.5d0*TEMR/THETA)*(1.d0+P01*TEMR)
      CT=1.d0+UP/DN*TEMR
         F=.66666667d0/THETA/dsqrt(THETA)
         X1=FERINV(F,1) ! non-relativistic result
     -     -    1.5d0*dlog(CT) ! Relativistic fit
         CMU1=TEMR*X1 ! Fit to chemical potential w/o mc^2
         CHEMFIT=CMU1+1.d0
!*  Fit to the DERIVATIVE of chemical potential:
         THETAlog=dlog(THETA)
         THETACOR=(THETA32+0.29)**.666667d0
         TB=dexp(-B1ICH*THETAlog)   ! = THETA^{-b}
         TB1=TB/THETA           ! = THETA^{-(b+1)}
         DENOM=(1.d0+AICH*TB)
         UP=(AICH*TB1+BICH*dsqrt(TB1))
        DMT=-1.5d0-
     -   (B1ICH+1.d0)*(AICH*TB1+BICH/2.d0*dsqrt(TB1))/DENOM+
     +   UP*B1ICH*AICH*TB/DENOM**2 ! d (mu/T) / d (ln theta)
         DTP=-THETA/TEMR*XSR**2/dsqrt(1.d0+XSR**2) !d(ln theta)/d ln p_F
!* INCLUDE RELATIVISTIC CORRECTION :
         DMT=DMT-
     -     1.5d0*THETA32*TEMR**2/
     /     ((1.d0+2.d0*TEMR)*THETACOR+TEMR**2)/(THETA32+.27)
         DDENR=3.d0*DENR/(DMT*DTP)/TEMR ! d n_e / d mu
         return
      end

!*  ========================   Fitting density   =====================  *
      subroutine DENRFIT(B,TEMR,CMU,DENR,DDENR)
!*                                                       Version 10.10.96
!* Input: Zeff, B (magn. field in 4.4*10^13 G),
!*        TEMR and CMU (temperature and chemical potential in mc^2)
!* Output: DENR - el.number density [rel.units], 
!*         DDENR - derivative of DENR with respect to CMU [rel.un.]
      implicit double precision (A-H), double precision (O-Z)
      save
      data PRECLN/7.d0/,PI/3.14159265d0/
      TWOB=2.*B
      DELTE=TEMR*PRECLN
      EMAX=CMU+DELTE
!*        Upper boundary of the integration
      if (CMU.lt.1.d0) EMAX=1.d0+DELTE
      if (B.gt.0.) then
         CNL0=(EMAX**2-1.d0)/TWOB
!*        Max Landau number
      else
         CNL0=10000.
      endif
      if (CNL0.gt.1000.) then ! non-magn. case
         TINV=1.d0/TEMR  ! inverse relativistic temperature
         CHI=(CMU-1.d0)/TEMR
         W1=BLINW(1,TINV,CHI)
         W0=BLINW(0,TINV,CHI)
         G1=(W1+TINV*W0)*dsqrt(2.d0*TINV)
         DENR=TEMR**3*G1/PI**2
         DUMMY=CHEMFIT(DENR,DDENR,TEMR) ! determine DDENR
         goto 10
      endif
      NL0=CNL0
      EN2=1.d0
      S0=0.
      S1=0.
      do 1 NL=0,NL0
      if (NL.eq.0) then
         WEIGHT=1.d0
      else
         WEIGHT=2.d0
         EN2=EN2+TWOB
      endif
      EN=sqrt(EN2)
      X=(CMU-EN)/TEMR
      Y=EN/TEMR
      S0=S0+FitFERMI(X,Y,Derivat)*WEIGHT
      S1=S1+Derivat*WEIGHT
    1 continue
      CN=B/2.d0/PI**2
      DENR=CN*S0*TEMR
      DDENR=CN*S1
   10 continue
      return
      end

!*   ===================  AUXILIARY SUBROUTINES  ====================   *
      function EXPINT(XI,L) ! = e^XI E_{L+1}(XI)
      implicit double precision (A-H), double precision (O-Z)
      save
      data GAMMA/.5772156649D0/,Nrep/21/
      if (XI.ge.1.) then ! continued fraction
         CL=L
         CI=Nrep
         C=0.
        do I=Nrep,1,-1
         C=CI/(XI+C)
         C=(CL+CI)/(1.d0+C)
   11    CI=CI-1.d0
        enddo
         Q0=1.d0/(XI+C)
      else ! power series
         PSI=-GAMMA
        do K=1,L
   21    PSI=PSI+1.d0/K ! Psi(L+1)
        enddo
         Q0=0.
         CMX=1.d0 ! (-XI)^M/M!
         CL=L
         CM=-1.d0
         do 22 M=0,Nrep
         CM=CM+1.d0
        if (M.ne.0) CMX=-CMX*XI/CM ! (-XI)^M/M!
        if (M.ne.L) then ! Actually DQ=-deltaQ0
           DQ=CMX/(CM-CL)
        else
           DQ=CMX*(dlog(XI+1.d-20)-PSI)
        endif
         Q0=Q0-DQ
   22    continue
         Q0=exp(XI)*Q0
      endif
   50 continue
      EXPINT=Q0
      return
      end

      function BLINW(K,TINV,CHI)
!* The integrals (Blinnikov et al.'96): W_k(\alpha,\chi)=
!* \int_0^\infty u^{k+1/2}sqrt{u+2\alpha}/(1+e^{u-\chi})du
!* DIVIDED by sqrt(2*TINV), - Fermi integrals if TEMR --> 0.
!*                                                       Version 29.04.98
      implicit double precision (A-H), double precision (O-Z)
      save
      dimension AC(5,0:2),AU(5,0:2),AX(5),AXI(5),AH(5),AV(5),
     & AA(5,0:2)
      data AC/.37045057 d0, .41258437 d0,
     &        9.777982 d-2, 5.3734153 d-3, 3.8746281 d-5, ! c_i^0
     &        .39603109 d0, .69468795 d0, 
     &        .22322760 d0, 1.5262934 d-2, 1.3081939 d-4, ! c_i^1
     &        .76934619 d0, 1.7891437 d0, 
     &        .70754974 d0, 5.6755672 d-2, 5.5571480 d-4/ ! c_i^2
      data AU/.43139881 d0, 1.7597537 d0, 
     &        4.1044654 d0, 7.7467038 d0, 13.457678 d0, ! u_i^0
     &        .81763176 d0, 2.4723339 d0, 
     &        5.1160061 d0, 9.0441465 d0, 15.049882 d0, ! u_i^1
     &        1.2558461 d0, 3.2070406 d0, 
     &        6.1239082 d0, 10.316126 d0, 16.597079 d0/ ! u_i^2
      data AX/7.265351 d-2, .2694608 d0, 
     &        .533122 d0, .7868801 d0, .9569313 d0/ ! x_i
      data AXI/.26356032 d0, 1.4134031 d0, 
     &         3.5964258 d0, 7.0858100 d0, 12.640801 d0/ ! \xi_i
      data AH/3.818735 d-2, .1256732 d0, 
     &        .1986308 d0, .1976334 d0, .1065420 d0/ ! H_i
      data AV/.29505869 d0, .32064856 d0, 7.3915570 d-2, 
     &        3.6087389 d-3, 2.3369894 d-5/ ! \bar{V}_i
      data PI/3.14159265 d0/
      data KRUN/0/
      if (K.lt.0.or.K.gt.2) stop'BLINW: K out of range'
      KRUN=KRUN+1
      if (KRUN.eq.1) then ! initialize
        do J=0,2
        do I=1,5
           AA(I,J)=dexp(-AU(I,J))
        enddo
        enddo
      endif
      ALT=2.d0*TINV
      SQA=dsqrt(ALT)
      W=0.
      WAUX=0.
!*   -----------------------   (end of preamble)   -------------------- 
        if (CHI.lt.-80.) stop'BLINW: CHI out of range' ! (23.08.97)
         ECHI=dexp(-CHI)
        do I=1,5
           W=W+AC(I,K)*dsqrt(AU(I,K)+ALT)/(AA(I,K)+ECHI)
        enddo
        if (CHI.lt..4) goto 60
!*   ----------------------------------------------------------------   *
        WAUX=W
        W=0.
        do I=1,5
           CHIX=CHI*AX(I)
           ECHI=dexp(CHI*(AX(I)-1.d0))
           DS1=AH(I)*CHIX**K*dsqrt(CHI*(CHIX+ALT))/(1.d0+ECHI)
           CHIX=AXI(I)+CHI
           DS2=AV(I)*CHIX**K*dsqrt(CHIX*(CHIX+ALT))
           W=W+CHI*DS1+DS2
        enddo
        W=W+(WAUX-W)*dexp(-(CHI-.4)*3.d0) ! for smoothness (29.04.98)
!*   ----------------------------------------------------------------   *
        if (CHI.lt.9.) goto 60
        WAUX=W
        W=0.
         P=dsqrt(2.d0*CHI/TINV)
         SQR=dsqrt(CHI*(ALT+CHI))
        if (P.gt.1.d-5) then
          F0=(CHI+TINV)*SQR-TINV**2*dlog(1.d0+(CHI+SQR)/TINV)
          F0=F0/2.d0
        else
          F0=SQA*dsqrt(CHI)**3/1.5d0
        endif
        if (K.eq.0) then
           W=F0
           goto 50
        endif
        if (P.gt.1.d-4) then
          F1=SQR**3/3.d0-TINV*F0
        else
          F1=SQA*dsqrt(CHI)**5/2.5d0
        endif
        if (K.eq.1) then
           W=F1
           goto 50
        endif
        if (P.gt..05) then
          F2=(CHI*SQR**3-5.d0*TINV*F1)/4.d0
        else
          F2=SQA*dsqrt(CHI)**7/3.5d0
        endif
        if (K.eq.2) then
           W=F2
           goto 50
        endif
   50    continue
         W=W+PI**2/6.d0*CHI**K*((K+1)*CHI+TINV*(2*K+1))/SQR
        W=W+(WAUX-W)*dexp(-.5d0*dsqrt(CHI-9.d0)) ! smoothness (29.04.98)
   60 BLINW=W/SQA
      return
      end
      
      function FERINV(F,N)   ! Inverse Fermi intergals
!*  X_q(f)=F^{-1}_q(f) : H.M.Antia 93 ApJS 84 101; relative error 0.01%
!*  q=N-1/2=-1/2,1/2,3/2,5/2 (N=0,1,2,3)                 Version 07.03.95
      implicit double precision (A-H), double precision (O-Z)
      save
      dimension A0(0:3),A1(0:3),A2(0:3),A3(0:3),
     *B0(0:3),B1(0:3),B2(0:3),B3(0:3),
     *C0(0:3),C1(0:3),D0(0:3),D1(0:3),D2(0:3)
      data KRUN/0/
      if (N.lt.0.or.N.gt.3) stop'FERINV: Invalid subscript'
      if (KRUN.eq.0) then   ! Initialize
         KRUN=1
         I=0  ! X_{-1/2}
         A0(I)=785.16685
         A1(I)=-140.34065
         A2(I)=13.257418
         A3(I)=1.
         B0(I)=1391.7278
         B1(I)=-804.63066
         B2(I)=158.54806
         B3(I)=-10.640712
         C0(I)=.0089742174
         C1(I)=-.10604768
         D0(I)=.035898124
         D1(I)=-.42520975
         D2(I)=3.6612154
         I=1  ! X_{1/2}
         A0(I)=44.593646
         A1(I)=11.288764
         A2(I)=1.
         A3(I)=0.
         B0(I)=39.519346
         B1(I)=-5.7517464
         B2(I)=.26594291
         B3(I)=0.
         C0(I)=34.873722
         C1(I)=-26.922515
         D0(I)=26.612832
         D1(I)=-20.452930
         D2(I)=11.808945
         I=2  ! X_{3/2}
         A0(I)=35.954549
         A1(I)=13.90891
         A2(I)=1.
         A3(I)=0.
         B0(I)=47.795853
         B1(I)=12.133628
         B2(I)=-.23975074
         B3(I)=0.
         C0(I)=-.98934493
         C1(I)=.090731169
         D0(I)=-.68577484
         D1(I)=.063338994
         D2(I)=-.1163584
         I=3  ! X_{5/2}
         A0(I)=213.89693
         A1(I)=35.399035
         A2(I)=1.
         A3(I)=0.
         B0(I)=710.85455
         B1(I)=98.73747
         B2(I)=1.0677555
         B3(I)=-.011827987
         C0(I)=-.51891788
         C1(I)=-.0091723019
         D0(I)=-.36278896
         D1(I)=-.0061502672
         D2(I)=-.03367354
      endif
      if (F.lt.4.) then
         T=F
         UP=A0(N)+T*(A1(N)+T*(A2(N)+T*A3(N)))
         DOWN=B0(N)+T*(B1(N)+T*(B2(N)+T*B3(N)))
         X=dlog(F*UP/DOWN)
      else
        if (N.eq.0) then
           T=1./F**2
        else
           T=dexp(-dlog(F)/(.5+N))
        endif
         UP=C0(N)+T*(C1(N)+T)
         DOWN=D0(N)+T*(D1(N)+T*D2(N))
         X=UP/DOWN/T
      endif
      FERINV=X
      return
      end

      function FitFERMI(X1,Y1,Derivat)
!*   FitFERMI=F(x1,y)=\int_0^\infty e^{t-x1}/(1+e^{t-x1})^2\sqrt{t(t+2y)}
!*   Derivat=dF/dx1
!*   Accuracy: 0.57% in F, 2% in Derivat.                Version 25.11.96
      implicit double precision (A-H), double precision (O-Z)
      save
      data PI/3.14159265d0/,
     &CA1/.103d0/,CA2/.043d0/,CB1/.0802d0/,CB2/.2944d0/
      Y=dmin1(Y1,1.d10)
      X=X1-1./(1.+.623*Y**1.6031)
      if (X.lt.-40.) then
         Derivat=0.
         FitFERMI=0.
         return
      endif
      if (X.lt.30.) then
         XI=dlog(1.+exp(X))
      else
         XI=X
      endif
!*   Correction from 25.11.96:   --------------
      if (XI.eq.0.) then
         FitFERMI=0.
         Derivat=0.
         return
      endif
!*   ------------------------------------------
      ETA=sqrt(XI+2.*Y)
      U0=sqrt(PI)/2.+(CA1+CA2*XI**2)*sqrt(XI)
      D0=1.+CB1*sqrt(XI)+CB2*XI+CA2*XI**3
      if (Y.lt.1.d4) then
         CA0=0.9422*Y**1.7262
         UP=1.+Y+XI+CA0*ETA*U0
         DOWN=1.+XI+CA0*D0
      else
         UP=ETA*U0
         DOWN=D0
      endif
      if (X1.lt.30.) then
         XI1=dlog(1.+exp(X1))
      else
         XI1=X1
      endif
      F=UP/DOWN*XI1
      DFU=sqrt(PI*XI)/4.+CA1*(XI+Y)+CA2*XI**2*(3.*XI+5.*Y)
      DFD=CB1/sqrt(XI)/2.+CB2+3.*CA2*XI**2
      if (Y.lt.1.d4) then
         DFUP=(1.+CA0/sqrt(XI)/ETA*DFU)/UP
         DFDOWN=(1.+CA0*DFD)/DOWN
      else
         DFUP=1./sqrt(XI)/ETA*DFU/UP
         DFDOWN=DFD/DOWN
      endif
      DFDX=(DFUP-DFDOWN)/(1.+exp(-X))+1./(1.+exp(-X1))/XI1
*        log derivative
      DFDX=F*DFDX
*        derivative dF/dX
      FitFERMI=F
      Derivat=DFDX
      return
      end

!*   ----------------------------------------------------

      subroutine TAUEESY(X,TEMR,TAUEE)
      implicit double precision (A-H), double precision (O-Z)
!* Relaxation time of electron-electron collisions
!*   according to Shternin & Yakovlev (2006),
!*   corrected in the nondegenerate regime so as to match Lampe (1968)
!* X - rel.parameter, TEMR - temperature [rel.un.]
!*                                                       Version 02.11.16
!* Stems from v.17.07.06. Diff.: use of asymptote at large Y
      E=dsqrt(1.d0+X**2)
      V=X/E ! =u in Eq.(12) of SY'06
      Y=.0963913/TEMR*X*dsqrt(V) ! =\theta; .0963913=2\sqrt{\alpha/\pi}
      if (Y.gt.1.d8) then ! Eq.(18)
         CIL=20.4013123/(V*Y**3) ! 20.4013123=\pi^5/15
         CIT=2.404*V/Y**2
         CITL=18.52*(V/Y**8)**(1.d0/3.d0)
      else ! Eqs.(21)-(23)
         C1=.123636+.016234*V**2
         C2=.0762+.05714*V**4
         A=12.2+25.2*V**3
         C=A*dexp(C1/C2)
         YV=Y*V
         CIL=dlog(1.d0+128.56/(37.1*Y+10.83*Y**2+Y**3))*
     *     (.1587-.02538/(1.+.0435*Y))/V
         CIT=V**3*dlog(1.d0+C/(A*YV+YV**2))*
     *     (2.404/C+(C2-2.404/C)/(1.d0+.1*YV))
         CILT=V*dlog(1.d0+C/(A*Y+10.83*YV**2+YV**(8.d0/3.d0)))*
     *   (18.52*V**2/C+(C2-18.52*V**2/C)/(1.d0+.1558*Y**(1.d0-.75d0*V)))
      endif
      FI=CIL+CIT+CILT
      FREQ=.00021381*X*Y*dsqrt(V)*FI ! .00021381=6\alpha^{3/2}/\pi^{5/2}
!* Correction for partial degeneracy:
      if (X.gt..001) then
         THETA=TEMR/(E-1.d0) ! degeneracy
      else
         THETA=2.d0*TEMR/X**2
      endif
      T=25.*THETA
      TAUEE=(1.d0+T+.4342*dsqrt(THETA)*T**2)/(1.d0+T**2)/FREQ
      return
      end

      subroutine CONDI(TEMR,DENRI,Zion,CMI,RKAPi)
!* ion thermal conductivities in the outer crust: the approximation
!*  of Chugunov & Haensel (2007, MNRAS 381, 1143)
!* Input: TEMR - temperature, DENRI - number density of ions
!*        Zion and CMI - ion charge and mass numbers
!* Output: RKAPi - ion thermal conduction
!* All quantities are in the relativistic units: \hbar = m_e = c = 1.
!*                                                       Version 06.10.08
!*                                                        Modif. 19.02.09
!*                                                         Corr. 10.09.19
      implicit double precision (A-H), double precision (O-Z)
      save
      data Uminus1/2.8/,Uminus2/13./,AUM/1822.9/,BOHR/137.036/,KAPi/1/
!* Dimensional quantities are in the relativistic units (m_e=\hbar=c=1)
!*        Uminus1,Uminus2 - dimensionless frequency moments of phonons
!*        AUM - atomic mass unit divided by the electron mass
!*        BOHR - radius of the first Bohr orbit in the rel.units
      parameter(PI=3.14159265d0,EPS=1.d-8,TINY=1.d-99)
!*   ----------------------   Preliminaries   -----------------------   *
      DENR=DENRI*Zion ! number density of electrons
      XSR=(3.d0*PI**2*DENR)**.33333333d0 ! x_r - density parameter
      PCL=XSR ! classical Fermi momentum; equality due to B=0.
      VCL=PCL/dsqrt(1.d0+PCL**2)
      SPHERION=(.75d0/PI/DENRI)**.33333333d0 ! Ion sphere radius
      QBZ=PCL*(2.d0/Zion)**.33333333d0 ! q_{BZ}
      GAMI=Zion**2/(BOHR*SPHERION*TEMR) ! corrected 06.10.08
      TRP=Zion/GAMI*dsqrt(CMI*AUM*SPHERION/3.d0/BOHR) ! =T/T_pi
      OMPI=TEMR/TRP ! ion plasma frequency
      if (GAMI.lt.TINY**.2) then
         print*,'CONDI: GAMI is too low: kappa>HUGE.'
         print*,'To continue, press <Enter>.'
      endif
!*   ---------------------- reduced ion heat capacity:
      call HLfit8(1.d0/TRP,F,U,CV,S,1)
!*   ---------------------- ion-electron:
      A0=1.683*sqrt(PCL/CMI/Zion) ! zero-vibr.param.(Baiko&Yakovlev95)
      WDW=A0*(.5d0*Uminus1*exp(-9.1*TRP)+Uminus2*TRP) ! DW fact.(B&Y'95)
      W=WDW ! because it is for outer crust
      F=.014+.03/(1.d0+dexp(.2/TRP)) ! CH'07, Eq.(45)
      Y=QBZ/(2.d0*PCL)
      WY2=W*Y**2
      if (W.gt.EPS) then
         EXPW=dexp(-W)
         EXPWY2=dexp(-WY2)
         CL=(EXPINT(WY2,0)*EXPWY2-EXPINT(W,0)*EXPW-
     -     VCL*(EXPWY2-EXPW)/W)/2.d0
      else ! small-W asymptote
         CL=dlog(1.d0/Y)-VCL*(1.d0-Y**2)/2.d0
      endif
      FLEN=8.32d5*SPHERION/(1.d0+PCL**2)/CL*F/Zion*
     *  dsqrt((XSR/1.00884)**3*CMI/Zion) ! Eq.(43)
      CS=OMPI/(3.d0*QBZ) ! sound speed
      RKAPie=CV*DENRI*CS*FLEN/3.d0 ! Eq.(42) of CH'07
!*   ---------------------- ion-ion:
      if (TRP.lt.1.d0/dlog(1.d0/TINY)) then
        if (KAPi.ne.-1) then ! 19.02.09
           write(*,'(''CONDI: T is too low: KAPii > HUGE.'',
     *       ''Press <Enter> to continue.''$)')
           read(*,'(A)')
           KAPi=-1
        endif
         RKAPi=RKAPie
      else
         RKAP0=OMPI*DENRI*SPHERION**2
         RKAPii=RKAP0*dsqrt(16.d0/GAMI**5/ ! corr.10.09.19 (1 --> 16)
     /     dlog(2.d0+.57735/dsqrt(GAMI)**3)**2+
     +     0.16+(GAMI/77.)**2*dexp(.6666667d0/TRP)) ! Eq.(26) of CH'07
!*   ---------------------- total ion:
         RKAPi=1.d0/(1.d0/RKAPii+1.d0/RKAPie)
      endif
      return
      end
      
      subroutine HLfit8(eta,F,U,CV,S,LATTICE)
!*                                                       Version 15.02.08
!* Baiko, Potekhin, & Yakovlev (2001). Stems from HLfit v.20.03.01
!* Fit to thermal part of the thermodynamic functions.
!* Zero-point lattice quantum energy 1.5u_1\eta INCLUDED (unlike HLfit)
!* Input: eta=Tp/T
!* Output: F and U (normalized to NkT),
!*   CV and S (normalized to Nk) in the HL model for bcc Coulomb lattice
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(EPS=1.d-5)
      if (LATTICE.eq.1) then ! bcc lattice
         CLM=-2.49389 ! 3*ln<\omega/\omega_p>
         U1=.5113875
         ALPHA=.265764
         BETA=.334547
         GAMMA=.932446
         A1=.1839
         A2=.593586
         A3=.0054814
         A4=5.01813d-4
         A6=3.9247d-7
         A8=5.8356d-11
         B0=261.66
         B2=7.07997
         B4=.0409484
         B5=.000397355
         B6=5.11148d-5
         B7=2.19749d-6
         C9=.004757014
         C11=.0047770935
      elseif (LATTICE.eq.2) then ! fcc lattice
         CLM=-2.45373
         U1=.513194
         ALPHA=.257591
         BETA=.365284
         GAMMA=.9167070
         A1=.0
         A2=.532535
         A3=.0
         A4=3.76545d-4
         A6=2.63013d-7
         A8=6.6318d-11
         B0=303.20
         B2=7.7255
         B4=.0439597
         B5=.000114295
         B6=5.63434d-5
         B7=1.36488d-6
         C9=.00492387
         C11=.00437506
      else
         stop'HLfit: unknown lattice type'
      endif
      if (eta.gt.1.d0/EPS) then ! asymptote of Eq.(13) of BPY'01
         U=3.d0/(C11*eta**3)
         F=-U/3.d0
         CV=4.d0*U
        goto 50
      elseif (eta.lt.EPS) then ! Eq.(17) of BPY'01
         F=3.d0*dlog(eta)+CLM-1.5d0*U1*eta+eta**2/24.d0 
         U=3.d0-1.5d0*U1*eta+eta**2/12.d0
         CV=3.d0-eta**2/12.d0
         goto 50
      endif
      B9=A6*C9
      B11=A8*C11
      UP=1.d0+A1*eta+A2*eta**2+A3*eta**3+A4*eta**4+A6*eta**6+A8*eta**8
      DN=B0+B2*eta**2+B4*eta**4+B5*eta**5+B6*eta**6+
     +  B7*eta**7+B9*eta**9+B11*eta**11
      EA=dexp(-ALPHA*eta)
      EB=dexp(-BETA*eta)
      EG=dexp(-GAMMA*eta)
      F=dlog(1.d0-EA)+dlog(1.d0-EB)+dlog(1.d0-EG)-UP/DN ! th.free en./NT
      UP1=A1+
     + 2.d0*A2*eta+3.d0*A3*eta**2+4.d0*A4*eta**3+
     + 6.d0*A6*eta**5+8.d0*A8*eta**7
      UP2=2.d0*A2+6.d0*A3*eta+12.d0*A4*eta**2+
     + 30.d0*A6*eta**4+56.d0*A8*eta**6
      DN1=2.d0*B2*eta+4.d0*B4*eta**3+5.d0*B5*eta**4+6.d0*B6*eta**5+
     +  7.d0*B7*eta**6+9.d0*B9*eta**8+11.d0*B11*eta**10
      DN2=2.d0*B2+12.d0*B4*eta**2+20.d0*B5*eta**3+30.d0*B6*eta**4+
     +  42.d0*B7*eta**5+72.d0*B9*eta**7+110.d0*B11*eta**9
      U=ALPHA*EA/(1.d0-EA)+BETA*EB/(1.d0-EB)+GAMMA*EG/(1.d0-EG)-
     -  (UP1*DN-DN1*UP)/DN**2 ! int.en./NT/eta
      CV=ALPHA**2*EA/(1.d0-EA)**2+BETA**2*EB/(1.d0-EB)**2+
     +  GAMMA**2*EG/(1.d0-EG)**2+
     +  ((UP2*DN-DN2*UP)*DN-2.d0*(UP1*DN-DN1*UP)*DN1)/DN**3 ! cV/eta^2
      U=U*eta
      CV=CV*eta**2
   50 continue
      S=U-F
!* Add zero-point lattice energy:
      E0=1.5d0*U1*eta
      U=U+E0
      F=F+E0
      return
      end
	
