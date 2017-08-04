C
c     here a subroutine 
c     to calculate travel-time tables for hyposat.f
c
C              JOHANNES SCHWEITZER
C
C     (developed from laufps.f (version 28. April 1997))
c
c     11. June  1997
c
c     Nov 14, 1997  indph corrected for multiples and surface focus.
c
c     March 1999    changes due to CRUST 5.1 model input
c
c     October 2000  including converted surface reflections and
c                   converted reflections from the discontinuities 
c                   Conrad and Moho
c
c     last changes/corrections 27. July 2001
c
c

      subroutine ttloc(ho,dis,czo1,nphas2,ttc,dtdd,dtdh,dpdh,dddp,
     +                 phcd,rmax,typctl,ierr,indph)

c
c     input:
c 
c             ho      = source depth in km (fixed or not, see fixho)
c
c             dis     = receiver distance in deg
c
c             czo1    = character to inform about fixed or not fixed 
c                       depth (D == no-fixed, other == fixed)
c
c             typctl  = verbosity level
c
c
c                       only P-onsets at receiver
c
c             indph   = 10000 for only direct P-onsets
c                     = 11000 for P and pP
c                     = 11100 for P, pP, and PP 
c                     = 11110 for P, pP, PP, and PbP and PmP 
c                     = 11111 for P, pP, PP, PbP, PmP, and sP and SmP
c
c                       only S-onsets at receiver
c
c             indph   = 20000 for only direct S-onsets
c                     = 22000 for S and sS
c                     = 22200 for S, sS, and SS 
c                     = 22220 for S, sS, SS, and SbS and SmS 
c                     = 22222 for S, sS, SS, SbS, SmS, and pS and PmS 
c
c                       P- and S-onsets at receiver
c
c              indph  = 30000 for P and S (direct only)
c                     = 33000 for P, S, pP, and sS
c                     = 33300 for P, S, pP, sS, PP, and SS
c                     = 33330 for P, S, pP, sS, PP, SS, PbP, PmP, SbS,
c                                and Sms
c                     = 33333 for P, S, pP, sS, PP, SS, PbP, PmP, SbS.
c                                 Sms, pS, sP, PmS, and SmP
c
c	 in common MODEL :
c 
c             v0(1,i) =  P velocity in layer i
c             v0(2,i) =  S velocity in layer i
c
c             h(i)    =  depth of layer i
c
c             elev    =  topograhic elevation at modelled point 
c
c             elat    =  latitude  to get the right CRUST 5.1 model 
c                        parameters
c             elon    =  longitude to get the right CRUST 5.1 model 
c                        parameters
c
c             elev2   =  topograhic elevation at second modelled point 
c
c             elat2   =  latitude  to get the right CRUST 5.1 model 
c                        parameters at second modelled point
c             elon2   =  longitude to get the right CRUST 5.1 model 
c                        parameters at second modelled point
c
c             imo     <  0 no new model reading
c                     =  1 reading of model from filloc
c                 (   =  2 CRUST 5.1 only used for travel-time 
c                          corrections  )
c                     =  3 using model from CRUST 5.1
c                     =  4 dito + CRUST 5.1 used for travel-time 
c                          corrections
c
c             jmod    =  number of layers
c
c             iread   =  0 data of CRUST 5.1 model must be read in
c                     =  1 data of CRUST 5.1 model are already read
c
c             filloc  =  filename for file with local velocity model
c
c             azo(i)  =  Conrad/Moho indicator
c
c             mtyp    =  CRUST 5.1 model type
c            
c
c     output 
c
c             nphas2  = number of found onsets at receiver distance
c
c             phcd    = array with names of found onsets
c
c             ttc     = travel times of onsets in [sec]
c
c             dtdd    = ray parameters of onsets in [sec/deg]
c
c             dtdh    = partial derivative d(travel time)/d(ho) in 
c                       [sec/km]
c
c             dddp    = partial derivative d(ray parameter)/d(dis) in 
c                       [(sec/deg**2]
c
c             dpdh    = partial derivative d(ray parameter)/d(ho) in 
c                       [(sec/deg)/km]
c
c             ierr    = 0     everything o.k.
c                       else  some error occurred.
c
c             rmax    =  maximum distance for which this model shall be
c                        used (read in from model file).
c
c	 in common MODEL :
c 
c             v0(1,i) =  P velocity in layer i
c             v0(2,i) =  S velocity in layer i
c
c             h(i)    =  depth of layer i
c
c             elev    =  topograhic elevation at modelled point 
c
c             elat    =  latitude  to get the right CRUST 5.1 model 
c                        parameters
c             elon    =  longitude to get the right CRUST 5.1 model 
c                        parameters
c
c             elev2   =  topograhic elevation at second modelled point 
c
c             elat2   =  latitude  to get the right CRUST 5.1 model 
c                        parameters at second modelled point
c             elon2   =  longitude to get the right CRUST 5.1 model 
c                        parameters at second modelled point
c
c
c             imo     <  0 no new model reading
c                     =  1 reading of model from filloc
c                 (   =  2 CRUST 5.1 only used for travel-time 
c                          corrections  )
c                     =  3 using model from CRUST 5.1
c                     =  4 dito + CRUST 5.1 used for travel-time 
c                          corrections
c
c             jmod    =  number of layers
c
c             iread   =  0 data of CRUST 5.1 model must be read in
c                     =  1 data of CRUST 5.1 model are already read
c
c             filloc  =  filename for file with local velocity model
c
c             azo(i)  =  Conrad/Moho indicator
c
c             mtyp    =  CRUST 5.1 model type
c
c             zmax    =  maximum depth for which this model can be
c                        used (read in from model file).
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      REAL*4    ho,dis,rmax
      Character*1 czo1
      INTEGER   nphas2,ierr,indph,typctl,trimle

      include 'ttimes.h'

      dimension dpdh(mphas),dddp(mphas)

c
c     ML        maximum number of allowed layers in model
c
c     NP        maximum number of calculated (defined) phases
c
c               (change also parameter ml and np in subroutine reflex )
c               (                      np in function phnum )
c               (                      ml in subroutine get_mod_mlm )
c               (                      ml in function crust )
c
c
      PARAMETER (ml=61,np=40)

      SAVE      

      DIMENSION Z(ml),H(ml),RR(2),TT(2),V0(2,ml),
     *          V(2,ml),G(2,ml),V2(2,ml),ttp(np,3,3),
     *          ppp(np,3,3),ion(np,3,3),PA1(2),VHQ(2),
     *          FA(2,ml),zo(3),pa(2),del(3),ndisc(ml),
     *          pm(3,3),tm(3,3),indx(np)

      real      tti(np)

      CHARACTER AZ(ml)*4,phase*8,phas1*8,phlist(np)*8,
     *          azo(ml)*4,line*34,mtyp*3,filloc*80,mpath*80

      COMMON    /MODEL/  v0,h,elev,elat,elon,zmax,elat2,elon2,elev2,
     *                   jmod,iread,imo,mtyp
      COMMON    /MODELC/ filloc,azo,mpath

      COMMON    /REF/ rr,tt,pa,fa,ttp,ppp,ion,IB,IQQ,IQL,jh1,
     *                PIM,AA,VHQ,PA1,PI,V,G,V2,del,phase,conv

      common    /phasel/ phlist

      integer   phnum, iqq

      logical   fixho, kp, ks, zof, surf, mul, surfc, conv

      DATA      PHLIST /'Pg','Pb','Pn','P','pPg','pPb','pPn','pP',
     *                  'PbP','PmP','PgPg','PbPb','PnPn','PP',
     *                  'pSg','pSb','pSn','pS','PbS','PmS',
     *                  'Sg','Sb','Sn','S','sSg','sSb','sSn','sS',
     *                  'SbS','SmS','SgSg','SbSb','SnSn','SS',
     *                  'sPg','sPb','sPn','sP','SbP','SmP'         /

      PI=4.d0*DATAN(1.d0)
      PIM=PI/180.d0
      RE=6371.d0
      AA=PIM*RE

      IB = 20
      IBN= IB*10

      if (typctl.gt.5) then
	 print *,'imo ',imo
      endif

      ierr = 0

      if (imo.le.2) then

         lfil = trimle(filloc)

         OPEN (UNIT=15,FILE=filloc(1:lfil))
C
C        Loop to get the models of P and S velocities
C

15       read(15,*) rmax

         I=0

50       I=I+1
         IF(I.GT.60) THEN
             write(*,'('' Model contains too much layers (> 60)'')')
	     ierr = 99
             close (15)
             GO TO 9000
         ENDIF

         READ(15,'(A)',err=55,end=56) line

         if (line(31:34).ne.'    ' .and. line(31:34).ne.'MOHO' .and.
     *	     line(31:34).ne.'CONR' ) go to 55

         READ(line,'(3F10.3,A4)',err=55,end=56)
     *             H(I),V0(1,i),V0(2,i),azo(I)

         if(h(i) + v0(1,i) + v0(2,i).eq.0.d0) go to 56
         go to 50

55       write(*,'('' Read ERROR for file: '',a)') filloc
         ierr = 99
         close (15)
         go to 9000

56       jmod=I-1
         close (15)

	 elev = 0.d0

         imoh = 600
         icon = 600
         ipd  = 0
         isd  = 0

         zoold = -999.d0
         zmax  = h(jmod)

	 go to 100

      endif

      if (imo.ge.3)  then

	 if(mtyp.ne.'MLM') go to 9000

	 itrue = 0
	 inum = 2
         call get_mod_mlm(itrue,inum,typctl,ierr)

	 if(ierr.ne.0) then
	    ierr = 99
	    go to 9000
	 endif

         imoh = 600
         icon = 600
         ipd  = 0
         isd  = 0

	 zoold = -999.d0
	 rmax  = 6.

	 go to 100

      else

	 print *,' No local/regional model defined! '
	 ierr = 99
	 go to 9000

      endif

100   continue

c
c     reset onset table
c
      if(ho.gt.sngl(zmax)) then
	 ierr = 99
	 print *,'Depth greater than maximum model depth'
	 go to 9000
      endif

      do 110 i=1,np
      do 110 j=1,3
      do 110 k=1,3
110   ion(i,j,k) = 0

      zof = .false.

      if (czo1.eq.'D') then
	fixho=.false.
	nzo  = 3
	jh1  = 2
	jh2  = 3
        zo(2) = dble(ho)
        zo(1) = zo(2) - 1.d0
	if(zo(1).lt.0.d0) zo(1)=0.d0
        zo(3) = zo(2) + 1.d0
      else
	fixho=.true.
	nzo   = 1
	jh1 = 1
	jh2 = 1
	zo(1) = dble(ho)
	if (zo(1).eq.zoold) zof=.true.
      endif

      del(2) = dble(dis)
      del(1) = del(2)-0.01d0
      if(del(1).lt.0.d0) del(1)=0.d0
      del(3) = del(2)+0.01d0

      do 8000 iql = 1,nzo

C
C     First work to establish the model
c
c     We will calculate at first P than S phases
c
C     (k-loop)
C

      DO 810 K=1,2

      if(k.eq.1) then
	kp=.true.
	ks=.false.
      else
	kp=.false.
	ks=.true.
      endif

      if (zof) go to 810

      ij = 0
      izo = 0

      DO 500 I=1,jmod

      i2  = i + 1
      ij = ij + 1

      D=V0(K,I)
      B=RE/(RE-H(I))

      if(zo(iql).eq.h(i).and.izo.eq.0) then
	 IQQ = IJ
	 izo = 1
      endif
      az(ij) = azo(i)

      V(K,IJ)=D*B
      IF(V(1,IJ).GE.10.D0 .AND. IPD.NE.0)  IPD  = I
      IF(V(2,IJ).GE.4.7D0 .AND. ISD.NE.0)  ISD  = I
      Z(IJ)=RE*DLOG(B)

      if(h(i2).gt.zo(iql).and.h(i).lt.zo(iql).and.izo.eq.0) then

         D=(zo(iql)-h(i))*(V0(K,I2)-V0(K,I))/(h(i2)-h(i)) + V0(K,I)
         B=RE/(RE-zo(iql))
	 ij = ij + 1
	 IQQ = IJ
	 az(ij) = ''
	 izo = 1
C
         V(K,IJ)=D*B
         IF(V(1,IJ).GE.10.D0 .AND. IPD.NE.0)  IPD  = I
         IF(V(2,IJ).GE.4.7D0 .AND. ISD.NE.0)  ISD  = I
         Z(IJ)=RE*DLOG(B)

      endif

500   continue

      j = IJ
      m = j - 1

      DO 800 I=1,M

      I2=I+1
      V2(K,I)=V(K,I2)

      IF(AZ(I).EQ.'CONR')  ICON = I
      IF(AZ(I).EQ.'MOHO')  IMOH = I


      IF(V2(K,I).EQ.V(K,I)) THEN
         V2(K,I)=1.0001d0*V(K,I)
         V(K,I2)=V2(K,I)
      ENDIF

      zdiff=Z(I2)-Z(I)
      ndisc(i) = 0
      IF(zdiff.EQ.0.d0)  then
	 zdiff=0.0001d0
	 z(i2)= z(i2) + zdiff
	 ndisc(i) = 1
      endif

      G(K,I)=(V2(K,I)-V(K,I))/zdiff

800   continue

810   continue

      DO 7500 K=1,2

      if(k.eq.1) then
        kp=.true.
        ks=.false.
      else
        kp=.false.
        ks=.true.
      endif

      indph1 = indph  / 10000
      indphx = indph  - indph1*10000
      indph2 = indph  / 1000
      indphx = indph  - indph2*1000
      indph3 = indphx / 100
      indphx = indphx - indph3*100
      indph4 = indphx / 10
      indph5 = indphx - indph4*10

      if(indph1.ne.3) then
        if(kp.and.indph1.ne.1) go to 7500
        if(ks.and.indph1.ne.2) go to 7500
      endif

      VHQ(K)=V(K,IQQ)*V(K,IQQ)
      PA1(K)=(RE-Z(IQQ))*PIM/V(K,IQQ)


      IF(IQQ.EQ.1)  GO TO 1000

C
C     direct waves ( if source deeper than 0. )
c
c     plus defining requested direct phases
C

      if(kp) phase(1:1)='P'
      if(ks) phase(1:1)='S'

      if(iqq.le.imoh) then
	 if(iqq.le.icon) phase(2:)='g     '
         if(iqq.gt.icon) phase(2:)='b     '
      else
         phase(2:)='n     '
      endif

      if((kp .and. iqq.ge.ipd .and. ipd.ne.0) .or.
     +   (ks .and. iqq.ge.isd .and. isd.ne.0) ) phase(2:)='      '

      CALL REFLEX(IQQ,K)

C
C     body waves 
C

1000  imul=1111

      mul   = .false.
      surf  = .false.
      surfc = .false.
      conv  = .false.

      IQ4=0
      MULT=0
      FFA=0.d0


1100  IQ5=1
      IQ6=0
      IF(IQQ.GT.1)  IQ5=IQQ
      IF(IQ4.GT.1)  IQ5=IQ4
      IF(imul.LT.M.AND.imul.GT.IQ6) IQ6=imul
      IQ5=MAX0(IQ6,IQ5)
      VMAX=V(K,IQ5)

      DO 1300 I=1,IQ5
      FA(1,I)=2.d0
      FA(2,I)=0.d0
      IF(I.GE.imul) FA(1,I)=FFA
      IF(I.LT.IQQ)  THEN
         FA(1,I)=FA(1,I)-1.d0
         GO TO 1300
      ENDIF
      IF(I.LT.IQ4.AND.surf) FA(1,I)=FA(1,I)+1.d0
      IF(I.LT.IQ4.AND.surfc) then
	 FA(2,I)=FA(2,I)+1.d0
         IF(k.eq.1 .and. VMAX.LT.V(2,I)) VMAX=V(2,I)
         IF(k.eq.2 .and. VMAX.LT.V(1,I)) VMAX=V(1,I)
      endif
      IF(VMAX.LT.V(K,I)) VMAX=V(K,I)
1300  continue
C
C
      DO 3000 I=IQ5,M

      ib2 = ib
      ib1 = 1
      ibm = 0

      FA(1,I)=2.d0
      IF(I.GE.imul) FA(1,I)=FFA

      if (i.eq.ndisc(i)) go to 3000

      if(kp) phase(1:1)='P'
      if(ks) phase(1:1)='S'

      if(i.lt.imoh) then
	 if(i.le.icon) phase(2:)='g     '
         if(i.gt.icon) phase(2:)='b     '
      else
         phase(2:)='n     '
      endif

      if((kp .and. i.ge.ipd .and. ipd.ne.0) .or.
     +   (ks .and. i.ge.isd .and. isd.ne.0) ) phase(2:)='      '

      D=V2(K,I)
      IF(D.LE.VMAX)  GO TO    3000
      IF(imul.LT.M.AND.I.LT.imul) GO TO 2600
      C=V(K,I)
      IF(C.LT.VMAX) C=VMAX
1350  G1=DBLE(IB2-1)
      B=(D-C)/G1

      DO 2500 I2=ib1,IB2
      G2=dble(i2-1)
      VV=C+G2*B
      R=0.D0
      T=0.D0
C
      DO 2000 KK=1,I

      E=V(K,KK)
      G1=E/VV
      P=DSQRT(DABS(1.D0-G1*G1))
      O=1.d0/G(K,KK)

      IF(KK.GE.I)  THEN
         F=VV
         Q=0.d0
      else
         F=V2(K,KK)
         G3=F/VV
         Q=DSQRT(DABS(1.D0-G3*G3))
      ENDIF

      R=R+FA(1,KK)*(P-Q)*O
      DT=FA(1,KK)*DLOG(F*(1.D0+P)/(E*(1.D0+Q)))*O
      T=T+DT

2000  CONTINUE

c
c     extension for converted surface reflections (i.e. pS or sP)
c
      if(surfc) then

	 if(k.eq.1) kc = 2
	 if(k.eq.2) kc = 1

         DO 2001 KK=1,IQ4

         if(FA(2,KK).lt.0.9d0) go to 2001

         E=V(KC,KK)
         G1=E/VV
         P=DSQRT(DABS(1.D0-G1*G1))
         O=1.d0/G(KC,KK)

         F=V2(KC,KK)
         G3=F/VV
         Q=DSQRT(DABS(1.D0-G3*G3))

         R=R+FA(2,KK)*(P-Q)*O
         DT=FA(2,KK)*DLOG(F*(1.D0+P)/(E*(1.D0+Q)))*O
         T=T+DT

2001     CONTINUE

      endif

      RR(2) = R*VV/AA
      TT(2) = T

      PA(2) = AA/VV

      phas1 = phase
      iphase = trimle(phase)

      if(mul) then
	 phas1 = phase(1:iphase)//phase(1:iphase)
      endif

      if(surf) then
	 if(kp) phas1 = 'p' // phase(1:iphase)
	 if(ks) phas1 = 's' // phase(1:iphase)
      else if(surfc) then
	 if(kp) phas1 = 's' // phase(1:iphase)
	 if(ks) phas1 = 'p' // phase(1:iphase)
      endif

      iph=phnum(phas1)

      do 2400 klr = 1,3

      if(iql.ne.jh1 .and. klr.ne.2) go to 2400

      del1 = del(klr)

      IF (RR(2).EQ.del1) THEN

	 if(ibm .eq. 0 ) then
	    ib2 = ibn
	    ib1 = (i2-1)*10
	    if(ib1.le.0) ib1 = 1
	    ibm = 1
	    go to 1350
	 endif
            
	 ion(iph,iql,klr) = ion(iph,iql,klr)+1

	 if(ion(iph,iql,klr).eq.1) then
	    ttp(iph,iql,klr)=TT(2)
	    ppp(iph,iql,klr)=PA(2)
	 else
	    if(TT(2).lt.ttp(iph,iql,klr)) then
	       ttp(iph,iql,klr)=TT(2)
	       ppp(iph,iql,klr)=PA(2)
	    endif
	 endif
         GO TO 2400
      ENDIF

      if(i2.le.1) go to 2400

      FCT1=DEL1-RR(1)
      FCT2=DEL1-RR(2)
      IF(FCT1*FCT2.LT.0.d0) THEN

	 if(ibm .eq. 0 ) then
	    ib2 = ibn
	    ib1 = (i2-1)*10
	    ibm = 1
	    go to 1350
	 endif
            
         FCT3=FCT1/(RR(2)-RR(1))
         TT1=FCT3*(TT(2)-TT(1))+TT(1)
         PA2=FCT3*(PA(2)-PA(1))+PA(1)

	 ion(iph,iql,klr) = ion(iph,iql,klr)+1

	 if(ion(iph,iql,klr).eq.1) then
	    ttp(iph,iql,klr)=TT1
	    ppp(iph,iql,klr)=PA2
	 else
	    if(TT1.lt.ttp(iph,iql,klr)) then
	       ttp(iph,iql,klr)=TT1
	       ppp(iph,iql,klr)=PA2
	    endif
	 endif
      ENDIF

2400  continue

      rr(1) = rr(2)
      tt(1) = tt(2)
      pa(1) = pa(2)

2500  continue
C
2600  VMAX=D
3000  CONTINUE

      if(mul) go to 3500

C
C     Now the surface reflections of the body waves will be
C     calculated (i.e. pP,sS...).
C
      IF(surf .or. iqq.eq.1) then
	 if(surf) then
	   IQQ=IQ4
	   IQ4=1
	   surf=.false.
	 endif
	 GO TO 3100
      endif

      if(indph2.ne.3) then
        if(kp.and.indph2.ne.1) go to 3100
        if(ks.and.indph2.ne.2) go to 3100
      endif

      surf = .true.
      IQ4=IQQ
      IQQ=1

      GO TO 1100

C
C     Now the converted surface reflections of the body waves 
C     will be calculated (i.e. pS,sP...).
C

3100  CONTINUE

      IF(surfc .or. iqq.eq.1) then
	 if(surfc) then
	   IQQ=IQ4
	   IQ4=1
	   surfc=.false.
	 endif
	 GO TO 3200
      endif

      if(indph5.ne.3) then
        if(kp.and.indph5.ne.1) go to 3200
        if(ks.and.indph5.ne.2) go to 3200
      endif

      surfc = .true.
      IQ4=IQQ
      IQQ=1

      GO TO 1100

C
C     Now the multiple phases will be done (e.g. PgPg, PnPn or SnSn...)
C     At the moment only single multiples can be calculated (i.e. e.g. 
C     no PPP or SnSnSn ...).
C

3200  continue

      if(indph3.ne.3) then
        if(kp.and.indph3.ne.1) go to 3500
        if(ks.and.indph3.ne.2) go to 3500
      endif

      imul = 1
      mult = 1
      mul = .true.

      FFA=dble(MULT*2+2)
      GO TO 1100
C
C     End of the body-phase and direct-wave loop
C

3500  CONTINUE

C
C     Reflections for P and S from the two possible layers:
C     the Conrad-discontinuity and the Mohorivicic-discontinuity.
C

      if(indph4.ne.3) then
        if(kp.and.indph4.ne.1) go to 6500
        if(ks.and.indph4.ne.2) go to 6500
      endif

      DO 6000 I=IQQ,J

      IF(I.NE.ICON .and. I.NE.IMOH) GO TO 6000

      if (i.eq.icon) then
	if(kp) phase = 'PbP'
        if(ks) phase = 'SbS'
      else if (i.eq.imoh) then
	if(kp) phase = 'PmP'
        if(ks) phase = 'SmS'
      endif

      I2=I
      CALL REFLEX(I2,K)

6000  CONTINUE


6500  continue
c
C
C     Converted Reflections for P and S from the two possible layers:
C     the Conrad-discontinuity and the Mohorivicic-discontinuity.
C     PbS, SbP, PmS, SmP
C

      if(indph5.ne.3) then
        if(kp.and.indph5.ne.1) go to 7500
        if(ks.and.indph5.ne.2) go to 7500
      endif

      conv = .true.

      DO 7000 I=IQQ,J

      IF(I.NE.ICON .and. I.NE.IMOH) GO TO 7000

      if (i.eq.icon) then
	if(kp) phase = 'SbP'
        if(ks) phase = 'PbS'
      else if (i.eq.imoh) then
	if(kp) phase = 'SmP'
        if(ks) phase = 'PmS'
      endif

      I2=I
      CALL REFLEX(I2,K)

7000  CONTINUE

      conv = .false.

7500  continue

c
c
c     Finally we have to do some interpolations
c
c

8000  continue

      nphas  = mphas + 1
      nphas1 = 0

      do 8800 i=1,np

      indx(i) = 0
      tti(i)  = 0.

      dtdh1 = 0.d0
      dtdh2 = 0.d0

      dpdh1 = 0.d0
      dpdh2 = 0.d0

      dpdd1 = 0.d0
      dpdd2 = 0.d0

      dh1   = 0.d0
      dh2   = 0.d0

      dd1   = 0.d0
      dd2   = 0.d0

      ni    = 0

      do 8700 j=1,nzo

      do 8500 k=1,3

      tm(j,k) = 0.d0
      pm(j,k) = 0.d0

      if (ion(i,j,k).eq.0) go to 8455

      tm(j,k) = ttp(i,j,k)
      pm(j,k) = ppp(i,j,k)

      if (fixho) go to 8450

      j1 = j-1

      if(k.eq.2) then
        if(j.eq.1) then
	  dtdh1 = tm(j,k)
	  dpdh1 = pm(j,k)
	  dh1   = zo(j) 
        else if(j.eq.2 .and. ion(i,j1,k).eq.0) then
	  dtdh1 = tm(j,k)
	  dpdh1 = pm(j,k)
	  dh1   = zo(j) 
        else if(j.eq.3 .and. ion(i,j,k).ne.0) then
	  dtdh2 = tm(j,k)
	  dpdh2 = pm(j,k)
	  dh2   = zo(j)
        else if(j.eq.3 .and. ion(i,j,k).eq.0) then
	  dtdh2 = tm(j1,k)
	  dpdh2 = pm(j1,k)
	  dh2   = zo(j1)
        endif
      endif

8450  k1 = k-1

      if(j.eq.jh1) then
        if(k.eq.1) then
	  dpdd1 = pm(j,k)
	  dd1   = del(k)
        else if(k.eq.2 .and. ion(i,j,k1).eq.0) then
	  dpdd1 = pm(j,k)
	  dd1   = del(k)
        else if(k.eq.3 .and. ion(i,j,k).ne.0) then
	  dpdd2 = pm(j,k)
	  dd2   = del(k)
        else if(k.eq.3 .and. ion(i,j,k).eq.0) then
	  dpdd2 = pm(j,k1)
	  dd2   = del(k1)
        endif
      endif

      if(j.eq.jh1 .and. k.eq.2) then

	ni = 1
        nphas1= nphas1 + 1
        nphas = nphas  - 1

        tti(nphas1) = SNGL(tm(j,k))
        dtdd(nphas) = SNGL(pm(j,k))
	phcd(nphas) = phlist(i)
	if(trimle(phcd(nphas)).ge.3) then
	   tti(nphas1) = tti(nphas1) + 0.00001
        endif

      endif

8455  if(j.eq.jh2 .and. k.eq.3 .and. ni.ne.0) then
	  
	dtdh(nphas) = 0.
	dpdh(nphas) = 0.d0
	dddp(nphas) = 0.d0

	if(.not.fixho) then
           dh3 = dh2-dh1
	   if(dh3.ne.0.d0) then
	      dtdh(nphas) = SNGL((dtdh2-dtdh1) / dh3)
	      dpdh(nphas) = (dpdh2-dpdh1) / dh3
	   endif
	endif

        dd3 = dd2-dd1
	if(dd3.ne.0.d0) then
	   dddp(nphas) = (dpdd2-dpdd1) / dd3
	endif

      endif

8500  continue

8700  continue
      if(typctl.gt.10 .and. nphas1.gt.0) then
         print *,i,nphas,nphas1,phcd(nphas),tti(nphas1),dtdd(nphas),
     *        dtdh(nphas),dpdh(nphas),dddp(nphas)
      endif
8800  continue

      call indexx(nphas1,tti,indx)

      do 8900 i=1,nphas1

      j       = indx(i)
      j2      = mphas + 1 - j

      ttc(i)  = tti(j)
      phcd(i) = phcd(j2)
      dtdd(i) = dtdd(j2)
      dtdh(i) = dtdh(j2)
      dddp(i) = dddp(j2)
      dpdh(i) = dpdh(j2)

c     if(typctl.gt.8) then
c        print *,i,j,j2,dis,phcd(i),ttc(i),dtdd(i),
c    *          dtdh(i),dpdh(i),dddp(i)
c     endif

8900  continue

      nphas2 = nphas1

      zoold =dble(ho)

9000  RETURN
      END
C
C
      SUBROUTINE  REFLEX(II,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ml=61,np=40)

      DIMENSION FA(2,ml),VHQ(2),PA1(2),V(2,ml),del(3),G(2,ml),
     *          V2(2,ml),RR(2),TT(2),pa(2),ttp(np,3,3),ppp(np,3,3),
     *          ion(np,3,3)

      CHARACTER phase*8

      integer phnum,ii,k

      logical conv

      COMMON    /REF/ rr,tt,pa,fa,ttp,ppp,ion,IB,IQQ,IQL,jh1,
     *                PIM,AA,VHQ,PA1,PI,V,G,V2,del,phase,conv

      iph=phnum(phase)

      L=II-1

      if(conv) then

	 VMAX=dmax1(V(1,IQQ),V(2,IQQ))

	 if(k.eq.1) kc = 2
	 if(k.eq.2) kc = 1

         DO  900  I=1,L

         FA(1,I)=1.D0
	 FA(2,I)=0.D0

         IF(V(k,I) .GT.VMAX)  VMAX=V(k,I)
         IF(V2(k,I).GT.VMAX)  VMAX=V2(k,I)

         IF(I.GE.IQQ)  then
            FA(2,I)=1.D0
            IF(V(kc,I) .GT.VMAX)  VMAX=V(kc,I)
            IF(V2(kc,I).GT.VMAX)  VMAX=V2(kc,I)
	 endif

900      CONTINUE

      else

         VMAX=V(K,IQQ)

         DO  1000  I=1,L
         FA(1,I)=2.D0
         IF(I.LT.IQQ)  FA(1,I)=1.D0
         IF(V(K,I) .GT.VMAX)  VMAX=V(K,I)
         IF(V2(K,I).GT.VMAX)  VMAX=V2(K,I)
1000     CONTINUE

      endif

      IBC=IB*30
      IF(IBC.GT.600) IBC=600

      RR(1)=0.d0
      RR(2)=0.d0

      B=PI/(2.d0*DBLE(IBC-1))

      DO  1500 I=1,IBC
      RVV=DSIN(B*DBLE(I-1))/VMAX
      T=0.D0
      R=0.D0

      DO  1100  KK=1,L
      E=V(K,KK)
      G1=E*RVV
      P=DSQRT(DABS(1.D0-G1*G1))
      O=1.d0/G(K,KK)
      F=V2(K,KK)
      G1=F*RVV
      Q=DSQRT(DABS(1.D0-G1*G1))
      R=R+FA(1,KK)*(P-Q)*O
      DT=FA(1,KK)*DLOG(F*(1.D0+P)/(E*(1.D0+Q)))*O
      T=T+DT
1100  CONTINUE

      if(conv) then

         DO  1110  KK=1,L
         if(FA(2,KK).lt.0.9d0) go to 1110
         E=V(KC,KK)
         G1=E*RVV
         P=DSQRT(DABS(1.D0-G1*G1))
         O=1.d0/G(KC,KK)
         F=V2(KC,KK)
         G1=F*RVV
         Q=DSQRT(DABS(1.D0-G1*G1))
         R=R+FA(2,KK)*(P-Q)*O
         DT=FA(2,KK)*DLOG(F*(1.D0+P)/(E*(1.D0+Q)))*O
         T=T+DT
1110     CONTINUE

       endif

      TT(2)=T

      IF(I.GT.1)  then
         RR(2)=R/RVV
         P=DSQRT(DABS(1.D0/(RVV*RVV*VHQ(K))-1.D0))
         IF(P.EQ.0.D0)  then
	   FI=0.5d0*PI
	 else
           FI=DATAN(1.D0/P)
	 endif
      else
         FI=0.D0
      endif

      IF(II.EQ.IQQ)  FI=pi-FI

      RR(2)=RR(2)/AA
      PA(2)=DSIN(FI)*PA1(K)

      DO 1400 KLR=1,3

      if(iql.ne.jh1 .and. klr.ne.2) go to 1400

      del1 = del(klr)

      IF (RR(2).EQ.del1) THEN

	 ion(iph,iql,klr) = ion(iph,iql,klr)+1

	 if(ion(iph,iql,klr).eq.1) then
	    ttp(iph,iql,klr)=TT(2)
	    ppp(iph,iql,klr)=PA(2)
	 else
	    if(TT(2).lt.ttp(iph,iql,klr)) then
	       ttp(iph,iql,klr)=TT(2)
	       ppp(iph,iql,klr)=PA(2)
	    endif
	 endif
         GO TO 1400
      ENDIF

      if(i.eq.1) go to 1400

      FCT1=del1-RR(1)
      FCT2=del1-RR(2)

      IF(FCT1*FCT2.LT.0.d0) THEN
         FCT3=FCT1/(RR(2)-RR(1))
         TT1=FCT3*(TT(2)-TT(1))+TT(1)
         PA2=FCT3*(PA(2)-PA(1))+PA(1)

	 ion(iph,iql,klr) = ion(iph,iql,klr)+1

	 if(ion(iph,iql,klr).eq.1) then
	    ttp(iph,iql,klr)=TT1
	    ppp(iph,iql,klr)=PA2
	 else
	    if(TT1.lt.ttp(iph,iql,klr)) then
	       ttp(iph,iql,klr)=TT1
	       ppp(iph,iql,klr)=PA2
	    endif
	 endif
      ENDIF

1400  CONTINUE

      tt(1) = tt(2)
      rr(1) = rr(2)
      pa(1) = pa(2)
 
1500  CONTINUE
C
      RETURN
      END
c
      function phnum(phase)
      parameter (np=40)
      CHARACTER phase*8,phlist(np)*8
      common  /phasel/ phlist
      integer phnum
      phnum = 999
      do 5 i = 1,np
      if(phase.eq.phlist(i)) then
	phnum = i
	return
      endif
5     continue
      return
      end
