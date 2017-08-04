C
      FUNCTION CRUST(PA,PHTYP,FA,TYPCTL)
c
c     CRUST is a funtion to calculate travel-time differences 
c     for phases passing a standard crust and the CRUST 5.1
c
C     JOHANNES SCHWEITZER
C
c     March 1999    
c
c     input:
c 
c             pa      = ray parameter in [s/deg]
c
c             phtyp   = gives the phase type (P or S)
c
c             fa      = factor how often this phase passes the crust
c                       [i.e. : at station fa = 1, at rfelctions point 
c                               (e.g.PP) fa = 2]
c
c             typctl  = verbosity level
c
c	 in common MODEL(G) parameter of standard model:
c 
c             imo     <  1 not used here
c                     =  1 not used here
c                     =  -2 CRUST 5.1 only used for travel-time 
c                           corrections at position elat,elon 
c                     =  3 not used here
c                     =  4 as (-2) CRUST 5.1 used for travel-time 
c                          corrections
c
c             elat    = latitude  to get the right CRUST 5.1 model 
c                       position for the station
c             elon    = longitude to get the right CRUST 5.1 model 
c                       position for the station
c
c             iread   =  0 data of CRUST 5.1 model must be read in
c                     =  1 data of CRUST 5.1 model are already read
c
c             mtyp    =  CRUST 5.1 model type
c            
c             v0g(1,i) =  P velocity in layer i
c             v0g(2,i) =  S velocity in layer i
c
c             hg(i)    =  depth of layer i
c
c             elevg    =  topograhic elevation at modelled point 
c
c             jmodg    =  number of layers
c
c             azog(i)  =  Conrad/Moho indicator (not used here)
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c     output:
c
c             crust   = travel-time correction for requested phase
c                       with respect to a standard model
c

      INTEGER   typctl

      PARAMETER (ml=61)

      DIMENSION Z(ml),H(ml),V0(2,ml),del(2),time(2),
     *          V(ml),G(ml),V2(ml),zg(ml),v0g(2,ml)

      CHARACTER azog(ml)*4,mtyp*3,phtyp*1

      COMMON    /MODEL/  v0,h,elev,elat,elon,zmax,elat2,elon2,elev2,
     *                   jmod,iread,imo,mtyp

      COMMON    /MODELG/ v0g,zg,elevg,jmodg,azog,zmaxg

      if (iabs(imo).ne.2 .and. imo.ne.4) go to 9998

      PI=4.d0*DATAN(1.d0)
      PIM=PI/180.d0
      re = 6371.d0
      AA=PIM*re

      mtyp = 'MLM'
      itrue = 1
      ierr = 0
      inum = 1
      call get_mod_mlm(itrue,inum,typctl,ierr)

      if(ierr.ne.0) then
	 go to 9000
      endif

c
c     reset onset table
c
      if(phtyp.eq.'P') iph = 1
      if(phtyp.eq.'S') iph = 2
C

      zma   = re - zmax - elev
      zmag  = re - zmaxg - elevg

      zmaxi = dmin1(zmag,zma)
      h(jmod)   = re - zmaxi + elev
      zg(jmodg) = re - zmaxi + elevg

      do 7000 k = 1,2

      if(k.eq.1) jmodk = jmod
      if(k.eq.2) jmodk = jmodg

      del(k)  = 0.d0
      time(k) = 0.d0

      vmax = 0.d0

      DO 500 I=1,jmodk

      if(k.eq.1) then
         V(I)=V0(iph,I)
	 Z(I)=H(I)
      endif 
      if(k.eq.2) then
         V(I)=V0g(iph,I)
	 Z(I)=ZG(I)
      endif

      if(vmax.lt.v(i)) vmax = v(i)

500   continue

      PA1=(re-Z(jmodk))*PIM
      RVV = pa / PA1
      if(rvv*vmax.gt.0.985d0) go to 9998

      m = jmodk - 1

      DO 800 I=1,M

      I2=I+1
      V2(I)=V(I2)

      IF(V2(I).EQ.V(I)) THEN
         V2(I)=1.0001d0*V(I)
         V(I2)=V2(I)
      ENDIF

      zdiff=Z(I2)-Z(I)
      IF(zdiff.EQ.0.d0)  then
          zdiff=0.0001d0
          z(i2)= z(i2) + zdiff
      endif

      G(I)=(V2(I)-V(I))/zdiff

800   continue

      T=0.D0
      R=0.D0

      DO  1100  KK=1,M
      E=V(KK)
      G1=E*RVV
      P=DSQRT(DABS(1.D0-G1*G1))
      O=1.d0/G(KK)
      F=V2(KK)
      G1=F*RVV
      Q=DSQRT(DABS(1.D0-G1*G1))
      R=R+FA*(P-Q)*O
      DT=FA*DLOG(F*(1.D0+P)/(E*(1.D0+Q)))*O
      T=T+DT
1100  CONTINUE

      del(k)  = r/(aa*rvv)
      time(k) = t

c     print *,'crust: ',k,del(k),time(k)

7000  continue

      ddel = del(1) - del(2)
      t1 = ddel * pa
      crust = time(1) - time(2) - t1

      if(typctl.ge.8) then
	 print *,'crust-correction: ',pa,phtyp,ddel,t1,crust
      endif
      go to 9999

9000  continue
      print *,'crust: no CRUST 5.1 correction possible'
      ierr = 0
      
9998  crust = 0.d0
9999  RETURN
      END
c
c     funtion crustc
c
c     driver to calculate the travel-time effect 
c     for reflections at the Earth's surface for 
c     different crustal velocity structures and 
c     topography and calculates station corrections
c     for all kind of phases.
c
c     we can correct for the following principle phases:
c
c                 pP, sP, sS, pS 
c                 PP, SS, PS, SP
c                 P'P', S'S', P'S', S'P'
c 
c
      function crustc(phase_t0,rayp0,ind,typctl)
c
c     input:
c              phase_t0 phase type (either P or S)
c
c              rayp0    ray parameter of phase in [s/deg]
c
c              ind      switch for reflection type
c
c                       = 1 station correction for all phases
c                       = 2 surphase reflection pP, sS, PP, SS
c                       = 3 converted surphase reflection sP, pS
c
c     calls FUNCTION CRUST
c

      implicit double precision (a-h,o-z)

      integer ind, typctl

      character*1 phase_t0, phase_r

      double precision rayp0, crustc, crust

      crustc = 0.d0

      rayp = rayp0

      phase_r = phase_t0

      if(ind.eq.1) then

        fmult   = 1.0d0
        crustc  = crust(rayp,phase_r,fmult,typctl)
        go to 9000

      else if(ind.eq.2) then

        fmult   = 2.d0
        crustc  = crust(rayp,phase_r,fmult,typctl)
        go to 9000

      else if (ind.eq.3) then

        fmult   = 1.d0
        phase_r = 'P'
        crustc1  = crust(rayp,phase_r,fmult,typctl)
        phase_r = 'S'
        crustc2  = crust(rayp,phase_r,fmult,typctl)

        crustc = crustc1 + crustc2
        go to 9000

      endif

9000  return
      end
