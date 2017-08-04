c
c     
c
      program HYPOSAT_4_4b

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character  version*20
      parameter (version='HYPOSAT Version 4.4b')

c
c     last update: 23 July 2003
c
c
c     This programm locates earthquakes with observed travel times,
c     backazimuth, and slowness values.
c
c     Different phases observed at one station are used to determine
c     travel time differences. These differences are also used to
c     invert for the hypocenter.
c
c     A preliminary epicenter will be located with the backazimuth 
c     values, or with other available informations.
c
c     A source time will be estimated with a Wadati-approach
c     from the travel time difference(s) assuming a constant 
c     v(p)/v(s)=sqrt(3.)
c
c     The final location is done with a least squares fit for
c     all four source parameters using the tau-spline-type tables
c     (i.e., IASP91, AK135, ...), or a local/regional model of
c     horizontal layers.
c
c     All available informations are used including standard 
c     deviations for observed data and 'a priori' informations
c     for the model parameter.
c
c     All calculations are done for the Earth as a sphere and
c     the travel times are corrected  for the ellipticity of the
c     Earth.
c
c--------------------------------------------------------------------
c
c                           Program History
c
c     First version progamed summer 1996 by     
c
c                                 Johannes Schweitzer, 
c                                 Institute of Geophysics
c                                 Ruhr-University Bochum
c                                 D-44780 BOCHUM
c                                 Germany
c
c
c     all changes and extension in the whole program package after 
c     1. July 1997:
c
c
c                                 Johannes Schweitzer
c    
c                                 NORSAR
c                                 P.O.Box 53
c                                 N-2027 KJELLER
c                                 Norway
c
c                                 johannes.schweitzer@norsar.no
c
c
c       major improvements, corrections or changes since 1997:
c
c                               12-15 February 1997
c                               czo = b  means: starting with 'f' and 
c                               ending with 'd'
c                               Flow of calculations for oscillating 
c                               solutions changed, especially for depth
c                               determinations.
c                               Correcting to-time for Wadati formula;
c                               included some maximum values for to.
c                               Handling of dtm changed.
c                               P1 and S1 included as phase name for
c                               (unknown) first P- and S-onsets. The
c                               program choose the right name depending
c                               on travel-time table and distance.
c
c                               13 March 1997
c                               Usage of PKiKP instead of PKPdf 
c                               whenever PKPdf does not exist (for 
c                               observations around the triplication). 
c                               Similar changes for P/Pdif and S/Sdif.
c                               Startsolution always at the closest 
c                               station, if no backazimuth information
c                               is available.
c
c                               23 April 1997
c                               some changes in hyposat_gmi to print 
c                               out the resolution, covariance, 
c                               correlation, and the information-density 
c                               matrix.
c
c               version 2.2     May 8, 1997
c                               Station corrections included with file 
c                               station.cor .
c                               Small bug to calculate dpdh removed.
c
c               version 3.0     June 2, 1997
c                               Local velocity model included.
c                               Checking if oscillating solution
c                               is running over 4 solutions.
c
c                               The new parameter SETCHECK included to 
c                               give a limit for stopping the iteration 
c                               process.
c
c               version 3.0a    June 5, 1997
c                               Changes in handling of TTLOC to reduce 
c                               runtime (indph as new parameter).
c                               June 11, 1997
c                               Further changes in handling of TTLOC to 
c                               reduce runtime.
c
c               version 3.0b    July 10, 1997 at NORSAR
c                               Switch diffflag included. If set to 
c                               .ne.0 , no travel-time differences will 
c                               be used.
c
c                               July 14, 1997
c                               Determined phase name only printed if 
c                               different from input phase name.
c
c                               July 18/23, 1997
c                               Handling of dtm changed again. Now 
c                               variable rms of resiudals of former 
c                               solution (only if last change is smaller 
c                               than 10*setcheck).
c                               Removing small bug in order of output 
c                               listing. Smaller changes in output-file.
c
c               version 3.0c    July 25-30, 1997
c
c                               Several mistyping removed. Phase-naming
c                               changed and logical errors removed.
c                               Smaller adjustments to relocate 
c                               REB-events.
c
c                               ERROR due to -180 deg / +180 deg border
c                               for calculating a mean starting solution 
c                               removed!
c
c                               Iteration also stops if the last change
c                               of hypocenter is smaller than 1/1000 of
c                               distance to the nearest station (as long
c                               as SETCHECK has the default value of 
c                               1.).
c                               All these ditances are measured in [km].
c                               
c               version 3.1     August 2-8, 1997
c                               SDLATG corrected for events on the 
c                               souther hemisphere.
c
c                               Handling of multiple-core phases 
c                               changed. Backazimuth handling of LR, LQ 
c                               fixed.
c
c                               Calculating standard deviations fixed, 
c                               if only two observations are available 
c                               for Wadati-curve.
c
c                               Because no ellipticity corrections are 
c                               available for events deeper than 700 km,
c                               we accept a smal error for events 
c                               between 700 km and 800 km.
c
c               version 3.2     September, 21 - October, 16 1997
c                               Plotting removed and the parameter file 
c                               decoupled from data file. Possibility to
c                               give a starting epicenter and its 
c                               standard errors included.
c                               Usage of travel-time differences 
c                               changed.
c                                
c               version 3.2a    October, 29 - November 28, 1997
c                               Startsolution for station nets 
c                               corrected. Cases with low information 
c                               (number of data very small) changed.
c
c                               Elevation correction corrected for 
c                               unknown phase-type. Inversion for 
c                               underestimated events changed. Handling 
c                               of unstable depth estimates changed. 
c                               Reducing the output of 'bad'-solutions.
c                               In parameter file HEIGHT changed to 
c                               ELEVATION.
c
c                               Missing backazimuth values corrected for
c                               LR phases in the REBs.
c
c                               Partial derivatives to define the source 
c                               depth must be at least 0.00001 [s/km] or
c                               0.00001 [(s/deg)/km].
c
c                               Comparing 'Wadati-source time' with 
c                               'final-source time' to get a better 
c                               depth in the case of a large mean travel-
c                               time residual.
c
c                               December 16, - December 18, 1997
c                               Calculating the partial derivatives for
c                               'tau-spline'-models numerically!
c
c                               Changed to new 'Flinn-Engdahl' regions.
c
c                               January 15, 1998
c                               DTKM adjusted for IASP91 crust.
c              
c               version 3.2b    June 17-22, 1998
c                               Smaller changes for one-phase 
c                               observation at only one array (station).
c
c               version 3.2c    September/October 1998
c                               Changes for local models (model can be 
c                               given with file name although no phases 
c                               are wanted), removing of smaller bugs in 
c                               calculating the starting source time.
c
c               version 3.2d    December 1998
c                               New parameters in 'hyposat-parameter'
c                               and surface wave 'Rg' included using a 
c                               constant group velocity.
c                               CSS 3.0 .site file format included.
c                               COMMON blocks for inversion calls 
c                               included.
c
c               version 4.0     March 1999
c                               CRUST 5.1 after Mooney, Laske, and 
c                               Masters (JGR, Jan. 1998) included.
c                               input-parameter 'touse' included.
c
c                               May 1999
c                               CSS 3.0 .site file format corrected.
c                               Correction for different structures at
c                               reflection points of the Earth's 
c                               surface.
c
c               version 4.1     September - January 2001
c                               Group velocities for LQ, LR, and Lg 
c                               included. Usage of (Lg-Pn) times for 
c                               source time estimates included. 
c
c                               Epicenter error ellipse inluded. Start 
c                               values for source time and its standard 
c                               deviations included. Travel-time table 
c                               for local models expanded (pS, sP, PmS, 
c                               SmP,...).
c
c                               Whenever needed, we changed from the
c                               standard Earth radius to the actual
c                               Earth radius as calculated with the
c                               standard ellipsoid parameter (see 
c                               FUNCTION RADLOC).
c
c                               HYPOSAT_CROSS corrected: error bars and
c                               non-crossing directions removed.
c                               RMS calculation also for slowness and
c                               backazimuth residuals.
c
c                               Calculation of magnitudes (Ms or mb)
c                               Attenuation model Gutenberg-Richter or
c                               Veith-Clawson for mb and for Ms the 
c                               IASPEI 1967 formula .
c
c               version 4.2     June 2001 (distributed)
c
c                               Possibility to use two different global 
c                               models during one inversion usage 
c                               included.
c
c                               Handling of oscillating solutions 
c                               changed and stabilized.
c
c               version 4.3     August/September 2001 (distributed)
c
c                               Calculation of weighted misfit for 
c                               all data  and of the azimuthal gap of 
c                               defining observations included.
c
c                               Ms calculation with Rezapour/Pearce
c                               (BSSA 88, 43-61) formula added.
c
c               version 4.3a    November 2001
c
c                               Calculation of azimuthal gap corrected,
c                               exponential output for large error
c                               ellipsis.
c
c               version 4.3b    January - March 2002
c
c                               Plane wave fit included to get better
c                               start solutions in the case of a distant
c                               event with respect to the network
c                               aperture.
c
c                               New version of libtau software (PKPdif!)
c                               included.
c
c               version 4.3c    March 2002
c               
c                               Re-compile and re-check of code for
c                               errors and inconsistencies.
c
c               version 4.3d    July 2002
c
c                               New input format for starting source 
c                               time. Usage of different global models
c                               extended to four!
c
c               version 4.4     October 2002 (distributed)
c
c                               Weighting of outlayers changed ( see
c                               DATMAX0 variable).
c                               Phase names changed as recommanded by
c                               IASPEI working group on phase names in
c                               location routine and libtau software,
c                               including ellipticity corrections.
c
c                               Some other smaller changes and 'bugs' 
c                               removed.
c
c               version 4.4a    March - May 2003 (distributed)
c                               
c                               New parameter to use azimuth only for
c                               the inital solution.
c                               Single array, single phase solutions
c                               improved and corrected.
c                               Small error for usage of different 
c                               models corrected.
c
c               version 4.4b    July  2003 (distribution exchanged)
c
c                               Bug for negative station elevation 
c                               correction fixed.
c
c
c     calls :   delazd, depi, dlsq, fetoh, fhtoe, findosci, get_station, 
c               hyposat_cross, ellip, ellip2, hyposat_gmi, indexx, 
c               tauget_mod, testphase, hyposat_geo, ttloc, get_mod_mlm,
c               ellcal, plane
c
c     functions: alpha1, alpha2, convlat, phase_type, trimle, crustc, 
c                ddmax, dirdel, dmean, q2, radloc, tauget_ray
c
c
c     PARAMETER settings for: mstat, mread, mvar, mosci0
c
c
c
c****6789012345678901234567890123456789012345678901234567890123456789012
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c     set variable type for functions
c
      integer          trimle
      double precision alpha1, alpha2, convlat, crustc, dirdel, q2, 
     +                 radloc, dmean, ddmax, getchi
      character        phase_type*1, file_check*80

c
c     mstat = maximum number of stations
c
      parameter  (mstat = 700)

      dimension stala(mstat),stalo(mstat),stael(mstat),del(mstat),
     +          azie(mstat),baz(mstat),stalae(mstat),stavp(mstat),
     +          stavs(mstat),istaph(mstat)
      character*5 sta(mstat),stat,stato,stationfile*80,statcorfile*80,
     +          stat1,outputfile*80,inputfile*80,magfile*80
c
c     mread = maximum number of phases
c
      parameter  (mread = 800, mr2=mread/2, mread2=mread*2)
      parameter  (mloc  = (mread/2 + 3)*mread)

      character phase(mread)*8,phid*8,used(mread)*5,text(mread)*125,
     +          phid2*8,text2(30)*125,phid1*8,phase_t*1,phidr0*8,
     +          string*130,touse(mread)*7,touse0*7,phidr*8,phipl*8
      dimension azi(mread),tt(mread),p(mread),azis(mread),tts(mread),
     +          ps(mread),period(mread),amplit(mread)
      real      arr(mread),epiaz(mread),dmazi,dmazi0,d1azi,d2azi
      integer   iev(mread),indx(mread),indx2(30)

      dimension elo(mloc),ela(mloc),elos(mloc),elas(mloc)
c
c     mvar = maximum varibales for lsq-fit or gmi
c
      parameter  (mvar = 4 )

      dimension a(mread2,mvar),var(mvar),dat(mread2),r(mvar),
     +          res(mread2),dats(mread2),rs(mvar),ax1(2),
     +          ax2(2),idtu(mread)

      logical   iellip

      common  /gmi/ a,var,dat,r,res,dats,rs,ax1,ax2,iellip

      integer   is(2)


c
c     mosci0 is the maximum number of solutions checked for
c     oscillating solutions.
c
      parameter (mosci0=10)
      dimension dzoos(mosci0),dtos(mosci0),dlaos(mosci0),
     +          dlo1os(mosci0),dloos(mosci0),dlo2os(mosci0),
     +          rzos(mosci0),rtos(mosci0),rlaos(mosci0),rloos(mosci0)

c
c     variables for travel-time calculations
c
      include 'ttimes.h'

      dimension dpdh(mphas),dddp(mphas)

      character phcd1(mphas)*8,phcd2(mphas)*8,art*16,modnam*20,
     +          modnam2*20,modn*20,modnam3*20,modnam4*20

      real*4 rzo,rdel,ecor,fla1,fla2,fla3,razi,rbaz,pa,rzo1,rzo2,
     +       dtdd2(mphas),rmax,rzoe,ttc1(mphas),dtdh1(mphas),
     +       rdel1,rdel2,drdel,dtdd1(mphas), rmcorr

c
c     common block with informations about the local/regional model
c     (used also for CRUST 5.1)
c
c     maxla = maximum number of layers for a local/regional velocity
c             model
c

      parameter (maxla = 61)
      DIMENSION V0(2,maxla),z(maxla),v0g(2,maxla),zg(maxla)
      character azo(maxla)*4,azog(maxla)*4,filloc*80,mtyp*3,mpath*80,
     +          imod*1
      integer   jmod,jmodg,iread,imo

      COMMON /MODEL/  v0,z,elev,elatc,elonc,zmax,elat2,elon2,elev2,
     +                jmod,iread,imo,mtyp

      COMMON    /MODELC/ filloc,azo,mpath
      COMMON    /MODELG/ v0g,zg,elevg,jmodg,azog,zmaxg
c
c     other variables
c
      integer yy,mon,dd,hh,idoy,ierr,typctl,mi,idum, isreg, regnum,
     +        typctlm
      character mm*4,name*48
      double precision lat,lon,kmdel
      real*4  elevs, sec, rlat, rlon

      character title*80, czo*1, region*80, czo1*1, magtypp*3, 
     +          magtyps*6

c
c     idtmax = number of different travel-time-difference definitions
c              used to calculate a start source time 
c
      parameter (idtmax = 4)
      dimension dt(idtmax,mr2),idtp(idtmax,mr2),idts(idtmax,mr2),
     +          idt(idtmax),to(idtmax),tos(idtmax),vpvs(idtmax),
     +          vpvss(idtmax),datho(mread),datla(mread),datlo(mread),
     +          ttt(mread)

      logical   iaspflag, zoflag , vlflag  , stcorfl, iloc, surf,
     +          diffflag, dtmflag, epistart, single , output, modout, 
     +          last, magflag, mod2flag, lastfix, direct, plflag,
     +          conr, tauget_ray, rayok, rayokf, mod3flag, mod4flag, 
     +          aziini

c
c     some constants and starting values
c

      pi      = 4.d0*datan(1.d0)
      deg2rad = pi / 180.d0
      rad2deg = 180.d0 / pi
      rearth  = 6371.d0
      grad1   = 1.d0/(deg2rad*rearth)
      eps     = q2(1.d0-1.d0/298.257d0)

      dtp0    = 600.d0
      dts0    = dtp0*2.d0
      dtm2    = dts0
      dtmin   = 9999.d0
      dtdt    = 10.d0

      ttray   = 0.d0
      ddel    = 0.d0

      dchang0 = 1.d0

      check    = 9999.d0
      disper   = 0.001d0
      rmso     = 9999.d0
      rmsold   = 9999.d0
      datmax0  = 9999.d0
      astmean0 = 0.d0
      nrms1    = 0

      iteras  = 0
      miteras = 0
      iterz   = 0
      iteraz  = 0
      ibad0   = 0
      idelw   = 0
      in0sw   = 0
      imodout = 0 

      insar   = 0

      nextiter1 = 0
      imaxiter  = 0
      ilastiter = 0
      mosci     = 4

      iaspflag = .false.
      vlflag   = .false.
      zoflag   = .false.
      stcorfl  = .false.
      diffflag = .true.
      iloc     = .false.
      epistart = .false.
      single   = .false.
      output   = .false.
      iellip   = .true.
      modout   = .false.
      magflag  = .false.
      mod2flag = .false.
      mod3flag = .false.
      mod4flag = .false.
      lastfix  = .true.
      plflag   = .true.
      dtmflag  = .false.
      rayok    = .false.
      rayokf   = .false.
      aziini   = .false.

c
c     set default values
c

      modnam = 'iasp91a'
      modnam2 = modnam
      modnam3 = modnam
      modnam4 = modnam
      rmax = 0.
      filloc = ''
      imo     = 1
      stationfile = 'stations.dat'
      outputfile  = 'hyposat-out'
      inputfile   = 'hyposat-in'
      mpath       = './'
      statcorfile = ''
      vpl = 99.d0
      vsl = 0.d0
      zo1 =  0.d0
      sdzo1 = 50.d0
      czo = 'F'
      maxiter = 80
      confl = 68.26895d0
      iellipi = 1
      dazim = 45.d0
      dpam  = 5.d0
      typctl = 4
      islow = 1
      setcheck1 = 0.d0
      setcheck  = 1.d0
      indph0 = 3333
      inddiff = 1
      epilat0 = -999.d0
      epilon0 = -999.d0
      sdlatg  = 10.d0
      sdlat   = 0.d0
      sdlon   = 10.d0
      tome0   = 0.d0
      stome0  = 120.d0

      vrg0    = 2.5d0
      vrg     = 0.d0

      vlg0    = 3.5d0
      vlg     = 0.d0

      vlr0    = 3.95d0
      vlr     = 0.d0

      vlq0    = 4.4d0
      vlq     = 0.d0

      magtypp = ''
      magtyps = ''

c
c     read in parameters from file
c

      open (unit=9,file='hyposat-parameter',status='old')

      do 1 jin = 1,100

      read (9,'(a)',end=2) string

      if(string(1:14).eq.'GLOBAL MODEL  ') then
          read (string(38:),'(a)') modnam
          modnam2 = modnam
          modnam3 = modnam
          modnam4 = modnam
          go to 1
      endif

      if(string(1:14).eq.'GLOBAL MODEL 1') then
          read (string(38:),'(a)') modnam
          modnam2 = modnam
          modnam3 = modnam
          modnam4 = modnam
          go to 1
      endif

      if(string(1:14).eq.'GLOBAL MODEL 2') then
          read (string(38:),'(a)') modnam2
          if(modnam2 .ne. '_' .and. modnam2.ne.' ') then
            mod2flag = .true.
          endif
          go to 1
      endif

      if(string(1:14).eq.'GLOBAL MODEL 3') then
          read (string(38:),'(a)') modnam3
          if(modnam3 .ne. '_' .and. modnam3.ne.' ') then
            mod3flag = .true.
          endif
          go to 1
      endif

      if(string(1:14).eq.'GLOBAL MODEL 4') then
          read (string(38:),'(a)') modnam4
          if(modnam4 .ne. '_' .and. modnam4.ne.' ') then
            mod4flag = .true.
          endif
          go to 1
      endif

      if(string(1:23).eq.'LOCAL OR REGIONAL MODEL') then
          read (string(38:),'(a)') filloc
	  go to 1
      endif

      if(string(1:24).eq.'OUTPUT OF REGIONAL MODEL') then

          read (string(38:),*) imodout
	  if(imodout.eq.1) then
	     modout = .true.
	  else
	     modout = .false.
	  endif
	  go to 1

      endif

      if(string(1:11).eq.'CRUST 5.1  ') then
          read (string(38:),*) imo
	  if(imo.le.0) imo=1
	  go to 1
      endif

      if(string(1:14).eq.'CRUST 5.1 PATH') then
          if(string(38:38).ne.'_' .and. string(38:38).ne.' ') 
     +       read (string(38:),'(a)') mpath
	  go to 1
      endif

      if(string(1:12).eq.'STATION FILE') then
          read (string(38:),'(a)') stationfile
	  stationfile = file_check(stationfile)
	  go to 1
      endif

      if(string(1:23).eq.'STATION CORRECTION FILE') then
          read (string(38:),'(a)') statcorfile
	  if (statcorfile(1:trimle(statcorfile)).ne.'' .and.
     +        statcorfile(1:trimle(statcorfile)).ne.'_'   ) then
	      statcorfile = file_check(statcorfile)
	      stcorfl = .true.
	  else
	      stcorfl = .false.
          endif
	  go to 1
      endif

      if(string(1:31).eq.'P-VELOCITY TO CORRECT ELEVATION') then
          read (string(38:),*) vpl
	  go to 1
      endif

      if(string(1:31).eq.'S-VELOCITY TO CORRECT ELEVATION') then
          read (string(38:),*) vsl
	  go to 1
      endif

      if(string(1:24).eq.'PLANE WAVE APPROXIMATION') then
          read (string(38:),*) iplane
	  if(iplane.ne.1) plflag = .false.
	  go to 1
      endif

      if(string(1:20).eq.'STARTING SOURCE TIME') then
c
c     time format: epochal time 
c     or in ASCII format: yyyy-doy:hh.mi.ss.sss
c     or in ASCII format: yyyy-mm-dd:hh.mi.ss.sss
c
	  if(string(42:42).ne.'-') then
             read (string(38:),*) tome0
	  else

	    mm   = '   '
	    idum = 0
	    idoy = 0
	    mon  = 0
	    dd   = 0

	    if(string(45:45).ne.'-') then
              read (string(38:51),'(i4,x,i3,2(x,i2))') yy,idoy,hh,mi
	      read (string(53:),*) sec
	    else
              read (string(38:53),'(i4,4(x,i2))') yy,mon,dd,hh,mi
	      read (string(55:),*) sec
	    endif

	    call fhtoe(tome0,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

	  endif
	  go to 1
      endif

      if(string(1:19).eq.'STARTING TIME ERROR') then
          read (string(38:),*) stome0
	  go to 1
      endif

      if(string(1:21).eq.'STARTING SOURCE DEPTH') then
          read (string(38:),*) zo1
	  go to 1
      endif

      if(string(1:20).eq.'STARTING DEPTH ERROR') then
          read (string(38:),*) sdzo1
	  if(sdzo1.eq.0.d0) sdzo1 = 50.d0
	  go to 1
      endif

      if(string(1:24).eq.'STARTING SOURCE LATITUDE') then
	  abc = -999.0d0
          read (string(38:),*) abc
          if(abc.ge.-90.d0 .and. abc.le.90) epilat0 = abc
	  go to 1
      endif

      if(string(1:23).eq.'STARTING LATITUDE ERROR') then
          read (string(38:),*) sdlatg
	  go to 1
      endif

      if(string(1:25).eq.'STARTING SOURCE LONGITUDE') then
	  abc = -999.0d0
          read (string(38:),*) abc
          if(abc.ge.-180.d0 .and. abc.le.180) epilon0 = abc
	  go to 1
      endif

      if(string(1:24).eq.'STARTING LONGITUDE ERROR') then
          read (string(38:),*) sdlon
	  go to 1
      endif

      if(string(1:10).eq.'DEPTH FLAG') then
          read (string(38:),'(a)') czo
	  go to 1
      endif

      if(string(1:23).eq.'MAXIMUM # OF ITERATIONS') then
          read (string(38:),*) maxiter
	  go to 1
      endif

      if(string(1:24).eq.'# TO SEARCH OSCILLATIONS') then
          read (string(38:),*) mosci
	  if(mosci .gt. mosci0) mosci = mosci0
	  if(mosci .lt. 1) mosci = 1
	  go to 1
      endif

      if(string(1:16).eq.'CONFIDENCE LEVEL') then

          read (string(38:),*) confl
	  if(confl.lt.68.26895d0)  confl = 68.26895d0
	  if(confl.gt.99.99d0) confl = 99.99d0

	  go to 1
      endif

      if(string(1:18).eq.'CONSTRAIN SOLUTION') then

          read (string(38:),*) ilastfix
	  if(ilastfix.ne.1) then
	     lastfix = .false.
	  else
	     lastfix = .true.
	  endif

	  go to 1
      endif

      if(string(1:23).eq.'EPICENTER ERROR ELLIPSE') then
          read (string(38:),*) iellipi

	  if (iellipi.eq.0) then
	    iellip = .false.
	  else if(iellipi.eq.1) then
	    iellip = .true.
	  endif

	  go to 1
      endif

      if(string(1:21).eq.'MAXIMUM AZIMUTH ERROR') then
          read (string(38:),*) dazim
	  go to 1
      endif

      if(string(1:29).eq.'AZIMUTH ONLY INITIAL SOLUTION') then
          read (string(38:),*) iaziini
	  if(iaziini.eq.1) aziini = .true.
	  go to 1
      endif


      if(string(1:22).eq.'MAXIMUM SLOWNESS ERROR') then
          read (string(38:),*) dpam
	  go to 1
      endif

      if(string(1:12).eq.'OUTPUT LEVEL') then
          read (string(38:),*) typctl
          typctlm = typctl
          if (typctl.gt.10) then
             typctl  = (typctl-10*(typctl/10))*2 - 3
	     if (typctl.lt.0)  typctl = 0
	     if (typctl.gt.10) typctl = 10
          endif

	  go to 1
      endif

      if(string(1:16).eq.'SLOWNESS [S/DEG]') then
          read (string(38:),*) islow
	  go to 1
      endif

      if(string(1:17).eq.'LOCATION ACCURACY') then
          read (string(38:),*) setcheck1
	  go to 1
      endif

      if(string(1:27).eq.'PHASE INDEX FOR LOCAL MODEL') then
          read (string(38:),*) indph0
	  if(indph0.gt.3333) then
	    print *,'Wrong input: PHASE INDEX FOR LOCAL MODEL'
            go to 9999
	  endif
	  go to 1
      endif

      if(string(1:34).eq.'FLAG USING TRAVEL-TIME DIFFERENCES') then
          read (string(38:),*) inddiff
	  go to 1
      endif

      if(string(1:15).eq.'INPUT FILE NAME') then
          if(string(38:38).ne.' ' .and. string(38:38) .ne.'_')
     +	                read (string(38:),*) inputfile
	  go to 1
      endif

      if(string(1:16).eq.'OUTPUT FILE NAME') then
          if(string(38:38).ne.' ' .and. string(38:38) .ne.'_')
     +	                read (string(38:),*) outputfile
	  go to 1
      endif

      if(string(1:13).eq.'OUTPUT SWITCH') then
          if(string(38:38).eq.'1') output=.true.
	  go to 1
      endif

      if(string(1:17).eq.'RG GROUP-VELOCITY') then
          read (string(38:),*) vrg
	  if(vrg.le.0.d0) vrg = vrg0
	  go to 1
      endif

      if(string(1:17).eq.'LG GROUP-VELOCITY') then
          read (string(38:),*) vlg
	  if(vlg.le.0.d0) vlg = vlg0
	  go to 1
      endif

      if(string(1:17).eq.'LR GROUP-VELOCITY') then
          read (string(38:),*) vlr
	  if(vlr.le.0.d0) vlr = vlr0
	  go to 1
      endif

      if(string(1:17).eq.'LQ GROUP-VELOCITY') then
          read (string(38:),*) vlq
	  if(vlq.le.0.d0) vlq = vlq0
	  go to 1
      endif

      if (string(1:21).eq.'MAGNITUDE CALCULATION') then
	  read (string(38:),*) imag
	  if (imag.eq.1) then
	    magflag = .true.
	  else
	    magflag = .false.
	  endif
	  go to 1
      endif

      if (string(1:19).eq.'P-ATTENUATION MODEL') then
	  read (string(38:),*) magtypp
	  go to 1
      endif

      if (string(1:19).eq.'S-ATTENUATION MODEL') then
	  read (string(38:),*) magtyps
	  go to 1
      endif

1     continue

2     close(9)

      fchi = dsqrt(getchi(1.0d0,confl))

      if(typctl.gt.4) then
         print *,'modnam   : ',modnam
         if(mod2flag) print *,'modnam 2: ',modnam2
         if(mod3flag) print *,'modnam 3: ',modnam3
         if(mod4flag) print *,'modnam 4: ',modnam4
         print *,'filloc = ',filloc
         print *,'imo = ',imo
	 print *,'CRUST 5.1 path',mpath
         print *,'stationfile = ',stationfile
         print *,'statcorfile = ',statcorfile 
         print *,'inputfile   = ',inputfile 
         print *,'outputfile  = ',outputfile 
         print *,'output switch ',output 
         print *,'vpl = ',vpl
         print *,'vsl = ',vsl
         print *,'vrg = ',vrg
         print *,'vlg = ',vlg
         print *,'vlq = ',vlq
         print *,'vlr = ',vlr
         print *,'zo1   = ',zo1
         print *,'sdzo1 = ',sdzo1
         print *,'czo   = ',czo
         print *,'epilat0 = ',epilat0
         print *,'sdlatg  = ',sdlatg
         print *,'epilon0 = ',epilon0
         print *,'sdlon   = ',sdlon
         print *,'maxiter = ',maxiter
         print *,'lastfix = ',lastfix
         print *,'confl   = ',confl  
         print *,'dazim = ',dazim 
         print *,'dpam  = ',dpam
         print *,'typctl = ',typctl
         print *,'islow = ',islow
         print *,'setcheck1 = ',setcheck1
         print *,'indph0 = ',indph0 
         print *,'inddiff = ',inddiff
	 print *,'Magnitude flags = ', magflag,magtypp,magtyps
	 print *,'Plane wave = ', plflag
      endif

      if(modnam(1:6).eq.'iasp91') iaspflag = .true.

      if(imo.gt.1) then

	 iread = 0
c
c     Default global crustal model is set to AK135
c
         mtyp = 'E6 '
         
	 if(modnam(1:2).eq.'jb')      mtyp = 'E1 '
	 if(modnam(1:4).eq.'prem')    mtyp = 'E2 '
	 if(modnam(1:7).eq.'iasp91 ') mtyp = 'E3 '
	 if(modnam(1:7).eq.'iasp91a') mtyp = 'E4 '
	 if(modnam(1:3).eq.'sp6')     mtyp = 'E5 '
	 if(modnam(1:5).eq.'ak135')   mtyp = 'E6 '

	 elatc = 0.d0
	 elonc = 0.d0

	 itrue = 1
	 inum  = 1
	 call get_mod_mlm(itrue,inum,typctl,ierr)
	 if (ierr.ne.0) then
	    rmax = 0.
	    write(*,'('' Cannot use CRUST 5.1 model '')')
	    imo = 1
	    ierr = 0
	 endif

	 jmodg = jmod
	 elevg = elev
	 zmaxg = zmax
	 do 21 i = 1,jmodg

	   v0g(1,i) = v0(1,i)
	   v0g(2,i) = v0(2,i)
	   zg(i)    = z(i)
	   azog(i)  = azo(i)

21       continue

	 if (imo.ge.3) then
	    mtyp = 'MLM'
	    rmax = 6.
	    iloc = .true.
	    filloc = 'CRUST 5.1'
         endif

      endif

      if(filloc(1:trimle(filloc)).ne.'' .and.
     +   indph0.ge.0                    .and.
     +   filloc(1:1).ne.'_'             .and.
     +   imo.le.2                            ) then
	 
	 filloc = file_check(filloc)
	 ierr   = 0
	 rzo    = 0.
	 rdel   = 0.
	 czo1   = ''
	 indph  = 0
	 elatc  = 0.d0
	 elonc  = 0.d0

	 call ttloc(rzo,rdel,czo1,nphas,ttc,dtdd,dtdh,dpdh,dddp,phcd,
     +              rmax,typctl,ierr,indph)

	 if(ierr.ne.0) then 
	    rmax = 0.
	    write(*,'('' Can only use global model: '',a)') modnam
	    ierr = 0
	 else

	    iloc = .true.
c
c           use local model as standard crustal model
c
	    jmodg = jmod
	    elevg = elev

	    do 22 i = 1,jmodg
	      v0g(1,i) = v0(1,i)
	      v0g(2,i) = v0(2,i)
	      zg(i)    = z(i)
	      azog(i)  = azo(i)
22          continue

	 endif
      endif

      dazim1 = dazim
      dpam1  = dpam

      if(vpl.ne.99.d0) then

        if(vpl.eq.0.d0) vpl=5.8d0
        if(vsl.eq.0.d0) vsl=vpl/dsqrt(3.d0)
	vlflag = .true.

	if (stcorfl) open(12,file=statcorfile)

      endif

      if(setcheck1.gt.0.d0) then
	 setcheck  = setcheck1
      else
	 setcheck1 = setcheck
      endif

      rminh    = 0.75d0*setcheck
      if(rminh.lt.0.1d0) rminh = 0.1d0
      rmint    = rminh/5.d0
      rming    = rminh/grad1

      setcheck2 = 15.d0*setcheck

      if(czo.eq.'f') czo='F'
      if(czo.eq.'d') czo='D'
      if(czo.eq.'b') czo='B'
      if(czo.ne.'F' .and. czo.ne.'D' .and. czo.ne.'B' ) czo='F'
      zo = zo1

      if(epilat0.ne.-999.d0 .and. epilon0.ne.-999.d0) epistart = .true.

      if(inddiff.eq.0) then
	 diffflag = .false.
      endif

      do 3 i=1,mstat
      sta(i) = '     '
3     continue

c
c     read in all available data
c

      open (unit=10,file=inputfile,status='old')

      title = version
      if(output) then
         open (unit=11,file=outputfile)
         write (11,'(a,/)') title(1:trimle(title))
      endif
      print *,'PROGRAM ',title(1:trimle(title))
      print *, ' '

      read (10,'(a)') title
      if(output) write (11,'(a,/)') title
      print *,'EVENT ',title

      ii = 0
      timemin = 9999999999999.d0

      stalam = 0.d0
      stalom = 0.d0

      terrm = 0.d0
      nobsst = 0

      sdpmean = 0.d0
      nsdp    = 0
      sdsmean = 0.d0
      nsds    = 0
      sdmeans = 0.d0

      do 12 i=1,mread+100

      read (10,'(a)',end=14) string

      if(string(1:1).eq.'*') go to 12
      if(string(1:1).eq.' ') go to 12

      ii = ii + 1

      if(ii.gt.mread) then
	print *,'Maximum number of input data reached: ',mread
	go to 9999
      endif

      ierr = 0

4     if(.not.magflag) then
         read (string,'(a5,1x,a8,1x,i4,4(1x,i2),1x,f6.3,1x,
     +               f5.3,1x,f6.2,3(1x,f5.2),1x,a6)',err=5) 
     +               stat,phase(ii),yy,mon,dd,hh,mi,sec,
     +               tts(ii),azi(ii),azis(ii),pin,ps(ii),
     +               touse0
      else
         read (string,'(a5,1x,a8,1x,i4,4(1x,i2),1x,f6.3,1x,
     +               f5.3,1x,f6.2,3(1x,f5.2),1x,a6,1x,f6.3,1x,
     +               f12.2)',err=5) 
     +               stat,phase(ii),yy,mon,dd,hh,mi,sec,
     +               tts(ii),azi(ii),azis(ii),pin,ps(ii),
     +               touse0,period(ii),amplit(ii)
      endif

      ierr = 0
      go to 6 

5     continue

      if(ierr.le.1) then

         indcs= 90
         indcs = index(string,'#')
         if(indcs.le.80 .and. indcs.gt.0) then
	    string(indcs:indcs)=' '
            go to 5
         endif
	 ierr = ierr + 1
         go to 4

      else

         print *,'Reading error with onset:',string
         ii = ii - 1
         go to 12

      endif

6     mm = '   '
      idoy = 0
      idum = 0
      timeo = 0.d0
      call fhtoe(timeo,idum,yy,mon,mm,dd,idoy,hh,mi,sec)
      tt(ii) = timeo

      touse(ii)=   'TASDR  '
      if(touse0.ne.'       ') then
	 if(touse0(1:1).eq.'t') touse0(1:1)='T'
	 if(touse0(1:1).ne.'T') touse0(1:1)=' '
	 if(touse0(2:2).eq.'a') touse0(2:2)='A'
	 if(touse0(2:2).ne.'A') touse0(2:2)=' '
	 if(touse0(3:3).eq.'s') touse0(3:3)='S'
	 if(touse0(3:3).ne.'S') touse0(3:3)=' '
	 if(touse0(4:4).eq.'d') touse0(4:4)='D'
	 if(touse0(4:4).ne.'D') touse0(4:4)=' '
	 if(touse0(5:5).eq.'r') touse0(5:5)='R'
	 if(touse0(5:5).ne.'R') touse0(5:5)=' '
	 if(touse0(6:6).ne.'2' .and. touse0(6:6).ne.'3' .and.
     +      touse0(6:6).ne.'4') touse0(6:6)=' '
	 touse(ii)=touse0
      endif

      if(touse(ii)(1:1).eq.'T' .or. touse(ii)(4:4).eq.'D') then
         if(tts(ii).le.0.d0) tts(ii) = 2.d0
         terrm  = terrm + tts(ii)
         nobsst = nobsst + 1
      endif

      if(phase(ii).eq.'PKP2')   phase(ii)="P'P'"
      if(phase(ii).eq.'PKPPKP') phase(ii)="P'P'"
      if(phase(ii).eq.'SKS2')   phase(ii)="S'S'"
      if(phase(ii).eq.'SKSSKS') phase(ii)="S'S'"
      if(phase(ii).eq.'PPP')    phase(ii)='P3'
      if(phase(ii).eq.'SSS')    phase(ii)='S3'
      if(phase(ii).eq.'PKhKP')  phase(ii)='PKPpre'

      incap = index(phase(ii),'diff')
      if(incap.ne.0) phase(ii)(incap+3:) = phase(ii)(incap+4:)
      incap = index(phase(ii),'DIFF')
      if(incap.ne.0) phase(ii)(incap:) = 'dif' // phase(ii)(incap+4:)

62    phid = phase(ii)
      incap =  index(phid,'N')
      if(incap.ne.0) then
	 phase(ii)(incap:incap) = 'n'
	 go to 62
      endif
      incap =  index(phid,'B')
      if(incap.ne.0) then
	 phase(ii)(incap:incap) = 'b'
	 go to 62
      endif
      incap =  index(phid,'A')
      if(incap.ne.0) then
	 phase(ii)(incap:incap) = 'a'
	 go to 62
      endif
      incap =  index(phid,'G')
      if(incap.ne.0) then
	 phase(ii)(incap:incap) = 'g'
	 go to 62
      endif
      incap =  index(phid,'C')
      if(incap.ne.0) then
	 phase(ii)(incap:incap) = 'c'
	 go to 62
      endif
      incap =  index(phid,'D')
      if(incap.ne.0) then
	 phase(ii)(incap:incap) = 'd'
	 go to 62
      endif

      if(phid(2:2).eq.' ') touse(ii)(5:5)=' '

      incap =  index(phid,'1')
      if(incap.ne.0 .and. phid(3:3).eq.' ') touse(ii)(5:5)=' '

      if(phid(4:4).eq.' ') then
         if(phid(2:2).eq.'b') touse(ii)(5:5)=' '
         if(phid(2:2).eq.'g') touse(ii)(5:5)=' '
         if(phid(2:2).eq.'n') touse(ii)(5:5)=' '
      endif

      if(phid(2:2).eq.'c') touse(ii)(5:5)=' '
      if(phid(2:2).eq.'d') touse(ii)(5:5)=' '
      if(phid(2:2).eq.'K') touse(ii)(5:5)=' '

      if(touse(ii)(1:1).ne.' ' .or. touse(ii)(4:4).ne.' ') then
         phase_t = phase_type(phid)
         if(phase_t.eq.'P') then
            nsdp    = nsdp + 1
            sdpmean = sdpmean + tts(ii)*tts(ii)
         endif
         if(phase_t.eq.'S') then
            nsds    = nsds + 1
            sdsmean = sdsmean + tts(ii)*tts(ii)
         endif
      endif

      if(azi(ii).eq.-1.d0) then
	 azi(ii) = -999.d0
	 azis(ii)=    0.d0
      endif

      if(azi(ii).eq.-999.d0) touse(ii)(2:2)=' '

c
c     The standard errors for azimuth or ray parameter are yet not 
c     given: set default values (30 or 40 [deg] and 5 [s/deg])!
c
      if(azi(ii).ne.-999.d0 .and. azis(ii).le.0.d0) then

	  azis(ii)= 30.d0

	  if(phase(ii)(1:2).eq.'LR'.or.phase(ii)(1:2).eq.'LQ' )
     +	     azis(ii)=40.d0

      endif

      do 10 j=1,mstat

      if(stat.eq.sta(j)) then

        iev(ii) = j
	go to 11

      else if (sta(j).eq.'     ') then

        istat = j
        call get_station(stationfile,stat,lat,lon,elevs,name,ierr)
	if(ierr.ne.0) then
	  print *,'Cannot find station: ',stat,' entry skipped'
	  ii = ii - 1
	  ierr = 0
	  go to 12
	endif

	if(vlflag .and. stcorfl) then
	  rewind(12)
7	  read(12,*,end=8,err=8) stat1,vpc,vsc
	  if(stat1.eq.stat) then
	    vp = vpc
	    vs = vsc
	    if(vs.le.0.d0) vs = vp / dsqrt(3.d0)
	    go to 9
	  endif
	  go to 7
        endif

8       if(vlflag) then
	  vp = vpl
	  vs = vsl
	endif

9	continue

 	if (typctl.gt.8) then
 	  print *,stat,lat,lon,elevs,name,vp,vs
 	endif

	sta(istat)   = stat
	stala(istat) = lat
	stalae(istat)= convlat(lat,1)
	stalo(istat) = lon
	stael(istat) = dble(elevs)/1000.d0
	stavp(istat) = vp
	stavs(istat) = vs
	istaph(istat) = 0
	stalam = stalam + stalae(istat)
	stalom = stalom + lon
	iev(ii) = istat

	drad = radloc(lat)/rearth

        go to 11

      endif

10    continue

11    if(timeo.lt.timemin) then
	 timemin = timeo
	 istatmin = iev(ii)
      endif

      if(pin.eq.-1.d0) then
	 pin    = -999.d0
	 ps(ii) =  0.d0
      endif

      if(islow.eq.0 .and. pin.ne.-999.d0) then
	 p(ii)  = drad/(grad1*pin)
	 ps(ii) = ps(ii)*p(ii)/pin
      else
	 p(ii) = pin
      endif

      if(p(ii).eq.-999.d0) touse(ii)(3:3) = ' '

      if(p(ii).ne.-999.d0 .and. ps(ii).le.0.d0) ps(ii)= 5.d0

      phase_t = phase(ii)(1:1)

      if(phase_t.eq.'P')  then
         if (istaph(iev(ii)).eq.0 .or. istaph(iev(ii)).eq.2) 
     +       istaph(iev(ii))=istaph(iev(ii)) + 1
      else if(phase_t.eq.'S') then
         if (istaph(iev(ii)).eq.0 .or. istaph(iev(ii)).eq.1)
     +       istaph(iev(ii))=istaph(iev(ii)) + 2
      endif

      if(phase(ii)(1:1).eq. 'P' .and. phase(ii)(3:3).eq. ' ')
     +             touse(ii)(7:7) = 'P'
      if(phase(ii)(1:3).eq. 'PKP' .and. phase(ii)(6:6).eq. ' ') 
     +             touse(ii)(7:7) = 'P'
      if(phase(ii)(1:4).eq. 'Pdif') touse(ii)(7:7) = 'P'

      if(typctl.gt.8) then
         print *,stat,tt(ii),phase(ii),azi(ii),azis(ii),p(ii),ps(ii),
     +           touse(ii)
      endif

12    continue

14    if(stcorfl) close (12)

      if(nobsst.gt.0) then
         terrm = terrm / dble(nobsst)
      else
	 terrm = 2.d0
      endif
      terrm = terrm * 10.d0 

      nobs  = ii
      nstat = istat
      stalam = stalam / dble(nstat)
      stalom = stalom / dble(nstat)

      if(nsdp.gt.0) then
	 sdmeans = sdpmean
         sdpmean = 2.d0*dsqrt(sdpmean / dble(nsdp))
         if(sdpmean.lt.1.d0) sdpmean = 1.d0
      else
	 sdpmean = 2.d0
      endif

      if(nsds.gt.0) then
	 sdmeans = sdmeans + sdsmean
         sdsmean = 2.d0*dsqrt(sdsmean / dble(nsds))
         if(sdsmean.lt.2.0d0) sdsmean = 2.0d0
      else
	 sdsmean = 4.0d0
      endif

      if(nsdp.gt.0 .or. nsds.gt.0 ) then
         sdmeans = dsqrt( sdmeans / dble(nsdp+nsds) )
      else
	 sdmeans = 5.d0
      endif

      close(10)

      if(nobs.eq.1) then

	 rzo1  = 0.

	 rdel1 = 3.
         call tauget_mod(rzo1,rdel1,nphas,phcd,ttc,dtdd,
     +                   dtdh,modnam)
	 pmoh  = dble(dtdd(1))

	 rdel1 = 105.
         call tauget_mod(rzo1,rdel1,nphas,phcd,ttc,dtdd,
     +                   dtdh,modnam)
	 pdif  = dble(dtdd(1))

	 single = .true.

	 if(typctl.gt.4) print *, 'Case: single array observation!'

      endif

c
c     Fix the time system at the earliest onset time for this event.
c

      do 15 i = 1,nobs
      tt(i) = tt(i)-timemin
15    continue

      if(epistart) then
      
	 elatm  = convlat(epilat0,1)
         elonm  = epilon0

         go to 65

      endif

c
c     At first, let us try to calculate an epicenter from all 
c     available azimuth observations.
c

      sela  = 0.d0
      svla  = 0.d0
      selo1 = 0.d0
      selo2 = 0.d0
      svlo  = 0.d0
      azim  = 0.d0
      azimr = 0.d0
      iazim = 0
      rpar  = 0.d0
      istater = 0

      jj = 0

      if (nobs.eq.1) then
         azim = azi(1)
	 if(p(1).ne.-999.d0) rpar = p(1)
         iazim = 1
         istataz = iev(1)
	 go to 51
      endif

      do 50 i=1,nobs-1

      if(touse(i)(2:2).eq.' ') go to 50
      azi1 = azi(i)
      if(index(phase(i),'pre').gt.0) go to 50
      if(index(phase(i),'KK').gt.0 .and. phase(1).ne.'S') 
     +         azi1 = alpha2(azi1-180.d0)
      if(index(phase(i),'P2K').gt.0) azi1 = alpha2(azi1-180.d0)
      if(index(phase(i),'P3K').gt.0) azi1 = alpha2(azi1-180.d0)
      if(index(phase(i),"P'P'").gt.0) azi1 = alpha2(azi1-180.d0)
      if(index(phase(i),"S'S'").gt.0) azi1 = alpha2(azi1-180.d0)

      azim = azim + azi1
      iazim = iazim + 1
      istataz = iev(i)

      slat1  = stala(iev(i))
      slat1e = stalae(iev(i))
      slon1  = stalo(iev(i))

      do 20 j=i+1,nobs

      if(touse(j)(2:2).eq.' ') go to 20
      azi2 = azi(j)
      if(index(phase(j),'pre').gt.0) go to 20
      if(index(phase(j),'KK').gt.0 .and. phase(1).ne.'S')       
     +         azi1 = alpha2(azi1-180.d0)
      if(index(phase(j),'P2K').gt.0) azi2 = alpha2(azi2-180.d0)
      if(index(phase(j),'P3K').gt.0) azi2 = alpha2(azi2-180.d0)
      if(index(phase(j),"P'P'").gt.0) azi2 = alpha2(azi2-180.d0)
      if(index(phase(j),"S'S'").gt.0) azi2 = alpha2(azi2-180.d0)

      if(i.eq.nobs-1 .and. j.eq.nobs) then
	 azim = azim + azi2
	 iazim = iazim + 1
	 istataz = iev(i)
      endif
      if(iev(i).eq.iev(j)) go to 20

      slat2  = stala(iev(j))
      slon2  = stalo(iev(j))

      jj = jj + 1

      if(jj.gt.mloc) then
	 print *,'Something wrong with number of locations!'
	 go to 9999
      endif

c     if(typctl.gt.7 ) then
c        print*,sta(iev(i)),'slat1,slon1,azi1,azis(i),slat2,',
c    +	        'slon2,azi2,azis(j)', slat1,slon1,azi1,
c    +          azis(i),sta(iev(j)),slat2,slon2,azi2,azis(j)
c     endif

c
c     Calculate distance and angles between the 2 stations
c
      call depi (slat1,slon1,slat2,slon2,del3,dk,ep2,ep1,d2km)
 
c     if(typctl.gt.7) then
c        print *,'station 1 (lat,lon,distance [km,deg],azimuth: ',
c    +           slat1,slon1,dk,del3,ep1
c        print *,'station 2 (lat,lon,azimuth; ',slat2,slon2,ep2
c     endif

c
c     Now a starting epicenter will be calculated
c

      ierr = 0
      call hyposat_cross(slat1e,slon1,azi1,azis(i),
     +               azi2,azis(j),del3,ep1,ep2,ela(jj),elas(jj),
     +               elo(jj),elos(jj),typctl,ierr)
      
      if(ierr.gt.0) then
	 jj = jj - 1
	 istater = istater + 1
	 ierr = ierr + 1
	 go to 20
      endif

      sela  = sela + ela(jj)/elas(jj)
      svla  = svla + 1.d0 / elas(jj)

      p1 = deg2rad*elo(jj)
      p2 = 1.d0/(deg2rad*elos(jj))

      selo1  = selo1 + dcos(p1)*p2
      selo2  = selo2 + dsin(p1)*p2
      svlo   = svlo  + p2

20    continue
50    continue

51    nloc = jj

      if(nloc.eq.0) then

	if(rpar .gt. 0.d0) then

	   phase_t = ' '
	   phidr0 = phase(1)
	   if(phidr0.eq.'P1') phidr0 = 'P'
	   if(phidr0.eq.'S1') phidr0 = 'S'

           rayokf = tauget_ray(phidr0,phase_t,rpar,modnam,zo,
     +                     ddel,ttray)

	   if(.not.rayokf) then

             ddel =  (14.d0 - rpar)*9.8d0 + 2.0d0
	     phidr0 = 'P'

	     if(rpar.lt.pdif) then
	        ddel = 150.d0
	        phidr0 = 'PKPdf'
	     endif

	     if(rpar.ge.10.0d0) then
	        ddel = 23.d0
	     endif

	     if(rpar.ge.13.0d0) then
	        ddel = 10.d0
	        phidr0 = 'Pn'
	     endif

	     if(rpar.ge.14.8d0) then
	        ddel = 2.d0
	        phidr0 = 'Pb'
	     endif

	     if(rpar.ge.17.0d0) then
	        ddel = 1.d0
	        phidr0 = 'Pg'
	     endif
	   
           endif
         
           if(typctl.gt.0 .and. .not.rayokf) then
  	      print *,'No distance found. Missing slowness values?'
  	      print *,'Distance set to ',ddel
	   else if(typctl.gt.4) then
              print *,'Starting distance from Station(net): ',ddel
	   endif

	endif
  
	if(iazim.eq.0 .or. istater.ne.0) then

c
c       Choose a point in the vicinity (1 deg) of the closest
c       station as starting solution:
c

	   istatd = istatmin
	   azim = 315.d0

	   if(typctl.gt.0 .and. .not.rayokf) then
  	      print *,'No epicenter found. Missing azimuth values?'
  	      print *,'Distance set to ',azim
	   endif

	else

	   istatd = istataz
	   azim   = azim/dble(iazim)
           azimr  = deg2rad*azim

	endif

	inddel = 1
	if (ddel.le.0.d0) ddel = 1.d0

	call delazd(stala(istatd),stalo(istatd),azim,ddel,
     +               inddel,elatmg,elonm)
	elatm = convlat(elatmg,1)

	go to 65

      endif

c
c     Now mean source coordinates have to be calculated.
c

      elatm  = sela / svla

      elonm  = rad2deg*datan2(selo2/svlo,selo1/svlo)

      if(nloc.eq.1) then

	sdlat  = elas(1)
	sdlon  = elos(1)

      else

        dla = 0.d0
        dlo = 0.d0
        do 60 i =1,nloc
        dla = dla + q2(ela(i)-elatm) / elas(i)
	p1 = alpha1(elo(i)-elonm)
        dlo = dlo + (p1*p1) / elos(i)
60      continue
        sdlat  = dsqrt(dla / svla)
        sdlon  = dsqrt(dlo / svlo)

      endif

c
c     The next step is to calculate a first source time. If already 
c     given with input parameter file, this given value will be used.
c
c     We are using the method of Wadati. If we have only one travel
c     time difference S-P we assume Vp/Vs = sqrt(3.). Otherwise
c     we calculate Vp/Vs as a constants for each specific phase type
c     (Pg,Pb,Pn,P,P1)-(Sg,Sb,Sn,S,S1).
c
c     Which travel-time differences between S and P onsets do we have?
c     These travel differences are also used to calculate a source 
c     distance from a station or an array.
c

65    dtkm  = 0.d0
      idtkm = 0

      do 66 i = 1,idtmax
      idt(i)=0
66    continue

      istatd = istatmin
      if(azimr.ne.0.d0) istatd = istataz

      if (nobs.eq.1) go to 71

      do 70 i = 1,nobs-1
      do 70 j = i+1,nobs

      if(iev(i).eq.iev(j)) then

	kmdel = 0.d0

        if((phase(i).eq.'P       ' .and. phase(j).eq.'S       ') .or.
     +	   (phase(i).eq.'P1      ' .and. phase(j).eq.'S1      ') .or.
     +     (phase(i).eq.'S       ' .and. phase(j).eq.'P       ') .or.
     +     (phase(i).eq.'S1      ' .and. phase(j).eq.'P1      ')) then

	  idt(1)        = idt(1) + 1
	  dt(1,idt(1))  = dabs(tt(j)-tt(i))

	  if(phase(i)(1:1).eq.'P') then
	     idtp(1,idt(1))= i
	     idts(1,idt(1))= j
	  else
	     idtp(1,idt(1))= j
	     idts(1,idt(1))= i
	  endif

	  if(iev(i).eq.istatd) then
	    kmdel = (dt(1,idt(1))/60.d0-1.d0)*1000.d0
	    if(kmdel.le.2000.d0) kmdel = dt(1,idt(1))*10.2d0
          endif

        else if(phase(i).eq.'Pg      '                           .and.
     +      (phase(j).eq.'Sg      ' .or. phase(j).eq.'Lg      ')) then

	  idt(2)        = idt(2) + 1
	  dt(2,idt(2))  = tt(j)-tt(i)
	  idtp(2,idt(2))= i
	  idts(2,idt(2))= j

	  if(iev(i).eq.istatd) kmdel = dt(2,idt(2))*8.58d0

        else if(phase(j).eq.'Pg      '                           .and.
     +      (phase(i).eq.'Sg      ' .or. phase(i).eq.'Lg      ')) then

	  idt(2)        = idt(2) + 1
	  dt(2,idt(2))  = tt(i)-tt(j)
	  idtp(2,idt(2))= j
	  idts(2,idt(2))= i

	  if(iev(i).eq.istatd) kmdel = dt(2,idt(2))*8.58d0

        else if((phase(i).eq.'Pn      '.and.phase(j).eq.'Sn      ') .or.
     +          (phase(j).eq.'Pn      '.and.phase(i).eq.'Sn      '))then

	  idt(3)        = idt(3) + 1
	  dt(3,idt(3))  = dabs(tt(j)-tt(i))

	  if(phase(i)(1:1).eq.'P') then
	     idtp(3,idt(3))= i
	     idts(3,idt(3))= j
	  else
	     idtp(3,idt(3))= j
	     idts(3,idt(3))= i
	  endif

	  if(iev(i).eq.istatd) kmdel = dt(3,idt(3))*10.2d0

        else if((phase(i).eq.'Pb      '.and.phase(j).eq.'Sb      ') .or.
     +          (phase(j).eq.'Pb      '.and.phase(i).eq.'Sb      '))then

	  idt(4)        = idt(4) + 1
	  dt(4,idt(4))  = dabs(tt(j)-tt(i))

	  if(phase(i)(1:1).eq.'P') then
	     idtp(4,idt(4))= i
	     idts(4,idt(4))= j
	  else
	     idtp(4,idt(4))= j
	     idts(4,idt(4))= i
	  endif

	  if(iev(i).eq.istatd) kmdel = dt(2,idt(4))*9.47d0

        else if(phase(i).eq.'Pn      '                          .and.
     +      (phase(j).eq.'Sg      '.or. phase(j).eq.'Lg      ')) then

	  if(iev(i).eq.istatd) kmdel = (tt(j)-tt(i))*6.02d0

        else if(phase(j).eq.'Pn      '                          .and.
     +      (phase(i).eq.'Sg      '.or. phase(i).eq.'Lg      ')) then

	  if(iev(i).eq.istatd) kmdel = (tt(i)-tt(j))*6.02d0

	endif

	if(kmdel.ne.0.d0) then
	   idtkm = idtkm + 1
	   dtkm  = dtkm + kmdel
	endif

      endif

70    continue

71    inet = 0

      if(idtkm.ne.0 .and. nloc.eq.0 .and. .not.epistart) then

         dtkm = dtkm / dble(idtkm)

         if(azimr.ne.0.d0) then

	    inddel = 2
	    call delazd(stala(istataz),stalo(istataz),azim,dtkm,
     +               inddel,elatmg,elonm)
            elatm = convlat(elatmg,1)
 
            if(typctl.gt.0) then
               print *,'Epicenter set from station ',
     +	             sta(istataz),': backazimuth',azim,' deg, delta',
     +               dtkm,' km' 
	    endif

            sdlatg = dtkm*grad1
	    sdlon  = dtkm*grad1

         else 

	    if(plflag) then

	       call plane(stala,stalo,iev,tt,nobs,azim,dazi,
     +                    ray,dray,phipl,touse,phase,jref,typctl)

	       if(jref.gt.0 .and. dazi.lt.90.d0 .and. dray.lt.4.d0) then

		   inddel = 2
		   call delazd(stala(iev(jref)),stalo(iev(jref)),azim,
     +                         dtkm,inddel,elatmg,elonm)
		   elatm = convlat(elatmg,1)

                   if(typctl.gt.0) then
                      print *,'Epicenter set from station ',
     +	                sta(iev(jref)),'after plane wave fit: ',
     +                  'backazimuth',azim,' deg, delta',dtkm,' km' 
	           endif

                   sdlatg = dazi*dtkm*grad1
	           sdlon  = dazi*dtkm*grad1

		   go to 72

                endif

	    endif

  	    if(dtkm.gt.120.d0 .and. nstat.gt.1) then
	    
	       elatm = stalam
	       elonm = stalom

	       inet = 1

               print *,'Epicenter set in center of station net '

	    endif

	 endif

      endif

72    elatmr  = deg2rad*elatm
      elatmg  = convlat(elatm,2)
      elatmgr = deg2rad*elatmg
      coelatm = 90.d0 - elatm
      coelatmr= deg2rad*coelatm
 
      if(sdlat.ne.0.d0) then
       sdlatg= sdlat / (eps*q2(dcos(elatmr))+q2(dsin(elatmr)))
      else
       sdlat = sdlatg*eps/(q2(dcos(elatmgr))+eps*q2(dsin(elatmgr)))
      endif

      if(typctl.gt.0) then
        if(nloc.gt.0) then
	   print*,'Mean epicenter calculated from ',nloc,
     +	          ' observation(s)'
        else
	   print*,'Set epicenter '
        endif
	print*,'(Mean) epicenter lat: ',elatmg,' +/- ',sdlatg
	print*,'(Mean) epicenter lon: ',elonm,' +/- ',sdlon
      endif

      if(output) then

        write(11,'(/''Parameters of starting solution ('',
     +              ''+/- 1 standard deviation):''/)') 

        if(nloc.gt.0) then

           write(11,'(''Mean epicenter calculated from'',i5,
     +              '' backazimuth observation pairs'')') nloc
           write(11,'(''Mean epicenter lat:'',f9.4,'' +/- '',f9.4,
     +              '' [deg]'')')        elatmg,sdlatg
           write(11,'(''Mean epicenter lon:'',f9.4,'' +/- '',f9.4,
     +              '' [deg]''/)')        elonm,sdlon

        else if(.not.epistart) then

  	   write(11,'(''No location from multiple azimuth '',
     +            ''observations found. Missing azimuth values?'')')

           if(inet.ne.0) then
              write(11,'(''Epicenter set in the center of station'',
     +	                  '' net '')')
	   else if (azimr.ne.0.d0 .and. dtkm.eq.0.d0) then
              write(11,'(''Epicenter set from station '',
     +	          a,'' with backazimuth'',f6.1,'' [deg], delta'',f7.2,
     +            '' [deg]'')') sta(istataz),azim,ddel
	   else 
              write(11,'(''Epicenter set from station '',
     +	          a,'' with backazimuth'',f6.1,'' [deg], delta'',f7.0,
     +            '' [km]'')') sta(istatd),azim,dtkm
	   endif

           write(11,'(''Epicenter lat:'',f9.4,'' [deg]'')')  elatmg
           write(11,'(''Epicenter lon:'',f9.4,'' [deg]''/)') elonm

        else if(epistart) then

           write(11,'(''Starting Epicenter set by input file'')')
	   write(11,'(''Epicenter lat:'',f9.4,'' [deg]'')')  elatmg
           write(11,'(''Epicenter lon:'',f9.4,'' [deg]''/)') elonm
	
        endif
      endif

      if(typctl.gt.7) then

        do 75 i=1,nstat

        call depi(stala(i),stalo(i),elatmg,elonm,del(i),dk,azie(i),
     +	          baz(i),d2km)
	if(azi(i).ne.-999.) then
          print *,sta(i),del(i),azie(i),baz(i),azi(i),azi(i)-baz(i)
	else
          print *,sta(i),del(i),azie(i),baz(i)
	endif

75      continue

      endif

      do 80 i=1,idtmax

      if(idt(i).eq.0) go to 80

      if(idt(i).eq.1) then

	vpvs(i) = dsqrt(3.d0)
	f1      = 1.d0/(vpvs(i)-1.d0)

	to(i)   = tt(idtp(i,1)) - dt(i,1)*f1
	tos(i)  = dsqrt( q2( (1.d0+f1)*tts(idtp(i,1)) ) +
     +	                 q2(       f1 *tts(idts(i,1)) )  )
	vpvss(i)= 0.5d0

      else if(idt(i).eq.2) then

	f1 = dt(i,2)-dt(i,1)
	f2 = 1.d0 / (tt(idtp(i,2))-tt(idtp(i,1)))

	am      = f1*f2
	am1     = 1.d0 /am
	vpvs(i) = am + 1.d0
	to(i)   = tt(idtp(i,1)) - dt(i,1)*am1

	f3 = am*f2
	f4 = am1*am1
	vpvss(i)= dsqrt ( q2(( f2+f3)*tts(idtp(i,1))) +
     +                    q2((-f2-f3)*tts(idtp(i,2))) +
     +                    q2(  f2    *tts(idts(i,1))) +
     +                    q2( -f2    *tts(idts(i,2))) )
	
	tos(i) = dsqrt( q2( (1.d0+am1+dt(i,1)*(f2+f3)*f4 )
     +	                                 *tts(idtp(i,1)) )   +
     +                  q2( (         dt(i,1)*(-f2-f3)*f4 )
     +                                   *tts(idtp(i,2)) )   +
     +                  q2( (    -am1+dt(i,1)*f2     *f4 )
     +                                   *tts(idtp(i,2)) )   +
     +                  q2( (        -dt(i,1)*f2     *f4 )
     +                                 *tts(idtp(i,2)) )   )

      else 

	do 77 j=1,idt(i)

        f1 = 1.d0
	if(tts(idtp(i,j)).ne.0.d0 .or. tts(idts(i,j)).ne.0.d0)
     +    f1 = 1.d0/dsqrt(q2(tts(idtp(i,j))) + q2(tts(idts(i,j))))

 	f2 = dt(i,j) * f1
 	f3 = tt(idtp(i,j)) * f1

	dat(j) = f2
	a(j,1) = f3
	a(j,2) = f1

c	if(typctl.gt.6) then
c          print *,f2,a(j,2),f3
c        endif

77      continue
        
        im = 2
	in = idt(i)
	call dlsq(in,im)

	if(typctl.gt.6) then
	  print *,r(1),var(1)
	  print *,r(2),var(2)
	  do 78 j=1,idt(i)
78        print *,dat(j),res(j)
        endif


	vpvs(i) = r(1) + 1.d0
 	vpvss(i)= var(1)
        to(i)   = -r(2) / r(1)
 	tos(i)  = dsqrt(q2(to(i)*var(2)/r(2)) +
     +	                q2(to(i)*var(1)/r(1))   )

      endif

      if(i.eq.1 .and. dabs(to(i)).gt.1300.d0) then
	 to(i)=dsign(1300.d0,to(i))
	 tos(i)=dabs(to(i))
      endif
      if(i.eq.2 .and. dabs(to(i)).gt.150.d0) then
	 to(i)=dsign(150.d0,to(i))
	 tos(i)=dabs(to(i))
      endif
      if(i.eq.3 .and. dabs(to(i)).gt.400.d0) then
	 to(i)=dsign(400.d0,to(i))
	 tos(i)=dabs(to(i))
      endif
      if(i.eq.4 .and. dabs(to(i)).gt.250.d0) then
	 to(i)=dsign(250.d0,to(i))
	 tos(i)=dabs(to(i))
      endif

      if(typctl.gt.0) then
         print *,'S-P Travel-time difference type ',i
         print *,'Source time from ',idt(i),' observation(s)'
         print *,'   to= ',to(i),  ' +/- ',tos(i)
         print *,'Vp/Vs= ',vpvs(i),' +/- ',vpvss(i)
      endif
      if(output) then
         write(11,'(''S-P Travel-time difference type'',i2,
     +          '' with'',i4,'' observation(s)'')') i,idt(i)
      endif

80    continue

c
c     now follows the statistics over all estimated to-values
c
      sto     = 0.d0
      stos    = 0.d0
      svpvs   = 0.d0
      svpvss  = 0.d0
      ito     = 0
 
      do 82 i = 1,idtmax
 
      if(idt(i).eq.0) go to 82

      ito = ito + 1

      sto   = sto  + to(i)/tos(i)
      stos  = stos + 1.d0 /tos(i)
 
      svpvs = svpvs  + vpvs(i)/vpvss(i)
      svpvss= svpvss + 1.d0   /vpvss(i)
 
82    continue
 
      dto   = 0.d0
      dvpvs = 0.d0
 
      if(ito.eq.0) then
	 if(nloc.ge.1) then

           call depi(stala(istatmin),stalo(istatmin),elatmg,elonm,
     +	       del(istatmin),dk,azie(istatmin),baz(istatmin),d2km)
	   rzo1 = 0.
	   rdel = real(del(istatmin))
	   call tauget_mod(rzo1,rdel,nphas,phcd1,ttc1,dtdd1,
     +	                 dtdh1,modnam)
           tom  = - dble(ttc1(1))

	 else

	   if(ttray.ne.0d0) then
	      tom  = -ttray
	   else
	      tom  = -dtp0/2.d0
	   endif

	 endif
	 sdto = dtp0
	 go to 85
      else
         tom   = sto / stos
         vpvsm = svpvs / svpvss
      endif
 
      if(ito.eq.1) then
	 sdto   = 1.d0 / stos
	 sdvpvs = 1.d0 / svpvss
	 go to 85
      endif

      do 83 i =1,idtmax

      if(idt(i).eq.0) go to 83
 
      dto   = dto   + q2(to(i)-tom) / tos(i)
      dvpvs = dvpvs + q2(vpvs(i)-vpvsm) / vpvss(i)
 
83    continue
 
      sdto    = dsqrt(dto / svpvs)
      sdvpvs  = dsqrt(dvpvs / svpvss)
 
85    continue

      tome = tom + timemin

      if(typctl.gt.0) then
        print*,'Mean source time: ',tome,' +/- ',sdto
        print*,'Mean       vp/vs: ',vpvsm,' +/- ',sdvpvs
      endif
      if(ito.gt.0) then
	 if(output) then
            write(11,'(''Mean source time:'',f15.3,'' +/- '',f7.3,
     +              '' [s]'')')  tome,sdto
            write(11,'(''Mean       vp/vs:'',f15.3,'' +/- '',f7.3)') 
     +               vpvsm,sdvpvs
	 endif
      else
	 if(output) then
            write(11,'(''Source time (set):'',f15.3,'' [s]'')') tome
	 endif
      endif

c
c     In any case, we use the starting source time and its standard
c     deviation from input file.
c     
      if (tome0 .ne. 0.d0) then
         tome = tome0

         if (stome0 .ne. 0.d0) then
	    sdto = stome0
         endif

         if(output) then
            write(11,'(/''Source time (from input):'',f15.3,'' +/- '',
     +                f7.3,'' [s]''/)') tome,sdto
         endif

      endif

c
c     Now new source parameters can be calculated by
c     several iterations using the GMI algorithm. 
c
c     For the first iteration we use as start solution the read in 
c     source depth, the source time to, and the epicenter coordinates 
c     elatmg and elonm.
c
      iter = 0

      if(czo.eq.'D') sdzo = sdzo1
      if(czo.eq.'F' .or. czo.eq.'B') sdzo = 1.d0

      rs(1) = sdto
      rs(2) = sdlat
      rs(3) = sdlon
      rs(4) = sdzo

      is(1) = 1
      is(2) = 0
      quot  = 1.d0

      ifail = 0
      nextiter = 0
      in = 0
      dtp = dtp0
      dts = dts0

      dchang = dchang0

c
c     we have a initial solution and will eventually not longer need the
c     azimuth information:
c
      if(aziini) then
	
	 do 99 i = 1,nobs
	   touse(i)(2:2) = '_'
99       continue

      endif
c
c     At first, we build the Jacobi-matrix
c     (loop 300 and 301)
c

100   continue

      last = .false.

      if (ibad0.gt.3) then
         print *,'Could not find a stable solution for these data'
	 if(output) then
            write(11,'(/''Could not invert these data!'')')
	 endif
	 go to 9999
      endif

      iter = iter + 1

      if((check.le.setcheck2 .or. nextiter1.eq.1) .and. 
     +    iteraz.eq.0 ) dtmflag = .true.

      if ((lastfix .and. check.le.setcheck2) .or. dtmflag) then

	 if(zo.le.0.1d0) zo=0.d0

	 f1 = dmax1(0.75d0,rmso)

         dtm    = f1 * 3.d0
         dazim0 = dmin1(15.d0,dazim1)
	 dpam0  = dmin1(2.d0,dpam1)

	 if (nrms1.gt.5) then
            dtm  = f1 * 2.d0
            dazim0 = dmin1(10.d0,dazim1)
	    dpam0  = dmin1(1.5d0,dpam1)
	 else if (nrms1.gt.10) then
            dtm    = f1 * 1.2d0
            dazim0 = dmin1(5.d0,dazim1)
	    dpam0  = dmin1(1.d0,dpam1)
	 endif

	 dtm0 = dtm * 2.d0

      else
         dtm  = dtm2
         dtm0 = dtm2
         dazim0 = dazim1
	 dpam0  = dpam1
	 dtmflag = .false.
      endif

      dtmp = dmin1(dtm,dtp)
      if(dtmp.lt.sdpmean) dtmp = sdpmean
      dtms = dmin1(2.d0*dtm,dts)
      if(dtms.lt.sdsmean) dtms = sdsmean
      dtm0 = dmin1(dtm0,dtm2)

      ifixaz = 0
      ifixto = 0

      rmsold = rmso
      iremo  = 0

101   continue

      if(dtmp.gt.1000.d0) dtmp = 1000.d0
      if(dtms.gt.1400.d0) dtms = 1400.d0
      if(dtm0.gt.1400.d0) dtm0 = 1400.d0
      if(dpam0.gt.15.d0)  dpam0 = 15.d0
      if(dazim0.gt.90.d0) dazim0 = 90.d0

      stato = '     '
      jj    = 0
      jdt   = 0
      jazi  = 0
      jpa   = 0
      rms1  = 0.d0
      nrms1 = 0
      datmax = 0.d0

      do 300 i = 1,nobs

      used(i)= '     '

      ttt(i) = 0.d0

      if (sta(iev(i)).ne.stato) then

	 stato = sta(iev(i))

	 imod2 = 1
	 if (touse(i)(6:6).eq.'2' .and. mod2flag) then
	   modn = modnam2
         else if (touse(i)(6:6).eq.'3' .and. mod3flag) then
	   modn = modnam3
         else if (touse(i)(6:6).eq.'4' .and. mod4flag) then
	   modn = modnam4
         else
	   if(iloc) imod2 = 0
	   modn = modnam
	 endif

         call depi(stala(iev(i)),stalo(iev(i)),elatmg,elonm,
     +	      del(iev(i)),dk,azie(iev(i)),baz(iev(i)),d2km)

c	 if(typctl.gt.8) then
c	   print *,'STATION EPI: ',stato,del(iev(i)),dk
c	 endif

	 rzo = sngl(zo)
	 rdel = sngl(del(iev(i)))
	 drdel = 0.1
	 if(single .and. (rdel.lt.20. .or. rdel.gt.100.)) 
     +      	 drdel=2.*drdel
	 rdel1= rdel-drdel
	 if(rdel1.lt.0.)   rdel1 = 0.
	 rdel2= rdel+drdel
	 if(rdel2.gt.180.) rdel2 = 180.

	 costalat  = 90.d0-stalae(iev(i))
	 costalatr = deg2rad*costalat

         fla1 = sngl(coelatm)
         fla2 = sngl(costalat)
	 fla3 = sngl(coelatmr)
	 razi = sngl(azie(iev(i)))
	 rbaz = sngl(baz(iev(i)))

         if(imod2.ne.0 .or. rdel.gt.rmax .or. zo.gt.zmax) then

	    nphas0 = 0
	    call tauget_mod(rzo,rdel,nphas0,phcd,ttc,dtdd,
     +	              dtdh,modn)


	    if(czo.eq.'D') then

	       rzo1 = rzo-1.
	       if(rzo1.lt.0.1) rzo1=0.1
	       rzo2 = rzo+1.

c
c              no event deeper than 799. km !
c
	       if(rzo2.ge.800.) then
	         rzo =799.
	         rzo1=798.
	         rzo2=799.
	       endif

               if(rzo-rzo1.gt.0.) then

	          nphas1 = 0
	          call tauget_mod(rzo1,rdel,nphas1,phcd1,ttc1,dtdd1,
     +	                    dtdh1,modn)

	       endif

               if(rzo2-rzo.gt.0.) then

	          nphas2 = 0
	          call tauget_mod(rzo2,rdel,nphas2,phcd2,ttc1,dtdd2,
     +	                    dtdh1,modn)

	       endif

               do 280 k=1,nphas0

	       dpdh(k) = 0.d0

	       if(rzo2-rzo1.eq.0.) go to 280

	       do 250 ki=1,nphas1
	       do 250 kj=1,nphas2
	       if(phcd1(ki).eq.phcd(k) .and. phcd(k).eq.phcd2(kj)) then

                  dpdh(k) = dble((dtdd2(kj)-dtdd1(ki))/(rzo2-rzo1))
	          go to 280
	       endif
250            continue

               if(rzo-rzo1.gt.0.) then
   	          do 260 ki=1,nphas1
	          if(phcd(k).eq.phcd1(ki)) then
	             dpdh(k) = dble((dtdd(k)-dtdd1(ki))/(rzo-rzo1))
	             go to 280
	          endif
260               continue
	       endif

   	       if(rzo2-rzo.gt.0.) then
	          do 270 kj=1,nphas2
	          if(phcd(k).eq.phcd2(kj)) then
	             dpdh(k) = dble((dtdd2(kj)-dtdd(k))/(rzo2-rzo))
	             go to 280
	          endif
270               continue
	       endif

280            continue

	    endif

	    if(rdel-rdel1.gt.0.)  then

	       nphas1 = 0
	       call tauget_mod(rzo,rdel1,nphas1,phcd1,ttc1,dtdd1,
     +                   dtdh1,modn)

	    endif

	    if(rdel-rdel2.lt.0.)  then

	       nphas2 = 0
	       call tauget_mod(rzo,rdel2,nphas2,phcd2,ttc1,dtdd2,
     +                   dtdh1,modn)

	    endif

            do 294 k=1,nphas0

	    dddp(k) = 0.d0

	    if (rdel2-rdel1.eq.0.) go to 294

	    do 291 ki=1,nphas1
	    do 291 kj=1,nphas2
	    if(phcd1(ki).eq.phcd(k) .and. phcd(k).eq.phcd2(kj)) then
	       dddp(k) = dble((dtdd2(kj)-dtdd1(ki))/(rdel2-rdel1))
	       go to 294
	    endif
291         continue

   	    if(rdel.ne.rdel1) then
	       do 292 ki=1,nphas1
	       if(phcd(k).eq.phcd1(ki)) then
	          dddp(k) = dble((dtdd(k)-dtdd1(ki))/(rdel-rdel1))
	          go to 294
	       endif
292            continue
	    endif

   	    if(rdel.ne.rdel2) then
	       do 293 kj=1,nphas2
	       if(phcd(k).eq.phcd2(kj)) then
	          dddp(k) = dble((dtdd2(kj)-dtdd(k))/(rdel2-rdel))
	          go to 294
	       endif
293            continue
	    endif

294         continue

	    nphas = nphas0

         else 

	    nphas = 0

	    ierr = 0
	    indph = istaph(iev(i))*10000 + indph0
	    elatc = elatmg
	    elonc = elonm

	    elat2 = stala(iev(i))
	    elon2 = stalo(iev(i))

	    call ttloc(rzo,rdel,czo,nphas,ttc,dtdd,dtdh,dpdh,dddp,
     +                 phcd,rmax,typctl,ierr,indph)

	    if(ierr.ne.0) then
	       print *, 'Error in getting travel-time tables for:'
	       print *, sta(iev(i)), ' in ',rdel,' deg --> not used!'
	       go to 300
	    endif

	 endif

	 f1 = dcos(coelatmr)
	 f3 = dsin(coelatmr)

	 f2 = dcos(costalatr)
	 f4 = dsin(costalatr)

	 f5 = deg2rad*alpha1(stalo(iev(i))-elonm)
	 f5 = dabs(f5)

	 f6 = dcos(f5)
         f8 = dsin(f5)

	 f7 = dsin(deg2rad*del(iev(i)))

	 if(baz(iev(i)).le.180.d0) then
            alpha = deg2rad*baz(iev(i))
	 else
            alpha = deg2rad*(360.d0-baz(iev(i)))
	 endif

   	 deldla =  (f4*f1*f6 - f3*f2) / f7
	 deldlo =  f4*dsin(alpha)

	 f9  = dcos(alpha) * f7*f7
	 f10 = dcos(deg2rad*del(iev(i)))

	 dazidla =  -f8*(f1*f7+f3*f10*deldla)/f9
  	 dazidlo =   f3*(f6*f7-f8*f10*deldlo)/f9

 	 if(baz(iev(i)).gt.180.d0) then
            dazidla = -dazidla
            deldlo  = -deldlo
 	 endif
	  
c         if(typctl.gt.8) then
c 	    print *,'baz ',baz(iev(i)),deldla,deldlo,dazidla,
c    +              dazidlo,(stalo(iev(i))-elonm),del(iev(i))
c
c         endif
      endif

      phid = phase(i)

      dpa = 0.d0

      if(phid(1:1).eq.'T') go to 298

      surf = .false.

      if(phid.eq.'Rg') then
	 if(dk.le.400.d0) then
	    surf = .true.
	    vsurf = vrg
	 else
	    phid = 'LR'
	 endif
      endif
		
      if(phid.eq.'Lg' .and. dk.le.3000.d0) then
	 surf = .true.
	 vsurf = vlg
      endif
		
      if(phid.eq.'LR') then
	 surf = .true.
	 vsurf = vlr
      endif
		
      if(phid.eq.'LQ') then
	 surf = .true.
	 vsurf = vlq
      endif

      if (surf) then
	 nphas       = nphas + 1
	 ttc(nphas)  = sngl(dk/vsurf)
	 dtdd(nphas) = sngl(d2km/vsurf)
	 dtdh(nphas) = 0.
	 phcd(nphas) = phid
	 dddp(nphas) = 0.d0
	 dpdh(nphas) = 0.d0
      endif
		
      icha = 0
      
      if((phid(2:3).eq.'n ' .or. phid(2:3).eq.'g ' .or. phid(2:3).eq.
     +   'b ') .and. rdel.gt.30. ) phid(2:3)='  '

295   continue

      if(single .and. insar.le.0 .and. iter.eq.1) then
	 phid = phidr0
	 phase(1) = phid
      endif
      
      nphass = nphas

      if(phase(i).eq.'P1') then
         nphass = 1
	 if (icha.ne.0) then
	    if (phid(1:1).eq.'P')   nphass = 2
            if (phid(1:2).eq.'PK' .or. phid(1:3).eq.'Pdi') nphass = 5
	 endif
      endif

      do 297 j = 1,nphass

c
c     Let's take the first phase from IASPEI-1991-type tables,
c     which fits the phase name of the onset
c
      phid1 = phcd(j)

      if(phid.eq.'P1' .and. phid1(1:1).eq.'P' .and. phid1(3:3).eq. 
     +   ' ') then

        phid = phid1

	if(p(i).lt.0.9d0*pdif .and.phid1(1:3).eq.'Pdi') then
	   if(del(iev(i)).gt.110.d0) then
	     phid = 'PKPdf'
	     go to 297
           endif
        endif

      endif

      if(phid.eq.'S1' .and. phid1(1:1).eq.'S' .and. phid1(3:3).eq.
     +   ' ') then

        phid = phid1

      endif

      if(phid.eq.'PKPdf' .and. nobs.gt.2) then
	 
	 if(phid1(1:5).eq.'PKiKP' .and. 
     +	               del(iev(i)).ge.90.d0) phid='PKiKP'

      endif

      if(icha.eq.0) then
 
          if(phid1(2:4).eq.'dif') then
             if(phid.eq.'P') phid='Pdif'
             if(phid.eq.'S') phid='Sdif'
          endif
 
      endif

      if(phid1.eq.phid) then

c     Excursion : any ellipticity correction for this phase?
c


	 pa   = dtdd(j)
	 dpa  = dble(pa)
	 dpaa = dabs(dpa)

	 if(single .and. rayokf) then

	    if(dabs(p(i)-dpaa) .ge. 0.01d0) go to 297

	 endif

	 rzoe = rzo

	 ierre = 0
	 ecor = 0.

	 if (.not.surf) then
c
c          Ellipticity corrections are yet not available for sources 
c          deeper than 700 km. Therefore we accept a small error
c          and set the source depth to 700. km
c
	   if(rzoe.gt.700.) rzoe=700.

	   if(iaspflag) then
             call ellip(fla1,fla2,razi,rbaz,
     +                  rdel,rzoe,phid1,pa,ecor,ierre)
	   else
             call ellip2(fla3,razi,rdel,rzoe,phid1,pa,ecor,ierre)
	   endif

	 endif

	 phase_t = phase_type(phid1)

         dtm = dtm0

         if(phase_t.eq.'P') dtm = dtmp
         if(phase_t.eq.'S') dtm = dtms
	 dtm = dmin1(dtm,datmax0)

	 th     = 0.d0

	 tcrust = 0.d0

	 if(vlflag .and. .not.surf) then

	   hsta = stael(iev(i))

	   vloc = 99999.d0
           if(phase_t.eq.'P') vloc = stavp(iev(i))
           if(phase_t.eq.'S') vloc = stavs(iev(i))

           if(vloc .ne. 99999.d0) then

              if(iabs(imo).eq.2 .or. imo.eq.4) then

                 elatc  = stalae(iev(i))
                 elonc  = stalo(iev(i))
                 indr = 1
                 tcrust = crustc(phase_t,dpaa,indr,typctl)

                 if(tcrust.ne.0.d0) hsta = hsta - elev

              endif

              if(hsta.ne.0.d0) then

                 radkm = deg2rad*radloc(stala(iev(i)))
                 phin = vloc*dpaa/radkm

                 if(phin.lt.1.d0) then

                    dl = hsta / dcos(dasin(phin))

                    if(dl.gt.0.d0) then
                      th = dl/vloc - 
     +                     dpaa*dsqrt(dl*dl-hsta*hsta)/radkm
                    else if(dl.lt.0.d0) then
                      th = dl/vloc + 
     +                     dpaa*dsqrt(dl*dl-hsta*hsta)/radkm
                    endif

                 endif
              endif
           endif
         endif

         if(typctl.gt.8) then
            print *,i,tt(i),tom,ttc(j),ecor,th
         endif

c
c        We have eventually to correct this phase for the local 
c        structure at the reflection point at the surface (if 
c        CRUST 5.1 is available).
c
	 trefl = 0.d0
	 if( .not.surf. and. touse(i)(5:5).eq.'R' .and.
     +       (rdel.gt.rmax .or. zo.gt.zmax .or. .not.iloc) ) then 

            phase_t = phid(1:1)

	    if(phid(1:1).eq.'p') phase_t = 'P'
	    if(phid(1:1).eq.'s') phase_t = 'S'
	     
	    fmult = 1.d0
	    del0  = dirdel(dpaa,zo,fmult,phase_t)
	    azi0  = azie(iev(i))

	    if(phid(1:1).ne.'p' .and. phid(1:1).ne.'s') then
	       if(dpa.ge.0.d0) then
	          del0 = (del(iev(i))-del0)/2.d0
               else
	          del0 = (360.d0 - del(iev(i)) - del0)/2.d0
		  azi0 = alpha2(azie(iev(i))-180.d0)
               endif
            endif
		
	    inddel = 1
	    call delazd(elatmg,elonm,azi0,del0,inddel,elatc,elonc)
    
c
c     correction for depth phases (e.g.: pP, sP...)
c 
	    if( (phid(1:1).eq.'p' .or. phid(1:1).eq.'s') 
     +          .and. zo.ge.zmax .and. phase_t.ne.' ') then

	      if ((phid(1:1).eq.'p' .and. phid(2:2).eq.'P') .or.
     +	          (phid(1:1).eq.'s' .and. phid(2:2).eq.'S')) then

		  indr = 2

	      else 
		
	          indr = 3

	      endif

	      trefl   = crustc(phase_t,dpaa,indr,typctl)
              used(i)(5:5) = 'R'

            endif

c
c     correction for surface reflections (e.g.: PnPn,...,PP,SS,P'P')
c
	    if( (phid(1:1).eq.phid(2:2) .or. 
     +	         phid(1:2).eq.phid(3:4)  ) .and. phase_t.ne.' ') then

	        indr = 2
	        trefl = crustc(phase_t,dpaa,indr,typctl)
	        used(i)(5:5) = 'R'
   
	    endif

c
c      correction for converted surface reflections (e.g.: PnSn,...)
c
	    conr = .false.
	    if( (phid(1:1).eq.'P' .or. phid(2:2).eq.'P') .and.
     +          (phid(1:1).eq.'S' .or. phid(2:2).eq.'S') .and.
     +          (phid(3:3).eq.' ' .or. phid(3:3).eq.'g' .or.
     +           phid(3:3).eq.'b' .or. phid(3:3).eq.'n')) then
		   conr=.true.
		   phidr = phid(2:)
	    endif

            if( (phid(1:1).eq.'P' .or. phid(3:3).eq.'P') .and.
     +          (phid(1:1).eq.'S' .or. phid(3:3).eq.'S') .and.
     +           phid(2:2).ne.'b' .and. phid(2:2).ne.'m' .and.
     +           phid(2:2).ne.'c' .and. phid(2:2).ne.'k' .and.
     +          (phid(2:2).eq.phid(4:4) .or. phid(2:2).eq.'g' .or.
     +           phid(2:2).eq.'n'                      )) then
		   conr=.true.
		   phidr = phid(3:)
	    endif

	    if(conr) then
		
		zor = 0.d0
		rayok = tauget_ray(phidr,phase_t,dpaa,modn,zor,
     +	                del0,ttray)

		azi0 = azie(iev(i))
		if(dpa.lt.0.d0) azi0 = alpha2(azie(iev(i))-180.d0)

	        inddel = 1
	        call delazd(stala(iev(i)),stalo(iev(i)),azi0,del0,
     +                      inddel,elatc,elonc)

		indr = 3
	        trefl   = crustc(phase_t,dpaa,indr,typctl)
                used(i)(5:5) = 'R'

	    endif

	    if (used(i)(5:5).ne.' ' .and.typctl.gt.6) then
	       print *,'dirdel: ',phid,' azi ',azi0,' del ',del0
     	       print *,'  lat ',elatc,' lon ',elonc,' trefl ',trefl
	    endif

	 endif

         ttt(i) = tom + dble(ttc(j) + ecor) + th + tcrust + trefl

	 dtt     = tt(i) - ttt(i)

c	 if(typctl.gt.8) then
c	    print*,'data input: onset time',tt(i),ttt(i),dtt
c        endif

	 if(dtt.gt.100.d0 .and. phase(i)(1:2).eq.'P1' .and. 
     +      phid(1:3).eq.'Pdi') then
	    phid='PKPdf'
            go to 297
         endif

	 if(dabs(dtt).le. dtm .and. touse(i)(1:1).eq.'T') then

           jj  = jj + 1
	   jdt = jdt + 1
	   rms1 = rms1 + dtt*dtt
	   nrms1= nrms1+ 1

	   dat(jj) = dtt
	   dats(jj)= tts(i)

	   if(dabs(dtt).gt.datmax) datmax=dabs(dtt/tts(i))

c
c          ierre .ne.0 > No ellipticity correction is available for this
c                       phase. We assume a larger data error!
c
	   if(ierre.ne.0) then
	      dats(jj)= dats(jj) + 1.d0
	      ierre = 0
	   endif

	   a(jj,1) = 1.d0
           a(jj,2) = dpa*deldla
           a(jj,3) = dpa*deldlo

	   a(jj,4) = dble(dtdh(j))
	   if(dabs(a(jj,4)).lt.1.d-5) a(jj,4)=0.d0

           used(i)(1:1) = 'T'

	   datla(i) = a(jj,2)
	   datlo(i) = a(jj,3) 
	   datho(i) = a(jj,4)

c          if(typctl.gt.5) then
c   	     print *,jj,' tt ',tt(i),a(jj,1),a(jj,2),a(jj,3),
c    +	             a(jj,4),dat(jj),used(i)
c          endif

	 else if (dabs(dtt).le.dtm+dtdt .and. touse(i)(4:4).eq.'D')then

c
c          phase can later eventually be used for a travel-time-
c          difference observation.
c

           used(i)(1:1) = 't'

	   datla(i) = dpa*deldla
	   datlo(i) = dpa*deldlo
	   datho(i) = dble(dtdh(j))

	 endif

	 if(touse(i)(3:3).eq.'S') then

	   ddpa = p(i) - dabs(dpa)
	   if(dabs(ddpa).lt.dpam0) then
	     jj  = jj + 1
	     jpa = jpa + 1
	     dat(jj)  = ddpa
	     dats(jj) = ps(i)
	     a(jj,1) = 0.d0
	     a(jj,2) = dddp(j)*deldla
	     a(jj,3) = dddp(j)*deldlo

	     partabl = dddp(j)

	     a(jj,4) = dpdh(j)
	     if(dabs(a(jj,4)).lt.1.d-5) a(jj,4)=0.d0

             used(i)(3:3) = 'S'

c            if(typctl.gt.5) then
c  	        print *,jj,' p ',p(i),a(jj,1),a(jj,2),a(jj,3),
c    +	               a(jj,4),dat(jj),used(i)
c             endif

	   endif

	 endif

	 go to 298

       endif

297   continue

c
c     Try it with another phase-name from the same phase-type.
c

      if(single) go to 298

      call testphase (phid,icha)

      if(icha.eq.999) go to 298

      go to 295

298   if(touse(i)(2:2).eq.'A') then

	 if(used(i)(1:1).eq.'T' .and. (dpa.lt.0.d0 
     +	                      .or.phase(i)(1:4).eq.'P3KP') ) then
	   ddazi = alpha1(azi(i) - alpha2(baz(iev(i))-180.d0))
	 else
	   ddazi = alpha1(azi(i) - baz(iev(i)))
	 endif

	 if(dabs(ddazi).lt.dazim0) then
	   jj = jj + 1
	   jazi = jazi + 1
 	   dat(jj)  = ddazi
 	   dats(jj) = azis(i)
 	   a(jj,1) = 0.d0
 	   a(jj,2) = dazidla
 	   a(jj,3) = dazidlo
 	   a(jj,4) = 0.d0
           used(i)(2:2) = 'A'

c          if(typctl.gt.5) then
c   	      print *,jj,' baz ',azi(i),a(jj,1),a(jj,2),a(jj,3),
c    +	             a(jj,4),dat(jj),used(i)
c          endif
	 endif

       endif

300   continue

      if(single) go to 302

      if(jazi .gt. jdt .and. jdt.le.3 .and. ifixaz.le.5) then
	 ifixaz = ifixaz + 1
	 dtmp   = dtmp * 1.5d0
	 dtms   = dtms * 1.5d0
	 dtm0   = dtmp + dtms
	 dazim0 = dazim0 * 1.5d0
	 dpam0  = dpam0  * 1.5d0
	 go to 101
      endif

      if (nrms1.ge.2 .and. jj.ge.3) then

	 rmso = dsqrt(rms1/dble(nrms1))

	 if(ilastiter.eq.1 .and. iremo.lt.2 .and. 
     +      nrms1.gt.2 .and. jj.gt.3) then

	    rms0 = 1.5d0 * rmso / sdmeans

 	    if(datmax.ge.rms0) then
	       datmax0 = rmso*2.d0
	       iremo = iremo+1
	       go to 101
            endif

 	 endif

	 iteraz = 0 


      else 
	 if(rmso.le.50.d0) then
	    rmso = rmso*10.d0
         else 
	    rmso = 9999.d0
         endif
	 iteraz = iteraz + 1
      endif

      if(.not.diffflag) go to 302

c
c     Add possible travel-time difference(s) as additional 
c     condition(s) to the equation system to be solved.
c

c
c     Travel-time differences can only be used in the case that we 
c     have more than 2 travel-time observations at one station.
c
      
      if(jdt.le.2) go to 302

      ndt = 0

      do 301 i = 1,nobs-1

      do 301 j = i+1,nobs

         if(sta(iev(i)).ne.sta(iev(j))) go to 301

         if(phase(i).eq.phase(j)) go to 301

         if(touse(i)(4:4).ne.'D' .or.  touse(j)(4:4).ne.'D' ) go to 301

         if((used(i)(1:1).ne.'T' .and. used(i)(1:1).ne.'t') ) go to 301

         if((used(j)(1:1).ne.'T' .and. used(j)(1:1).ne.'t') ) go to 301


	 dtt = (tt(j) - ttt(j)) - (tt(i) - ttt(i))

 	 if(dabs(dtt).le.dtm0 .and. dabs(dtt).gt.0.d0) then

	   jj = jj + 1

	   ndt = ndt + 1

	   idtu(ndt) = j+i

	   dat(jj) = dtt
	   dats(jj)= dsqrt(q2(tts(i)) + q2(tts(j)))

	   a(jj,1) = 0.d0
           a(jj,2) = datla(j) - datla(i)
           a(jj,3) = datlo(j) - datlo(i)

	   a(jj,4) = datho(j) - datho(i)
	   if(dabs(a(jj,4)).lt.1.d-5) a(jj,4)=0.d0

           used(i)(4:4) = 'D'
           used(j)(4:4) = 'D'

c          if(typctl.gt.5) then
c   	      print *,jj,' dt ',a(jj,1),a(jj,2),a(jj,3),
c    +	             a(jj,4),dat(jj),used(i),used(j)
c          endif

 	 endif

301   continue

c
c     Everything is ready for a 'final' inversion
c
c     hyposat_gmi will do it
c

302   if(czo.eq.'D') then

        if(nobs.lt.4 .or. jdt.lt.4) then
           if(zo/(deg2rad*radloc(stala(istatmin))).lt.0.25d0 .and. 
     +        iter.gt.2) then
             if( zoflag ) go to 9998
             zo = 0.1d0
             czo = 'B'
	     if(typctl.gt.0) then
	        print *,'No depth resolution, fixed at',zo,' [km]'
	     endif

	     if(output) then
	        write(11,'(/''No resolution for depth, fixed at''
     +	              ,f7.2)') zo
	     endif

             go to 101
           endif
        endif

        im = 4

      else if (czo.eq.'F' .or. czo.eq.'B') then
        im = 3
      endif

      in = jj

      if(in.gt.mread2 .and. .not.single) then

	 print *, 'Inversion matrix: ',in,' (data) ',mread2,
     +            ' (wrong dimension of Jacobian)'
	 go to 9999

      endif

      if(in.le.1 .and. .not.single) then

	 if(in0sw .ge. 6) then
	    print*,'Inversion failed'
	    print*,'No data fit within limits of start model!'
	    go to 9999
	 endif

         if(plflag .and. in0sw.eq.5) then

           call plane(stala,stalo,iev,tt,nobs,azim,dazi,
     +                ray,dray,phipl,touse,phase,jref,typctl)

           if(jref.gt.0 .and. dazi.lt.90.d0 .and. dray.lt.4.d0) then

	     phase_t = ' '
	     itray = 1
	     rayok = .false.
3021         rayok = tauget_ray(phipl,phase_t,ray,modn,zo,
     +	                ddel,ttray)

             if(rayok) then

               inddel = 1
               call delazd(stala(iev(jref)),stalo(iev(jref)),azim,
     +                     ddel,inddel,elatmg,elonm)
               elatm = convlat(elatmg,1)

               if(typctl.gt.0) then
                  print *,'Epicenter set from station ',
     +              sta(iev(jref)),'after plane wave fit: backazimuth',
     +              azim,' deg, delta',ddel,' deg'
               endif

               sdlatg = 45.d0/ray
               sdlon  = 90.d0/ray

	       tome = tt(jref) + timemin - ttray

               if(output) then
                  write(11,'(''Epicenter set from station '',a8,
     +              ''after plane wave fit: backazimuth'',f7.2,
     +              '' deg, delta'',f7.2,'' deg'')') sta(iev(jref)),
     +              azim,ddel
                  write(11,'(''Epicenter lat:'',f9.4,'' [deg]'')')  
     +  	    elatmg,sdlatg
		  write(11,'(''Epicenter lon:'',f9.4,'' [deg]''/)') 
     + 	            elonm,sdlon
		  write(11,'(''Source time set to: '',f15.2)') tome
               endif

               in0sw = in0sw + 1
               go to 101

	     else
	       
	       if(itray.lt.2) then

		  itray = itray + 1
		  ray = ray - dray
		  go to 3021

	       endif

             endif

           endif

         endif

         in0sw = in0sw + 1

	 dtmp   = dtmp * 2.d0
	 dtms   = dtms * 2.d0
	 dtm0   = dtms
	 dazim0 = dazim0 * 2.0d0
	 dpam0  = dpam0 * 2.0d0

	 if(typctl.ge.4) then
	    print *, 'Inversion matrix error: in = ',in,' (too less'
	    print *, 'data to invert!). Time boundaries changed'
	 endif

	 go to 101

      endif

c
c     special case: only one, single array observation
c
      if(single .and. in.eq.3) then

	dtmflag = .true.

	ddel = dat(2) / partabl

 	if((p(1).lt.pdif .or. p(1).gt.9.d0) .and. 
     +                   dabs(ddel).gt.1.d0    )    ddel=ddel/2.d0
 	if(p(1).gt.9.d0 .and. ddel.gt.1.d0    )     ddel=1.d0
 	if(p(1).gt.9.d0 .and. ddel.lt.-1.d0    )    ddel=-1.d0

	deln = del(1) + ddel 

	inddel = 1

        call delazd(stala(1),stalo(1),azi(1),deln,inddel,elat1,elon1)

        elatm1 = convlat(elat1,1)

	r(1) = dat(1)
	r(2) = elatm1 - elatm
	r(3) = elon1 - elonm
	r(4) = 0.d0

	var(1) = dats(1)

	dvar = dabs(dats(2) / partabl)

	ddel = ddel + dvar
	call delazd(stala(1),stalo(1),azi(1),ddel,inddel,elat1,elon1)
	ddel = ddel - 2.d0*dvar
	call delazd(stala(1),stalo(1),azi(1),ddel,inddel,elat2,elon2)
	var1 = (elat2-elat1)/2.d0
	var2 = (elon2-elon1)/2.d0

	ddel = ddel + dvar
	aziv1 = azi(1) + dats(3)
	call delazd(stala(1),stalo(1),aziv1,ddel,inddel,elat1,elon1)
	aziv1 = azi(1) - 2.d0*dats(3)
	call delazd(stala(1),stalo(1),aziv1,ddel,inddel,elat2,elon2)

	var(2) = dsqrt(q2(elat2-elat1)/4.d0 + var1*var1)
	var(3) = dsqrt(q2(elon2-elon1)/4.d0 + var2*var2)

	var(4) = 0.d0

	res(1) = 0.d0
	res(2) = 0.d0
	res(3) = 0.d0
	res(4) = 0.d0

	go to 304

      else if (single .and. in.lt.3) then

           if(insar.le.15) then

	      ddel = deln

	      if(insar.le.0) then 

		 phase(1) = 'P'
		 ddel = 50.d0

              else

                if(p(1).le.pdif) then

		   if(insar.ge.7) go to 3028

	           ddel = 148.d0

	           if(phid(1:2).eq.'P ')   phase(1) = 'PKPab'
	           if(phid(1:3).eq.'Pdif') phase(1) = 'PKPab'
	           if(phid.eq.'PKPab')     phase(1) = 'PKPbc'
	           if(phid.eq.'PKPdif')    phase(1) = 'PKPbc'
	           if(phid.eq.'PKPbc')     phase(1) = 'PKPdf'
	           if(phid.eq.'PKPdf')     phase(1) = 'PKiKP'
              
                else if(p(1).gt.9.d0) then

		   if(insar.ge.4) go to 3028

	           if(phid(1:2).eq.'P ' .or. phid(1:2).eq.'P1')  then
		      phase(1) = 'Pn'
		      ddel = 10.d0
	           else if(phid.eq.'Pn') then
		      phase(1) = 'Pb'
		      ddel = 2.d0
	           else if(phid.eq.'Pb') then
		      phase(1) = 'Pg'
		      ddel = 1.d0
	           endif

	        endif

	      endif

	      insar = insar + 1
	      inddel = 1
	      call delazd(stala(1),stalo(1),azi(1),ddel,inddel,
     +	      elatmg,elonm)
              elatm = convlat(elatmg,1)

	      go to 100

	   endif

3028	   print*,'Single phase case, but no inversion possible'
	   if(output) then
	      write(11,*)'Single phase case, but no inversion ',
     +	                 'possible'
	   endif

	   if(p(1).ge. pmoh) then
	      print*, 'due to missing direct',
     +	             ' crustal phase matching theoretical ray',
     +               ' parameter in chosen model.'
	      if(output) then
	         write(11,*) 'due to missing direct',
     +	             ' crustal phase matching theoretical ray',
     +               ' parameter in chosen model.'
	      endif
	   else if(p(1).gt. 10.1d0 .and. p(1).lt. pmoh) then
	      print*, 'due to missing direct upper mantle',
     +	             ' crustal phase matching theoretical ray',
     +               ' parameter in chosen model.'
	      if(output) then
	         write(11,*) 'due to missing direct upper mantle',
     +	             ' crustal phase matching theoretical ray',
     +               ' parameter in chosen model.'
	      endif
	   else
	      print *, 'due to core phase travel-time-curve ',
     +	             'triplication!'
	      if(output) then
	         write(11,*) 'due to core phase travel-time-curve ',
     +	             'triplication!'
	      endif
	   endif
	   go to 9999
	   
      endif

      in0sw = 0

      if (iellipi.eq.1 .and. .not.iellip) iellip=.true.

303   call hyposat_gmi(in,im,nq,ierr,typctlm)

      if(ierr.ne.0) then
	 if (ifail.eq.0 .and. is(1).eq.3) then
	    ifail = 1
	    quot = quot*2.d0
	    go to 303
	 endif
	 if(output) then
	    write(11,*) 'GMI failed: no new solution!'
	 endif
	 print*,'GMI failed: no new solution!'
	 go to 9999
      endif

304   if(typctl.gt.4) then
	  print *, 'in,im,nq',in,im,nq
	  print *,r(1),var(1)
	  print *,r(2),var(2)
	  print *,r(3),var(3)
	  print *,r(4),var(4)
	  do 305 j=1,jj
	  f1 = dat(j)-res(j)
305       print *,j,dat(j),res(j),f1
      endif

c
c     estimating the new hypocenter
c
         
c
c     the new source depth
c
      if(im.eq.4) then
        ar4 = dabs(r(4))
        if(ar4.gt.200.d0 ) then
	   r40 = r(4)*0.1875d0
        else if(ar4.le.200.d0 .and. ar4.gt.100.d0 ) then
           r40 = r(4)*0.375d0
        else if(ar4.le.100.d0 .and. ar4.gt.40.d0) then
	   r40 = r(4)*0.75d0
        else
	   r40 = r(4)
        endif

	zo = zo + r40*dchang

        if(var(4).gt.1.d-5) then
	   rs(4) = var(4)
	else
	   rs(4) = sdzo
	endif

        if(zo.lt.0.1d0) then
	   r40 = 0.1d0 - zo
	   rs(4) = dsqrt(rs(4)*rs(4) + q2(r(4)-r40))
	   zo = 0.1d0
	endif

	if(zo.ge.800.d0) then
	   r40 = 799.d0 - zo
	   rs(4) = dsqrt(rs(4)*rs(4) + q2(r(4)-r40))
	   zo = 799.d0
	endif

        if(rs(4).gt.250.d0) rs(4) = 250.d0

      endif

c
c     the new source time
c
      rsi = 1.d0
      ar1 = dabs(r(1))
      if(r(1).lt.0.d0) rsi = -1.d0
      if(ar1.le.60.d0 .or. single) then
         tome = tome   + r(1) * dchang
      else if(ar1.le.120.d0.and.ar1.gt.60.d0) then
	 tome = tome + r(1)/2.d0
      else
	 tome = tome + rsi*120.d0
      endif

      tom    = tome - timemin

      if(var(1).gt. 1.d-3) rs(1) = var(1)
      if(rs(1).gt.250.d0) rs(1) = 250.d0

c
c     save the old epicenter solution 
c
      elatmo = elatmg
      elonmo = elonm

c
c     the new source latitude
c
      rsi = 1.d0
      ar2 = dabs(r(2))
      if(r(2).lt.0.d0) rsi = -1.d0
      if(ar2.le.5.d0) then
         elatm = elatm + r(2)*dchang
      else if(ar2.le.10.d0.and.ar2.gt.5.d0) then
         elatm = elatm + r(2)/2.d0
      else
         elatm = elatm + rsi*10.d0
      endif

      ilon = 0
      if(elatm.gt. 90.d0) then
	 elatm = 180.d0 - elatm
	 ilon  = 1
      else if(elatm.lt.-90.d0) then
	 elatm = -(elatm + 180.d0)
	 ilon  = 1
      endif

      elatmr  = deg2rad*elatm
      elatmg  = convlat(elatm,2)
      sdlatg  = var(2) /( eps*q2(dcos(elatmr))
     +                            +q2(dsin(elatmr)) )

      coelatm = 90.d0 - elatm
      coelatmr= deg2rad*coelatm

      if(var(2).gt. 1.d-5) rs(2) = var(2)
      if(rs(2) .gt.180.d0)  rs(2) = 180.d0

c
c     the new source longitude
c
      rsi = 1.d0
      ar3 = dabs(r(3))
      if(r(3).lt.0.d0) rsi = -1.d0
      if(ar3.lt.5.d0) then
         elonm   = elonm  + r(3)*dchang
      else if(ar3.le.10.d0.and.ar3.gt.5.d0) then
         elonm   = elonm  + r(3)/2.d0
      else
         elonm   = elonm  + rsi*10.d0
      endif

      elonm = alpha1(elonm)

      if(ilon.eq.1) elonm = alpha1(elonm+180.d0)

      if(var(3).gt. 1.d-5) rs(3) = var(3)
      if(rs(3).gt.180.d0) rs(3) = 180.d0

      last = .false.

      if(iter.eq.maxiter) then
	 dtmflag = .true.
	 if(zo.le.0.1d0) zo = 0.d0
      else if(iter.gt.maxiter) then
	 print *,'Location stopped, more than ',maxiter,
     +           ' iterations'
	 if(output) then
	    write(11,'(/''Location stopped, more than '',i4,
     +              '' iterations'')') maxiter
	 endif

	 go to 399
      endif

      if(iter.le.mosci) then

        dtos(iter)  = tom
        dlaos(iter) = elatm
        dloos(iter) = elonm
        dlo1os(iter)= dcos(deg2rad*elonm)
        dlo2os(iter)= dsin(deg2rad*elonm)
        dzoos(iter) = zo
        rtos(iter)  = dabs(r(1)) + var(1)
        rlaos(iter) = dabs(r(2)) + var(2)
        rloos(iter) = dabs(r(3)) + var(3)
        rzos(iter)  = dabs(r(4)) + var(4)

      endif
 
c
c     We will check if the new solution is close to the former solution.
c     If CHECK [km] is smaller than SETCHECK [km] we will stop.
c

c     The change in the horizontal components

      dk = 0.d0
      call depi (elatmg,elonm,elatmo,elonmo,del3,dk,ep2,ep1,d2km)

c     The change in source time is compensated eventually by a
c     change in depth

      dtokm = 0.d0

      if(im.eq.4) then

         dtokm = r(1)*6.d0
         if (zo.gt.20.d0) then
	    dtokm = r(1)*7.d0
         else if (zo.gt.30.d0) then
	    dtokm = r(1)*8.d0
         else if (zo.gt.260.d0) then
	    dtokm = r(1)*9.d0
         else if (zo.gt.450.d0)  then
   	    dtokm = r(1)*10.d0
         else if (zo.gt.660.d0)  then
	    dtokm = r(1)*11.d0
         endif

         if(var(4).gt.1.d-5) dtokm = dtokm - r(4)

      endif

      check = dsqrt(dtokm*dtokm + dk*dk)

      call depi (elatmg,elonm,stala(istatmin),stalo(istatmin),del3,
     +	            dk,ep2,ep1,d2km)

      ilastiter = 0

      if (check.le.setcheck .or. (check.le.disper*dk .and. 
     +      (iteraz.ge.5 .or. imaxiter.ge.5 .or. 
     +       dble(maxiter)/dble(iter).lt.1.3d0   )) ) then


	 if((dtmflag .or. .not.lastfix .or. in.le.im) .and. 
     +       rmsold/rmso.lt.1.2d0) then
	    go to 400
         else

	    if(zo.le.0.1d0) zo=0.d0

	    if(var(4).eq.0.d0 .and. czo.eq.'D') then

	       czo = 'B'
	       if(typctl.gt.0) then
	          print *,'No depth resolution, fixed at',zo,' [km]'
	       endif

	       if(output) then
	          write(11,'(/''No resolution for depth, fixed at''
     +	                ,f7.2)') zo
	       endif

	    endif

	    setcheck2 = check*1.1d0

	    ilastiter = 1
	    dchang = dchang0

	    go to 390

         endif
      endif

      dtmflag = .false.
      direct  = .false.

      if(iter.gt.mosci) then

        if( ifixto.le.5  .and. (var(1).eq.0.d0 .or. 
     +	   (jdt.lt.in/3 .and. jdt.lt.nstat/2)) ) then

	   ifixto = ifixto + 1

	   rzo = sngl(zo)
	   rdel = sngl(del3)
	   call tauget_mod(rzo,rdel,nphas,phcd,ttc,dtdd,
     +	                   dtdh,modnam)

	   tom  = - dble(ttc(1))
	   tome = timemin + tom

	   rs(1) = rs(1) + (-tom)

 	   dtp   = dtp  * 1.2d0
 	   dts   = dts  * 1.2d0
	   dtm2  = dtp + dts

	   setcheck = setcheck*1.2d0
	   setcheck2= setcheck*10.d0

	   go to 101

        endif

        nextiter1 = 0

c
c       check for oscillation in the solutions
c

	mosci2 = 0

        do 370 i = mosci,1,-1

        if( dabs(dzoos(i)-zo).lt.rminh            .and.
     +	    dabs(dtos(i)-tom).lt.rmint            .and.
     +      dabs(dlaos(i)-elatm).lt.rming         .and.
     +      dabs(alpha1(dloos(i)-elonm)).lt.rming ) then
	    mosci2 = i
	    go to 371
        endif

370     continue

371     if(mosci2.ne.0 .and. mosci2.ne.mosci) then
c
c       we have to calculate a new starting solution from all
c       oscillating solutions!
c

            nextiter = nextiter + 1
            if(nextiter.gt.8) then
    	      if(output) then
   	         write(11,'(/''Oscillating solution: after'',i3,
     +	                '' iterations stopped!'')') iter
	      endif
	      go to 399
            endif

            nextiter1 = 1

	    if (czo.eq.'D') then

               zo  = dmean(dzoos,mosci,mosci2) 
	       rs(4) = ddmax(rzos,mosci,mosci2)
	       if(rs(4).le.1.d-5) rs(4) = sdzo

	       if(zo.le.0.1d0 .or. nextiter.gt.mosci*2/3) then
		 if(iterz.eq.1 .and. zoflag) go to 9998
	         iterz = iterz + 1
		 czo = 'B'
		 var(4) = 0.d0
		 rs(4)  = 1.d0
	         if(typctl.gt.0) then
	            print *,'No depth resolution, fixed at',zo,' [km]'
	         endif

	         call findosci(zmin,zmax1,dzoos,mosci,mosci2,1)

		 if(output) then

	            write(11,'(/''No resolution for depth, '',
     +                          ''oscillating between:'',f6.1,
     +                      '' and'',f6.1,'' [km]''/
     +               ''depth fixed at:'',f6.1,'' km'')') zmin,zmax1,zo
		 endif

               endif

            endif

	    tom = dmean(dtos,mosci,mosci2)
	    tome = tom  + timemin

	    rs(1) = ddmax(rtos,mosci,mosci2)
	    if(rs(1).lt.1.d-3) rs(1) = sdto

            elatm = dmean(dlaos,mosci,mosci2)
	    elatmr= deg2rad*elatm
            elatmg  = convlat(elatm,2)
            coelatm = 90.d0 - elatm
            coelatmr= deg2rad*coelatm
	    rs(2) = ddmax(rlaos,mosci,mosci2)
	    if(rs(2).lt.1.d-5) rs(2) = sdlat

	    elonm1= dmean(dlo1os,mosci,mosci2)
	    elonm2= dmean(dlo2os,mosci,mosci2)
	    elonm = rad2deg*datan2(elonm2,elonm1)

	    rs(3) = ddmax(rloos,mosci,mosci2)
	    if(rs(3).lt.1.d-5) rs(3) = sdlon

            if (nrms1.gt.5) then

 	       dtp = dtp  * 0.9d0
 	       dtp = dmin1(dtp,rmso*2.d0)
	       if(dtp.lt.3.d0) dtp = 3.d0


 	       dts = dts  * 0.9d0
 	       dts = dmin1(dts,rmso*4.d0)
	       if(dts.lt.6.d0) dts = 6.d0

	       dtm2  = dtp + dts

	    endif

	    if(nobs.gt.1) then

               dazim1 = dazim1*0.9d0
               if(dazim1.lt.3.d0) dazim1 = 3.d0

               dpam1  = dpam1 *0.9d0
               if(dpam1.lt.1.5d0) dpam1 = 1.0d0

	    endif

            rminh = rminh*1.2d0
            rming = rming*1.2d0
            rmint = rmint*1.2d0

	    disper = disper*1.1d0
	    if(disper.gt.0.05d0) disper=0.05d0

	    setcheck = setcheck*1.5d0
	    setcheck2= setcheck*10.d0

 	    dchang = dchang * 0.8d0

	    direct = .true.

        endif
	 
        if(var(4).gt.4.d0*zo .and. czo.eq.'D' .and. var(4).gt.30.d0 
     +     .and. iter.gt.5 ) then

 	   if(zoflag .and. iterz.ge.6) go to 9998

	   czo = 'B'
	   rs(4) = 1.d0
	   var(4) = 0.d0
           zo  = dmean(dzoos,mosci,mosci2)

	   iterz = iterz + 1

	   if(typctl.gt.0) then
	      print *,'Bad resolution for depth; depth fixed at',zo
	   endif

	   if(output) then
	      write(11,'(/''Bad resolution for depth; depth fixed'')')
	   endif

	   direct = .true.

	   go to 380

        endif

        if(czo.eq.'D' .and. nextiter1.eq.0 .and. iterz.le.3) then

c
c         check for oscillation in the solutions for the focal depth
c
	  mosci3 = 0

          do 377 i = mosci,1,-1

          if(dabs(dzoos(i)-zo).lt.rminh) then
	    mosci3 = i
	    go to 378
          endif

377       continue

378       if((mosci3.ne.0 .and. mosci3.ne.mosci) .or. 
     +       (mosci3.eq.mosci .and. zo.eq.0.1d0)           ) then
c
c         we have to calculate a new startvalue for the depth and fix it 
c

            rs(4)  = 1.d0
	    var(4) = 0.d0
            zo  = dmean(dzoos,mosci,mosci3)
	    czo = 'B'

	    iterz = iterz + 1

            if(typctl.gt.0) then
              print *,'Oscillating solution, ',
     +                'depth fixed at: ',zo,' km'
            endif

	    call findosci(zmin,zmax1,dzoos,mosci,mosci3,1)

	    if(output) then
               write(11,'(/''Oscillating solution between:'',f6.1,
     +                 '' and'',f6.1,'' [km]''/
     +             ''depth fixed at:'',f6.1,'' km'')') zmin,zmax1,zo
	    endif

	    direct = .true.

          endif

        endif

380     continue
	
	do 381 i = 1,mosci-1

	i2 = i + 1

        dtos(i)  = dtos(i2)
        dlaos(i) = dlaos(i2)
        dloos(i) = dloos(i2)
        dlo1os(i)= dlo1os(i2)
        dlo2os(i)= dlo2os(i2)
        dzoos(i) = dzoos(i2)
        rtos(i)  = rtos(i2)
        rlaos(i) = rlaos(i2)
        rloos(i) = rloos(i2)
        rzos(i)  = rzos(i2)

381     continue

        dtos(mosci)  = tom
        dlaos(mosci) = elatm
        dloos(mosci) = elonm
        dlo1os(mosci)= dcos(deg2rad*elonm)
        dlo2os(mosci)= dsin(deg2rad*elonm)
        dzoos(mosci) = zo
        rtos(mosci)  = dabs(r(1)) + var(1)
        rlaos(mosci) = dabs(r(2)) + var(2)
        rloos(mosci) = dabs(r(3)) + var(3)
        rzos(mosci)  = dabs(r(4)) + var(4)

	if(direct) go to 100

      endif

390   if(typctl.ge.4) then
         print*,'Iteration: ',iter,'   # of def.: ',in
         print*,'New source time  : ',tome,' +/- ',var(1)
         print*,'New epicenter lat: ',elatmg,' +/- ',sdlatg
         print*,'New epicenter lon: ',elonm,' +/- ',var(3)
         print*,'New source depth : ',zo,' +/- ',var(4)
      endif

c

      if(iter.gt.nint(maxiter*0.75) .and. imaxiter.lt.5 .and.
     +   ilastiter.eq.0) then

c
c     If we were coming close to the end of all iterations,
c     let's try it with a mean solution of the last 
c     4 solutions as new start solution.
c

         imaxiter = imaxiter + 1
	 maxiter  = maxiter + nint(maxiter*0.25/imaxiter)

         if (czo.eq.'D') then
 
            zo  = dmean(dzoos,mosci,1)
            rs(4) = ddmax(rzos,mosci,1)
            if(rs(4).le.1.d-5) rs(4) = sdzo

         endif
 
         tom = dmean(dtos,mosci,1)
	 tome = tom  + timemin
         rs(1) = ddmax(rtos,mosci,1)
         if(rs(1).lt.1.d-3) rs(1) = sdto
 
         elatm = dmean(dlaos,mosci,1)
         elatmr= deg2rad*elatm
	 elatmg  = convlat(elatm,2)
         coelatm = 90.d0 - elatm
         coelatmr= deg2rad*coelatm
         rs(2) = ddmax(rlaos,mosci,1)
         if(rs(2).eq.1.d-5) rs(2) = sdlat
 
	 elonm1= dmean(dlo1os,mosci,1)
	 elonm2= dmean(dlo2os,mosci,1)
	 elonm = rad2deg*datan2(elonm2,elonm1)

         rs(3) = ddmax(rloos,mosci,1)
         if(rs(3).eq.1.d-5) rs(3) = sdlon

         dazim1 = dazim1*2.d0
         if(dazim1.gt.90.d0) dazim1 = 90.d0

	 dpam1  = dpam1 *2.d0
	 if(dpam1.gt.15.d0) dpam1 = 15.d0

 	 dtp   = dtp * 2.d0
 	 dts   = dts * 2.d0
	 dtm2  = dtp + dts

	 setcheck = setcheck*1.5d0
	 setcheck2= setcheck*15.d0
	 disper   = disper*2.d0
	 if(disper.gt.0.05d0) disper=0.05d0

         rminh = rminh*1.5d0
         rming = rming*1.5d0
         rmint = rmint*1.5d0

      endif

      go to 100

399   continue

      last = .true.

      call findosci(tomin,tomax,dtos,mosci,mosci2,1)
      call findosci(dlamin,dlamax,dlaos,mosci,mosci2,1)
      call findosci(dlomin,dlomax,dloos,mosci,mosci2,2)
      call findosci(zmin,zmax1,dzoos,mosci,mosci2,1)

      if(output) then
         write(11,'(''Rel. source time between'',f8.2,'' and'',
     +           f8.2,'' [s]'')') tomin,tomax
      endif

      flamin = convlat(dlamin,2)
      flamax = convlat(dlamax,2)
      if(output) then
         write(11,'(''Latitude         between'',f8.2,'' and'',
     +  	 f8.2,'' [deg]'')') flamin,flamax

         write(11,'(''Longitude        between'',f8.2,'' and'',
     +  	 f8.2,'' [deg]'')') dlomin,dlomax

         if(czo.eq.'d') then
            write(11,'(''Depth            between'',f7.1,''  and'',
     +              f7.1,''  [km]''/)') zmin,zmax1
         endif

   	 write(11,'(//''Following the last (must not be the '',
     +                  ''best!) solution:''/)') 

      endif

400   continue
      if(output) then
         write(11,'(/''Iterations        :'',i5)') iter
      endif

      ibad = 0

      if(rmso .gt. 50.d0) ibad = ibad + 1

      if(output) then
         write(11,'(''Number of defining:'',i5)') in
         if(iloc) then
	    if(mod2flag .or. mod3flag .or. mod4flag) then
                write(11,'(''First reference models  : '',a,'' and '',
     +                a)') filloc(1:trimle(filloc)),modnam
	        if(mod2flag) then
	           write(11,'(''Second reference model  : '',a)') 
     +                  modnam2
		endif
	        if(mod3flag) then
	           write(11,'(''Third reference model   : '',a)') 
     +                  modnam3
		endif
	        if(mod4flag) then
	           write(11,'(''Fourth reference model  : '',a)') 
     +                  modnam4
		endif
	    else
                write(11,'(''Reference models  : '',a,'' and '',a)') 
     +                filloc(1:trimle(filloc)),modnam
	    endif
         else
	    if(mod2flag .or.  mod3flag .or. mod4flag) then
                write(11,'(''First reference model   : '',a)') modnam
	        if(mod2flag) then
	           write(11,'(''Second reference model  : '',a)') 
     +                  modnam2
		endif
	        if(mod3flag) then
	           write(11,'(''Third reference model   : '',a)') 
     +                  modnam3
		endif
	        if(mod4flag) then
	           write(11,'(''Fourth reference model  : '',a)') 
     +                  modnam4
		endif
	    else
                write(11,'(''Reference model   : '',a)') modnam
            endif
         endif
      endif

      call fetoh(tome,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

      if (fchi.ne.1.d0) then
	 sdlatg = sdlatg * fchi
	 var(3) = var(3) * fchi
	 var(4) = var(4) * fchi
	 var(1) = var(1) * fchi
      endif

      if(output) then

         write(11,'(/''The new source parameters''/)')

	 write(11,'(''Confidence level of given uncertainties:'',
     +         f7.2,'' %''/)') confl

         write(11,'(''Source time  :'',i5,4i3.2,f7.3,'' +/- '',
     +           f8.3,'' [s]'')')  yy,mon,dd,hh,mi,sec,var(1)
         write(11,'(''        or'',14x,f14.3,'' +/- '',f8.3,
     +           '' [s]'')') tome,var(1)
      endif

      if(var(1).gt.zo/6.d0 .and. czo.eq.'D') ibad = ibad + 1

      if(sdlatg.lt.45.d0) then
	 if(output) then
            write(11,'(''Epicenter lat:'',14x,f10.4,'' +/- '',f8.4,
     +              '' [deg]'')')  elatmg,sdlatg
	 endif
      else
	 ibad = ibad + 1
	 if(output) then
	    write(11,'(''Epicenter lat:'',14x,f10.4,'' +/- '',f8.4,
     +              '' (!!) [deg]'')') elatmg,sdlatg
	 endif
      endif
      if(var(3).lt.90.d0) then
	 if(output) then
            write(11,'(''Epicenter lon:'',14x,f10.4,'' +/- '',f8.4,
     +              '' [deg]'')')  elonm,var(3)
	 endif
      else
	 ibad = ibad + 1
	 if(output) then
            write(11,'(''Epicenter lon:'',14x,f10.4,'' +/- '',f8.4,
     +              '' (!!) [deg]'')') elonm,var(3)
	 endif
      endif
      if(czo.eq.'D') then

	 if(output) then
            write(11,'(''Source depth :'',15x,f7.2,''   +/- '',f6.2,
     +            ''   [km]''/)') zo,var(4)
	 endif

	 if(var(4).gt.zo*5.d0) ibad = ibad + 1

      else if(czo.eq.'F' .or. czo.eq.'B') then
	 if(output) then
            write(11,'(''Source depth :'',15x,f7.2,''   [km] Fixed''
     +	          /)') zo
	 endif
      endif

      if(nq.lt.3 .and.single) iellip =.false.

      if (iellip .and. .not.last) then

	 call ellcal(elatmg,ax1,ax2,fchi,elmax,elmin,eazi,earea)

         if(output) then

	    if(earea.lt.10000.d0) then
	      write(11,'( ''Epicenter error ellipse:''/ 
     +           ''Major axes: '',f8.2,'' [km]  Minor axes: '',
     +           f8.2,'' [km]'',/''Azimuth:'',f11.1,''  [deg] '',
     +           ''Area: '',f14.2,'' [km**2]''/)')    
     +           elmax,elmin,eazi,earea
	    else
	      write(11,'( ''Epicenter error ellipse:''/ 
     +           ''Major axes: '',f8.2,'' [km]  Minor axes: '',
     +           f8.2,'' [km]'',/''Azimuth:'',f11.1,''  [deg] '',
     +           ''Area: '',2p e14.4,'' [km**2]''/)')    
     +           elmax,elmin,eazi,earea
	    endif

         endif

      else if (iellipi.eq.1 .and. .not.iellip) then

         if(output) then

	    write(11,'( ''Epicenter error ellipse calculation '',
     +         ''not possible (too less parameter resolution)''/)')

         endif

      endif

      if(ibad.gt.2 .and. .not.last) go to 405

      rlat = sngl(elatmg)
      rlon = sngl(elonm)
      call hyposat_geo( rlat,rlon, isreg, regnum, region , ierr )
      
      if(output) then
         write(11,'(''Flinn-Engdahl Region ('',i4,'' ): '',a/)')  
     +            isreg, region(1:trimle(region))
      endif
c
c     let us now calculate the final residuals and print them out
c
405   stmean    = 0.d0
      strmean   = 0.d0
      samean    = 0.d0
      sarmean   = 0.d0
      rmsazi    = 0.d0
      spmean    = 0.d0
      sprmean   = 0.d0
      rmsp      = 0.d0
      rms       = 0.d0
      rmsisc    = 0.d0
      wisc      = 0.d0

c
c     misfit parameters for all input data!
c
      tmisf     = 0.d0
      tmisfl    = 0.d0
      ntmisf    = 0
      dmisf     = 0.d0
      dmisfl    = 0.d0
      ndmisf    = 0
      amisf     = 0.d0
      amisfl    = 0.d0
      namisf    = 0
      pmisf     = 0.d0
      pmisfl    = 0.d0
      npmisf    = 0
      wmisf     = 0.d0
      wmisfl    = 0.d0
      nwmisf    = 0

      nobst     = 0
      nobsa     = 0
      nobsp     = 0
      stato     = '     '

      if(magflag) then
         imsm     = 0
         dmsm     = 0.d0
         imbm     = 0
         dmbm     = 0.d0
      endif 

      do 450 i = 1,nobs

      epiaz(i) = -999.

      if (sta(iev(i)).ne.stato) then

	 stato = sta(iev(i))

         call depi(stala(iev(i)),stalo(iev(i)),elatmg,elonm,
     +	      del(iev(i)),dk,azie(iev(i)),baz(iev(i)),d2km)
	 rzo = sngl(zo)
	 rdel = sngl(del(iev(i)))

         fla1 = sngl(90.d0-elatm)
         fla2 = sngl(90.d0-stalae(iev(i)))
	 fla3 = sngl(deg2rad*(90.d0-elatm))
	 razi = sngl(azie(iev(i)))
	 rbaz = sngl(baz(iev(i)))

	 nphas = 0

	 imod2 = 1
	 if (mod2flag .and. touse(i)(6:6).eq.'2') then
	   modn = modnam2
	 else if (mod3flag .and. touse(i)(6:6).eq.'3') then
	   modn = modnam3
	 else if (mod4flag .and. touse(i)(6:6).eq.'4') then
	   modn = modnam4
	 else
	   if(iloc) imod2 = 0
	   modn = modnam
	 endif

         if(imod2.ne.0 .or. rdel.gt.rmax .or. zo.gt.zmax) then

	   call tauget_mod(rzo,rdel,nphas,phcd,ttc,dtdd,
     +	                 dtdh,modn)

         else

	   ierr = 0
	   indph = istaph(iev(i))*10000 + indph0
	   elatc = elatmg
	   elonc = elonm

	   elat2 = stala(iev(i))
	   elon2 = stalo(iev(i))

	   call ttloc(rzo,rdel,czo,nphas,ttc,dtdd,dtdh,dpdh,dddp,
     +                phcd,rmax,typctl,ierr,indph)

         endif

      endif

      phid = phase(i)
      ttobs  = timemin + tt(i)

      dpa   = 0.d0
      phid1 = phase(i)
      ttres  = -9999.d0
      pares  = -999.d0

      call fetoh(ttobs,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

      if(phid(1:1).eq.'T') go to 430

      surf = .false.

      if(phid.eq.'Rg') then
         if(dk.le.400.d0) then
            surf = .true.
            vsurf = vrg
         else
            phid = 'LR'
         endif
      endif

      if(phid.eq.'Lg' .and. dk.le.3000.d0) then
	 surf = .true.
	 vsurf = vlg
      endif
		
      if(phid.eq.'LR') then
	 surf = .true.
	 vsurf = vlr
      endif
		
      if(phid.eq.'LQ') then
	 surf = .true.
	 vsurf = vlq
      endif

      if (surf) then
	 nphas       = nphas + 1
	 ttc(nphas)  = sngl(dk/vsurf)
	 dtdd(nphas) = sngl(d2km/vsurf)
	 dtdh(nphas) = 0.
	 phcd(nphas) = phid
	 dddp(nphas) = 0.d0
	 dpdh(nphas) = 0.d0
      endif
		
      icha = 0
      imin = 0

      if((phid(2:3).eq.'n ' .or. phid(2:3).eq.'g ' .or. phid(2:3).eq.
     +   'b ') .and. rdel.gt.30. ) phid(2:3)='  '

410   continue

      if(phase(i).eq.'P1') then
	 nphass = 1
	 if (icha.ne.0) then
	    if (phid(1:1).eq.'P')   nphass = 2
            if (phid(1:2).eq.'PK' .or. phid(1:3).eq.'Pdi') nphass = 5
	 endif
      else
         nphass = nphas
      endif

      dpa2= 0.d0

      do 420 j = 1,nphass

      dpa = 0.d0
c
c     Let's take the first phase from IASPEI-1991-type tables,
c     which fits the phase name of the onset
c
      phid1 = phcd(j)

      if(phid.eq.'P1' .and. phid1(1:1).eq.'P' .and. phid1(3:3).eq.
     +   ' ') then

        phid = phid1

	if(p(i).lt.pdif .and.phid1(1:3).eq.'Pdi') then
	   if(del(iev(i)).gt.110.d0) then
	     phid = 'PKPdf'
	     go to 420
           endif
        endif

      endif

      if(phid.eq.'S1' .and. phid1(1:1).eq.'S' .and. phid1(3:3).eq.
     +   ' ') then

        phid = phid1

      endif

      if(phid.eq.'PKPdf') then
	 
	 if(phid1(1:5).eq.'PKiKP' .and.
     +	               del(iev(i)).ge.90.d0) phid='PKiKP'

      endif

      if(icha.eq.0) then
	 
	  if(phid1(2:4).eq.'dif') then
             if(phid.eq.'P') phid='Pdif'
             if(phid.eq.'S') phid='Sdif'
          endif

      endif

      if(phid1.eq.phid .or. imin.eq.1) then

c
c     checking : any ellipticity correction for this phase?
c

         phid1 = phcd(j)
	 pa    = dtdd(j)
	 dpa   = dble(pa)
	 dpaa  = dabs(dpa)

	 if(single .and. rayokf) then

	    if(dabs(p(i)-dpaa) .ge. 0.01d0) go to 420

	 endif

	 rzoe = rzo

	 ecor = 0.

	 if (.not.surf) then
c
c           Ellipticity corrections are yet no available for sources 
c           deeper than 700 km. Therefore we accept a small error
c           and set the source depth to 700. km
c
	    if(rzoe.gt.700.) rzoe=700.

            if(iaspflag) then
                call ellip(fla1,fla2,razi,rbaz,
     +                     rdel,rzoe,phid1,pa,ecor,ierr)
            else
                call ellip2(fla3,razi,rdel,rzoe,phid1,pa,ecor,ierr)
            endif

         endif

	 phase_t = phase_type(phid1)

	 th = 0.d0

         tcrust = 0.d0

	 if(vlflag .and. .not.surf) then

	   hsta = stael(iev(i))

	   vloc = 99999.d0
           if(phase_t.eq.'P') vloc = stavp(iev(i))
           if(phase_t.eq.'S') vloc = stavs(iev(i))

	   if(vloc .ne. 99999.d0) then

	      if(iabs(imo).eq.2 .or. imo.eq.4) then
	      
	         elatc  = stalae(iev(i))
	         elonc  = stalo(iev(i))
	         indr = 1
	         tcrust = crustc(phase_t,dpaa,indr,typctl)

	         if(tcrust.ne.0.d0) hsta = hsta - elev

	      endif

	      if(hsta.ne.0.d0) then

	         radkm = deg2rad*radloc(stala(iev(i)))
	         phin = vloc*dpaa/radkm

	         if(phin.lt.1.d0) then 

	            dl = hsta / dcos(dasin(phin))

                    if(dl.gt.0.d0) then
                      th = dl/vloc - 
     +		           dpaa*dsqrt(dl*dl-hsta*hsta)/radkm
		    else if(dl.lt.0.d0) then
                      th = dl/vloc + 
     +		           dpaa*dsqrt(dl*dl-hsta*hsta)/radkm
                    endif

	         endif
	      endif
	   endif
	 endif

c
c        We have eventually to correct this phase for the local 
c        structure at the reflection point at the surface (if 
c        CRUST 5.1 is available).
c
	 trefl = 0.d0
	 if( .not.surf .and. touse(i)(5:5).eq.'R' .and.
     +       (rdel.gt.rmax .or. zo.gt.zmax .or. .not.iloc) ) then 

            phase_t = phid1(1:1)
	    if(phid1(1:1).eq.'p') phase_t = 'P'
	    if(phid1(1:1).eq.'s') phase_t = 'S'
	     
	    fmult = 1.d0
	    del0  = dirdel(dpaa,zo,fmult,phase_t)
	    azi0 = azie(iev(i))

	    if(phid1(1:1).ne.'p' .and. phid1(1:1).ne.'s') then
	       if(dpa.ge.0.d0) then
	          del0 = (del(iev(i))-del0)/2.d0
               else
	          del0 = (360.d0 - del(iev(i)) - del0)/2.d0
	          azi0 = alpha2(azie(iev(i)) + 180.d0)
               endif
            endif
		
	    inddel = 1
	    call delazd(elatmg,elonm,azi0,del0,inddel,elatc,elonc)
    
c
c     correction for depth phases (e.g.: pP, sP...)
c
            if( (phid(1:1).eq.'p' .or. phid(1:1).eq.'s')
     +          .and. zo.ge.zmax) then

              if ((phid(1:1).eq.'p' .and. phid(2:2).eq.'P') .or.
     +            (phid(1:1).eq.'s' .and. phid(2:2).eq.'S')) then

                  indr = 2

              else

                  indr = 3

              endif

              trefl   = crustc(phase_t,dpaa,indr,typctl)
              used(i)(5:5) = 'R'

            endif

c
c     correction for surface reflections (e.g.: PnPn,...,PP,SS,P'P')
c
            if( phid(1:1).eq.phid(2:2) .or.
     +          phid(1:2).eq.phid(3:4)  ) then

                indr = 2
                trefl = crustc(phase_t,dpaa,indr,typctl)
                used(i)(5:5) = 'R'

            endif

c
c      correction for converted surface reflections (e.g.: PnSn,...)
c
            conr = .false.
            if( (phid(1:1).eq.'P' .or. phid(2:2).eq.'P') .and.
     +          (phid(1:1).eq.'S' .or. phid(2:2).eq.'S') .and.
     +          (phid(3:3).eq.' ' .or. phid(3:3).eq.'g' .or.
     +           phid(3:3).eq.'b' .or. phid(3:3).eq.'n')) then
                   conr=.true.
                   phidr = phid(2:)
            endif

            if( (phid(1:1).eq.'P' .or. phid(3:3).eq.'P') .and.
     +          (phid(1:1).eq.'S' .or. phid(3:3).eq.'S') .and.
     +           phid(2:2).ne.'b' .and. phid(2:2).ne.'m' .and.
     +           phid(2:2).ne.'c' .and. phid(2:2).ne.'k' .and.
     +          (phid(2:2).eq.phid(4:4) .or. phid(2:2).eq.'g' .or.
     +           phid(2:2).eq.'n'                      )) then
                   conr=.true.
                   phidr = phid(3:)
            endif

            if(conr) then

		zor = 0.d0
		rayok = tauget_ray(phidr,phase_t,dpaa,modn,zor,
     +	                del0,ttray)

                azi0 = azie(iev(i))
                if(dpa.lt.0.d0) azi0 = alpha2(azie(iev(i))-180.d0)

                inddel = 1
                call delazd(stala(iev(i)),stalo(iev(i)),azi0,del0,
     +                      inddel,elatc,elonc)

                indr = 3
                trefl   = crustc(phase_t,dpaa,indr,typctl)
                used(i)(5:5) = 'R'

            endif

            if (trefl.ne.0.d0.and.typctl.gt.6) then
               print *,'dirdel: ',phid,' azi ',azi0,' del ',del0
               print *,'  lat ',elatc,' lon ',elonc,' trefl ',trefl
            endif

         endif

         ttt(i) = tome + dble(ttc(j) + ecor) + th + tcrust + trefl

	 ttres  = ttobs - ttt(i)

	 if(ttres.gt.100.d0 .and. phase(i)(1:2).eq.'P1' .and. 
     +      phid(1:3).eq.'Pdi' .and. used(i)(1:1).eq.'T') then
	    phid='PKPdf'
            go to 420
         endif

	 if(imin.eq.0) then

c	    if(typctl.ge.4) then
c	       print *,'i,ttt,ttobs,ttres,ecor,th,tcrust,used,tts'
c	       print *,i,ttt(i),ttobs,ttres,ecor,th,tcrust,used(i),tts(i)
c	    endif

	    if(used(i)(1:1).eq.'T') then
	       stmean  = stmean + ttres
	       strmean = strmean + dabs(ttres)
	       rms     = rms    + ttres*ttres
	       rmsisc  = rmsisc + q2(ttres)*tts(i)
	       wisc    = wisc   + tts(i)
	       epiaz(i)= razi
	       nobst   = nobst + 1
            endif

	    pares = p(i)-dabs(dpa)
	    if(used(i)(3:3).eq.'S') then
	       spmean  = spmean + pares
	       sprmean = sprmean + dabs(pares)
	       rmsp    = rmsp + q2(pares)
	       epiaz(i)= razi
	       nobsp   = nobsp + 1
            endif

	    go to 430

         else if(imin.eq.1.and.dabs(ttres).lt.dabs(dtmin)) then

	    phid2 = phid1
	    dtmin = ttres
	    pares2 = p(i)-dabs(dpa)
	    dpa2  = dpa

	 endif

      endif

420   continue

      if (imin.eq.1) go to 425

c
c     Try it with another phase-name from the same phase-type.
c
 
      if(single) go to 410

      call testphase (phid,icha)

      if(icha.ne.999) go to 410

      if(imin.eq.0) then
         imin  = 1
         dtmin = 9999.d0
         go to 410
      endif

425   if(dabs(dtmin).le.15.d0 .or. used(i)(1:1).eq.'T' .or.
     +   (used(i)(1:1).eq.'t' .and. used(i)(4:4).eq.'D')) then
	 ttres  = dtmin
	 pares  = pares2
	 phid1  = phid2
	 dpa    = dpa2
      else
         ttres  = -9999.d0
         pares  = -999.d0
	 dpa    = 0.d0
         phid1  = ''
      endif

430   continue

      arr(i) = rdel
      if(used(i)(1:1).eq.'t') used(i)(1:1) = ' '

      if(dpa.lt.0.d0 .or.phase(i)(1:4).eq.'P3KP') then
	 azires = alpha1(azi(i) - alpha2(baz(iev(i))-180.d0))
      else
         azires = alpha1(azi(i) - baz(iev(i)))
      endif

      if(touse(i)(1:1).eq.'T' .and. ttres.ne.-9999.d0) then
	 fmis    = ttres/tts(i)
	 tmisfl  = tmisfl + dabs(fmis)
 	 tmisf   = wmisf + q2(fmis)
	 ntmisf  = ntmisf + 1
      endif

      if(touse(i)(3:3).eq.'S' .and. pares.ne.-999.d0) then
	 fmis    = pares/ps(i)
 	 pmisfl  = pmisfl + dabs(fmis)
 	 pmisf   = pmisf + q2(fmis)
	 npmisf  = npmisf + 1
      endif

      if(touse(i)(2:2).eq.'A') then
	 fmis    = azires/azis(i)
 	 amisfl  = amisfl + dabs(fmis)
 	 amisf   = amisf + q2(fmis)
	 namisf  = namisf + 1
      endif

      if(used(i)(2:2).eq.'A') then
         samean  = samean + azires
         sarmean = sarmean + dabs(azires)
	 rmsazi  = rmsazi + q2(azires)
	 epiaz(i)= razi
         nobsa   = nobsa + 1
      endif

      imod = ' '
      if(touse(i)(6:6).eq.'2' .and. mod2flag) imod = '2'
      if(touse(i)(6:6).eq.'3' .and. mod3flag) imod = '3'
      if(touse(i)(6:6).eq.'4' .and. mod4flag) imod = '4'

      if(phase(i).eq.phid1) then

        write(text(i),'(a5,f8.3,f7.2,1x,a8,8x,2i3.2,f7.3,f8.3,
     +                f7.2,f8.2,f6.2,f7.2,1x,a5,a1)') 
     +                stato,rdel,razi,phase(i),hh,mi,sec,
     +                ttres,azi(i),azires,p(i),pares,used(i),imod

      else

        write(text(i),'(a5,f8.3,f7.2,1x,2a8,2i3.2,f7.3,f8.3,
     +                f7.2,f8.2,f6.2,f7.2,1x,a5,a1)') 
     +                stato,rdel,razi,phase(i),phid1,hh,mi,sec,
     +                ttres,azi(i),azires,p(i),pares,used(i),imod

      endif

      if(ttres.eq.-9999.d0) text(i)(51:58)='        '
      if(ttres.lt.-999.d0)  text(i)(51:58)='        '
      if(ttres.gt.9999.d0)  text(i)(51:58)='        '
      if(azi(i).le.-999.d0) text(i)(59:73)='               '
      if(azires.lt.-360.d0) text(i)(66:73)='        '
      if(azires.gt.360.d0)  text(i)(66:73)='        '
      if(p(i).le.-999.d0)   text(i)(74:86)='             '
      if(pares.le.-999.d0)  text(i)(80:86)='       '

      if (magflag) then
	if(amplit(i).gt.0.d0 .and. period(i).gt.0.d0) then

c
c     the standard IASPEI (1967) formula:
c     Ms = log (A / T )  + 1.66 log (DELTA) + 0.3 (A in nanometer)
c
c     or the Rezapour/Pearce(BSSA 88, 43-61) formula:
c
c     Ms = log (A / T ) + 1/3 log (DELTA) + 1/2 log (sin(DELA) +
c          0.0046 DELTA + 2.370  (A in nanometer)
c

           dmag = 0.d0

	   if(phid1.eq.'LR' ) then

	      if(magtyps(1:6).eq.'IASPEI') then
	        dmag = dlog10(amplit(i)/period(i)) +
     +              dlog10(del(iev(i)))*1.66d0  + .3d0
	      else if(magtyps(1:3).eq.'R-P') then
		d1 = del(iev(i))
	        dmag = dlog10(amplit(i)/period(i)) +
     +              dlog10(d1)/3.d0 + dlog10(d1*deg2rad)/2.d0
     +              + 0.0046d0*d1 + 2.370d0
	      else
		print *,' Ms attenuation model not defined!'
	      endif

	      dmsm = dmsm + dmag
	      imsm = imsm + 1

	   else if(phid1.eq.'Lg' .or. phid1.eq.'Sg ') then
c
c             we will use ML attenuation file
c
c             magfile = 'MLCORR.DAT'
c
c             not yet implemented !!
c

	   else if(ttres.le.8.d0) then

	      magfile = ''
	      if((phid1(1:3).eq.'Pn ' .or. phid1(1:2).eq.'P ' .or.
     +              phid1(1:4).eq.'Pdif') .and. rdel.le.110. ) then

		if(magtypp.eq.'G-R') magfile = 'MB_G-R.DAT'
		if(magtypp.eq.'V-C') magfile = 'MB_V-C.DAT'

              else if(phid1(1:5).eq.'PKPbc' .and. magtypp.eq.'V-C' .and.
     +                rdel.le.150. ) then

                magfile = 'MB_V-C.DAT'

	      else if(phid1(1:5).eq.'PKPdf' .and. magtypp.eq.'V-C') then

		magfile = 'MB_V-C.DAT'

	      endif

              if(magfile.ne.'') then

	         magfile = file_check(magfile)
	         irc = 0
	         call magfact(magfile, rdel, rzo, rmcorr, irc)
	         if(irc.eq.0) then

		    dmag = dlog10(amplit(i)/period(i)) + dble(rmcorr)

	            dmbm = dmbm + dmag
	            imbm = imbm + 1

	         endif

	      endif

	   endif

	   if(dmag.ne.0.d0) then
	      write (text(i)(94:),'(1x,f12.2,1x,f6.3,1x,f4.2)') 
     +               amplit(i),period(i),dmag
	   else
	      write (text(i)(94:),'(1x,f12.2,1x,f6.3)') 
     +               amplit(i),period(i)
	   endif
	endif
      endif

      if(typctl.ge.10) then
	 print*,i,text(i)
      endif

450   continue

      if(nobst.gt.0) then
	 stmean  = stmean/dble(nobst)
	 strmean = strmean/dble(nobst)
	 rms     = dsqrt(rms/dble(nobst))
	 rmsisc  = dsqrt(rmsisc/wisc)
      endif 
      astmean = dabs(stmean)

      if(.not.single .and. astmean.ge.terrm) then
	ibad = ibad + NINT(real(astmean/terrm))
      endif

      if(ibad.ge.3) then
	 ibad0 = ibad0 + 1
      else if(ibad.eq.0) then
	 ibad0 = 0
      endif

      if(.not.last .and. astmean.gt.10.d0 .and. dabs(samean).lt.5.d0 
     +     .and.     nstat.le.5     .and. nobsa.gt.0           .and. 
     +     idelw.le.4     .and. czo.ne.'D'  .and.
     +     astmean.gt.rms .and. astmean.gt.strmean) then

	 idelw = idelw + 1
	 
	 miteras = miteras + iter
	 iter = 0
	 nextiter = 0
	 nextiter1= 0
	 dchang   = dchang0
	 iteraz = 0
	 rmso   = 9999.d0
	 datmax0  = 9999.d0
	 dazim1 = dazim0*2.0d0
	 dpam1  = dpam0*2.0d0
	 dtmin  = 9999.d0
	 dtp    = dtp0
	 dts    = dts0
	 dtm2   = dts
         check  = 9999.d0
	 setcheck = setcheck1
         setcheck2 = 15.d0*setcheck

	 if (nobs.eq.1) go to 100

	 inddel = 1
	 ddel = del(istatmin) / 2.d0
	 call delazd(stala(istatmin),stalo(istatmin),baz(istatmin),
     +               ddel,inddel,elatmg,elonm)
	 
         elatm = convlat(elatmg,1)
	 elatmr= deg2rad*elatm
         coelatm = 90.d0 - elatm
         coelatmr= deg2rad*coelatm

	 tom  = tom  + stmean
	 tome = tom  + timemin
	 rs(1) = rs(1) + astmean

	 if(output) then
          write(11,'(''Mean travel-time residual ('',f8.3,'') too '',
     +          ''large for source-time uncertainty ('',f8.3,'')'')') 
     +          stmean,var(1)
          write(11,'(''Epicenter set (lat,lon):'',2f9.3)') elatmg,elonm
	 endif
	 go to 100

      endif

      if(.not.last .and.(astmean.gt.2.0d0*var(1) .or. astmean.gt.5.d0)
     +    .and. iteras.le.1  .and. .not.single .and.
     +    dabs(astmean-astmean0) .gt. 2.d0 ) then

	 iteras = iteras + 1
	 miteras = miteras + iter
	 iter = 0
	 nextiter = 0
	 iteraz = 0
	 dchang   = dchang0

	 dazim1 = dazim0 * 1.5d0
	 dpam1  = dpam0  * 1.5d0
	 dtmin  = 9999.d0
	 dtm2   = dtm0   * 2.0d0
	 check  = 9999.d0

	 dtzo   = stmean

    	 if(iteras.le.1 .and. (czo.eq.'B' .or. czo.eq.'F') .and. 
     +      astmean.gt.5.d0) then

	     dzo = 0.d0
	     dtzo =  stmean/1.8d0

	     dzo = dtzo*8.d0
	     zo  = zo + dzo

	     if(zo.lt.0.1d0)  zo = 0.1d0
	     if(zo.gt.700.d0) zo = 700.d0
	     rs(4) = rs(4) + dabs(dzo)

	     if(typctl.gt.0) then
	        print *,'New depth set at ',zo,' [km]'
	     endif

	     if(output) then
	        write(11,'(/''New depth set at '',f6.2)') zo
	     endif

	 endif

	 tom  = tom  + dtzo
	 tome = tom  + timemin
c	 rs(1) = rs(1) + dabs(dtzo)

	 if(output) then
            write(11,'(''Mean travel-time residual ('',f8.3,'') too '',
     +            ''large for source-time uncertainty ('',f8.3,'')'')')
     +            stmean,var(1)
	 endif

	 astmean0 = astmean

	 go to 100

      endif

      if (ibad.ge.3 .and. .not.last) go to 455

      call indexx(nobs,arr,indx)

      ntext = 0
      stato = sta(iev(indx(1)))

      if(output) then
	if(.not.magflag .or. (imsm.eq.0 .and. imbm.eq.0)) then
	  write(11,'('' Stat  Delta   Azi   Phase   [used]    Onset '',
     +          ''time    Res     Baz     Res   Rayp   Res  Used''/)')
	else

	  if(imsm.gt.0) then
	     dmsm = dmsm / dble(imsm)
	  endif

	  if(imbm.gt.0) then
	     dmbm = dmbm / dble(imbm)
	  endif

	  if(imsm.gt.0 .and.imbm.eq.0) then
	     write(11,'(''Magnitude: '',f4.2,'' (Ms, '',a,'')''/)') 
     +               dmsm,magtyps(1:trimle(magtyps))
	  else if (imsm.eq.0 .and.imbm.gt.0) then
	     write(11,'(''Magnitude: '',f4.2,'' (mb, '',a,'')''/)')
     +               dmbm,magtypp

	  else if (imsm.gt.0 .and.imbm.gt.0) then
             write(11,'(''Magnitudes: '',f4.2,'' (mb, '',a,'') '',
     +                f4.2,'' (Ms, '',a,'')''/)') dmbm,magtypp,dmsm,
     +                magtyps(1:trimle(magtyps))
	  endif

	  write(11,'('' Stat  Delta   Azi   Phase   [used]    Onset '',
     +          ''time    Res     Baz     Res   Rayp   Res  Used '',
     +          ''   Amplitude   Period  MAG''/)')

        endif
      endif

      do 452 i=1,nobs

      if(sta(iev(indx(i))).ne.stato .or. i.eq.nobs) then

	ntext0 = ntext

	if(i.eq.nobs) then 

	  ntext        = ntext + 1
	  arr(ntext)   = sngl(tt(indx(i)))
	  text2(ntext) = text(indx(i))

	  ntext0       = ntext
          if(sta(iev(indx(i))).ne.stato) ntext0 = ntext - 1

	endif

        call indexx(ntext0,arr,indx2)

        if(output) then

         do 451 j=1,ntext0
         write(11,'(a)') text2(indx2(j))(1:trimle(text2(indx2(j))))
451      continue
	endif

	if(i.eq.nobs) then
           if(ntext0.eq.ntext-1) then
	       if(output) write(11,'(a)') 
     +		  text(indx(i))(1:trimle(text(indx(i))))
	   endif
	   goto 452
	endif

	stato = sta(iev(indx(i)))
	ntext = 1

      else

	ntext = ntext + 1

      endif

      arr(ntext)   = sngl(tt(indx(i)))
      text2(ntext) = text(indx(i))

452   continue
      
455   continue

      if(nobsa.gt.0) then
	 samean  = samean/dble(nobsa)
	 sarmean = sarmean/dble(nobsa)
	 rmsazi  = dsqrt(rmsazi/dble(nobsa))
      endif
      if(nobsp.gt.0) then
	 spmean  = spmean/dble(nobsp)
	 sprmean = sprmean/dble(nobsp)
	 rmsp    = dsqrt(rmsp/dble(nobsp))
      endif

      if (ibad.ge.3 .and. .not.last) go to 470

c
      if(ndt.eq.0) go to 466
      if(output) then
        write(11,'(/''Defining travel-time differences:''/)')
	write(11,'('' Stat  Delta  Phases'',11x,''Observed   Res''/)')
      endif

      i2 = 0
      sdmean  = 0.d0
      sdrmean = 0.d0
      rmsdt   = 0.d0

      do 461 i = 1,nobs-1

      do 460 j = i+1,nobs

	 if(sta(iev(i)).ne.sta(iev(j))) go to 460

	 if(phase(i).eq.phase(j)) go to 460

         if(idtu(i2+1).ne.i+j .or. used(j)(4:4).ne.'D' .or.
     +      used(i)(4:4).ne.'D' ) go to 460

         i2    = i2 + 1
 
         arr(i2) = sngl(del(iev(i)))
         dtth  = ttt(j) - ttt(i)
         dtobs = tt(j) - tt(i)
         dtres = dtobs - dtth

	 sdmean  = sdmean  + dtres
	 sdrmean = sdrmean + dabs(dtres)
	 rmsdt   = rmsdt   + q2(dtres)

	 fmis   = dtres/dsqrt(q2(tts(i))+q2(tts(j)))
	 dmisfl = dmisfl + dabs(fmis)
	 dmisf  = dmisf + q2(fmis)
	 ndmisf = ndmisf + 1
 
	 art = phase(j)(1:trimle(phase(j)))//' - '
     +         //phase(i)(1:trimle(phase(i)))

c	 if(typctl.gt.5) then
c	    print *,sta(iev(j)),del(iev(j)),
c     +	       art(1:trimle(art)),dtobs,dtres
c	 endif

	 if(output) then
	    write(text(i2),'(a5,f8.3,1x,a16,f9.3,f8.3)') 
     +	          sta(iev(j)),del(iev(i)),art,dtobs,dtres
	 endif

460   continue

461   continue

      sdmean  = sdmean / dble(ndt)
      sdrmean = sdrmean / dble(ndt)
      rmsdt   = dsqrt(rmsdt  / dble(ndt))

      call indexx(i2,arr,indx)

      if(output) then
         do 463 i=1,i2
         write(11,'(a)') text(indx(i))(1:trimle(text(indx(i))))
463      continue
      endif

466   continue

      miteras = miteras + iter
      if(output) then
         if(miteras.gt.iter) then
            write(11,'(/''Total number of iterations: '',i5)') miteras
         endif
      endif

c
c     We will now calculate the maximum azimuthal gap for as
c     defining used observations.
c

      if(nstat.gt.1) then

         call indexx(nobs,epiaz,indx)

         dmazi0 = 0.
         d1azi  = 0.
         d2azi  = 360.

         im1 = 0
         im2 = 0

         do 468 i = 1,nobs

         if (epiaz(indx(i)).lt.0.) go to 468

         if(im1.eq.0) im1 = i
         im2 = i

468      continue

         if(im2.gt.im1) then

           do 469 i =im1+1,im2

           dmazi = epiaz(indx(i)) - epiaz(indx(i-1))
           if(dmazi.gt.dmazi0) then
              d1azi  = epiaz(indx(i-1))
              d2azi  = epiaz(indx(i))
              dmazi0 = dmazi
           endif

           if(i.eq.im2) then

              dmazi = epiaz(indx(im1)) + 360. - epiaz(indx(i))
              if(dmazi.gt.dmazi0) then
                 d1azi  = epiaz(indx(i))
                 d2azi  = epiaz(indx(im1))
                 dmazi0 = dmazi
              endif

           endif

469        continue

           if(output) then
              write(11,'(/''Maximum azimuthal gap of defining '',
     +        ''observations: '',f5.1,'' [deg] -> '',f5.1,'' [deg]'',
     +        '' = '',f5.1,'' [deg]'')') d1azi,d2azi,dmazi0

           endif

         endif

      endif

c
c     output of mean errors
c
      if(output) then
         write(11,'(/''Residuals of defining data'',10x,
     +               ''RMS    MEAN-ERROR      MEAN'')')

	 if(nobst.eq.1) 
     +   write(11,'(i6,'' onset time              : '',f8.3,x,2(3x,
     +	 f8.3),''  [s]'')') nobst,rms,strmean,stmean

	 if(nobst.gt.1) 
     +   write(11,'(i6,'' onset times             : '',f8.3,x,2(3x,
     +	 f8.3),''  [s]'')') nobst,rms,strmean,stmean

	 if(nobsa.eq.1) 
     +   write(11,'(i6,'' azimuth value           : '',f8.3,x,2(3x,
     +	 f8.3),''  [deg]'')') nobsa,rmsazi,sarmean,samean

	 if(nobsa.gt.1) 
     +   write(11,'(i6,'' azimuth values          : '',f8.3,x,2(3x,
     +	 f8.3),''  [deg]'')') nobsa,rmsazi,sarmean,samean

         if(nobsp.eq.1)
     +   write(11,'(i6,'' ray parameter           : '',f8.3,x,2(3x,
     +	 f8.3),''  [s/deg]'')') nobsp,rmsp,sprmean,spmean

         if(nobsp.gt.1)
     +   write(11,'(i6,'' ray parameters          : '',f8.3,x,2(3x,
     +	 f8.3),''  [s/deg]'')') nobsp,rmsp,sprmean,spmean

	 if(ndt.eq.1) 
     +   write(11,'(i6,'' travel-time difference  : '',f8.3,x,2(3x,
     +	 f8.3),''  [s]'')') ndt,rmsdt,sdrmean,sdmean

	 if(ndt.gt.1) 
     +   write(11,'(i6,'' travel-time differences : '',f8.3,x,2(3x,
     +	 f8.3),''  [s]'')') ndt,rmsdt,sdrmean,sdmean

	 write(11,'(/''Weighted RMS of onset times (ISC type): '',
     +              f8.3,'' [s]''/)') rmsisc

         write(11,'(''Weighted misfit of input data'',9x,
     +              ''L1      L2'')')

	 if(ntmisf.ge.1) then
	   tmisf1 = dsqrt(tmisf/ntmisf)
	   tmisfl1 = tmisfl/ntmisf
           write(11,'(i6,'' onset times             :'',2(x,f8.3))')
     +     ntmisf,tmisfl1,tmisf1
	 endif

	 if(namisf.ge.1) then
	   amisf1 = dsqrt(amisf/namisf)
	   amisfl1 = amisfl/namisf
           write(11,'(i6,'' azimuth values          :'',2(x,f8.3))')
     +     namisf,amisfl1,amisf1
	 endif

	 if(npmisf.ge.1) then
	   pmisf1 = dsqrt(pmisf/npmisf)
	   pmisfl1 = pmisfl/npmisf
           write(11,'(i6,'' ray parameters          :'',2(x,f8.3))')
     +     npmisf,pmisfl1,pmisf1
	 endif

	 if(ndmisf.ge.1) then
	   dmisf1 = dsqrt(dmisf/ndmisf)
	   dmisfl1 = dmisfl/ndmisf
           write(11,'(i6,'' travel-time differences :'',2(x,f8.3))')
     +     ndmisf,dmisfl1,dmisf1
	 endif

	 nwmisf = ntmisf + namisf + npmisf + ndmisf
	 wmisf = dsqrt((tmisf + amisf + pmisf + dmisf)/nwmisf)
	 wmisfl = (tmisfl + amisfl + pmisfl + dmisfl)/nwmisf
	 write(11,'(i6,'' misfit over all         :'',2(x,f8.3))')
     +   nwmisf,wmisfl,wmisf

      endif
c
c     output of one line with all calculated source parameters and
c     quality parameters
c

      if(output) then
         write(11,'(/''T0'',25x,''LAT'',6x,''LON'',7x,''Z'',5x,
     +      ''VPVS'',3x,''DLAT'',5x,''DLON'',6x,''DZ'',7x,''DT0'',4x,
     +      ''DVPVS DEF'',4x,''RMS'' )' )
      endif
      write(6,'(/''T0'',25x,''LAT'',6x,''LON'',7x,''Z'',5x,
     +      ''VPVS'',3x,''DLAT'',5x,''DLON'',6x,''DZ'',7x,''DT0'',4x,
     +      ''DVPVS DEF'',4x,''RMS'' )' )
      
      call fetoh(tome,idum,yy,mon,mm,dd,idoy,hh,mi,sec)
 
      isec1 = nint(sec*1000)
      isec  = isec1/1000
      isec2 = isec1-isec*1000

      if(czo.eq.'D') then

        if(output) then
           write(11,'(i4,''-'',i2.2,''-'',i2.2,3i3.2,''.'',i3.3,2f9.3,
     +           f8.2,f7.2,2f9.4,f8.2,f9.3,f7.2,i5,f9.3)') 
     +           yy,mon,dd,hh,mi,isec,isec2,elatmg,elonm
     +           ,zo,vpvsm,sdlatg,var(3),var(4),var(1),sdvpvs,in,rms
	endif
        write(6,'(i4,''-'',i2.2,''-'',i2.2,3i3.2,''.'',i3.3,2f9.3,
     +           f8.2,f7.2,2f9.4,f8.2,f9.3,f7.2,i5,f9.3)') 
     +           yy,mon,dd,hh,mi,isec,isec2,elatmg,elonm
     +           ,zo,vpvsm,sdlatg,var(3),var(4),var(1),sdvpvs,in,rms

      else if(czo.eq.'F' .or. czo.eq.'B') then

        if(output) then
           write(11,'(i4,''-'',i2.2,''-'',i2.2,3i3.2,''.'',i3.3,2f9.3,
     +           f8.2,f7.2,2f9.4,''  Fixed '',f9.3,f7.2,i5,f9.3)') 
     +           yy,mon,dd,hh,mi,isec,isec2,elatmg,elonm
     +           ,zo,vpvsm,sdlatg,var(3),var(1),sdvpvs,in,rms

	endif
        write(6,'(i4,''-'',i2.2,''-'',i2.2,3i3.2,''.'',i3.3,2f9.3,
     +           f8.2,f7.2,2f9.4,''  Fixed '',f9.3,f7.2,i5,f9.3)') 
     +           yy,mon,dd,hh,mi,isec,isec2,elatmg,elonm
     +           ,zo,vpvsm,sdlatg,var(3),var(1),sdvpvs,in,rms
      endif

      if( output .and. iloc .and. modout) then

	if(imo.ge.3) then

	   write(11,'(/,''CRUST 5.1 model for source-region (max. '',
     +            ''delta ='',f6.2,'' deg):'',/,''    DEPTH   '',
     +            ''    VP        VS    DISCON'')') rmax
	
	   itrue = 0
	   elatc = elatmg
	   elonc = elonm
	   inum  = 1
	   call get_mod_mlm(itrue,inum,typctl,ierr)
	   write(11,'(3f10.3,a6)') (z(i),v0(1,i),v0(2,i),
     +           azo(i),i=1,jmod)

	else

	   write(11,'(/,''Local model '',a,'' used (max. delta ='',
     +            f6.2,'' deg):'',/,''    DEPTH       VP'',
     +            ''        VS'',''    DISCON'')') 
     +            filloc(1:trimle(filloc)),rmax
	   write(11,'(3f10.3,a6)') (z(i),v0(1,i),v0(2,i),
     +           azo(i),i=1,jmod)

	endif

      endif
 
c
c     Now let us try to locate the event better without a fixed depth.
c

c     print *,'czo ',czo,' iterz ',iterz,' zoflag ',zoflag

470   if(czo.eq.'B' .and. ((.not. zoflag) .or. iterz.lt.2)) then
	 zoflag=.true.
	 czo ='D'
	 if(zo.le.0.1d0) zo = 33.d0
	 rs(4)  = sdzo
	 rs(1)  = dsqrt(rs(1)*rs(1)+25.d0)

	 if(idelw.gt.1) idelw = idelw - 2

	 miteras = miteras + iter
	 iter = 0
	 nextiter = 0
	 nextiter1= 0
	 dchang   = dchang0
	 iteras = 0
	 iteraz = 0
	 rmso   = 9999.d0
	 rmsold   = 9999.d0
	 datmax0  = 9999.d0
	 dazim1 = dazim0*2.0d0
	 dpam1  = dpam0*2.0d0
	 dtmin  = 9999.d0
	 dtp    = dtp0
	 dts    = dts0
	 dtm2   = dts
         check  = 9999.d0
	 setcheck = setcheck1
         setcheck2 = 15.d0*setcheck
	 astmean0 = 0.d0

	 go to 100
      endif

      go to 9999

9998  continue
      if(output) then
	 write(11,'(/''Could not find a better or non-oscillating'',
     +            '' solution without fixed depth'')')
      endif
      write(6,'(/''Could not find a better or non-oscillating'',
     +            '' solution without fixed depth'')')

9999  continue

      close(11)

      stop

      end
