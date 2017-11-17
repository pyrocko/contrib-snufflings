c
      program HYPOMOD_1_1b

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character  version*20
      parameter (version='HYPOMOD Version 1.1b')

c
c     Version 1.0  September 2001
c
c                       Version 1.0a October 2001
c
c                       Version 1.0b November 2001
c
c                       Version 1.0c February 2002 - new libtau software
c
c                       Version 1.0d July 2002 - Four different global 
c                                    models, new format for source time
c
c                       Version 1.1  October 2002 - all changes needed
c                                    to be compatible with HYPOSAT_4_4
c
c                       Version 1.1a May 2003 - all changes needed
c                                    to be compatible with HYPOSAT_4_4a
c
c                       Version 1.1b July 2003 - all changes needed
c                                    to be compatible with HYPOSAT_4_4b
c
c     last corrections:  24 July 2003
c
c
c     This programm calculates for a given seismic event all residuals
c     observed travel times, backazimuth, and slowness values.
c
c     All input files  and the output file are identical to hyposat. 
c     See HYPOSAT manual for details. However some features are just 
c     ignored because we do not invert any data!
c
c
c                                 Johannes Schweitzer
c    
c                                 NORSAR
c                                 P.O.Box 51
c                                 N-2027 KJELLER
c                                 Norway
c
c                                 johannes@norsar.no
c
c
c
c     calls :   delazd, depi, dlsq, fetoh, fhtoe, get_station, 
c               ellip, ellip2, indexx, 
c               tauget_mod, testphase, hyposat_geo, ttloc, get_mod_mlm,
c
c     functions: alpha1, alpha2, convlat, phase_type, trimle, crustc, 
c                ddmax, dirdel, dmean, q2, radloc, tauget_ray
c
c
c     PARAMETER settings for: mstat, mread
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
     +                 radloc
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
      parameter  (mread = 800)

      character phase(mread)*8,phid*8,text(mread)*125,
     +          phid2*8,text2(30)*125,phid1*8,phase_t*1,phase0*8,
     +          string*130,touse(mread)*6,touse0*6,phidr*8
      dimension azi(mread),tt(mread),p(mread),azis(mread),tts(mread),
     +          ps(mread),period(mread),amplit(mread)
      real      arr(mread),epiaz(mread),dmazi,dmazi0,d1azi,d2azi
      integer   iev(mread),indx(mread),indx2(30)

c
c     variables for travel-time calculations
c
      include 'ttimes.h'

      dimension dpdh(mphas),dddp(mphas)

      character art*16,modnam*20,modnam2*20,modn*20,modnam3*20,
     +          modnam4*20

      real*4 rzo,rdel,ecor,fla1,fla2,fla3,razi,rbaz,pa,rzo1,
     +       rmax,rzoe,rdel1,rmcorr

c
c     common block with informations about the local/regional model
c     (used also for CRUST 5.1)
c
c     maxla = maximum number of layers for a local/regional velocity
c             model
c

      parameter (maxla = 61)
      DIMENSION V0(2,maxla),z(maxla),v0g(2,maxla),zg(maxla)
      character azo(maxla)*4,azog(maxla)*4,filloc*80,mtyp*3,mpath*80
      integer   jmod,jmodg,iread,imo

      COMMON /MODEL/  v0,z,elev,elatc,elonc,zmax,elat2,elon2,elev2,
     +                jmod,iread,imo,mtyp

      COMMON    /MODELC/ filloc,azo,mpath
      COMMON    /MODELG/ v0g,zg,elevg,jmodg,azog,zmaxg
c
c     other variables
c
      integer yy,mon,dd,hh,idoy,ierr,typctl,mi,idum, isreg, regnum
      character mm*4,name*48
      double precision lat,lon
      real*4  elevs, sec, rlat, rlon

      character title*80, region*80, czo1*1, magtypp*3,
     +          magtyps*6

c
      dimension ttt(mread)

      logical   iaspflag, vlflag  , stcorfl, iloc, surf,
     +          output, modout,
     +          magflag, mod2flag, conr, tauget_ray, rayok,
     +          mod3flag,mod4flag

c
c     some constants and starting values
c

      pi      = 4.d0*datan(1.d0)
      deg2rad = pi / 180.d0
      rearth  = 6371.d0
      grad1   = 1.d0/(deg2rad*rearth)

      dtmin   = 9999.d0

      imodout = 0 

      iaspflag = .false.
      vlflag   = .false.
      stcorfl  = .false.
      iloc     = .false.
      output   = .false.
      modout   = .false.
      magflag  = .false.
      mod2flag = .false.
      rayok    = .false.

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
      outputfile  = 'hypomod-out'
      inputfile   = 'hyposat-in'
      mpath       = './'
      statcorfile = ''
      vpl = 99.d0
      vsl = 0.d0
      zo1 =  0.d0
      typctl = 4
      islow = 1
      indph0 = 3333
      inddiff = 1
      epilat0 = -999.d0
      epilon0 = -999.d0
      tome0   = 0.d0

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

      if(string(1:21).eq.'STARTING SOURCE DEPTH') then
          read (string(38:),*) zo1
	  go to 1
      endif

      if(string(1:24).eq.'STARTING SOURCE LATITUDE') then
	  abc = -999.0d0
          read (string(38:),*) abc
          if(abc.ge.-90.d0 .and. abc.le.90) epilat0 = abc
	  go to 1
      endif

      if(string(1:25).eq.'STARTING SOURCE LONGITUDE') then
	  abc = -999.0d0
          read (string(38:),*) abc
          if(abc.ge.-180.d0 .and. abc.le.180) epilon0 = abc
	  go to 1
      endif

      if(string(1:12).eq.'OUTPUT LEVEL') then
          read (string(38:),*) typctl
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

      if(typctl.gt.4) then
         print *,'modnam ',modnam
         if(mod2flag) print *,'modnam 2',modnam2
         if(mod3flag) print *,'modnam 3',modnam3
         if(mod4flag) print *,'modnam 4',modnam4
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
         print *,'epilat0 = ',epilat0
         print *,'epilon0 = ',epilon0
         print *,'typctl = ',typctl
         print *,'islow = ',islow
         print *,'indph0 = ',indph0 
         print *,'inddiff = ',inddiff
	 print *,'Magnitude flag = ', magflag
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

      if(vpl.ne.99.d0) then

        if(vpl.eq.0.d0) vpl=5.8d0
        if(vsl.eq.0.d0) vsl=vpl/dsqrt(3.d0)
	vlflag = .true.

	if (stcorfl) open(12,file=statcorfile)

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

      stalam = 0.d0
      stalom = 0.d0

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

      touse(ii)=   'TASDR '
      if(touse0.ne.'      ') then

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

      if(tts(ii).le.0.d0) tts(ii) = 2.d0

      if(azi(ii).eq.-1.d0) azi(ii) = -999.d0

      if(azi(ii).eq.-999.d0) touse(ii)(2:2)=' '

c

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

11    continue

      if(pin.eq.-1.d0) then
	 pin    = -999.d0
      endif

      if(islow.eq.0 .and. pin.ne.-999.d0) then
	 p(ii)  = drad/(grad1*pin)
      else
	 p(ii) = pin
      endif

      if(p(ii).eq.-999.d0) touse(ii)(3:3) = ' '

      phase_t = phase(ii)(1:1)

      if(phase_t.eq.'P' .or. phase_t.eq.'p')  then
         if (istaph(iev(ii)).eq.0 .or. istaph(iev(ii)).eq.2) 
     +       istaph(iev(ii))=istaph(iev(ii)) + 1
      else if(phase_t.eq.'S' .or. phase_t.eq.'s') then
         if (istaph(iev(ii)).eq.0 .or. istaph(iev(ii)).eq.1)
     +       istaph(iev(ii))=istaph(iev(ii)) + 2
      endif

      if(typctl.gt.8) then
         print *,stat,tt(ii),phase(ii),azi(ii),p(ii)
      endif

12    continue

14    if(stcorfl) close (12)

      nobs  = ii

      close(10)

      if(nobs.eq.1) then
	 rzo1  = 0.
	 rdel1 = 105.
         call tauget_mod(rzo1,rdel1,nphas,phcd,ttc,dtdd,
     +                   dtdh,modnam)
	 pdif   = dble(dtdd(1))
	 phase0 = phase(1)
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

50    continue

      tome = tome0

      zo = zo1

      elonm  = epilon0

      elatmg = epilat0
      elatm  = convlat(elatmg,1)

      if(iloc) then
         if(mod2flag .or. mod3flag .or. mod4flag) then
             write(11,'(''First reference models  : '',a,'' and '',
     +             a)') filloc(1:trimle(filloc)),modnam
             if(mod2flag) then
                write(11,'(''Second reference model  : '',a)')
     +               modnam2
             endif
             if(mod3flag) then
                write(11,'(''Third reference model   : '',a)')
     +               modnam3
             endif
             if(mod4flag) then
                write(11,'(''Fourth reference model  : '',a)')
     +               modnam4
             endif
         else
             write(11,'(''Reference models  : '',a,'' and '',a)')
     +             filloc(1:trimle(filloc)),modnam
         endif
      else
         if(mod2flag .or.  mod3flag .or. mod4flag) then
             write(11,'(''First reference model   : '',a)') modnam
             if(mod2flag) then
                write(11,'(''Second reference model  : '',a)')
     +               modnam2
             endif
             if(mod3flag) then
                write(11,'(''Third reference model   : '',a)')
     +               modnam3
             endif
             if(mod4flag) then
                write(11,'(''Fourth reference model  : '',a)')
     +               modnam4
             endif
         else
             write(11,'(''Reference model   : '',a)') modnam
         endif
      endif

      call fetoh(tome,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

      if(output) then

         write(11,'(/''The source parameters''/)')

         write(11,'(''Source time  :'',i5,4i3.2,f7.3
     +           )')  yy,mon,dd,hh,mi,sec
         write(11,'(''        or'',14x,f14.3
     +           )') tome

         write(11,'(''Epicenter lat:'',14x,f10.4,
     +           '' [deg]'')')  elatmg
         write(11,'(''Epicenter lon:'',14x,f10.4,
     +           '' [deg]'')')  elonm
         write(11,'(''Source depth :'',15x,f7.2,
     +         ''   [km]''/)') zo
      endif

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
      stmean    = 0.d0
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

	 imod2 = 0

	 if (mod2flag .and. touse(i)(6:6).eq.'2') then
	   modn = modnam2
	 else if (mod3flag .and. touse(i)(6:6).eq.'3') then
	   modn = modnam3
	 else if (mod4flag .and. touse(i)(6:6).eq.'4') then
	   modn = modnam4
         else
	   modn = modnam
	   if(iloc) imod2 = 1
	 endif

         if(imod2.eq.0 .or. rdel.gt.rmax .or. zo.gt.zmax) then

	   call tauget_mod(rzo,rdel,nphas,phcd,ttc,dtdd,
     +	                 dtdh,modn)

         else

	   ierr = 0
	   indph = istaph(iev(i))*10000 + indph0
	   elatc = elatmg
	   elonc = elonm

	   elat2 = stala(iev(i))
	   elon2 = stalo(iev(i))

	   call ttloc(rzo,rdel,czo1,nphas,ttc,dtdd,dtdh,dpdh,dddp,
     +                phcd,rmax,typctl,ierr,indph)

         endif

      endif

      phid = phase(i)
      ttobs  = tt(i) 

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

      if(nobs.eq.1 .and. phid.eq.phase0) then
	 if(p(i).lt.pdif)  phid='PKPab'
	 if(p(i).ge.pdif)  phid='P'
	 if(p(i).ge.20.2d0) phid='S1'
	 if(phase0.ne.phase(i)) then
	    phid     = phase(i)
	    phase(i) = phase0
         endif
      endif
      
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
	 
	  if(phid1(2:5).eq.'dif') then
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
     +                     dpaa*dsqrt(dl*dl-hsta*hsta)/radkm
                    else if(dl.lt.0.d0) then
                      th = dl/vloc + 
     +                     dpaa*dsqrt(dl*dl-hsta*hsta)/radkm
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

            endif

c
c     correction for surface reflections (e.g.: PnPn,...,PP,SS,P'P')
c
            if( phid(1:1).eq.phid(2:2) .or.
     +          phid(1:2).eq.phid(3:4)  ) then

                indr = 2
                trefl = crustc(phase_t,dpaa,indr,typctl)

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
     +                  del0,ttray)

                azi0 = azie(iev(i))
                if(dpa.lt.0.d0) azi0 = alpha2(azie(iev(i))-180.d0)

                inddel = 1
                call delazd(stala(iev(i)),stalo(iev(i)),azi0,del0,
     +                      inddel,elatc,elonc)

                indr = 3
                trefl   = crustc(phase_t,dpaa,indr,typctl)

            endif

            if (trefl.ne.0.d0.and.typctl.gt.6) then
               print *,'dirdel: ',phid,' azi ',azi0,' del ',del0
               print *,'  lat ',elatc,' lon ',elonc,' trefl ',trefl
            endif

         endif

         ttt(i) = tome + dble(ttc(j) + ecor) + th + tcrust + trefl

	 ttres  = ttobs - ttt(i)

	 if(ttres.gt.100.d0 .and. phase(i)(1:2).eq.'P1' .and. 
     +      phid(1:3).eq.'Pdi' .and. touse(i)(1:1).eq.'T') then
	    phid='PKPdf'
            go to 420
         endif

	 if(imin.eq.0) then

	    if(touse(i)(1:1).eq.'T') then
	       stmean  = stmean + ttres
	       strmean = strmean + dabs(ttres)
	       rms     = rms    + ttres*ttres
	       rmsisc  = rmsisc + q2(ttres)*tts(i)
	       wisc    = wisc   + tts(i)
	       epiaz(i)= razi
	       nobst   = nobst + 1
            endif

	    pares = p(i)-dabs(dpa)
	    if(touse(i)(3:3).eq.'S') then
	       spmean  = spmean + pares
	       sprmean = sprmean + dabs(pares)
	       rmsp    = rmsp + q2(pares)
	       epiaz(i)= razi
	       nobsp   = nobsp + 1
            endif

            ttti = ttt(i) 

	    go to 430

         else if(imin.eq.1.and.dabs(ttres).lt.dabs(dtmin)) then

	    if(touse(i)(1:1).eq.'T') touse(i)(1:1)= ' '
	    phid2 = phid1
	    dtmin = ttres
	    pares2 = p(i)-dabs(dpa)
	    dpa2  = dpa

            ttti = ttt(i) 

	 endif

      endif

420   continue

      if (imin.eq.1) go to 425

c
c     Try it with another phase-name from the same phase-type.
c
 
      if(nobs.eq.1) go to 410

      call testphase (phid,icha)

      if(icha.ne.999) go to 410

      if(imin.eq.0) then
         imin  = 1
         dtmin = 9999.d0
         go to 410
      endif

425   if(dabs(dtmin).le.15.d0 .or. touse(i)(1:1).eq.'T' .or.
     +   (touse(i)(4:4).eq.'D')) then
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
         samean  = samean + azires
         sarmean = sarmean + dabs(azires)
	 rmsazi  = rmsazi + q2(azires)
	 epiaz(i)= razi
         nobsa   = nobsa + 1

      endif

      if(phase(i).eq.phid1) then

        write(text(i),'(a5,f8.3,f7.2,1x,a8,8x,2i3.2,f7.3,f8.3,
     +                f7.2,f8.2,f6.2,f7.2,1x,a6)')
     +                stato,rdel,razi,phase(i),hh,mi,sec,
     +                ttres,azi(i),azires,p(i),pares,touse(i)

      else

        write(text(i),'(a5,f8.3,f7.2,1x,2a8,2i3.2,f7.3,f8.3,
     +                f7.2,f8.2,f6.2,f7.2,1x,a6)') 
     +                stato,rdel,razi,phase(i),phid1,hh,mi,sec,
     +                ttres,azi(i),azires,p(i),pares,touse(i)

      endif

      phase(i) = phid1
      ttt(i)   = ttti

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
c     the standard IASPEI formula:
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
     +              phid1(1:5).eq.'Pdif') .and. rdel.le.110. ) then

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

450   continue

      if(nobst.gt.0) then
	 stmean  = stmean/dble(nobst)
	 strmean = strmean/dble(nobst)
	 rms     = dsqrt(rms/dble(nobst))
	 rmsisc  = dsqrt(rmsisc/wisc)
      endif 

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
     +          ''time    Res     Baz     Res   Rayp   Res  Used    '',
     +          ''Amplitude   Period  MAG''/)')
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

      if(inddiff.ne.1) go to 466

      i2 = 0
      ndt = 0
      sdmean  = 0.d0
      sdrmean = 0.d0
      rmsdt   = 0.d0

      do 461 i = 1,nobs-1

      if(phase(i)(1:2).eq.'LR'.or.phase(i)(1:2).eq.'LQ') go to 461

      do 460 j = i+1,nobs

         if(phase(j)(1:2).eq.'LR'.or.phase(j)(1:2).eq.'LQ') go to 460

	 if( ( sta(iev(i)).ne.sta(iev(j))       ) .or. 
     +       ( phase(j)(1:trimle(phase(j))).eq. 
     +         phase(i)(1:trimle(phase(i)))     ) ) go to 460

         i2    = i2 + 1
	 ndt   = ndt + 1
 
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

      if (ndt.le.0) go to 466

      sdmean  = sdmean / dble(ndt)
      sdrmean = sdrmean / dble(ndt)
      rmsdt   = dsqrt(rmsdt  / dble(ndt))

      call indexx(i2,arr,indx)

      if(output) then

         write(11,'(/''Travel-time differences:''/)')
	 write(11,'('' Stat  Delta  Phases'',11x,''Observed   Res''/)')

         do 463 i=1,i2
         write(11,'(a)') text(indx(i))(1:trimle(text(indx(i))))
463      continue
      endif

466   continue
c
c     We will now calculate the maximum azimuthal gap for as
c     defining used observations.
c

      if(istat.gt.1) then

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
         write(11,'(/''Residuals of data'',19x,
     +               ''RMS    MEAN-ERROR      MEAN'')')

	 if(nobst.eq.1) 
     +   write(11,'(i6,'' onset time              : '',f8.3,x,3(3x,
     +	 f8.3),''  [s]'')') nobst,rms,strmean,stmean

	 if(nobst.gt.1) 
     +   write(11,'(i6,'' onset times             : '',f8.3,x,2(3x,
     +	 f8.3),''  [s]'')') nobst,rms,strmean,stmean

	 if(nobsa.eq.1) 
     +   write(11,'(i6,'' backazimuth value       : '',f8.3,x,2(3x,
     +	 f8.3),''  [deg]'')') nobsa,rmsazi,sarmean,samean

	 if(nobsa.gt.1) 
     +   write(11,'(i6,'' backazimuth values      : '',f8.3,x,2(3x,
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
     +              f8.3,'' [s]'')') rmsisc

         write(11,'(/''Weighted misfit of input data'',9x,
     +              ''L1       L2'')')

	 if(ntmisf.ge.1) then
	   tmisf1 = dsqrt(tmisf/ntmisf)
	   tmisfl1 = tmisfl/ntmisf
           write(11,'(i6,'' onset times             :'',2(x,f8.3))')
     +     ntmisf,tmisfl1,tmisf1
	 endif

	 if(namisf.ge.1) then
	   amisf1 = dsqrt(amisf/namisf)
	   amisfl1 = amisfl/namisf
           write(11,'(i6,'' backazimuth values      :'',2(x,f8.3))')
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

9999  continue

      close(11)

      stop

      end
