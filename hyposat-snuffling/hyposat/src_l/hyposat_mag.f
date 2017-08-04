      subroutine magfact (inpfile, dist, depth, corr, irc)
      
      include 'magpar.h'
      

c.======================================================================
c.    PURPOSE                                                           
c     Common_for_magnitude_calculation_tables                       GL<<
c.----------------------------------------------------------------------
c.    PARAMETER                                                         
c                                                                       
c===> Mb calculation of P-waves
c
c---> Distance detpth correction
c
c..   mbntbd      - Number of distance samples in tables                      
c..   mbntbz      - Number of depth samples in tables
c
c     i = 1, mbntbd
c
c..   mbtbd(i)    - Angular distance (degrees)
c
c     j = 1, mbntbz
c
c..   mbtbz(i)    - Depth (km)
c
c..   mbtbqf(i,j) - Distance depth correction for mb after 
c                   Veith Clawsson
c
c---> Help array
c     
c     mtbam (i,j) - Help array
c
c.======================================================================
                                                                        
      integer*4 mbntbd, mbntbz 
                                                                  
      real*4    mbtbd(maxtbd) 
                                                                  
      real*4    mbtbz(maxtbz)
                                                                  
      real*4    mbtbqf(maxtbd, maxtbz)
      real*4    mtbam(maxtbd, maxtbz)                             

      character*(*) inpfile
      real          dist, depth, corr
      integer*4     irc
      
c.======================================================================
c.    PURPOSE
c     Get_correction_factors_for_mb_calculations                     <<
c
c     Most of the routines used are built from subroutines made
c     Hans Israelsson at CSS.
c.----------------------------------------------------------------------
c.    KEYWORDS
c.----------------------------------------------------------------------
c.    INPUT
c     
c     
c..   inpfile   - Input file with correction values. E.g.,
c                 'iasp91.pfact'
c                 'iasp91.qfvc'
c                 
c..   dist      - Distance in degrees
c..   depth     - Depth in km
c     
c.    OUTPUT
c     
c..   corr      - Magnitude correction
c..   irc       - Return code (0=> O.K.) (-1=> Fatal error)
c                 (1=>returned null value)
c     
c.----------------------------------------------------------------------
c.    PROGRAMMER    Tormod Kvaerna
c.    CREATION_DATE 040293
c.    MADE_AT  NTNF/NORSAR
c     Pb. 51
c     N-2007 Kjeller
c     
c.    MODIFICATION
c.    CORRECTION
c.======================================================================
      
c---- 
c     Internal declarations
c---- 

      integer*4    iderrdz, luerr
      logical*4    first
      data      first/.true./
c---- 
c     Set print flag
c---- 
      irc = 0
c---- 
c     Read table data the first time the routine is called
c---- 
      luerr = 6
      if (first) then
         first  = .false.
         call rdtab2 (inpfile, maxtbd, maxtbz, luerr, mbntbd, mbntbz,
     $        mbtbd, mbtbz, mbtbqf, mtbam, ierr)
         if (ierr .ne. 0) then
            irc = -1
            return
         endif
      endif

      call qfcal0 (dist,depth,
     &     mbntbd,mbntbz,
     &     mbtbd,mbtbz,mbtbqf,
     &     corr,iderrdz)

      return

      end
      subroutine rdtab2 (filnam,maxtbd,maxtbz,luerr,
     .                   ntbd,ntbz,tbd,tbz,tbtt,tbam,ierr)
c
c Read travel-time tables for a given wave from file FILNAM.
c
c Input
c
c   filnam :  name of file to be read
c   maxtbd :  maximum number of distance samples in tables
c   maxtbz :  maximum number of depth samples in tables
c   luerr  :  logical unit number for error output
c
c Output
c
c   ntbd      :  number of distance samples in tables
c   ntbz      :  number of depth samples in tables
c     for i = 1..ntbd :
c   tbd(i)    :  angular distance of (i,j)'th sample
c     for j = 1..ntbz :
c   tbz(j)    :  depth of (i,j)'th sample
c     for i = 1..ntbd, j = 1..ntbz :
c   tbtt(i,j) :  T-T (i,j)'th sample (period for R wave)
c   tbam(i,j) :  amplitude of (i,j)'th sample
c   ierr      :  error flag; 0 no error,
c                            1 file won't open,
c                            2 unexpected end-of-file
c
c Note: if file FILNAM will not open, 
c       then NTBD and NTBZ are returned zero.
c
      character*(*) filnam
      character*80  string
      character*1   toka
      integer*4     iunit, ierr, trimle
      real*4        tbd(maxtbd),tbz(maxtbz)
      real*4        tbtt(maxtbd,maxtbz), tbam(maxtbd,maxtbz)
c
c Open file
c

      luerr = 6
      k = trimle(filnam)
      iunit = 72

      open (iunit,file=filnam,status='old',iostat=ios)
      if (ios.ne.0) then
	write (luerr,'(3a)') '? File ',filnam(1:k),' will not open'
        ntbd = 0
        ntbz = 0
	ierr = 1
        return
      end if
c
c Read whether amplitude and auxiliary information is included
c  and the title of the table.
c
      read (iunit,'(a1)',end=90) toka
c
c Read depth sampling
c
      read (iunit,*,end=90) ntbzx
      ntbz = min(maxtbz,ntbzx)
      if (ntbzx.gt.maxtbz) then
        write (luerr,'(2a)') '? Too many depth samples in file ',
     &                       filnam(1:k)
        write (luerr,'(2(a,i4))') '  Number in file:',ntbzx,
     &                            '  Number kept:',maxtbz
      end if
      read (iunit,*,end=90) (tbz(i),i=1,ntbz),(dum,i=ntbz+1,ntbzx)
c
c Read distance sampling
c
      read (iunit,*,end=90) ntbdx
      ntbd = min(maxtbd,ntbdx)
      if (ntbdx.gt.maxtbd) then
        write (luerr,'(2a)') '? Too many distance samples in file ',
     &                       filnam(1:k)
        write (luerr,'(2(a,i4))') '  Number in file:',ntbdx,
     &                            '  Number kept:',maxtbd
      end if
      read (iunit,*,end=90) (tbd(i),i=1,ntbd),(dum,i=ntbd+1,ntbdx)
c
c Read tables
c
      if (toka.eq.'y' .or. toka.eq.'Y') then
        do 100 j = 1,ntbz
          read (iunit,'(a)',end=90) string
          read (iunit,*,end=90) (tbtt(i,j),tbam(i,j),dum,dum, i=1,ntbd),
     &                       (dum,dum,dum,dum, i=ntbd+1,ntbdx)
100     continue
      else
        do 110 j = 1,ntbz
          read (iunit,'(a)',end=90) string
          read (iunit,*,end=90) (tbtt(i,j), i=1,ntbd),
     &                       (dum, i=ntbd+1,ntbdx)
110     continue
      endif
c
      close (iunit)
      ierr = 0
      return
c
c Error
c
90    write (luerr,'(3a)') '? Unexpected end-of-file on file ',
     $     filnam(1:k)
      close (iunit)
      ierr = 2
      return
      end
      subroutine qfcal0 (delta,zfoc,
     &                   ntbd,ntbz,tbd,tbz,tbqf,
     &                   dcalx,iterr)
c
c Computes q-correction for magnitude calculation from
c q-tables as a function of distance and depth.
c
c Correction is obtained from the interpolation
c  of tables.  Values are obtained by extrapolation
c  for points off the ends of the curves.  A point in a hole in
c  the curve is considered bad.  No values are returned in this case.
c
c  NOTE: This can be used for interpolation in other two-dimensional
c  tables not necessarily a function of distance and depth.
c
c Calls BRACK and HOLIN2 and directly.
c Calls BRACK, HOLINT, QUAINT, FIXHOL and HERMIT indirectly.
c
c Input
c
c   zfoc :  event focal depth (km below sea level)
c   delta : distance from event to station (deg)
c   ntbd : number of distance samples in  tables.
c   ntbz : number of depth samples in  tables.
c   tbd(i) : distance samples for tables (deg).
c   tbz(j) : depth samples in tables (km).
c   qftt(i,j) :  correction tables (magnitude units).
c
c Output
c
c   dcalx : calculated correction
c   iterr: error code 
c        =  0, no problem, normal interpolation
c        = 11, distance-depth point (x0,z0) in hole in T-T curve
c        = 12, x0 < x(1)
c        = 13, x0 > x(max)
c        = 14, z0 < z(1)
c        = 15, z0 > z(max)
c        = 16, x0 < x(1) and z0 < z(1)
c        = 17, x0 > x(max) and z0 < z(1)
c        = 18, x0 < x(1) and z0 > z(max)
c        = 19, x0 > x(max) and z0 > z(max)
c        [NOTE: if any of these codes is negative (e.g. iderr = -17),
c               the datum was used to compute event location]
c
c  AUTHOR 
c     Hans Israelsson, March 6, 1989,
c     updating of original subroutine by SAIC, Lajolla
c

      include 'magpar.h'
      dimension tbd(maxtbd), tbz(maxtbz)
      dimension tbqf(maxtbd,maxtbz)

c Find relevant range of table depths
 
      call brack (ntbz,tbz,zfoc, ileft)
      jz = max(1,ileft-1)
      nz = min(ntbz,ileft+2) - jz + 1
 
c HOLIN2 does bivariate interpolation
 
      call holin2 (ntbd,nz,tbd,tbz(jz),tbqf(1,jz),
     &             maxtbd,-1.0,delta,zfoc, dcalx,dtddel,dtdz,
     &             dcross,iext,jext,ibad)

c Interpolation point in hole in curve. Value no good.
	    if (ibad.ne.0) then
	      iterr = 11

c Interpolation point less than first distance point in curve.
	    else if (iext.lt.0 .and. jext.eq.0) then
              iterr = 12

c Interpolation point greater than last distance point in curve.
	    else if (iext.gt.0 .and. jext.eq.0) then
              iterr = 13

c Interpolation point less than first depth point in curve.
	    else if (iext.eq.0 .and. jext.lt.0) then
              iterr = 14

c Interpolation point greater than last depth point in curve.
	    else if (iext.eq.0 .and. jext.gt.0) then
              iterr = 15

c Interpolation point less than first distance point in curve
c  and less than first depth point in curve.
	    else if (iext.lt.0 .and. jext.lt.0) then
              iterr = 16

c Interpolation point greater than last distance point in curve
c  and less than first depth point in curve.
	    else if (iext.gt.0 .and. jext.lt.0) then
              iterr = 17

c Interpolation point less than first distance point in curve
c  and greater than first depth point in curve.
	    else if (iext.lt.0 .and. jext.gt.0) then
              iterr = 18

c Interpolation point greater than last distance point in curve
c  and greater than first depth point in curve.
	    else if (iext.gt.0 .and. jext.gt.0) then
              iterr = 19

c Reset error code to 0 if valid table interpolation.
            else 
	      iterr = 0
            end if

      return
      end

c NAME
c	brack -- Bracket an array of interpolative values.

c FILE
c	brack.f

c SYNOPSIS
c	Using bi-section, brack an array of interpolative values by 
c	performing a binary search.

c DESCRIPTION
c	Subroutine.  Perform a binary search to find those elements of 
c	array x() that bracket x0.  Given the array x(i), i = 1,.,N, in 
c	non-decreasing order, and given the number x0, this routine finds 
c	ileft from 0..n, such that (pretend x(0) = -infinity, 
c	x(n+1) = +infinity):

c		x(ileft) <= x0 <= x(ileft+1)
c		x(ileft) < x(ileft+1)

c	Note that x() may contain duplicate values, but ileft will still 
c	point to a non-zero interval.

c	---- On entry ----
c	n:	Dimension of input vector (array), x()
c	x(n):	One-dimensional input array of values to be bracketed
c	x0:	Value being compared against

c	---- On return ----
c	ileft:	Left bracketed indice

c DIAGNOSTICS
c

c FILES
c

c NOTES
c

c SEE ALSO
c

c AUTHOR
c


      subroutine brack (n, x, x0, ileft)
 


c     ---- On entry ----

      integer*4 n
      real*4    x(n), x0

c     ---- On return ----

      integer*4 ileft

c     ---- Internal variables ----

      integer*4 i, imid, iright


c     Initialize

      ileft  = 0
      iright = n + 1

 1000 imid = (ileft+iright)/2
      if (imid.eq.ileft) then
         return
      else if (x0.lt.x(imid)) then
         iright = imid
         goto 1000
      else if (x0.gt.x(imid)) then
         ileft = imid
         goto 1000
      end if 
 
c     Special case: The point x(imid) found to equal x0.  Find bracket 
c     [x(ileft),x(ileft+1)], such that x(ileft+1) > x(ileft).

      do 1010 i = imid+1, n
      if (x(i).gt.x0) then
         ileft = i-1
         return
      end if
 1010 continue

      do 1020 i = imid-1, 1, -1
      if (x(i).lt.x0) then
         ileft = i
         return
      end if
 1020 continue

      ileft = 0

      return
      end


c NAME
c	fixhol -- Deal with bad values (holes) in function.

c FILE
c	fixhol.f

c SYNOPSIS
c	When "holes" in a given function are found, look for "good" samples.

c DESCRIPTION
c	Subroutine.  Fix up the sampled function f(x) to allow for "holes"
c	with bad values.  A bad function value at a sample point x(i) is 
c	assumed to occur when the function sample f(i) has the value fbad.  
c	For interpolation purposes, the function is assumed then to be fbad 
c	in the intervals surrounding x(i) up to the next "good" samples, 
c	where a jump discontinuity is assumed to occur.

c	Given the original function samples -- x(i), f(i), i = 1, m -- this 
c	routine creates a new set of samples -- xs(i), fs(j), j = 1, ms --
c	in which the intervals of bad values and discontinuities between
c	good and bad values are explicitly sampled.  To create a 
c	discontinuity at a given point, the point is included as two 
c	samples with different function values.

c	---- Example ----
c	x =  0.0  2.0   3.0  7.0  8.5  12.0  14.5  18.0  19.0  20.5 21.5  22.0
c	f =  2.3  1.1  fbad  7.6  4.5  fbad  fbad  12.1  fbad   6.2  4.3  fbad

c	xs =  0.0  2.0   2.0  7.0  7.0   8.5   8.5  20.5 20.5  21.5  21.5
c	fs =  2.3  1.1  fbad fbad  7.6   4.5  fbad  fbad  6.2   4.3  fbad

c	---- Indexing ----
c	i = 1, m;	j = 1, ms;

c	---- On entry ----
c	m:	Number of x() samples
c	x(i):	Sample values of x(); must be ordered
c	f(i):	Value of function at x(i)
c	fbad:	Function value signifying that f(x) is not well-defined (bad)

c	---- On return ----
c	ms:	Number of new x() samples; May be as large as 1 + (4*m)/3
c	xs(j):	New sample values of x()
c	fs(j):	New value of function at x(j)

c DIAGNOSTICS
c

c FILES
c

c NOTES
c

c SEE ALSO
c

c AUTHOR
c


      subroutine fixhol (m, x, f, fbad, ms, xs, fs)


c     ---- On entry ----

      integer*4 m
      real*4    f(m), fbad, x(m)

c     ---- On return ----

      integer*4 ms
      real*4    fs(*), xs(*)

c     ---- Internal variables ----

      integer*4 i


c     Trivial case
 
      if (m.le.0) then
         ms = 0
         return
      end if

c     Set up first point

      ms = 1
      xs(1) = x(1)
      fs(1) = f(1)

c     Do the rest

      do 1010 i = 2, m
         if (f(i).ne.fbad) then
            if (fs(ms).ne.fbad) then
               if (x(i).eq.xs(ms)) then
                  if (f(i).eq.fs(ms)) go to 1010
               end if
            else
               if (ms.gt.1) ms = ms + 1
               xs(ms) = x(i)
               fs(ms) = fbad
            end if
            ms     = ms + 1
            xs(ms) = x(i)
            fs(ms) = f(i)
                else
            if (fs(ms).ne.fbad) then
               if (ms.gt.1) then
                  if (fs(ms-1).eq.fbad) then
                     ms = max(1,ms-2)
                     go to 1010
                  end if
               end if
               ms     = ms + 1
               xs(ms) = xs(ms-1)
               fs(ms) = fbad
            end if
         end if
 1000    if (ms.gt.2) then
            if (xs(ms).eq.xs(ms-2)) then
               fs(ms-1) = fs(ms)
               ms       = ms - 1
               go to 1000
            end if
         end if
 1010 continue

      return
      end


c NAME
c	hermit -- Two-point Hermite cubic interpolation routine.

c FILE
c	hermit.f

c SYNOPSIS
c	A simple two-point Hermitian cubic interpolation routine.

c DESCRIPTION
c	Subroutine.  Perform a Hermite cubic interpolation of function y(x) 
c	bewteen two sample points.

c	---- On entry ----
c	x1, x2:		Sample values of independent variable
c	y1, y2:		Values of function at x1 and x2, respectively
c	yp1, yp2:	Values of derivative of function at x1 and x2
c	x0:		Value of independent variable for interpolation

c	---- On return ----
c	y0:		Interpolated value of function at x0
c	yp0:		Interpolated value of derivative at x0

c DIAGNOSTICS
c

c FILES
c

c NOTES
c

c SEE ALSO
c

c AUTHOR
c


      subroutine hermit (x1, x2, y1, y2, yp1, yp2, x0, y0, yp0)



c     ---- On entry ----

      real*4   x0, x1, x2, y1, y2, yp1, yp2

c     ---- On return ----

      real*4   y0, yp0

c     ---- Internal variables ----

      real*4   a, b, c, d, df, dx, f1, f2, fp1, fp2, sfp, t


      dx = x2 - x1
      t  = (x0 - x1) / dx
      if (t.le.0.5) then
         f1  = y1
         f2  = y2
         fp1 = yp1
         fp2 = yp2
      else
         t   = 1.0 - t
         dx  = -dx
         f1  = y2
         f2  = y1
         fp1 = yp2
         fp2 = yp1
      end if

      fp1 = fp1*dx
      fp2 = fp2*dx
      df  = f2 - f1
      sfp = fp1 + fp2
      a   = f1
      b   = fp1
      c   = 3.0*df - sfp - fp1
      d   = -2.0*df + sfp

      y0  = ((d*t + c)*t + b)*t + a
      yp0 = ( (3.0*d*t + 2.0*c)*t + b ) / dx

      return
      end


c NAME
c	holint -- Monotone, quadratic interpolation routine.

c FILE
c	holint.f

c SYNOPSIS
c	Perform a monotone, quadratic interpolation, even when holes in the
c	data are present.

c DESCRIPTION
c	Subroutine.  Monotone, quadratic interpolation of function f(x) 
c	which might have holes (bad values).   Bad function samples are 
c	given the value fbad.

c	---- Indexing ----
c	i = 1, n;

c	---- On entry ----
c	n:	Number of function samples
c	x(i):	Sample values of independent variable; Must be ordered: 
c		x(i) >= x(i-1)
c	f(i):	Value of function at x(i)
c	fbad:	Function value denoting bad sample
c	x0:	Value of independent variable for interpolation

c	---- On return ----
c	f0:	Interpolated value of function at x0
c	fp0:	Interpolated value of derivative at x0
c	iext:	Flag indicating whether extrapolation has occurred;
c		  =  0,	No extrapolation
c		  = -1,	Yes, x0 < x(1)
c		  = +1,	Yes, x0 > x(N)
c	ibad:	Flag indicating whether interopolation point is in a hole;
c		  =  0,	No
c		  =  1,	Yes

c	---- Subroutines called ----
c	Local
c		- Calls brack, fixhol and quaint, directly
c		- Calls brack and hermit, indirectly
c

c DIAGNOSTICS
c

c FILES
c

c NOTES
c	- f(x) may be discontinuous.  A discontinuity is presumed to occur
c	  when x(i) repeats for consecutive i.
c	- If x0 is out of range (iext = -1 or +1), then f0 and fp0 are 
c	  defined through linear extrapolation of function.
c	- If f(i) = fbad, then the function is assumed to be fbad in the
c	  intervals on either side of x(i), and then jump discontinuously
c	  to the next good sample value.  See subroutine fixhol for details
c	  on how holes are defined.
c	- If x0 is in a hole (ibad = 1), then f0 is returned fbad and fp0 
c	  is zero.

c SEE ALSO
c

c AUTHOR
c


      subroutine holint (n, x, f, fbad, x0, f0, fp0, iext, ibad)



c     ---- On entry ----

      integer*4 n
      real*4    f(n), fbad, x(n), x0

c     ---- On return ----

      integer*4 ibad, iext
      real*4    f0, fp0

c     ---- Internal variables ----

      integer*4 ileft, imax, imin, nh, nuse
      real*4    fh(6), xh(6)


c     Find four relevant samples and then construct a version of this
c     piece of the function with holes fixed

      call brack (n, x, x0, ileft)
      imin = max(1,ileft-1)
      imax = min(n,ileft+2)
      nuse = imax - imin + 1
      call fixhol (nuse, x(imin), f(imin), fbad, nh, xh, fh)

c     Interpolate fixed function

      if (nh.le.1) then
         f0  = fh(1)
         fp0 = 0.0
      else
         call quaint (nh, xh, fh, x0, f0, fp0, iext)
      endif

c     Now check if interpolation point is in a hole

      if (f0.eq.fbad .and. fp0.eq.0.0) then
         ibad = 1
      else
         ibad = 0
      end if

      return
      end


c NAME
c	quaint -- Monotone, quadratic interpolation with linear derivatives.

c FILE
c	quaint.f

c SYNOPSIS
c	Constrain derivative to be linear during monotone, quadratic 
c	interpolation.

c DESCRIPTION
c	Subroutine.  Perform monotone, quadratic interpolation of function 
c	f(x).  The interpolating function between two points is montone in 
c	value and linear in derivative.

c	---- Indexing ----
c	i = 1, n;

c	---- On entry ----
c	n:	Number of function samples
c	x(i):	Sample values of independent variable; must be ordered: 
c		x(i) >= x(i-1)
c	f(i):	Value of function at x(i)
c	x0:	Value of independent variable for interpolation

c	---- On return ----
c	f0:	Interpolated value of function at x0
c	fp0:	Interpolated value of derivative at x0
c	iext:	Flag indicating whether extrapolation has occurred;
c		  =  0,	No extrapolation
c		  = -1,	Yes, x0 < x(1)
c		  = +1,	Yes, x0 > x(N)

c	---- Subroutines called ----
c	Local
c		- Calls brack and hermit

c DIAGNOSTICS
c

c FILES
c

c NOTES
c	- f(x) may be discontinuous.  A discontinuity is presumed to occur
c	  when x(i) repeats for consecutive i.
c	- If x0 is out of range (iext = -1 or +1), then f0 and fp0 are 
c	  defined through linear extrapolation of function.

c SEE ALSO
c

c AUTHOR
c


      subroutine quaint (n, x, f, x0, f0, fp0, iext)



c     ---- On entry ----

      integer*4 n
      real*4    f(n), x(n), x0

c     ---- On return ----

      integer*4 iext
      real*4    f0, fp0

c     ---- Internal variables ----

      integer*4 i1, i2, i3, i4, ileft
      real*4    f1, f2, f3, f4, fac, fp2, fp3, fpdev, fpdev2, fpdev3
      real*4    h12, h23, h34, s12, s23, s34, x1, x2, x3, x4


c     Binary search for samples bounding x0

      call brack (n, x, x0, ileft)

c     x0 < x(1)

      if (ileft.lt.1) then
         if (x(2).gt.x(1)) then
            fp0 = (f(2) - f(1)) / (x(2) - x(1))
         else
            fp0 = 0.0
         end if
         f0   = f(1) + fp0*(x0-x(1))
         iext = -1
         return
      end if

c     x0 > x(n)

      if (ileft.ge.n) then
         if (x(n).gt.x(n-1)) then
            fp0 = (f(n) - f(n-1)) / (x(n) - x(n-1))
         else 
            fp0 = 0.0
         end if
         f0 = f(n) + fp0*(x0-x(n))
         iext = +1
         return
      end if

c     Normal case

c     Define points 1..4, such that x1 <= x2 <= x0 <= x3 <= x4 ----
c     If necessary, make x1 = x2 or x3 = x4

      i1 = max(1,ileft-1)
      i2 = ileft
      i3 = ileft + 1
      i4 = min(n,ileft+2)

      x1 = x(i1)
      x2 = x(i2)
      x3 = x(i3)
      x4 = x(i4)

      f1 = f(i1)
      f2 = f(i2)
      f3 = f(i3)
      f4 = f(i4)

c     Find widths of three intervals
c     Note 'brack' guarantees x(ileft) < x(ileft+1), and thus h23 > 0

      h12 = x2 - x1
      h23 = x3 - x2
      h34 = x4 - x3

c     Set finite-difference derivative in center interval

      s23 = (f3 - f2) / h23

c     Assign a function derivative to point 2; call it fp2.  The derivative 
c     of the parabola fitting points 1, 2 and 3 c (evaluated at x2) is used,  
c     howvever, if h12 is zero, s23 is used.

      if (h12.gt.0.0) then
         s12 = (f2 - f1) / h12
         fp2 = (s23*h12 + s12*h23) / (h12 + h23)
      else
         fp2 = s23
      end if

c     Assign a function derivative to point 3; call it fp3.  The derivative 
c     of the parabola fitting points 2, 3 and 4 (evaluated at x3) is used,  
c     howvever, if h34 is zero, s23 is used.

      if (h34.gt.0.0) then
         s34 = (f4 - f3) / h34
         fp3 = (s23*h34 + s34*h23) / (h34 + h23)
      else
         fp3 = s23
      end if

c     Adjust fp2 and fp3 such that they average to s23, but neither gets 
c     farther from s23

      fpdev2 = s23 - fp2
      fpdev3 = fp3 - s23
      if (fpdev2*fpdev3.le.0.0) then
         fpdev = 0.0
      else if (fpdev2.lt.0.0) then
         fpdev = -min(-fpdev2,-fpdev3)
      else
         fpdev = min(fpdev2,fpdev3)
      end if

c     Adjust derivatives such that Hermite cubic interpolant is monotonic

      if (s23.ne.0.0) then
         fac = abs(fpdev/s23)
         if (fac.gt.1.0) fpdev = fpdev/fac
      end if
      fp2 = s23 - fpdev
      fp3 = s23 + fpdev

c     Now do a straight Hermite cubic interpolation bewteen points 2 and 3

      call hermit (x2, x3, f2, f3, fp2, fp3, x0, f0, fp0)

      iext = 0

      return
      end

      subroutine holin2 (m,n,x,y,f,ldf,fbad,x0,y0,
     .                   f0,fx0,fy0,fxy0,iext,jext,ibad)
c
c Monotone, quadratic interpolation of function f(x,y) which
c  might have holes (bad values).
c Bad function samples are given the value FBAD.
c
c Calls BRACK, HOLINT and QUAINT directly.
c Calls BRACK, FIXHOL and HERMIT indirectly.
c 
c Input
c
c   M   :   number of x samples.
c   N   :   number of y samples.
c   X(i), i=1..M :  sample values of x.
c   Y(j), j=1..N :  sample values of y.
c   F(i,j), i=1..M, j=1..N :  value of function at (X(i),Y(j)).
c   LDF :  leading dimension of array F.
c   FBAD:  function value denoting bad sample.
c   X0  :  value of x for interpolation.
c   Y0  :  value of y for interpolation.
c
c Output
c
c   F0   :  interpolated value of function at (X0,Y0).
c   FX0  :  interpolated value of x-derivative of function at (X0,Y0).
c   FY0  :  interpolated value of y-derivative of function at (X0,Y0).
c   FXY0 :  interpolated value of x-y-derivative of function at (X0,Y0).
c   IEXT :  error flag;  0, no error;  -1, X0 < X(1);  1, X0 > X(M).
c   JEXT :  error flag;  0, no error;  -1, Y0 < Y(1);  1, Y0 > Y(N).
c   IBAD :  flag indicating whether interopolation point is in a hole;
c           0, no ... 1, yes
c
c Notes
c
c   1. If IBAD = 1, F0 is set to FBAD and the derivatrives to zero.
c   2. If X0 or Y0 is out of range (IEXT != 0 or JEXT != 0) then F0, FX0,
c       FY0 and FXY0 are defined through linear extrapolation of function.
c
      dimension x(m),y(n),f(ldf,n)
      dimension f0s(4),fx0s(4)
c
c Bracket X0.  Find 4 relevant x samples, or as many as needed.
c
      call brack (m,x,x0, ileft)
      imin = max(1,ileft-1)
      imax = min(m,ileft+2)
      muse = imax - imin + 1
c
c Do the same for y.
c
      call brack (n,y,y0, jleft)
      jmin = max(1,jleft-1)
      jmax = min(n,jleft+2)
      nuse = jmax - jmin + 1
c
c Now interpolate to (X0,Y(j)), j=jmin,jmax).
c
      do 100 j = jmin,jmax
      js = j - jmin + 1
100   call holint (muse,x(imin),f(imin,j),fbad,x0,
     .             f0s(js),fx0s(js),iext,ibad)
c
c Now interpolate to (X0,Y0).
c
      call holint (nuse,y(jmin),f0s,fbad,y0, f0,fy0,jext,ibad)
      if (ibad.gt.0) then
	fx0 = 0.0
	fxy0 = 0.0
      else
        call quaint (nuse,y(jmin),fx0s,y0, fx0,fxy0,jext)
      end if
c
c Find minimum interpolated x-derivative.
c
      return
      end
