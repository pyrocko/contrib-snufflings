      subroutine hyposat_cross(slat1,slon1,azi1,azi1s,
     +			   azi2,azi2s,del3,ep1,ep2,
     +			   elat,elats,elon,elons,typctl,ierr)
c
c****6789012345678901234567890123456789012345678901234567890123456789012
c
c     subroutine hyposat_cross locates an epicenter of a seismic event
c     with the knowledge of 2 stations and the backazimuth
c     (i.g. angle from the station to the event measured from
c     north). 
c
c     input:  slat1, slon1   geocentric coordinates of station 1
c
c             azi1, azi2     backazimuth at station 1 and station 2
c
c             azi1s, azi2s   corresponding standard deviations
c
c             del3           distance between the 2 stations
c
c             ep1, ep2       the 2 angles from one station to the other
c
c             typctl       verbosity level
c
c     output: elat, elon   geocentric coordinates of the event
c
c             elats, elons corresponding standard deviations
c
c
c version:  19. August 1996,  johannes schweitzer
c
c   last change:   Nov 10, 1997 (correct for source at lat +/- 90 deg)
c                  Nov 22, 2000 non-crossing directions removed
c                               error bars for lat/lon corrected.
c

      implicit double precision (a-h,o-z)

      integer typctl,ierr

c
c     Exclude combinations which give no crossing point because
c     azimuth directions point in different hemispheres.
c
      if(azi1.eq.ep1 .or. azi1.eq.alpha2(ep1-180.d0)) go to 9999
      if(ep1.eq.alpha2(azi1-180.d0)) go to 9999

      if(azi2.eq.ep2 .or. azi2.eq.alpha2(ep2-180.d0)) go to 9999
      if(ep2.eq.alpha2(azi2-180.d0)) go to 9999

c
c     definition of several constants:
c
      ierr    = 0

      pi      = 4.d0*datan(1.d0) 
      pi2     = 2.d0*pi
      pih     = pi / 2.d0
      pi9     = pi*0.9444444444d0
      deg2rad = pi/180.d0
      rad2deg = 180.d0/pi
c
c     convert angles to radian
c 
      dla1    = pih - deg2rad*slat1
      slon1r  = deg2rad*slon1

      azi1r   = deg2rad*azi1
      azi2r   = deg2rad*azi2

      del3    = deg2rad*del3
      d3c     = dcos(del3)
      ep1r    = deg2rad*ep1
      ep2r    = deg2rad*ep2

      if(azi1r.ge.ep1r) then
	ga1 = azi1r-ep1r
	ri1 = 1.d0
      else 
	ga1 = ep1r-azi1r
	ri1 = -1.d0
      endif

      if(ga1.lt.0.d0) ga1=ga1 + pi2

      if(ga1.ge.pi)   then
	 ga1=pi2 - ga1
	 ri1=-ri1
      endif

      if(azi2r.ge.ep2r) then
	ga2 = azi2r-ep2r
	ri2 = 1.d0
      else 
	ga2 = ep2r-azi2r
	ri2 = -1.d0
      endif
      if(ga2.lt.0.d0) ga2=ga2 + pi2
      if(ga2.ge.pi)   then
	ga2=pi2 - ga2
	ri2=-ri2
      endif

      if(azi1r.le.pi) then
        al1  = azi1r
  	ra1  = 1.d0
      else
        al1  = pi2 - azi1r
  	ra1  = -1.d0
      endif

      g1s = dsin(ga1)
      g1c = dcos(ga1)

      g2s = dsin(ga2)
      g2c = dcos(ga2)

      g3c = g1s*g2s*d3c-g1c*g2c
      ga3 = f2(g3c,2)
      g3s = dsin(ga3)

      d1 = f2((g2c+g1c*g3c)/(g1s*g3s),2)

      d2 = f2((g1c+g2c*g3c)/(g2s*g3s),2)

c      if(typctl.gt.8) then
c	 print 1000,'ga1,ga2,ga3,al1 :',ga1*rad2deg,ga2*rad2deg,
c     +	          ga3*rad2deg,al1*rad2deg,d1,d2
c1000  format (a,6f9.4)
c      endif
      
      if(d1.gt.pi9 .or. d2.gt.pi9) then
	 ierr = 1
         if(typctl.gt.4) then
	    print *,'distance d1 or d2 .gt. 170 deg'
	 endif
	 return
      endif

      if(d1.lt.0.d0 .or. d2.lt.0.d0) then
	 ierr = 2
         if(typctl.gt.4) then
	    print *,'distance d1 or d2 .lt. 0 deg'
	 endif
	 return
      endif

      if(typctl.gt.5) then
	 print *,'event distance [deg] from station 1: ',d1*rad2deg
	 print *,'event distance [deg] from station 2: ',d2*rad2deg
      endif

      p1   = dcos(dla1)
      p2   = dcos(d1)

      p3   = dsin(dla1)
      p4   = dsin(d1)

      p5   = dcos(al1)
      p6   = dsin(al1)

c
c     calculating event cos of lat and lon and
c     preparing the output of the event data
c

      elax  = p1*p2 + p3*p4*p5
      elatr = f2(elax,2)

      p7    = dsin (elatr)
c    
c     event latitude
c
      elat = (pih - elatr)*rad2deg

      if(dabs(elat).lt.90.d0) then
         elox  = (p2 - p1*elax) / (p3 * p7)
      else
	 elox  = 1.d0
      endif

      eloxr = f2(elox,2)

      if(azi1r .lt. pi) then
        elonr = slon1r + eloxr
      else if(azi1r .eq. pi) then
	elonr = slon1r
      else if(azi1r .gt. pi) then
	elonr = slon1r - eloxr
      endif

      if(elonr.lt.-pi) elonr = elonr + pi2
      if(elonr.gt.pi)  elonr = elonr - pi2
c    
c     event longitude
c
      elon = elonr*rad2deg

c
c     finally we have to calculate the error bars
c

      dga3daz1  = -(g1s*g2c + g1c*g2s*d3c)/g3s
      dga3daz2  = -(g1c*g2s + g1s*g2c*d3c)/g3s

      dd1daz1   = ri1*(g3s*(g3c + g1c*g2c) + dga3daz1*g1s*(g1c+g2c*g3c))
      dd1daz1   = dd1daz1 / (p4*g1s*g1s*g3s*g3s)

      dd1daz2   = ri2*(g1s*g2s*g3s + dga3daz2*g1s*(g1c+g3c*g2c))
      dd1daz2   = dd1daz2 / (p4*g1s*g1s*g3s*g3s)

      if(p7.ne.0.d0) then

         delatdaz1 = ((p1*p4 - p2*p3*p5)*dd1daz1 + p3*p4*p6*ra1) / p7
         delatdaz2 = ((p1*p4 - p2*p3*p5)*dd1daz2               ) / p7

         f3 = (p2 * elax - p1)/p7
      else
	 delatdaz1 = 1.d0
	 delatdaz2 = 1.d0
	 f3 = 0.d0
      endif

      elats = dsqrt((delatdaz1*azi1s)**2.d0 + (delatdaz2*azi2s)**2.d0)

      f4 = dsin(eloxr)*p3*p7

      if(dabs(f4).gt.0.00005d0) then
         delondaz1 = (p4*dd1daz1 + delatdaz1*f3) / f4
         delondaz2 = (p4*dd1daz2 + delatdaz2*f3) / f4
      else
	 delondaz1 = 1.d0
	 delondaz2 = 1.d0
      endif

      elons = dsqrt((delondaz1*azi1s)**2.d0 + (delondaz2*azi2s)**2.d0)

      if(dabs(elat).eq.90.d0) then
	 elons = 180.d0
      endif

      if(elons.gt.180.d0) then
 	 ierr = -3
	 elons=180.d0
      endif
      if(elats.gt.90.d0) then
 	 ierr = -ierr
	 elats=90.d0 
      endif

      if (typctl.gt.4) then
	 print*,'event latitude : ',elat,' +/- ',elats
	 print*,'event longitude: ',elon,' +/- ',elons
	 if(typctl.gt.4) print*,'ierr = ',ierr
      endif

      return

9999  if (typctl.gt.8) then
	  print*,'no crossing point posible'
      endif
      ierr = 10

      return
      end
