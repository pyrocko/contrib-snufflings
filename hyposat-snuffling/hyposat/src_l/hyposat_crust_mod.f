c
c  subroutine version of a SCRIPS program called 'getCNpoint'
c  to get the parameters of the crust as descriped in:
c  Mooney, Laske and Masters, Crust 5.1: a global
c  crustal model at 5x5 degrees, JGR, January 1998
c
c  Because we don't need any densities to locate an event,
c  this parameter was thrown out from the original data and program
c  files.
c
c     NORSAR, Johannes Schweitzer, March 1999
c
c
c layer one and two flipped, after the read statement!
c layer 1: water
c layer 2: ice
c

      subroutine get_mod_mlm(itrue,inum,typctl,ierr)

      save

      integer          typctl,ierr,trimle,itrue,inum
      integer          jmod,iread,imo

      parameter (ml=61)

      double precision rlat,rlon,elev,v0(2,ml),z(ml),
     *                 x(4),y(4),u(4),vp(2,ml),vs(2,ml),
     *                 dz(2,ml),bilinear,dzw(2),dcolo,dcola,
     *                 zmax,ele(2),elon,elat,elev2,dum,dum2,fac
      character        mtyp*3,azo(ml)*4,filloc*80,mpath*80,file*120


      COMMON /MODEL/  v0,z,elev,rlat,rlon,zmax,elat,elon,elev2,
     *                jmod,iread,imo,mtyp

      COMMON /MODELC/ filloc,azo,mpath

      parameter(ityp=145)
      dimension fvel(ityp,8),fvels(ityp,8)
      dimension fthi(ityp,8)
      dimension ampvp(8,72,36),
     +          ampvs(8,72,36),ampthi(7,72,36),
     +          ampele(72,36)
      character*2 ctype(ityp),line*506,cdum*1
      character*2 types(72),atype(72,36)
      character*12 names(7)
      character*80 file_check
      data names/'water','ice','soft sed.','hard sed.',
     +         'upper crust','middle crust','lower crust'/

      ierr = 0

      if(iread.eq.0) then 

	 file = mpath(1:trimle(mpath)) // 'CNtype_key'
         open(65,file=file,err=10,status='old')
	 go to 11

10       file = file_check('CNtype_key')
         open(65,file=file)

11       file = mpath(1:trimle(mpath)) // 'CNtype'
         open(66,file=file,err=12,status='old')
	 go to 13

12       file = file_check('CNtype')
         open(66,file=file)

13	 file = mpath(1:trimle(mpath)) // 'CNelevatio'
         open(67,file=file,err=14,status='old')
	 go to 15

14       file = file_check('CNelevatio')
         open(67,file=file)

15       continue

c... read in key for crust types
c...............................
         read(65,'(///a)')cdum
         print*,' ... reading global crustal model file ...'
         do 101 i=1,ityp
            read(65,'(a)')ctype(i)
            read(65,*)(fvel(i,l),l=1,8)
            read(65,*)(fvels(i,l),l=1,8)
            read(65,*)(fthi(i,l),l=1,7)
c	    print *,i,ctype(i)
c flip layers
            aux=fvel(i,1)
            fvel(i,1)=fvel(i,2)
            fvel(i,2)=aux
            aux=fvels(i,1)
            fvels(i,1)=fvels(i,2)
            fvels(i,2)=aux
            aux=fthi(i,1)
            fthi(i,1)=fthi(i,2)
            fthi(i,2)=aux
 101     continue

c... read CNtype file
c...............................
         read(66,*)cdum
         read(67,'(a)')line
         do 45 j=1,36
            read(67,*)ilat,(ampele(i,j),i=1,72)
            read(66,'(a)')line       
            read(line,*)ilat,types
c           print*,ilat
            do 40 i=1,72
               do 30 l=1,ityp
               if(types(i).eq.ctype(l))then
                 atype(i,j)=ctype(l)
                 do 20 k=1,8
                 ampvp(k,i,j)=fvel(l,k)
                 ampvs(k,i,j)=fvels(l,k)
 20              continue
                 do 21 k=1,7
 21              ampthi(k,i,j)=fthi(l,k)
                 goto 40
               endif
 30            continue
               print*,' crust type code not found: ',types(i)
               print*,' latitude: ',ilat,' long index: ',i
 40         continue
 45      continue

	 iread = 1
	 close (65)
	 close (66)
	 close (67)
      endif

*-------------------
c     now look up coordinates

      do 80 k = 1,inum

      if(mtyp.ne.'MLM') then

	 if(k.gt.1) go to 900
         do 200 i = 1,ityp

	 if (mtyp(1:2).eq.ctype(i)) then
	    imod = i
	    elev = 0.d0
	    go to 202
	 endif
200      continue
         print *,'Could not find standard model type ',mtyp
	 go to 999

      else

	 if(k.eq.1) then
	    dcola = 90.d0 - rlat
            dcolo = rlon+180.d0
         else if(k.eq.2) then
	    dcola = 90.d0 - elat
            dcolo = elon+180.d0
         endif

	 y(1) = dnint(dcola/5.d0)*5.d0 - 2.5d0
	 y(2) = y(1)
	 y(3) = y(1) + 5.d0
	 y(4) = y(3)

	 x(1) = dnint(dcolo/5.d0)*5.d0 - 2.5d0
	 x(2) = x(1) + 5.d0
	 x(3) = x(2)
	 x(4) = x(1)

         ilat1 = int(real(y(1)/5.d0))+1
	 if(ilat1.le.0) ilat1 = 1
	 if(ilat1.gt.36) ilat1 = 36
	 ilat2 = ilat1 + 1
	 if(ilat2.gt.36) ilat2 = 36
	 ilat3 = ilat2
	 ilat4 = ilat1

	 ilon1 = int(real(x(1)/5.d0))+1
	 if(ilon1.le.0) ilon1 = 72 + ilon1
	 if(ilon1.gt.72) ilon1 = ilon1 - 71
	 ilon2 = ilon1 + 1
	 if(ilon2.gt.72) ilon2 = ilon2 - 71
	 ilon3 = ilon1
	 ilon4 = ilon2

         u(1) = dble(ampele(ilon1,ilat1))
         u(2) = dble(ampele(ilon2,ilat2))
         u(3) = dble(ampele(ilon3,ilat3))
         u(4) = dble(ampele(ilon4,ilat4))

	 indbi = 0
	 ele(k) = bilinear(x,y,u,dcolo,dcola,indbi)/1000.d0
	 indbi = 1

      endif

202   jmod = 1

      iconr = 0
      do 70 i=1,8

      if (mtyp.eq.'MLM') then
         u(1)  = dble(ampvp(i,ilon1,ilat1))
         u(2)  = dble(ampvp(i,ilon2,ilat2))
         u(3)  = dble(ampvp(i,ilon3,ilat3))
         u(4)  = dble(ampvp(i,ilon4,ilat4))
	 vp(k,i) = bilinear(x,y,u,dcolo,dcola,indbi)

         u(1)  = dble(ampvs(i,ilon1,ilat1))
         u(2)  = dble(ampvs(i,ilon2,ilat2))
         u(3)  = dble(ampvs(i,ilon3,ilat3))
         u(4)  = dble(ampvs(i,ilon4,ilat4))
	 vs(k,i) = bilinear(x,y,u,dcolo,dcola,indbi)

	 if(i.le.7) then
            u(1)  = dble(ampthi(i,ilon1,ilat1))
            u(2)  = dble(ampthi(i,ilon2,ilat2))
            u(3)  = dble(ampthi(i,ilon3,ilat3))
            u(4)  = dble(ampthi(i,ilon4,ilat4))
	    dz(k,i) = bilinear(x,y,u,dcolo,dcola,indbi)
	 else
	    dz(k,i) = 20.d0
	 endif

      else

	 vp(k,i) = dble(fvel(imod,i))
	 vs(k,i) = dble(fvels(imod,i))
	 if(i.le.7) then
	    dz(k,i) = dble(fthi(imod,i))
	 else
	    dz(k,i) = 20.d0
	 endif

      endif
c
c     negative elevation: water layer at top!
c
      if(i.eq.1) then
         if(ele(k).lt.0.d0.or.dz(k,i).ne.0.d0) then
	    dzw(k)=-ele(k)
	    dz(k,i) = 0.d0
         endif
      endif

      if(k.eq.1) elev  = ele(k)
      if(k.eq.2) elev2 = ele(k)

      if(typctl.gt.8) then
         if(i.eq.8) then
          if(mtyp.eq.'MLM') then
             print *,' Mooney, Laske and Masters, CRUST 5.1:'
	     if(k.eq.1) then
             print 792,'type, latitude, longitude, elevation: ',
     +         atype(ilon4,ilat4),rlat,rlon,elev,inum
          else if(k.eq.2) then
             print 792,'type, latitude, longitude, elevation: ',
     +         atype(ilon4,ilat4),elat,elon,elev2,inum
          endif
792   format(a,2x,a2,2x,3f11.4,i3)
         else
            print *,'Standard Earth model : ',mtyp
         endif
         print 793,'Mantle below Moho: ave. vp, vs:  ',
     +        vp(k,i),vs(k,i)
 793  format(a,2x,2f11.4)
      endif

	 if (i.eq.1) then
            print *,' ' 
            print *,' 7-layer crustal model (thickness, vp,vs)'
	 endif
         if(i.lt.8) print 794, dz(k,i),vp(k,i),vs(k,i),names(i)
794      format (3f10.4,2x,a12)
      endif

      if(k.eq.inum) then

         if(inum.eq.1) then

            if(dz(k,i).ne.0.d0) then
	    
	       jmod1 = jmod +1
	       if(jmod.eq.1) then
	          z(jmod)    = 0.d0
	          if(itrue.eq.0) then
	             z(jmod1)   = dz(k,i) - ele(k)
	          else
	             z(jmod1)   = dz(k,i)
	          endif
               else
	          z(jmod)    = z(jmod-1)
	          z(jmod1)   = z(jmod)+dz(k,i)
               endif

	       v0(1,jmod) = vp(k,i)
	       v0(1,jmod1)= vp(k,i)

	       v0(2,jmod) = vs(k,i)
	       v0(2,jmod1)= vs(k,i)

	       azo(jmod)  = ' '
	       azo(jmod1) = ' '

	       if(v0(1,jmod).ge.6.3d0 .and.iconr.eq.0) then
	          azo(jmod-1) = 'CONR'
	          iconr = 1
	       endif

	       if(i.eq.8) then
	          azo(jmod-1) = 'MOHO'
	          jmod = jmod1
	       else
	          jmod = jmod1+1
	       endif

            endif

         else if(inum.eq.2) then

            if(dz(1,i).ne.0.d0 .or. dz(2,i).ne.0.d0) then

	       jmod1 = jmod +1
	       if(jmod.eq.1) then
	          z(jmod)    = 0.d0
	          if(itrue.eq.0) then
	             z(jmod1)   = (dz(1,i)-ele(1)+dz(2,i)-ele(2))/2.d0
	          else
	             z(jmod1)   = (dz(1,i) + dz(2,i) ) / 2.d0
	          endif
               else
	          z(jmod)    = z(jmod-1)
	          z(jmod1)   = z(jmod)+(dz(1,i)+dz(2,i))/2.d0
               endif

	       fac  = 0.d0
	       dum  = 0.d0
	       dum2 = 0.d0

               do 59 i2 =1,2
	          if(dz(i2,i).ne.0.d0) then
		     fac  = fac  + 1.d0
		     dum  = dum  + vp(i2,i)
		     dum2 = dum2 + vs(i2,i)
	          endif
59             continue

	       v0(1,jmod) = dum  /fac
	       v0(1,jmod1)= v0(1,jmod)
	       v0(2,jmod) = dum2 /fac
	       v0(2,jmod1)= v0(2,jmod)

	       azo(jmod)  = ' '
	       azo(jmod1) = ' '

	       if(v0(1,jmod).ge.6.3d0 .and.iconr.eq.0) then
	          azo(jmod-1) = 'CONR'
	          iconr = 1
	       endif

	       if(i.eq.8) then
	          azo(jmod-1) = 'MOHO'
	          jmod = jmod1
	       else
	          jmod = jmod1+1
	       endif

            endif

         endif

      endif
	    
70    continue

80    continue

      zmax = z(jmod)

      if(typctl.gt.8) then
	 do 90 i=1,jmod
	 print*,i,z(i),v0(1,i),v0(2,i),azo(i)
90       continue
      endif

900   return

999   print *,'Something wrong with CRUST 5.1 - input files'
      ierr = 99
      return

      end
