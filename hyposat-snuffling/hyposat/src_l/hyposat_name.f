      function phase_type(phid0)
c
c           last corrections :  16 February 1997
c                               correcting range of characters to check
c
c                               November 3, 1997
c                               corrected for AB and AC-branches
c
c                               October 8, 2002
c                               corrected for P1 and S1 phase names
c                               and new IASPEI phase names like P4 ...
c
      character phid0*8, phase_type*1, phid*8
      integer icph, trimle

	 phid = phid0
	 phase_type=' '

         icph = trimle(phid)
	 if(icph.le.1) go to 10

 	 if(phid(icph:icph).eq.'1') then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif
	 if(phid(icph:icph).eq.'2') then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif
	 if(phid(icph:icph).eq.'3') then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif
	 if(phid(icph:icph).eq.'4') then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif
	 if(phid(icph:icph).eq.'5') then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif
	 if(phid(icph:icph).eq.'6') then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif
	 if(phid(icph:icph).eq.'7') then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif
	 if(phid(icph:icph).eq.'8') then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif
	 if(phid(icph:icph).eq.'9') then
      	    icph=icph-1
	 endif

5	 if(icph.le.1) go to 10

 	 if(icph.ge.4) then
	    if(phid(icph-2:icph).eq.'dif'.or.phid(icph-2:icph).eq.
     +	       'DIF') then
               icph=icph-3
	       go to 5
	    endif
	    if(phid(icph-2:icph).eq.'pre'.or.phid(icph-2:icph).eq.
     +	       'PRE') then
               icph=icph-3
	       go to 5
	    endif
         endif

	 if(icph.ge.3) then

	    if(phid(icph-1:icph).eq.'ab'.or.
     +         phid(icph-1:icph).eq.'AB') then
      	       icph=icph-2
	       go to 5
	    endif
	    if(phid(icph-1:icph).eq.'ac'.or.
     +         phid(icph-1:icph).eq.'AC') then
      	       icph=icph-2
	       go to 5
	    endif
	    if(phid(icph-1:icph).eq.'bc'.or.
     +         phid(icph-1:icph).eq.'BC') then
      	       icph=icph-2
	       go to 5
	    endif
	    if(phid(icph-1:icph).eq.'df'.or.
     +         phid(icph-1:icph).eq.'DF') then
      	       icph=icph-2
	       go to 5
	    endif

	 endif

	 if(phid(icph:icph).eq.'n'.or.phid(icph:icph).eq.'N') then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif
	 if(phid(icph:icph).eq.'g'.or.phid(icph:icph).eq.'G') then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif
	 if(phid(icph:icph).eq.'b'.or.phid(icph:icph).eq.'B') then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif

	 if(phid(icph:icph).eq."'") then
      	    icph=icph-1
	    if(icph.le.1) go to 10
	 endif
	 if(phid(icph:icph).eq.'*') then
	    icph=icph-1
	 endif

10       if(icph.le.1) icph=1

	 if(phid(icph:icph).eq.'P') phase_type='P'
	 if(phid(icph:icph).eq.'S') phase_type='S'

      return
      end
c
c     subroutine testphase (phid0,icha)
c
c     some changes added for surface-reflections
c
c     3 August 1997, Johannes Schweitzer, NORSAR
c
c     13 February 2002: PKPab + PKPdif added
c  
      subroutine testphase (phid0,icha)
      integer icha
      character phid*8,phid0*8,s*1
      logical flag

      flag = .false.
      phid = phid0
      s    = ''

      if(phid0(1:1).eq.'p' .or. phid0(1:1).eq.'s') then
	 flag = .true.
	 s    = phid0(1:1)
	 phid = phid0(2:)
      endif
      if(phid(1:1).eq.'P' .and. icha.ge.18) go to 50
      if(phid(1:1).eq.'S' .and. icha.ge.15) go to 50

      if(phid(1:4).eq.'PKP ') then
	 phid='PKPdf'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:5).eq.'PKPdf') then
	 phid='PKPdif'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:6).eq.'PKPdif') then
	 phid='PKPbc'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:5).eq.'PKPbc') then
	 phid='PKPab'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:5).eq.'PKPab') then
	 phid='Pdif'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:4).eq.'Pdif') then
	 phid='P'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'P') then
	 phid='Pn'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'PmP') then
	 phid='Pn'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'Pn') then
	 phid='Pb'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'PbP') then
	 phid='Pb'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'Pb') then
	 phid='Pg'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'Pg') then
	 phid='PKPdf'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'P1') then
	 phid='Pg'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:4).eq."P'P'ab") then
	 phid="P'P'bc"
         icha = icha+1
	 goto 100
      endif
      if(phid(1:4).eq."P'P'bc") then
	 phid="P'P'df"
         icha = icha+1
	 goto 100
      endif
      if(phid(1:4).eq.'PKKP  ') then
	 phid='PKKPab'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:4).eq.'PKKPab') then
	 phid='PKKPbc'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:4).eq.'PKKPbc') then
	 phid='PKKPdf'
         icha = icha+1
	 goto 100
      endif

      if(phid(1:4).eq.'SKS ') then
         phid='SKSdf'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:5).eq.'SKSdf') then
         phid='SKSac'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:5).eq.'SKSac') then
         phid='Sdif'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:4).eq.'Sdif') then
         phid='S'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'S') then
         phid='Sn'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'SmS') then
	 phid='Sn'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'Sn') then
         phid='Sb'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'SbS') then
	 phid='Sb'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'Sb') then
         phid='Sg'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'S1') then
         phid='Sg'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'Lg') then
         phid='Sg'
         icha = icha+1
	 goto 100
      endif
      if(phid.eq.'Rg') then
         phid='SKSac'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:4).eq."S'S'ac") then
	 phid="S'S'df"
         icha = icha+1
	 goto 100
      endif
      if(phid(1:4).eq.'SKKSac') then
	 phid='SKKSdf'
         icha = icha+1
	 goto 100
      endif
      if(phid(1:4).eq.'Sg') then
	 phid='Lg'
         icha = icha+1
	 goto 100
      endif

50    icha = 999
      return

100   if(flag) then
	phid0= s // phid(1:7)
      else
	phid0 = phid
      endif

      return
      end
