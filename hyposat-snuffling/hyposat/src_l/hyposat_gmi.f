      subroutine hyposat_gmi(n,m,kq,ierr,typctl)

C
C
C     hyposat_gmi developed as a subroutine to do
C     generalized matrix inversion with single value decomposition
C     
C
C      last changes :  23. July    1998       JOHANNES SCHWEITZER
c
c                      17. October 2000 vector of error ellipse included
c
c                      10. October 2002 special case for n = m fixed
c
c
c     input:
c
c     a           matrix of partial derivatives
c
c     d           data vector
c
c     dst         vector of 'a priori' standard deviations of data
c
c     amst        vector with 'a priori' standard deviations of model
c                 parameter
c
c     n           number of data
c
c     m           number of model parameter
c
c     typctl      verbosity level
c
c           = 9, 10  the matrix of partial derivatives before and
c                 after weighting is printed
c
c
c           special output in an output file 'hyposat_gmi.out'
c
c           = 11  the resolution matrix is printed in the output file
c
c           = 12  the covariance matrix is printed in the output file
c
c           = 13  the correlation matrix is printed in the output file
c
c          => 13  all three are printed in the output file
c
c          => 20  additionally the diagonal elements of the 
c                 information-density matrix is printed
c
c          => 30  additionally the whole information-density matrix is 
c                 printed
c
c     output:
c
c     xg          solution vector of model parameter
c
c     xst         'a posteriori' standard deviations of model parameter
c
c     dg          vector of calculated data
c
c     ierr        if svd failed ierr = maximum of svd-iterations
c
c     a           contains the correlation matrix
c
c     Several other results are calculated internaly, but not yet given
c     back as output. See code and change respectivly!
c
      parameter (mm=4,nn=1600)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NN,MM),AT(NN,NN),SV(MM),U(NN,NN),V(MM,MM),XST(MM),
     *          D(NN),DST(NN),R(MM,MM),AMST(MM),DG(NN),
     *          EWMOD(MM),XG(MM),ax1(2),ax2(2)

      logical iellip

      common  /gmi/ a,xst,d,xg,dg,dst,amst,ax1,ax2,iellip
      common  /gm2/ u,sv,v

      integer ierr, typctl, unit

      double precision q2

      ierr = 0

      quot = 1.d-4

      if(typctl.gt.10) then
         unit = 35
         open (unit,file='hyposat_gmi.out')
      endif

C
C 'a priori' data variance, zero setting of xg,
C  and 'a priori' model variance
C
      DO 10 J=1,M

      xg(j) = 0.d0
      xst(j)= 0.d0

      DO 10 I=1,N

      U(I,J)=AMST(J)*A(I,J)/DST(I)

c     if(typctl.gt.8 .and. typctl.lt.11) then
c        print *,i,j,u(i,j),AMST(J),DST(I)
c     endif

10    continue

C
C single value decomposition (SVD)
C
      CALL SVDCMP(n,m,ierr)

      if(ierr.ge.60) then
	print *,'SVD: too much iterations ( ',ierr,' )'
	return
      endif

C
C modifying the singular values
C
      kq = m

      DO 50 K=1,m

      IF (SV(K).LT.QUOT) then
         SV(K)    = 1.d0
	 if (k.eq.kq) kq = kq-1
	 EWMOD(K)=0.d0
      else
         EWMOD(K)=1.d0/SV(K)
      endif

50    CONTINUE

C
C calucating the generalised inverse of A , stored in array AT
c (and zero setting of field dg )
C
        DO 110 I=1,N

	dg(i) = 0.d0

        DO 110 J=1,M

        AT(J,I)=0.d0

        DO 110 K=1,M
          AT(J,I)=AT(J,I)+V(J,K)*EWMOD(K)*U(I,K)
110     CONTINUE

C
C weighting with 'a priory' data variances
C and       with 'a priory' model variances
C
        DO 120 I=1,N
        DO 120 J=1,M
          AT(J,I)=AMST(J)*AT(J,I)/DST(I)
120     CONTINUE

C
C calculating the solution vector
C
        DO 130 J=1,M
        DO 130 I=1,N
          XG(J)=XG(J)+AT(J,I)*D(I)
130     continue

C
C calculating the resolution matrix
C
        DO 145 J2=1,M
        DO 145 J1=1,M

        R(J1,J2)=0.d0

        DO 140 K=1,m
          R(J1,J2)=R(J1,J2)+V(J1,K)*V(J2,K)
140     CONTINUE

        R(J1,J2)=AMST(J1)*R(J1,J2)/AMST(J2)

145     CONTINUE

 	if (typctl.eq.11 .or. typctl.gt.13) then

	   write(unit,'(/'' The Resolution Matrix:''/)')
	   do 150 j1=1,m
	   write(unit,'(5f12.9)') (r(j1,j2),j2=1,m)
150        continue

 	endif

C
C residual variances (mean value stored in rv)
C

	rv = 1.d0

        if(n.gt.m) then


          DO 160 J=1,M
          DO 160 I=1,N
            DG(I)=DG(I)+A(I,J)*XG(J)
160       CONTINUE

          RV0=0.d0

          DO 165 I=1,N
          EPS=(DG(I)-D(I))/DST(I)
          rv0=RV0+q2(EPS)
165       CONTINUE

	  if(rv0.gt.1.d-8) RV=RV0/DBLE(N-M)

	endif

C
C covariance matrix (temporarly stored in AT)
C
c	   write(unit,'(/'' The Matrix V:''/)')
c           do 169 j1=1,m
c 	     write(unit,'(4f12.8,2x,2f13.8)') 
c     +           (v(j1,j2),j2=1,m),sv(j1),EWMOD(j1)
c169        continue


        DO 175 J2=1,M
        DO 175 J1=1,M

          AT(J1,J2)=0.d0

          DO 170 K=1,m
          AT(J1,J2) = AT(J1,J2) +
     +                   V(J1,K)*q2(EWMOD(K))*V(J2,K)
170       continue

        AT(J1,J2)=AMST(J1)*AT(J1,J2)*AMST(J2)*RV

175     CONTINUE

	if (typctl.eq.12 .or. typctl.gt.13) then

	   write(unit,'(/'' The Covariance Matrix:''/)')
	   do 176 j1=1,m
	     write(unit,'(6f12.8)') 
     +           (AT(j1,j2),j2=1,m)
176        continue

	endif

C
C vector with standard deviations of parameters modelled
C
        DO 180 J=1,M
        XST(J)=DSQRT(AT(J,J))
180     continue

C
C correlation matrix  ( stored in A )
C
        DO 190 J2=1,M
        DO 190 J1=1,M
	a(j1,j2) = 0.d0
	if(XST(J1).eq.0.d0 .or. XST(J2).eq.0.d0) go to 190
        A(J1,J2)=AT(J1,J2)/(XST(J1)*XST(J2))
190     continue


	if (typctl.ge.13) then

	   write(unit,'(/'' The Correlation Matrix:''/)')
	   do 195 j1=1,m
	   write(unit,'(5f12.8)') (a(j1,j2),j2=1,m)
195        continue

	endif

c
c       Now we have to calculate the vector of the principal axes 
c       of the error ellipse of the epicenter. This only in case 
c       that we have enough enough resolution for at least 3 
c       parameters (to, lat, and lon), depth is optional.
c

	if (iellip .and. kq.ge.3 .and. rv.gt.0.d0) then

	 slat1 = 0.d0
	 slon1 = 0.d0

	 slat2 = 0.d0
	 slon2 = 0.d0

	 do 196 j = 1,m

	   if(ewmod(j).le.0.d0) go to 198

	   if(j.ne.3) then

 	      slat1 = slat1 + q2(ewmod(j)*v(2,j))
 	      slon1 = slon1 + q2(ewmod(j)*v(3,j))

	   endif

	   if(j.ne.2) then

 	      slat2 = slat2 + q2(ewmod(j)*v(2,j))
 	      slon2 = slon2 + q2(ewmod(j)*v(3,j))

	   endif

196      continue

	 f1    = sv(2)*xst(2)
	 slat1 = f1*dsqrt(slat1)
	 slat2 = f1*dsqrt(slat2)

         f1    = sv(3)*xst(3)
	 slon1 = f1*dsqrt(slon1)
	 slon2 = f1*dsqrt(slon2)

	 AX1(1) = dsign(slat1,v(2,2))
	 AX1(2) = dsign(slon1,v(3,2))

	 AX2(1) = dsign(slat2,v(2,3))
	 AX2(2) = dsign(slon2,v(3,3))

	 go to 199

        endif

198     iellip = .false.
C
C calculating the information density matrix (stored in AT )
C
199     continue
 	if (typctl.ge.20) then
           DO 210 I2=1,N
           DO 210 I1=1,N
           AT(I1,I2)=0.d0
           DO 200 K=1,m
           AT(I1,I2)=AT(I1,I2)+U(I1,K)*U(I2,K)
200        CONTINUE
           AT(I1,I2)=DST(I1)*AT(I1,I2)/DST(I2)
210        CONTINUE


	   write(unit,'(/'' The Diagonal Elements of the '',
     +	                ''Information-Density Matrix:''/)')
	   do 215 j1=1,n
	   write(unit,'(f12.8)') at(j1,j1)
215        continue

 	endif


 	if (typctl.ge.30) then

	   write(unit,'(/'' The Information-Density Matrix:''/)')
	   do 216 j1=1,n
	   write(unit,'(10f12.8)') (at(j1,j2),j2=1,n)
216        continue

 	endif

C
        if(typctl.gt.10) then
          close (unit)
	endif

        return
        END
