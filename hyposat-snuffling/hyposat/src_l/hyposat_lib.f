CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE indexx(n,arr,indx)
C  (C) Copr. 1986-92 Numerical Recipes Softwar
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SVDCMP(M,N,INDERR)
C  (C) Copr. 1986-92 Numerical Recipes Software
      PARAMETER (NMAX=1600,mm=4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RV1(NMAX)
      common  /gm2/ a(nmax,nmax),w(mm),v(mm,mm)
      G=0.0d0
      SCALE=0.0d0
      ANORM=0.0d0
      INDERR = 0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0d0
        S=0.0d0
        SCALE=0.0d0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+DABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0d0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-DSIGN(DSQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0d0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0d0
        S=0.0d0
        SCALE=0.0d0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+DABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0d0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-DSIGN(DSQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0d0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=DMAX1(ANORM,(DABS(W(I))+DABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0d0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0d0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0d0
            V(J,I)=0.0d0
31        CONTINUE
        ENDIF
        V(I,I)=1.0d0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0d0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0d0) THEN
          G=1.0d0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0d0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0d0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0d0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,60
          DO 41 L=K,1,-1
            NM=L-1
            IF ((DABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((DABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0d0
          S=1.0d0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((DABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=dpythag(f,g)
              W(I)=H
              H=1.0d0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0d0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.60) THEN
             INDERR=60
             RETURN
          ENDIF
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0d0*H*Y)
          G=dpythag(f,1.0d0)
          F=((X-Z)*(X+Z)+H*((Y/(F+DSIGN(G,F)))-H))/X
          C=1.0d0
          S=1.0d0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=dpythag(F,H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=dpythag(F,H)
            W(J)=Z
            IF (Z.NE.0.0d0) THEN
              Z=1.0d0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0d0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION dpythag(a,b)
C  (C) Copr. 1986-92 Numerical Recipes Software
      DOUBLE PRECISION a,b,dpythag
      DOUBLE PRECISION absa,absb,q2
      absa=dabs(a)
      absb=dabs(b)
      if(absa.gt.absb)then
        dpythag=absa*dsqrt(1.0d0+q2(absb/absa))
      else
        if(absb.eq.0.0d0)then
          dpythag=0.0d0
        else
          dpythag=absb*dsqrt(1.0d0+q2(absa/absb))
        endif
      endif
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function alpha1(a)
      double precision a,alpha1
      if(a.gt. 180.d0)  then
         alpha1 = a  - 360.d0
      else if(a.lt.-180.d0) then
         alpha1 = a  + 360.d0
      else
         alpha1 = a
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function alpha2(a)
      double precision a,alpha2
      if(a.gt. 360.d0)  then
         alpha2 = a  - 360.d0
      else if(a.lt.0.d0) then
         alpha2 = a  + 360.d0
      else
         alpha2 = a
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c       function to count the number of characters of a string
c       (from left).
c
      function trimle(t)
      integer trimle
      character t*(*)
           do 1 trimle=len(t),1,-1
              if (t(trimle:trimle).ne.' ') return
1          continue
      trimle=1
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        SUBROUTINE DLSQ(N,M) 
C 
C       method: least squares fit
c
c       here double precision version
C 
C       linear equation sytem: G*A=D (n equations for m unknowns, n>=M)
C       kernel matrix G (n times m), model A (m times 1),
C       data D (n times 1)
C       VAR = STANDARD DEVIATIONS (sqrt of variances) of A 
C       RES = RESIDUALS = G*A-D
C 
C       code from Gerhard Mueller, Inst. of Met. and Geophysics,
C                                  University Frankfurt/Main
C 
C       NN and MM must have the actual values from the calling routine.
C 
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)   

	include 'cgmi.h'

        DIMENSION GT(MM,NN),Q(MM,MM),QI(MM,MM),AL(MM,MM),ALI(MM,MM),  
     *  ALIT(MM,MM),GTD(MM),DTH(NN)

C 
C       MATRIX Q UND VEKTOR GTD   
C 
        if(m.gt.mm .or. n.gt.nn) then
           print *,' Dimensions in DLSQ-routine are too small!!'
	   stop
	endif
        DO  50  I=1,M 
	var(i) = 0.d0
        DO  50  J=1,N 
50      GT(I,J)=G(J,I)
        DO  100  I=1,M
        DO  100  J=1,M
	al(i,j) = 0.d0
        CC=0.d0
        DO  80  K=1,N 
80      CC=CC+GT(I,K)*G(K,J)  
100     Q(I,J)=CC 
        DO  200  I=1,M
        CC=0.d0
        DO  150  K=1,N
150     CC=CC+GT(I,K)*D(K)
200     GTD(I)=CC 
C 
C       MATRIX AL 
C 
        DO  2000  I=1,M   
        S=0.d0
        IF(I.EQ.1)  GO TO 700 
        DO  500  K=1,I-1  
500     S=S+q2(AL(I,K))
700     AL(I,I)=SQRT(Q(I,I)-S)
        IF(I.EQ.M)  GO TO 2000
        DO  1000  J=I+1,M 
        S=0.d0
        IF(I.EQ.1)  GO TO 1000
        DO  800  K=1,I-1  
800     S=S+AL(I,K)*AL(J,K)   
1000    AL(J,I)=(Q(I,J)-S)/AL(I,I)
2000    CONTINUE  
C 
C       MATRIX ALI
C 
        DO  4000  I=1,M   
        ALI(I,I)=1.d0/AL(I,I)   
        IF(I.EQ.1)  GO TO 3200
        DO  3000  J=1,I-1 
        S=0.d0
        DO  2400  K=J,I-1 
2400    S=S+AL(I,K)*ALI(K,J)  
3000    ALI(I,J)=-S/AL(I,I)   
3200    DO  3500  J=I+1,M 
3500    ALI(I,J)=0.d0 
4000    CONTINUE  
C 
C       MATRIX QI 
C 
        DO  4500  I=1,M   
        DO  4500  J=1,M   
4500    ALIT(I,J)=ALI(J,I)
        DO  4800  I=1,M   
        DO  4800  J=1,M   
        CC=0.d0
        DO  4600  K=1,M   
4600    CC=CC+ALIT(I,K)*ALI(K,J)  
4800    QI(I,J)=CC
C 
C       MODEL A  
C 
        DO  5000  I=1,M   
        CC=0.d0
        DO  4900  K=1,M   
4900    CC=CC+QI(I,K)*GTD(K)  
5000    A(I)=CC   
C 
C       RESIDUALS AND MODEL STANDARD DEVIATIONS  
C 
        DO  5500  I=1,N   
        CC=0.d0
        DO  5300  K=1,M   
5300    CC=CC+G(I,K)*A(K) 
5500    DTH(I)=CC 
        S=0.d0
        DO  6000  I=1,N   
        RES(I)=D(I)-DTH(I)
6000    S=S+q2(RES(I))
        if(n.gt.m) then
	  E=S/(N-M) 
          DO  7000  I=1,M   
7000      VAR(I)=DSQRT(QI(I,I)*E)
	endif
        RETURN
        END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine findosci(amin,amax,a,n1,n2,ind)
c
c       fixed for longitudes: Nov 5, 1997 JS
c
	implicit double precision (a-h,o-z)
	dimension a(*)
	integer ind,n1,n2

        deg2rad = datan(1.d0)/45.d0

	amin = 9.d99
	amax = -9.d99
	bmin = 9.d99
	bmax = -9.d99
	cmin = 9.d99
	cmax = -9.d99

	do 10 i=n1,n2,-1

	if(ind.eq.1) then

	   if(a(i).lt.amin) amin = a(i)
	   if(a(i).gt.amax) amax = a(i)

	else if (ind.eq.2) then
c
c       We have to handle the longitude values especially!
c

	   p1 = deg2rad*a(i)
	   p2 = dcos(p1)
	   p3 = dsin(p1)

	   if(p2.lt.bmin) bmin = p2
	   if(p2.gt.bmax) bmax = p2

	   if(p3.lt.cmin) cmin = p3
	   if(p3.gt.cmax) cmax = p3

	endif

10      continue

	if (ind.eq.2) then

           amin = alpha1(datan2(cmin,bmax)/deg2rad)
           amax = alpha1(datan2(cmax,bmin)/deg2rad)

	endif
	      
        return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        function convlat(rlat,ind)
c
c       Funtion to convert geographic latitude into geocentric
c       and back.
c       
c       ind = 1 geographic to geocentric latitude
c
c           = 2 geocentric to geographic latitude
c
c       Johannes Schweitzer, NORSAR, October 2, 1997
c
        implicit double precision (a-h,o-z)
	integer ind

	pi      = 4.d0*datan(1.d0)
	deg2rad = pi / 180.d0
	rad2deg = 180.d0 / pi

c
c     after Stacey (1992), Physcis of the Earth
c
        rada = 6378.136d0
        radb = 6356.751d0

        flatt = (rada - radb)/rada

	eps = q2(1.d0-flatt)

	elat = rlat*deg2rad

	if(ind.eq.1) then

	  convlat = rad2deg * datan(eps*dtan(elat))

	else if(ind.eq.2) then

	  convlat = rad2deg * datan(dtan(elat)/eps)

	else 

	  convlat = rlat

	endif

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       function bilinear to interpolate the value of a
c       point inside a given area of 4 cornerpoints and
c       the values at these cornerpoints respectively.
c 
c       Johannes Schweitzer, March 1999
c
c       assuming a cartesian grid is spanned by:
c
c              x(4),y(4)  .       x(3),y(3)
c                         .
c               ....... x1,y1 ...........
c                         .
c              x(1),y(1)  .       x(2),y(2)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function bilinear(x,y,z,x1,y1,indbi)
	implicit double precision (a-h,o-z)

	dimension x(4), y(4), z(4)

	save

	if(indbi.ne.0) go to 10

	dx = x(2)-x(1)
	if(dx.ne.0.d0) then
	  t = (x1-x(1))/dx
	else
	  t = 1.d0
	endif

	dy = y(4)-y(1)
	if(dy.ne.0.d0) then
	  u = (y1-y(1))/dy
	else
	  u  = 1.d0
	endif


10	bilinear = (1.d0 - t)*(1.d0 - u) *z(1) +
     +                     t *(1.d0 - u) *z(2) +
     +                     t *        u  *z(3) +
     +             (1.d0 - t)*        u  *z(4)
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function getchi(dfree0,confi0)
c
c     uses subroutines and functions of the Numerical Recipes 
c     Software to calculate the chi**2. values for a specific
c     given confidence level and degree of freedom.
c

      double precision dfree,dfree0,confi,getchi,confi0

      double precision rtbis,xacc,xa,xe

      confi = 1.0d0 - (confi0 / 100.d0)

      dfree = dfree0 / 2.0d0

      call zbrak(xa,xe,dfree,confi)

      xacc=0.5d-6*(xa+xe)

      getchi = 2.d0 * rtbis(xa,xe,xacc,dfree,confi)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE zbrak(xa,xe,dfree,confi)
C  (C) Copr. 1986-92 Numerical Recipes Software
c
c     here adopted for calculating chi**2 for a specific
c     confidence level.
c
c     Johannes Schweitzer, NORSAR, October 2000
c
      DOUBLE PRECISION confi,xa,xe,dfree,gammq
      INTEGER i,n
      DOUBLE PRECISION dx,fc,fp,x,x1,x2
      n = 50
      x1= .5d0
      x2= 15.0d0
      x = x1
      dx=(x2-x1)/n
      fp= confi - gammq(dfree,x)

      do 11 i=1,n
        x=x+dx
        fc= confi - gammq(dfree,x)
        if(fc*fp.lt.0.d0) then
          xa=x-dx
          xe=x
	  go to 1
        endif
        fp=fc
11    continue
1     continue
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION rtbis(x1,x2,xacc,dfree,confi)
C  (C) Copr. 1986-92 Numerical Recipes Software
c
c     Here adopted to get chi**2 for a specific confidence level
c     Johannes Schweitzer, October 2000, NORSAR
c
      DOUBLE PRECISION rtbis,x1,x2,xacc,gammq,dfree,confi
      DOUBLE PRECISION dx,f,fmid,xmid
      fmid= confi - gammq(dfree,x2)
      f=    confi - gammq(dfree,x1)
      if(f.lt.0.d0)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,100
        dx=dx*.5d0
        xmid=rtbis+dx
        fmid= confi - gammq(dfree,xmid)
        if(fmid.le.0.d0)rtbis=xmid
        if(dabs(dx).lt.xacc .or. fmid.eq.0.d0) return
11    continue
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION gammq(a,x)
      DOUBLE PRECISION a,gammq,x
C  (C) Copr. 1986-92 Numerical Recipes Software
CU    USES gcf,gser
      DOUBLE PRECISION gammcf,gamser,gln
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammq=1.d0-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION gammln(xx)
C  (C) Copr. 1986-92 Numerical Recipes Software
      DOUBLE PRECISION gammln,xx
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+dlog(stp*ser/x)
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gser(gamser,a,x,gln)
C  (C) Copr. 1986-92 Numerical Recipes Software
      INTEGER n, itmax
      DOUBLE PRECISION a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
CU    USES gammln
      DOUBLE PRECISION ap,del,sum,gammln

      gln=gammln(a)
      if(x.le.0.d0)then
        if(x.lt.0.d0) print *,'GAMMQ: x < 0 in gser'
        gamser=0.d0
        return
      endif
      ap=a
      sum=1.d0/a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.d0
        del=del*x/ap
        sum=sum+del
        if(dabs(del).lt.dabs(sum)*EPS)goto 1
11    continue
      print *,'GAMMQ: a too large, ITMAX too small in gser'
1     gamser=sum*dexp(-x+a*dlog(x)-gln)
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gcf(gammcf,a,x,gln)
C  (C) Copr. 1986-92 Numerical Recipes Software
      INTEGER ITMAX, I
      DOUBLE PRECISION a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
CU    USES gammln
      DOUBLE PRECISION an,b,c,d,del,h,gammln

      gln=gammln(a)
      b=x+1.d0-a
      c=1.d0/FPMIN
      d=1.d0/b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.d0
        d=an*d+b
        if(dabs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(dabs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=d*c
        h=h*del
        if(dabs(del-1.d0).lt.EPS)goto 1
11    continue
      print *, 'GAMMQ: a too large, ITMAX too small in gcf'
1     gammcf=dexp(-x+a*dlog(x)-gln)*h
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      function q2(a)
      double precision a,q2
      q2 = a*a
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function dmean(a,n1,n2)
      double precision a(*),dmean, b
      integer n1,n2
      if(n1.ne.n2) then
         b = 0.d0
         do 10 i = n1,n2,-1
10       b = b + a(i)
         dmean = b /dble(n1-n2+1)
      else
	 dmean = 0.d0
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function ddmax(a,n1,n2)
      double precision a(*),ddmax
      integer n1,n2
      ddmax = -9.d99
      do 10 i = n1,n2,-1
      if(a(i).gt.ddmax) ddmax = a(i)
10    continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine plane(stala,stalo,iev,data,ndat,azi,dazi,
     +                 ray,dray,phas,touse,phase,jref,typctl)
C
C     programmer: Johannes Schweitzer NORSAR
C                 January 2002
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      dimension stala(*),stalo(*),data(*),iev(*)
      character*(*) touse(*),phase(*)

      character phas*8

      integer ndat,typctl

      include 'cgmi.h'

c
      pi=4.0d0*datan(1.0d0)
      pi2 = 2.d0*pi
      rad2phi=180.d0/pi
      phi2rad=pi/180.d0

      azi  = 0.d0
      dazi = 180.d0
      ray  = 0.d0
      dray = 100.d0

      j= 0
      d2k0 = 0.d0
      is   = 1
      jref = 0

51    DO 100 i=is,ndat

      if(touse(i)(7:7).ne.'P') go to 100

      if(j.eq.0) then
	 phas = phase(i)
	 stla = stala(iev(i))
	 stlo = stalo(iev(i))
	 dat0 = data(i)
	 jref = i
	 go to 90
      endif
      if(phase(i).ne.phas) go to 100
      call depi(stala(iev(i)),stalo(iev(i)),stla,stlo,dg,dk,
     +          azi0,baz,d2k)
      if(dk.le.0.d0) go to 100
      g(j,1) = dk*dsin(phi2rad*azi0)
      g(j,2) = dk*dcos(phi2rad*azi0)
      d(j) = data(i) - dat0
      d2k0 = d2k0 + d2k
90    j = j + 1
100   continue

      n=j-1

      if(n.lt.3) then

	j  = 0
	is = jref + 1
	jref = 0

	if (is.ge.ndat-2 .or. n.lt.0) then
	   print *,'Tried, but cannot find enough data for a ',
     +	    'plane wave approximation'
	   go to 900
	endif
	go to 51

      endif

      m=2
      d2k0 = d2k0 / dble(n)
c
      if(typctl .ge. 8) then
	write(6,
     +	'(//''The kernel for the plane wave inversion'')')
	 do 150 i=1,n
           write(6,'(i3,2f10.1)') i,g(i,1),g(i,2)
150      continue
      endif
c
c     the least squares fit
c
      call dlsq(n,m)
c

      phi = datan2(a(1),a(2))
      sph = dsin(phi)
      cph = dcos(phi)

      v  = sph / a(1)

      a1 = sph*v
      a2 = cph*cph/a(2)

      dazi  = rad2phi*dsqrt((a2*var(1))**2.d0 +
     +                      (a1*var(2))**2.d0 )

      dray  = d2k0*dsqrt((var(1)*sph)**2.d0 +
     +                   (var(2)*cph)**2.d0 )

      phi = phi + pi
      if(phi.lt.0.d0)  phi = phi + pi2
      if(phi.ge.pi2)   phi = phi - pi2

      azi = phi*rad2phi
      ray = d2k0 / v

c
      if(typctl.gt.5 ) then

         write (6,'('' PLANE: Plane wave fit of '',i3,
     +       '' measured onsets:'')') n+1

         write (6,'(''ray parameter'',f7.2,
     +      ''+/-'',f7.2,''sec/deg azimuth:'',f7.2,''+/-'',f7.2)')
     +      ray,dray,azi,dazi

      end if

      if(typctl.gt.8 ) then
	write(6,
     +	   '(//''Model parameter and standard deviations'')')
	   do 200 i=1,m
              write(6,'(i3,3f8.4)') i,a(i),var(i)
200        continue

        write(6,'(//''Observed data, calculated data, and residuals'')')
        do 250 i=1,n
           dt=d(i)-res(i)
           write(6,'(i3,3f10.4)') i,d(i),dt,res(i)
250     continue
      end if
c
900   continue
      return
      end
