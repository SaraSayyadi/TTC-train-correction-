      PROGRAM TTCLAM
C**************************************************************************
C*                                                                        *
C*  This program carries out a Total Terrain Correction on                *
C*  gridded data using the line mass approximation for distant            *
C*  prisms, and Nagy's Equation for those close to the gravity station.   *
C*                                                                        *
C*  Before use, the spheroid constants and the UTM projection constants   *
C*                must be defined. (Not needed in this version, mtg)      *
C*                                                                        *
                    *
C**************************************************************************
C*
C	Output for 25 m, 200 m and 1000 m grid
C
C
C
C**************************************************************************
      DOUBLE PRECISION PDLONG,PDX,PDY,SMIN,XXS,YYS,PHI,PLAT,PLONG
     +,xlo,xhi,ylo,yhi,zlo,zhi,dx,dy
     +,x1,x2,y1,y2,z,tc,FA
      integer*2 nnx,nny
      INTEGER STANUM
      CHARACTER*30 ANAFN,BNAFN,CNAFN
      character*4 aiii
      character*6 stna
      DIMENSION ZARR(5001,5001),PLATS(5001,5001),PLONGS(5001,5001)
      COMMON/STA/XSTA,YSTA,ZSTA
      COMMON/ALL/G,GRISEP
      COMMON/CEQN9/X1,X2,Y1,Y2,Z,TC
      COMMON/CPRIS/PLAT,PLONG,PZ
      COMMON/CFAR/SCRX,SCRY,SCRZ
      COMMON/LATCON/F,SMAJ,SMIN,ECC2,SF,DLONO,RLAT,RLON
      common/grinm/fxmax,fxmin,fymax,fymin,gxmax,gxmin,gymax,gymin

      G=6.673E-11
       PHI=3.141592654D0
c
       write (*,*) 'ATH: •essi £tg fa leiŒr‚ttir £t ¡ 50 km '
C    Density (kg/m3)

      !RHO=1200.


C input and output files
C
        WRITE (*,*) ' Grid file (elevation) '
        READ (*,'(A30)') ANAFN
        WRITE (*,*) ' Data file (gravity stations) '
        READ (*,'(A30)') BNAFN
        WRITE (*,*) ' Outfile (terrain corrected gravity) '
        READ (*,'(A30)') CNAFN
        OPEN (UNIT=3,FILE=anafn,status='OLD')
        OPEN (UNIT=9,FILE=bnafn,STATUS='OLD')
        OPEN (UNIT=4,FILE=CNAFN,STATUS='NEW')
c        open (unit=10,file='log.log',status='old')
C read elevation data
C essi innlestur er fyrir skr  £r GRID
C
        !WRITE (*,*)' glacier surface=1, bedrock=2 '
        !READ (*,*)ITYPE
        READ (3,*) AIII
        READ (3,*) nnx,nny
        READ (3,*) xlo,xhi
        READ (3,*) ylo,yhi
        READ (3,*) zlo,zhi
C  Breytir n£ yfir ¡ venjuleg Lambertshnit
        numx=nnx
        numy=nny
        xmin=xlo
        xmax=xhi
        ymin=ylo
        ymax=yhi
C
      DO 1000 I=1,NUMY
           READ (3,*) (zarr(I,J), J=1,NUMX)
1000  CONTINUE

C  Kannar hvort blank-gildi leynist i grid-skra
	write (*,*) 'ATH inniheldur serstakt prof hvort gildi '
	write (*,*) 'laegri en 0 seu i griddi.  Setur tha punkta '
	write (*,*) '= 0.  Forritid raedur ekki vid sjo enntha '
	write (*,*) '22. okt. 1997 '
        DO 2222 I=1,numy
           DO 2222 J=1,NUMX

             IF (ZARR(I,J) .GT. 1.E30) THEN
               write (*,*) ' grid skra inniheldur blank gildi!'
	           goto 9900
             ENDIF
2222    CONTINUE
C
      GRISEP=(YMAX-YMIN)/(NUMY-1)
      DO 2100 I=1,NUMY
        DO 2200 J=1,NUMX
                PDX=XMAX-(J-1)*GRISEP
                PDY=YMIN+(I-1)*GRISEP
                PLATS(I,J)=PDY
                PLONGS(I,J)=PDX
2200    CONTINUE
2100  CONTINUE
c
c
      do 2900 nn=1,1000
c
C Read Gravity Data here
c
C 2020-06-14 Changed read from file to format 5000:
      READ(9,3000,END=9900)
     + STNA,XSTA,YSTA,ZSTA,FA,T1,T1sea,T2,T2sea,TCORR,ZDIFF
3000  FORMAT(A6,2F8.0,1f7.1,7F7.7)
      !write (*,*)FA
5000  FORMAT(A6,2F8.0,1f7.1,F9.2,6F9.7)
c      WRITE(*,5000)STNA,XSTA,YSTA,ZSTA,FA,T1,T1sea,T2,T2sea,TCORR,ZDIFF
C

      XSTA=XSTA
      TCT=0.
      TCL=0.
      TCLsea=0.
      TCTsea=0.
      XXS=XSTA
      YYS=YSTA

c
c	finnur nu hnit punkta i grid-skra sem nota a i utreikn.
c
      call gridnm(grisep,xsta,ysta,xmax,xmin,ymax,ymin,indi)
!C

!C
      NX0=(XMAX-FXMAX)/GRISEP+2
      NXF=NX0+(FXMAX-FXMIN)/GRISEP-1
      NY0=(FYMIN-YMIN)/GRISEP+2
      NYF=NY0+(FYMAX-FYMIN)/GRISEP-1
C 2020-06-14 - insert write statement below for testing:
        !write(*,*) NX0, NXF, NY0, NYF
!c
c:  2020-06-20 - correct setup in next two lines - shorter do-loop only used for testing
      DO 4100 I=NY0,NYF
        DO 4000 J=NX0,NXF
C       DO 4100 I=NYF-350,NYF-350
C        DO 4000 J=NXF-350,NXF-350

           TCA=0.
C 14.6.2020 changed XMin-(J-1)*GRISEP to XMin+(J-1)*GRISEP below
           PX=XMin+(J-1)*GRISEP
           PY=YMIN+(I-1)*GRISEP
           PZ=ZARR(I,J)
           !write (*,*) Px,Py,Pz, zsta
           PLAT=PLATS(I,J)
           PLONG=PLONGS(I,J)

!C  Determine which method of calculating TC is required

           DX=PX-XSTA
           DY=PY-YSTA
!c    write (*,*) ' lina 209-211 '
           DX=SQRT(DX*DX)
           DY=SQRT(DY*DY)
           DD=SQRT(DX*DX+DY*DY)
 !c          write (10,*)Pz
!C  Initially only terrain corrected up to 50km from station
      IF(DD.GT.50000.)GOTO 4000

       IF (INDI.EQ.1) THEN
       !IF (DX.LT.1600. .AND. DY.LT.1600.) THEN

          IF (DX.LT. 2.5 .AND. DY.LT. 2.5 .AND. ITYPE.EQ.1) THEN
          ZDIFF=ZSTA-PZ
          PZ=ZSTA
          ENDIF
	!write (*,*) ' lina 222 '
         RHO=1200.
        If (pz > 0) then
        Z=SQRT((PZ-ZSTA)**2)
        RHO=RHO
        CALL NAGYin(PX,PY,DX,DY,TCA)
        TC1=TCA
        Z=ZSTA
        CALL NAGYin(PX,PY,DX,DY,TCA)
         TCA=Rho*(TCA-TC1)
         !print *, Rho
        Elseif (pz <0) then
c 2020-06-20  - TESTING temporarily setting: if pz < 0 then pz = 0)
c           pz = 0
          pZ=abs(PZ)+(2*ZSTA)
          Z=PZ-ZSTA
          RHO=RHO-1025
          !print*,RHO
          CALL NAGYin(PX,PY,DX,DY,TCAsea)
          TC1sea=TCAsea
          Z=ZSTA
          CALL NAGYin(PX,PY,DX,DY,TCAsea)
          TCAsea=Rho*(TCAsea-TC1sea)
c  --- writes out TCAsea in mGal
c          write (*,*)TCAsea
c  --- test ends
c          print*,RHO
      endif

         !print*,TCA
      !end if
      ELSE IF (INDI.EQ.2) THEN

        !IF (DX.LT. 100.AND. DY.LT. 100 .AND. ITYPE.EQ.2) THEN
          !ZDIFF=ZSTA-PZ

          !PZ=ZSTA
      ! ENDIF

       RHO=1200.
         IF (pz > 0 ) then
         Z=SQRT((PZ-ZSTA)**2)
         RHO=RHO
         CALL NAGYout(PX,PY,DX,DY,TCA)
         TC1=TCA
         Z=ZSTA
         CALL NAGYout(PX,PY,DX,DY,TCA)
         TCA=RHO*(TCA-TC1)

         !print *,TCA
       elseif (pz < 0) then
          Z=abs(PZ)+(2*ZSTA)
          Z=PZ-ZSTA
          RHO=RHO-1025
          !print*,rho
          CALL NAGYout(PX,PY,DX,DY,TCAsea)
          TC1sea=TCAsea
          Z=ZSTA
          CALL NAGYout(PX,PY,DX,DY,TCAsea)
          TCAsea=Rho*(TCAsea-TC1sea)

       endif

      ELSE
      IF (PX.LT.GXMAX.AND.PX.GT.GXMIN.AND.PY.LT.GYMAX.AND.
     1        PY.GT.GYMIN) GOTO 4000

       ENDIF
            !write (*,*) Px,Py,Pz, zsta
       TCL=TCL+TCA

       TCLsea=TCLsea+TCAsea
      !print *, TCLsea
c      write (*,*) TcLsea
4000    CONTINUE
4100  CONTINUE

C  Convert ms-2 into milligals
      !print*, TCL
      TCTsea=TCLsea*100000.
      TCT=TCL*100000.
      TCORR=TCORR+TCT+TCTsea

      IF (INDI .EQ. 1 ) THEN
       T1=T1+TCT

       T1sea=T1sea+TCTsea

          ELSEIF (INDI .EQ. 2 ) THEN
            T2=T2+TCT

             T2sea=T2sea+TCTsea
       ENDIF
 !c     write(*,*)T1sea
      ! XSTA=XSTA
        !print*,TCT
        !print *,T2

         WRITE(4,5000)
     + STNA,XSTA,YSTA,ZSTA,FA,T1,T1sea,T2,T2sea,TCORR,ZDIFF
        WRITE(*,5000)
     + STNA,XSTA,YSTA,ZSTA,FA,T1,T1sea,T2,T2sea,TCORR,ZDIFF


2900  continue
9900  CONTINUE
      CLOSE (3)
      CLOSE (4)
      CLOSE (9)
      STOP
      END











c***********************************************************************
c                                                                      *
c	GRIDNM calculates the coordinates of points in grid used in    *
c	calculation of correction for gravity points with coordinates  *
c	XSTA,YSTA.						       *
c								       *
c***********************************************************************
C
      SUBROUTINE GRIDNM(GRISEP,XSTA,YSTA,XMAX,XMIN,YMAX,YMIN,INDI)
      common/grinm/fxmax,fxmin,fymax,fymin,gxmax,gxmin,gymax,gymin

       IF ((GRISEP-5).LE. 0.1) THEN
        FXMAX=INT((XSTA+1500.+100.)/200.)*200.+100.
        FYMAX=INT((YSTA+1500.+100.)/200.)*200.+100.
        FXMIN=FXMAX-2*1500.
        FYMIN=FYMAX-2*1500.
        !write (*,*)fxmin,fxmax,fymin,fymax
        INDI=1
      ELSEIF ((GRISEP-200) .LE. 0.1) THEN
        FXMAX=INT((XSTA+50000.+500.)/1000.)*1000.+500.
        FYMAX=INT((YSTA+50000.+500.)/1000.)*1000.+500.
        FXMIN=FXMAX-2*50000.
        FYMIN=FYMAX-2*50000.
        INDI=2
      ENDIF
       IF (INDI.EQ.2) THEN
         GXMAX=INT((XSTA+1500.+100.)/200.)*200.+100
         GYMAX=INT((YSTA+1500.+100.)/200.)*200.+100.
         GXMIN=GXMAX-2*1500.
         GYMIN=GYMAX-2*1500.
      ! ELSEIF (INDI.EQ.3) THEN
         !GXMAX=INT((XSTA+5000.+2500.)/1000.)*1000.+200.
         !GYMAX=INT((YSTA+5000.+500.)/1000.)*1000.+200.
         !GXMIN=GXMAX-2*6000.
         !GYMIN=GYMAX-2*6000.
      ! ELSEIF (INDI.EQ.3) THEN
        !GXMAX=INT((XSTA+50000.+2500.)/5000.)*5000.+500.
        !GYMAX=INT((YSTA+50000.+2500.)/5000.)*5000.-500.
        !GXMIN=GXMAX-2*50000.
        !GYMIN=GYMAX-2*50000.
       ENDIF
       !WRITE (10,*)FXMIN,FXMAX,FYMIN,FYMAX
c      WRITE (10,*)GXMIN,GXMAX,GYMIN,GYMAX
c      write (10,*) grisep
c      WRITE (*,*)INDI
      RETURN
      END





C**********************************************************************
C                                                                     *
C     NAGY is for prisms for innerloop                 *
C                                                                     *
C**********************************************************************

      SUBROUTINE NAGYin(X,Y,DX,DY,TCN)
      double precision x1,x2,y2,z,tc,rho,grisep,g
      double precision dx,dy,p,q,r,s,pp,qq,rr,ss,y1
      COMMON/STA/XSTA,YSTA,ZSTA
      COMMON/CEQN9/X1,X2,Y1,Y2,Z,TC
c
c	  Checks for size of grid elements for calculation
c   Uses parameter INDI to determina grid spacing
      PP=SQRT((XSTA-(X-2.5))**2)
      QQ=SQRT((XSTA-(X+2.5))**2)
      RR=SQRT((YSTA-(Y-2.5))**2)
      SS=SQRT((YSTA-(Y+2.5))**2)
c      write (*,*)pp,qq,rr,ss
      P=MIN(PP,QQ)
      Q=MAX(PP,QQ)
      R=MIN(RR,SS)
      S=MAX(RR,SS)
      TC=0.
      IF(DX.LT.2.5.AND.DY.LT.2.5)THEN
       X1=0.
       X2=P
       Y1=0.
       Y2=R
       CALL EQN9
       X2=Q
       CALL EQN9
       X2=P
       Y2=S
       CALL EQN9
       X2=Q
       CALL EQN9
      ELSEIF(DX.LT.2.5)THEN
       X1=0.
       X2=P
       Y1=R
       Y2=S
       CALL EQN9
       X2=Q
       CALL EQN9
      ELSE IF(DY.LT.2.5)THEN
       X1=P
       X2=Q
       Y1=0
       Y2=R
       CALL EQN9
       Y2=S
       CALL EQN9
      ELSE
       X1=P
       X2=Q
       Y1=R
       Y2=S
       !write (*,*)x1,x2,y1,y2
       CALL EQN9
      END IF
      TCN=TC
      RETURN
      END








C**********************************************************************
C                                                                     *
C     NAGY is for prisms for which the line mass approximation        *
C         is not valid.                                               *
C                                                                     *
C**********************************************************************


      SUBROUTINE NAGYout(X,Y,DX,DY,TCN)
      double precision x1,x2,y2,z,tc,rho,grisep,g
      double precision dx,dy,p,q,r,s,pp,qq,rr,ss,y1
      COMMON/STA/XSTA,YSTA,ZSTA
      COMMON/CEQN9/X1,X2,Y1,Y2,Z,TC

      PP=SQRT((XSTA-(X-100))**2)
      QQ=SQRT((XSTA-(X+100))**2)
      RR=SQRT((YSTA-(Y-100))**2)
      SS=SQRT((YSTA-(Y+100))**2)
c      write (*,*)pp,qq,rr,ss
      P=MIN(PP,QQ)
      Q=MAX(PP,QQ)
      R=MIN(RR,SS)
      S=MAX(RR,SS)
      TC=0.
      IF(DX.LT.100.AND.DY.LT.100)THEN
       X1=0.
       X2=P
       Y1=0.
       Y2=R
       CALL EQN9
       X2=Q
       CALL EQN9
       X2=P
       Y2=S
       CALL EQN9
       X2=Q
       CALL EQN9
      ELSEIF(DX.LT.100)THEN
       X1=0.
       X2=P
       Y1=R
       Y2=S
       CALL EQN9
       X2=Q
       CALL EQN9
      ELSE IF(DY.LT.100)THEN
       X1=P
       X2=Q
       Y1=0
       Y2=R
       CALL EQN9
       Y2=S
       CALL EQN9
      ELSE
       X1=P
       X2=Q
       Y1=R
       Y2=S
       !write (*,*)x1,x2,y1,y2
       CALL EQN9
      END IF
      TCN=TC
      RETURN
      END

C***************************************************************************
C                                                                          *
C   EQN9 is Nagy's equation (9)                                            *
C                                                                          *
C***************************************************************************

      SUBROUTINE EQN9
      double precision st1,st2,st3,st4,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      double precision st5,tc,x1,x2,y1,y2,z,Rho
      COMMON/ALL/G,GRISEP
      COMMON/CEQN9/X1,X2,Y1,Y2,Z,TC
      ST1=X2**2+Y2**2
      ST2=X2**2+Y1**2
      ST3=X1**2+Y2**2
      ST4=X1**2+Y1**2
      ST5=Z**2
c      write (10,'(5f12.1)')st1,st2,st3,st4,st5
      IF(ST1.EQ.0) THEN
      T1=0.
      T5=0.
      ELSE
	!write (*,*) ' lina 480 - eqn9 '
      T1=LOG((Y2+SQRT(ST1))/(Y2+SQRT(ST1+ST5)))
      T5=LOG((X2+SQRT(ST1))/(X2+SQRT(ST1+ST5)))
      END IF
      IF(ST2.EQ.0.)THEN
      T2=0.
      T7=0.
      ELSE
      T2=LOG((Y1+SQRT(ST2))/(Y1+SQRT(ST2+ST5)))
      T7=LOG((X2+SQRT(ST2))/(X2+SQRT(ST2+ST5)))
      END IF
      IF(ST3.EQ.0)THEN
      T3=0.
      T6=0.
      ELSE
      T3=LOG((Y2+SQRT(ST3))/(Y2+SQRT(ST3+ST5)))
      T6=LOG((X1+SQRT(ST3))/(X1+SQRT(ST3+ST5)))
      END IF
      IF(ST4.EQ.0.)THEN
      T4=0.
      T8=0.
      ELSE
      T4=LOG((Y1+SQRT(ST4))/(Y1+SQRT(ST4+ST5)))
      T8=LOG((X1+SQRT(ST4))/(X1+SQRT(ST4+ST5)))
      END IF
      IF(Y2.EQ.0.AND.ST5.EQ.0.)THEN
      T9=0.
      T10=0.
      ELSE
c      write (10,'(4f10.1,2f12.1)')x1,x2,y1,y2,st1,st5
c	write (*,*) ' lina 510 - eqn9 '
      T9=((Y2**2+ST5+Y2*SQRT(ST1+ST5))/((Y2+SQRT(ST1+ST5))*
     +SQRT(Y2**2+ST5)))
c      write (10,'(f12.8)')t9
      if ((1.-t9) .lt. 1.e-9) then
        t9=asin(1.0)
        else
         t9=asin(t9)
      endif
	!write (*,*) ' lina 519 - eqn9 '
      T10=((Y2**2+ST5+Y2*SQRT(ST3+ST5))/((Y2+SQRT(ST3+ST5))*
     +SQRT(Y2**2+ST5)))
      if ((1.-t10) .lt. 1.e-9) then
        t10=asin(1.0)
        !print *, t10
        else
         t10=asin(t10)
      endif
      END IF
      IF(Y1.EQ.0.AND.ST5.EQ.0.)THEN
      T11=0.
      T12=0.
      ELSE
      !RHo=1200




c	write (*,*) ' lina 532 - eqn9 '
      T11=((Y1**2+ST5+Y1*SQRT(ST2+ST5))/((Y1+SQRT(ST2+ST5))*
     +SQRT(Y1**2+ST5)))
      if ((1.-t11) .lt. 1.e-9) then
        t11=asin(1.0)
        else
         t11=asin(t11)
      endif
c	write (*,*) ' lina 540 - eqn9 '
      T12=((Y1**2+ST5+Y1*SQRT(ST4+ST5))/((Y1+SQRT(ST4+ST5))*
     +SQRT(Y1**2+ST5)))
      if ((1.-t12) .lt. 1.e-9) then
        t12=asin(1.0)
        else
         t12=asin(t12)
      endif
      END IF
      !print*,RHO
      TC=TC+(G*(X2*(T1-T2)-X1*(T3-T4)+Y2*(T5-T6)
     +-Y1*(T7-T8)+Z*(T9-T10-T11+T12)))
c      write (10,'(3f12.8)')t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12
      RETURN
      END






