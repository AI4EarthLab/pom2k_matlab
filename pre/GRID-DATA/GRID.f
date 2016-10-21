      PROGRAM GRID
C
C 1. Generate horizontal and vertical grid
C 2. Read bottom topog. & interpolate to grid - subroutine BATH
C 3. Read temp. & sal.  & interpolate to grid - subroutine TAND
C 4. Read wind velocity & interpolate to grid - subroutine WIND
C 5. Write grid & initial conditions for model - WRITE(40)
C 6. Write data for MATLAB plot - WRITE(43,44,45,46)
C----------------------------------------------------------------------
C NOTE: This is a simple version modified to run as is without any
C       external data and no in-program plotting. All f90 formats
C       were also removed (tested on Linux g77 compiler)
C
C                                              T.E (March, 2003)
C----------------------------------------------------------------------
C
C**********************************************************************!
C                                                   VERSION(7-25-88)   !
C                  GENERAL CIRCULATION MODEL                           !
C          ORTHOGONAL CURVILINEAR COORDINATE GRID GENERATOR            !
C                            ECOM3D                                    !
C                                                                      !
C                                                                      !
C     THIS COMPUTER CODE COMPUTES AN ORTHOGONAL CURVILINEAR COORDINATE !
C     GRID FOR USE WITH THE THREE DIMENSIONAL, TIME DEPENDENT,         !
C     PRIMITIVE EQUATION, CIRCULATION MODEL DEVELOPED BY GEORGE MELLOR !
C     AND ME. A GRID AND INTERPOLATED BOTTOM TOPOGRAPHY ARE PROVIDED.  !
C                                                                      !
C     THE GRID GENERATOR WAS WRITTEN BY JOHN WILKIN OF WOODS HOLE      !
C     OCEANOGRAPHIC INSTITUTION, WOODS HOLE, MASS 02543.               !
C                                                                      !
C                                                                      !
C     FOR DETAILS OF THE GOVERNING EQUATIONS AND SOLUTION              !
C     TECHNIQUES THE INTERESTED READER IS REFERRED TO:                 !
C                                                                      !
C     WILKIN, J.L. A COMPUTER PROGRAM FOR GENERATING TWO-DIMENSIONAL   !
C     ORTHOGONAL CURVILINEAR COORDINATE GRIDS. UNPUBLISHED REPORT,1987 !
C                                                                      !
C     IVES, D.C. AND R.M.ZACHARIAS, CONFORMAL MAPPING AND ORTHOGONAL   !
C     GRID GENERATION, PAPER NO. 87-2057, AIAA/SAE/ASME/ASEE 23RD      !
C     JOINT PROPULSION CONFERENCE, SAN DIEGO, CALIFORNIA, JUNE 1987.   !
C                                                                      !
C     BLUMBERG, A.F. AND G.L. MELLOR, DIAGNOSTIC AND PROGNOSTIC        !
C     NUMERICAL CIRCULATION STUDIES OF THE SOUTH ATLANTIC BIGHT        !
C     J. GEOPHYS. RES., 88, 4579-4592, 1983.                           !
C                                                                      !
C     BLUMBERG, A.F. AND G.L. MELLOR, A DESCRIPTION OF A THREE         !
C     COASTAL OCEAN CIRCULATION MODEL, THREE DIMENSIONAL SHELF         !
C     MODELS, COASTAL AND ESTUARINE SCIENCES, 5, N.HEAPS, ED.,         !
C     AMERICAN GEOPHYSICAL UNION, 1987.                                !
C                                                                      !
C     AND                                                              !
C                                                                      !
C     BLUMBERG, A.F. AND H.J. HERRING, CIRCULATION MODELING USING      !
C     ORTHOGONAL CURVILINEAR COORDINATES, THREE DIMENSIONAL MODELS     !
C     OF MARINE AND ESTUARINE DYNAMICS, J.C.J.NIHOUL AND B.M.JAMART    !
C     ED., ELSEVIER OCEANOGRAPHY SERIES, 45,1987.                      !
C                                                                      !
C                                                                      !
C                                                                      !
C     PLEASE DIRECT CRITICISMS, SUGGESTIONS AND REPORTS OF ERRORS TO ME!
C                                                                      !
C                                             ALAN BLUMBERG            !
C                                               HYDROQUAL              !
C                                               Mahwah,NJ              !
C                                                7/25/88               !
C                                                                      !
C    This version has been changed and extended by George Mellor       !
C     and Tal Ezer                                                     !
C**********************************************************************!
C
C 
C     DOUBLE PRECISION VERSION
C     USES COMPLEX
C
C --- read model size (IM,JM,KB) ----
       INCLUDE 'gridcom'
C
C*********************************************************

      PARAMETER (     L2=IM-1   , M2=JM-1    )
      PARAMETER (     MM=M2/2   , LM=L2/2    )

      PARAMETER (    KEP=10,NWRK=(KEP-2)*(2**(KEP+1))+KEP+5*L2+6*M2+49)
      PARAMETER (    LIJ=IM*JM   )
C     COMPLEX*16     Z(M2+L2+M2+L2)
      COMPLEX        Z(M2+L2+M2+L2)
      DIMENSION      XB(M2+L2+M2+L2),YB(M2+L2+M2+L2),WRK(M2+L2),
     1               WXI(L2+1),WETA(M2+1),XINT(M2+L2),YINT(M2+L2),
     2               X(0:L2,0:M2),Y(0:L2,0:M2),RHS(0:L2,0:M2),EWRK(NWRK)             
      COMMON / XIEJ / SXI(0:L2),SETA(0:M2)
C
      REAL KM,KH
      COMMON/OUT40/AZ(KB),ZZ(KB),DZ(KB),DZZ(KB),ALON(IM,JM),ALAT(IM,JM),
     1    DX(IM,JM),DY(IM,JM),H(IM,JM),FSM(IM,JM),DUM(IM,JM),
     2    DVM(IM,JM),ART(IM,JM),ARU(IM,JM),ARV(IM,JM),COR(IM,JM),
     3    T(IM,JM,KB),S(IM,JM,KB),TMEAN(IM,JM,KB),SMEAN(IM,JM,KB),
     4    RMEAN(IM,JM,KB),TBE(JM,KB),TBW(JM,KB),TBS(IM,KB),
     5    SBE(JM,KB),SBW(JM,KB),SBS(IM,KB)
CC   
      DIMENSION       UW(IM,JM),VW(IM,JM),WUSURF(IM,JM),WVSURF(IM,JM)
      DIMENSION       LBL(20),DT(IM,JM),DS(IM,JM),ANG(IM,JM)
      DIMENSION       ERRPLT(LM,MM),IVAR(LM,MM),ZD(LM,MM)
      CHARACTER*50    TITLE   
      EQUIVALENCE     (RHS,ERRPLT)
      EQUIVALENCE     (X,ALON),(Y,ALAT)
      EXTERNAL        COFX,COFY
      DATA PI/3.141593/,RAD/.01745329/,RE/6371.E3/,GRAV/9.807/
C
C     OPEN(6,FILE='printout',form='formatted')
C
C --- read parameter file (IGRID,IBATH,ITAND,IWIND) from rungrid
       INCLUDE 'params'
C
C --- SPECIFIED GRID INSTEAD OF CURVILINEAR GRID GENERATION
C                               (ignore gridborder)
      IF(IGRID.EQ.0) THEN
        ASOUTH=32.
        ANORTH=39.
        AWEST=-77.
        AEAST=-71.
C -- for matlab plot only
         WRITE(45,'(2F8.2)') AWEST,ASOUTH
         WRITE(45,'(2F8.2)') AWEST,ANORTH
         WRITE(45,'(2F8.2)') AEAST,ASOUTH
         WRITE(45,'(2F8.2)') AEAST,ANORTH
       DO I=1,IM
       DO J=1,JM
        ALON(I,J)=AWEST+I*(AEAST-AWEST)/IM
        ALAT(I,J)=ASOUTH+J*(ANORTH-ASOUTH)/JM
       ENDDO
       ENDDO
        GO TO 1100
      ENDIF
C
C     INITIALIZE VECTOR Z (COMPLEX) WITH CONTOUR OF PHYSICAL BOUNDARY
      N1 = M2
      N2 = M2+L2
      N3 = M2+L2+M2
      N4 = M2+L2+M2+L2
      N  = M2+L2+M2+L2
C   
C    ---------------------------------------------------------------
C    | input example:                 
C    |    
C    |                         boundary #4
C    |    
C    |                18   17   16   15   14   13
C    |                                           
C    |                 1                       12
C    |                                             
C    |  boundary #1    2                       11   boundary #3
C    |                                           
C    |     ^^          3                       10
C    |     ||                                    
C    |     ||          4    5    6    7    8    9   
C    |                         boundary #2
C    |                    I direction ===>
C    |    
C    ---------------------------------------------------------------
C
C==========================
C      Read input option and latitude of original point
C
       IOPT=2
       ORCOR=0.
C==========================
C       OPTION 1 :
C       GENERATE BOUNDARY FROM MATH. EXPRESSION
C
       IF(IOPT.EQ.1)THEN
C
      DO 1 I=1,N1
 1       CALL Z1(FLOAT(I)/FLOAT(N1),XB(I),YB(I),Z(I))
      DO 2 I=N1+1,N2
 2       CALL Z2(FLOAT(I-N1)/FLOAT(N2-N1),XB(I),YB(I),Z(I))
      DO 3 I=N2+1,N3
 3       CALL Z3(FLOAT(I-N2)/FLOAT(N3-N2),XB(I),YB(I),Z(I))
      DO 4 I=N3+1,N4
 4       CALL Z4(FLOAT(I-N3)/FLOAT(N4-N3),XB(I),YB(I),Z(I))
      ELSE
C=========================
C       OPTION 2 :
C       READ BOUNDARY FROM SUBROUTINE BORDER
C
      CALL BORDER(IM,JM,X,Y)
C
c     CALL PRXY (' X BF MERCATOR',0.0,X,L2+1,1   ,M2+1,1   ,0.)    
c     CALL PRXY (' Y BF MERCATOR',0.0,Y,L2+1,1   ,M2+1,1   ,0.)   
C  Stretch latitude as in Mercator projection
       DO I=0,L2
       DO J=0,M2
      Y(I,J)=alog(tan(0.5*rad*Y(I,J)+0.25*pi))/rad
       END DO
       END DO
C
C     INSERT BOUNDARY VALUES OF GRID FOR DISPLAY
      DO 90 I=1,N1
         XB(I)=X(0,N1-I)
         YB(I)=Y(0,N1-I)
 90      Z(I)=CMPLX(XB(I),YB(I))         
      DO 91 I=N1+1,N2
         XB(I)=X(I-N1,0)
         YB(I)=Y(I-N1,0)
 91      Z(I)=CMPLX(XB(I),YB(I))          
      DO 92 I=N2+1,N3
         XB(I)=X(L2,I-N2)        
         YB(I)=Y(L2,I-N2)        
 92      Z(I)=CMPLX(XB(I),YB(I))                  
      DO 93 I=N3+1,N4
         XB(I)=X(N4-I,M2)       
         YB(I)=Y(N4-I,M2)       
 93      Z(I)=CMPLX(XB(I),YB(I))                    
c     CALL PRXY (' X BF MAPPING',0.0,X,L2+1,1   ,M2+1,1   ,0.)    
c     CALL PRXY (' Y BF MAPPING',0.0,Y,L2+1,1   ,M2+1,1   ,0.)   
       ENDIF 
C
C     MAP PHYSICAL BOUNDARY TO A RECTANGLE
      DO 11 K=1,9    
         CALL RECT(Z,N,N1,N2,N3,N4)
C
C        CALCULATE DEPARTURE OF CONTOUR FROM RECTANGULAR
         ERROR = 0.
         DO 6 I=1,N1
 6          ERROR = ERROR + ABS(DBLE(Z(I))-DBLE(Z(1)))
         DO 7 I=N1+1,N2
 7          ERROR = ERROR + ABS(AIMAG(Z(I))-AIMAG(Z(N1+1)))
         DO 8 I=N2+1,N3
 8          ERROR = ERROR + ABS(DBLE(Z(I))-DBLE(Z(N2+1)))
         DO 9 I=N3+1,N4
 9          ERROR = ERROR + ABS(AIMAG(Z(I))-AIMAG(Z(N3+1)))
         ERROR = ERROR/FLOAT(N4)
         WRITE(6,10)K,ERROR
 10      FORMAT(' RECTANGULARITY ERROR IN MAPPED CONTOUR AT',
     1          ' ITERATION ',I2,' IS ',1PE10.4)
 11   CONTINUE
CC    WRITE(6,'(1x,2F10.3)') (Z(I),I=1,N) 
C
C     CUBIC SPLINE INTERPOLATION OF MAPPING ON BOUNDARIES 3 AND 4 TO
C     MATCH DISTRIBUTION OF POINTS WITH THOSE ON BOUNDARIES 1 AND 2
C
C     BOUNDARY 3
      DO 13 I=1,N3-N2
         WETA(I) = AIMAG(Z(N2+I))
         XINT(I) = XB(N2+I)
 13      YINT(I) = YB(N2+I)
      CALL SPLINE(WETA,XINT,N3-N2,1.E+35,1.E+35,WRK)
      DO 14 I=1,N1-1
         CALL SPLINT(WETA,XINT,WRK,N3-N2,AIMAG(Z(I)),XB(N3-I))
 14   CONTINUE   
      CALL SPLINE(WETA,YINT,N3-N2,1.E+35,1.E+35,WRK)
      DO 15 I=1,N1-1
         CALL SPLINT(WETA,YINT,WRK,N3-N2,AIMAG(Z(I)),YB(N3-I))
 15   CONTINUE
C
C     BOUNDARY 4
      DO 16 I=1,N4-N3
         WXI(I)  = REAL(Z(N4+1-I))
         XINT(I) = XB(N4+1-I)
 16      YINT(I) = YB(N4+1-I)
      CALL SPLINE(WXI,XINT,N4-N3,1.E+35,1.E+35,WRK)
      DO 17 I=1,N2-N1-1
         CALL SPLINT(WXI,XINT,WRK,N4-N3,REAL(Z(N1+I)),XB(N4-I))
 17   CONTINUE 
      CALL SPLINE(WXI,YINT,N4-N3,1.E+35,1.E+35,WRK)
      DO 18 I=1,N2-N1-1
         CALL SPLINT(WXI,YINT,WRK,N4-N3,REAL(Z(N1+I)),YB(N4-I))
 18   CONTINUE
C
C     STORE DISTRIBUTION OF XI,ETA POINTS ALONG BOUNDARIES (USED TO
C     COMPUTE COEFFICIENTS OF ELLIPTIC EQUATION)
      DO 20 I=0,L2
 20      SXI(I) = REAL(Z(N1+I))
      DO 21 I=0,M2-1
 21      SETA(I) = AIMAG(Z(N1-I))
         SETA(M2) = AIMAG(Z(N4))
C
C     SET BOUNDARY VALUES OF THE GRID
      DO 30 I=1,N1
         X(0,N1-I)  = XB(I)
 30      Y(0,N1-I)  = YB(I)
      DO 31 I=N1+1,N2
         X(I-N1,0)  = XB(I)
 31      Y(I-N1,0)  = YB(I)
      DO 32 I=N2+1,N3
         X(L2,I-N2) = XB(I)
 32      Y(L2,I-N2) = YB(I)
      DO 33 I=N3+1,N4
         X(N4-I,M2) = XB(I)
 33      Y(N4-I,M2) = YB(I)
c     CALL PRXY (' X BF SEPELI ',0.0,X,L2+1,1   ,M2+1,1   ,0.)    
c     CALL PRXY (' Y BF SEPELI ',0.0,Y,L2+1,1   ,M2+1,1   ,0.)   
C
C     SET RIGHT HAND SIDE OF ELLIPTIC EQUATION TO ZERO
      DO 40 J=0,M2
         DO 40 I=0,L2
 40         RHS(I,J) = 0.
C
C     SOLVE ELLIPTIC EQUATION TO FILL IN GRID
      EWRK(1)=NWRK
      CALL SEPELI(0,2,0.,FLOAT(L2),L2,1,WRK,WRK,WRK,WRK,
     1                0.,FLOAT(M2),M2,1,WRK,WRK,WRK,WRK,
     2                COFX,COFY,RHS,X,L2+1,EWRK,PERTRB,IERR)
      IF(IERR.NE.0)CALL CRASH(1,IERR)
      EWRK(1)=NWRK
      CALL SEPELI(0,2,0.,FLOAT(L2),L2,1,WRK,WRK,WRK,WRK,
     1                0.,FLOAT(M2),M2,1,WRK,WRK,WRK,WRK,
     2                COFX,COFY,RHS,Y,L2+1,EWRK,PERTRB,IERR)
      IF(IERR.NE.0)CALL CRASH(2,IERR)
C
C-- Restore stretched coordinate to latitude ---------
C
       DO I=0,L2
       DO J=0,M2
      Y(I,J)=(-0.5*PI+2.*atan(exp(Y(I,J)*rad)))/rad
       END DO
       END DO
C
 1100  CONTINUE
C
C --- SET ARTIFICIAL BATHYMETRY (SEAMOUNT)
C
C
      DO I=1,IM
      DO J=1,JM
      H(I,J)= 4500.*(1.-0.9*EXP(-((ALON(I,J)-ALON(IM/2,J))**2+
     &    (ALAT(I,J)-ALAT(I,JM/2))**2)/ 4.0**2) )
       IF(IBATH.EQ.1) H(I,J)= 0.
      END DO
      END DO
C
C*******************************************************************
C     READ REAL BATHYMETRY AND INTERPOLATE INTO GRID
C*******************************************************************
C
      IF(IBATH.EQ.1) CALL BATH
C
C                  MIN DEPTH (M)
      HMIN=10.
C
c     CALL PRXY(' H BEFORE ',TIME,H,IM,1,JM,1,0.)
      DO 50 I=1,IM
      DO 50 J=1,JM
       FSM(I,J)=1.
       DUM(I,J)=1.
       DVM(I,J)=1.
C
      IF(H(I,J).LT.HMIN)THEN
      H(I,J)=1.
      FSM(I,J)=0.
      DUM(I,J)=0.
      DVM(I,J)=0.
      ENDIF
   50 CONTINUE
      DO 52 J=1,JM-1
      DO 52 I=1,IM
      IF(FSM(I,J).EQ.0..AND.FSM(I,J+1).NE.0.) DVM(I,J+1)=0.
 52   CONTINUE
      DO 53 J=1,JM
      DO 53 I=1,IM-1
      IF(FSM(I,J).EQ.0..AND.FSM(I+1,J).NE.0.)DUM(I+1,J)=0.
 53   CONTINUE
C
C     SUBJECTIVE TOPOG. SMOOTHER - DH/H < 2*SLMIN
C     (smaller SLPMIN gives more smoothing)
      SLMIN=0.2
C     SLMIN=0.1
cc    CALL SLPMIN(H,IM,JM,FSM,SLMIN)
C Laplasian filter
C     CALL SMOOTH(H,FSM,IM,JM,1)
      DO 54 J=1,JM
      DO 54 I=1,IM
      IF(FSM(I,J).EQ.0.) H(I,J)=1.
   54 CONTINUE
C              
C               CALC. DX DY AND  MAX TIME STEP
C              
      DO 150 J=2,JM-1
      DO 150 I=2,IM-1
      DX(I,J)=0.5*RAD*RE*SQRT(((ALON(I+1,J)-ALON(I-1,J))
     1   *COS(ALAT(I,J)*RAD))**2+(ALAT(I+1,J)-ALAT(I-1,J))**2)
      DY(I,J)=0.5*RAD*RE*SQRT(((ALON(I,J+1)-ALON(I,J-1))
     1   *COS(ALAT(I,J)*RAD))**2+(ALAT(I,J+1)-ALAT(I,J-1))**2)
  150 CONTINUE
      DO 152 I=1,IM
      DX(I,1)=DX(I,2)
      DY(I,1)=DY(I,2)
      DX(I,JM)=DX(I,JM-1)
  152 DY(I,JM)=DY(I,JM-1)
      DO 153 J=1,JM
      DX(1,J)=DX(2,J)
      DY(1,J)=DY(2,J)
      DX(IM,J)=DX(IM-1,J)
  153 DY(IM,J)=DY(IM-1,J)
      DO 154 J=1,JM
      DO 154 I=1,IM
      UMAX=1.
      CMAX=2*SQRT(GRAV*H(I,J))+UMAX
      DXY=1./SQRT(1./(DX(I,J)**2) + 1./(DY(I,J)**2))
      DT(I,J)=DXY/CMAX
  154 CONTINUE
C
      DO 155 J=1,JM
      DO 155 I=1,IM
      ART(I,J)=DX(I,J)*DY(I,J)
      ARU(I,J)=.25*(DX(I,J)+DX(I-1,J))*(DY(I,J)+DY(I-1,J))
      ARV(I,J)=.25*(DX(I,J)+DX(I,J-1))*(DY(I,J)+DY(I,J-1))
C
      COR(I,J)=2.*7.29E-5*SIN(ALAT(I,J)*RAD)
C
  155 CONTINUE
C
      TIME=0.
      CALL PRXY(' H ',TIME,H,IM,ISKP,JM,JSKP,0.)
      CALL PRXY('LON',TIME,ALON,IM,ISKP,JM,JSKP,0.)
      CALL PRXY('LAT',TIME,ALAT,IM,ISKP,JM,JSKP,0.)
      CALL PRXY('FSM',TIME,FSM,IM,ISKP,JM,JSKP,0.)
C     CALL PRXY('DUM',TIME,DUM,IM,ISKP,JM,JSKP,0.)
C     CALL PRXY('DVM',TIME,DVM,IM,ISKP,JM,JSKP,0.)
C     CALL PRXY('ART',TIME,ART,IM,ISKP,JM,JSKP,0.)
C     CALL PRXY('ARU',TIME,ARU,IM,ISKP,JM,JSKP,0.)
C     CALL PRXY('ARV',TIME,ARV,IM,ISKP,JM,JSKP,0.)
C     CALL PRXY('COR',TIME,COR,IM,ISKP,JM,JSKP,1.E-7)
      CALL PRXY(' DX   ',TIME,DX,IM,ISKP,JM,JSKP,0.)
      CALL PRXY(' DY   ',TIME,DY,IM,ISKP,JM,JSKP,0.)
      CALL PRXY(' CFL-DT ',TIME,DT,IM,ISKP,JM,JSKP,0.)
C              
C*******************************************************************
C     GENERATE VERCTCAL GRID,
C*******************************************************************
C
      CALL DEPTH(AZ,ZZ,DZ,DZZ,KB)
      WRITE(6,'(''     K      Z         ZZ        DZ       DZZ '')')
      WRITE(6,'(1X,I5,4F10.4)') (K,AZ(K),ZZ(K),DZ(K),DZZ(K),K=1,KB)
C
C --- SET ARTIFICIAL TEMP., SAL. & DENSITY
C
       IF(ITAND.EQ.0) THEN
      DO I=1,IM
      DO J=1,JM
      DO K=1,KB
       S(I,J,K)=35.
       T(I,J,K)=5.+15.*EXP(ZZ(K)*H(I,J)/1000.)
       SMEAN(I,J,K)=35.
       TMEAN(I,J,K)=5.+15.*EXP(ZZ(K)*H(I,J)/1000.)
      END DO
      END DO
      END DO
      CALL DENS(S,T,RMEAN,FSM,ZZ,H,IM,JM,KB)
       ELSE
C
C*******************************************************************
C     READ TEMP. AND SAL. FROM FILE AND INTERPOLATE INTO GRID
C*******************************************************************
C 
      CALL TANDS
       ENDIF
C
C --- boundary fields TBS,TBN,TBE,TBW should now be calc. in pom
C
      CALL PRXYZ ('T       ',TIME,T ,6,IM,ISKP,JM,JSKP,KB,0.01)
      CALL PRXYZ ('S       ',TIME,S ,6,IM,ISKP,JM,JSKP,KB,0.01)
      CALL PRXYZ ('RMEAN-1.0',TIME,RMEAN,6,IM,ISKP,JM,JSKP,KB,1.E-5)
C
C
C
C     CHECK ORTHOGONALITY BY EVALUATING DX/DXI*DX/DETA+DY/DXI*DY/DETA
C     EVERYWHERE.  NORMALIZE WITH RESPECT TO LOCAL GRID CELL AREA.
C     STORE RESULT IN ERRPLT.
C     DO 80 J=2,JM-1
C        DO 80 I=2,JM-1
C80         ERRPLT(I,J) = ((X(I+1,J)-X(I-1,J))*
C    1                     (X(I,J+1)-X(I,J-1))+
C    2                     (Y(I+1,J)-Y(I-1,J))*
C    3                     (Y(I,J+1)-Y(I,J-1)))/
C    4                     ( DX(I,J)*DY(I,J))*100.
C
C     WRITE(6,108)
C 108 FORMAT(//' ORTHOGONALARITY ERROR ANALYSIS IN PERCENT'//)
C     CALL PRXY('ORTH.ERROR ',TIME,ERRPLT,IM,ISKP,JM,JSKP,0.)
C
      do i=1,im
      do j=1,jm
C --- for IWIND=0 one may specify wind (m/s) here instead of data 
       UW(i,j)=0.
       VW(i,j)=0.
       WUSURF(i,j)=0.
       WVSURF(i,j)=0.
C --- calc. tilting angle of curvilinear grid
       if(j.eq.jm) then
      DLON=(ALON(I,J)-ALON(I,J-1))*COS(ALAT(I,J)*RAD)
      DLAT=ALAT(I,J)-ALAT(I,J-1)
       else
      DLON=(ALON(I,J+1)-ALON(I,J))*COS(ALAT(I,J)*RAD)
      DLAT=ALAT(I,J+1)-ALAT(I,J)
       endif
      DLNT=(DLON**2+DLAT**2)**.5
      ANG(I,J)=ASIN(DLON/DLNT)
      end do
      end do
C
C
C
C*******************************************************************
C     READ WIND VELOCITY FROM FILE AND INTERPOLATE INTO GRID
C*******************************************************************
C
C  for monthly data use do loop example
C     do mon=1,12
C
      IF(IWIND.EQ.1) CALL WIND(UW,VW,IM,JM)
C
C --- RHOair/RHOwater
      RWA=1.226e-3/1.025
C
      DO I=1,IM
      DO J=1,JM
C --- transfer coordinates to model grid direction
       BETA=ANG(I,J)
       WUSURF(I,J)=UW(I,J)*COS(BETA)-VW(I,J)*SIN(BETA)
       WVSURF(I,J)=UW(I,J)*SIN(BETA)+VW(I,J)*COS(BETA)
C --- calculate wind stress in model units from wind velocity
C --- note: sign opposite to wind vector
       SPD=SQRT(UW(I,J)**2+VW(I,J)**2)
       CD=(7.5E-4+6.7E-5*SPD)
         WUSURF(I,J)=-1.*RWA*CD*SPD*WUSURF(I,J)
         WVSURF(I,J)=-1.*RWA*CD*SPD*WVSURF(I,J)
      END DO
      END DO
C
      CALL PRXY('ANG(RAD)',0.,ANG  ,IM,ISKP,JM,JSKP,0.001)
      CALL PRXY(' U-VEL(X,Y)',0., UW   ,IM,ISKP,JM,JSKP,0.01)
      CALL PRXY(' V-VEL(X,Y)',0., VW   ,IM,ISKP,JM,JSKP,0.01)
      CALL PRXY('WUSURF ON GRID',0.,WUSURF,IM,ISKP,JM,JSKP,1.E-5)
      CALL PRXY('WVSURF ON GRID',0.,WVSURF,IM,ISKP,JM,JSKP,1.E-5)
C
C     WRITE WIND STRESS FOR MODEL FORCING
C
cc    WRITE(30) WUSURF,WVSURF
C     end do
c      
C*******************************************************************
C     WRITE INITIAL CONDITION FILE FOR POM
C*******************************************************************
C
C     (NOTE: variables names are as in pom98.f)
C
C     WRITE(40) KB,AZ,ZZ,DZ,DZZ,IM,JM,ALON,ALAT,DX,DY,H,FSM,
C    1    DUM,DVM,ART,ARU,ARV,COR,T,S,TMEAN,SMEAN,RMEAN
C
C     IC FOR POM2K: write only necessary fields and calc. others
C      in subroutine "file2ic" in pom2k.f.   (T.E. Dec2004)
C      
C **** FORMATTED OUTPUT ****   var name in pom2k.f
C--- 1D ---
      WRITE(40,'('' Z  '')')
      WRITE(40,'(8E12.5)') AZ       ! Z
      WRITE(40,'('' ZZ '')')
      WRITE(40,'(8E12.5)') ZZ
      WRITE(40,'('' DZ '')')
      WRITE(40,'(8E12.5)') DZ
      WRITE(40,'('' DZZ'')')
      WRITE(40,'(8E12.5)') DZZ
C--- 2D ---
      WRITE(40,'(''ALON '')')
      WRITE(40,'(8E12.5)') ALON     ! east_e
      WRITE(40,'(''ALAT '')')
      WRITE(40,'(8E12.5)') ALAT     ! north_e
      WRITE(40,'('' H   '')')
      WRITE(40,'(8E12.5)') H
C--- 3D ---
      WRITE(40,'('' T   '')')
      WRITE(40,'(8E12.5)') T        ! tclim
      WRITE(40,'('' S   '')')
      WRITE(40,'(8E12.5)') S        ! sclim
      WRITE(40,'(''RMEAN'')')
      WRITE(40,'(8E12.5)') RMEAN
C--- CONSTANT WIND STRESS 
      WRITE(40,'(''WUSUR'')')
      WRITE(40,'(8E12.5)') WUSURF 
      WRITE(40,'(''WVSUR'')')
      WRITE(40,'(8E12.5)') WVSURF 
c      
C*******************************************************************
C     WRITE GRID & IC FOR PLOTTING WITH MATLAB
C     (files 43, 44, 45, 46, 47)
C*******************************************************************
C
       WRITE(43,'(3I4,16F8.4)') IM,JM,KB,(AZ(K),K=1,KB) 
      DO 221 I=1,IM 
      DO 221 J=1,JM 
       WRITE(44,'(2I4,2F9.3,F8.2,16F7.2)') 
     1  I,J,ALON(I,J),ALAT(I,J),H(I,J),(T(I,J,K),K=1,KB)
       WRITE(46,'(2I4,4F10.4)') 
     1  I,J,ALON(I,J),ALAT(I,J),UW(I,J),VW(I,J)
       WRITE(47,'(2I4,4F10.4)') 
     1  I,J,ALON(I,J),ALAT(I,J),DX(I,J)*0.001,DY(I,J)*0.001
  221 CONTINUE
C
C-----------------------------------------------------------C
C  
      STOP
      END
C
C ------------------------------------------------------------------
      SUBROUTINE TANDS 
C
c     This program reads Temperature & Salinity of the N. Atl. at
c      1x1 DEG. GRID ( LEVITUS ANNUAL CLIM.)
      INCLUDE 'gridcom'
      PARAMETER( KS=33)
      COMMON/OUT40/Z(KB),ZZ(KB),DZ(KB),DZZ(KB),ALON(IM,JM),ALAT(IM,JM),
     1    DX(IM,JM),DY(IM,JM),H(IM,JM),FSM(IM,JM),DUM(IM,JM),
     2    DVM(IM,JM),ART(IM,JM),ARU(IM,JM),ARV(IM,JM),COR(IM,JM),
     3    T(IM,JM,KB),S(IM,JM,KB),TMEAN(IM,JM,KB),SMEAN(IM,JM,KB),
     4    RMEAN(IM,JM,KB),TBE(JM,KB),TBW(JM,KB),TBS(IM,KB),
     5    SBE(JM,KB),SBW(JM,KB),SBS(IM,KB)
      DIMENSION ZS(KS),TB(IM,JM,KS),SB(IM,JM,KS)
      DIMENSION TAA(KS),SAA(KS),AREA(KS)
C LEVITUS VERTICAL LEVELS
      DATA ZS /0.,10.,20.,30.,50.,75.,100.,125.,150.,200.,
     2  250.,300.,400.,500.,600.,700.,800.,900.,1000.,1100.,
     3 1200.,1300.,1400.,1500.,1750.,2000.,2500.,3000.,3500.,4000.,
     4 4500.,5000.,5500./
      data n1,n2,n3,n4,n5/1,3,6,11,15/
C  
C
      DO 10 I=1,IM
      DO 10 J=1,JM
      DO 10 K=1,KB
       S(I,J,K)=0.
       T(I,J,K)=0.
       SB(I,J,K)=0.
       TB(I,J,K)=0.
 10   CONTINUE

C
C   MAPLEV  reads and spreads 1 deg X 1 deg and 33 level data on 
C   the curvilinear 33 level grid. 
C   ZTOSIG places the 33 level data on the KB level sigma coordinate
C   grid.
C
C levitus
      CALL MAPLEV(11,ZS,TB,SB,IM,JM,KS)
C
      WRITE(6,'(''  DEPTH         AREA      TAVE      SAVE'')')
      DO 20 K=1,KS
       TAA(K)=0.
       SAA(K)=0.
       AREA(K)=0.
 20   CONTINUE
      DO 75 K=1,KS
      P=1.025*ZS(k)
      DO 70 J=1,JM
      DO 70 I=1,IM
      IF(TB(I,J,K).EQ.0.) GOTO 70
c
c --- convert to POTENTIAL TEMPERATURE
c
      TB(I,J,K)=POTEM(TB(I,J,K),SB(I,J,K),P)
c     IF((ZS(K)+1.001).LT.H(I,J)) THEN
        TAA(K)=TAA(K)+TB(I,J,K)*ART(I,J)
        SAA(K)=SAA(K)+SB(I,J,K)*ART(I,J)
        AREA(K)=AREA(K)+ART(I,J)
c     ENDIF
   70 CONTINUE
c
c --- calc. area averaged values & save them in TB,SB
c
      TAA(K)=TAA(K)/AREA(K)
      SAA(K)=SAA(K)/AREA(K)
      WRITE(6,'(1X,F10.3,E12.3,2F10.3)') ZS(K),AREA(K),TAA(K),SAA(K)
   75 CONTINUE
c
c --- T,S are the 3D IC fields interp. to SIG grid
c
      CALL ZTOSIG(ZS,TB,ZZ,H,T,IM,JM,KS,KB)
      CALL ZTOSIG(ZS,SB,ZZ,H,S,IM,JM,KS,KB)
C
      DO 80 K=1,KS
       DO 80 I=1,IM
       DO 80 J=1,JM
      TB(I,J,K)=TAA(K)      
      SB(I,J,K)=SAA(K)      
   80 CONTINUE
c
c --- calc. area avr. density RMEAN (func(z) on sig. layers)
c
      CALL ZTOSIG(ZS,TB,ZZ,H,TMEAN,IM,JM,KS,KB)
      CALL ZTOSIG(ZS,SB,ZZ,H,SMEAN,IM,JM,KS,KB)
      CALL DENS(SMEAN,TMEAN,RMEAN,FSM,ZZ,H,IM,JM,KB)
C
       DO 99 I=1,IM
       DO 99 J=1,JM
       DO 99 K=1,KB
      TMEAN(I,J,K)=T(I,J,K)   
      SMEAN(I,J,K)=S(I,J,K)   
 99    CONTINUE
C
      RETURN
      END
C
C ------------------------------------------------------------------
      SUBROUTINE MAPLEV(NU,ZS,TB,SB,IMM,JMM,KS)
C read levitus N. Atlantic clim (1x1 deg. 79W-70W,31N-40N) 
C
       INCLUDE 'gridcom'
C --- data grid
      PARAMETER(IS=10,JS=10)
      COMMON/OUT40/Z(KB),ZZ(KB),DZ(KB),DZZ(KB),ALON(IM,JM),ALAT(IM,JM),
     1    DX(IM,JM),DY(IM,JM),H(IM,JM),FSM(IM,JM),DUM(IM,JM),
     2    DVM(IM,JM),ART(IM,JM),ARU(IM,JM),ARV(IM,JM),COR(IM,JM),
     3    T(IM,JM,KB),S(IM,JM,KB),TMEAN(IM,JM,KB),SMEAN(IM,JM,KB),
     4    RMEAN(IM,JM,KB),TBE(JM,KB),TBW(JM,KB),TBS(IM,KB),
     5    SBE(JM,KB),SBW(JM,KB),SBS(IM,KB)
      DIMENSION ZS(KS),TB(IMM,JMM,KS),SB(IMM,JMM,KS),TPS(IM,JM)
      DIMENSION BLON(IS,JS),BLAT(IS,JS),TS(IS,JS),FTS(IS,JS)
      DIMENSION TTB(IS,JS,KS),SSB(IS,JS,KS)
      CHARACTER*10 TTSS(2)
      CHARACTER*6 TIT
      DATA TTSS /'TEMP. DATA','SAL. DATA '/
      data k1,k2,k3,k4,k5/1,8,16,24,33/
C ---- south-western most data point of input data
      wlon=-79.
      slat=31.
C
C
      WRITE(6,'(//)')
C
C read LON/LAT of data 
C
      READ(NU,'(A6)') TIT
       DO 101 J=JS,1,-1
 101  READ(NU,'(10F7.2)') (BLON(I,J),I=1,IS)
      READ(NU,'(A6)') TIT
       DO 102 J=JS,1,-1
 102  READ(NU,'(10F7.2)') (BLAT(I,J),I=1,IS)
C
C read T,S data on 33 levels (land=-100) 
C
       DO 104 K=1,KS
      READ(NU,'(A6)') TIT
       DO 103 J=JS,1,-1
 103  READ(NU,'(10F7.2)') (TTB(I,J,K),I=1,IS)
      READ(NU,'(A6)') TIT
       DO 104 J=JS,1,-1
 104  READ(NU,'(10F7.2)') (SSB(I,J,K),I=1,IS)
C
      DO 10 ITS=1,2
      DO 10 K=1,KS
        DO I=1,IS
        DO J=1,JS
      IF(ITS.EQ.1)TS(I,J)=TTB(I,J,K)
      IF(ITS.EQ.2)TS(I,J)=SSB(I,J,K)
        END DO
        END DO
C
C put latitudal mean value on land points, 
C to eliminate bound. prob. when interp. to model grid.
C 
      DO 60 J=1,JS
       TSAV=0.
       ISAV=0
      DO 58 I=1,IS
        IF(TS(I,J).GT.0.) THEN
       ISAV=ISAV+1
       TSAV=TSAV+TS(I,J)
        ENDIF
 58   CONTINUE
       IF(ISAV.NE.0) THEN 
        TSAV=TSAV/ISAV
       ELSE
c --- if no data available at this latitude, put value of deepest point
        TSAV=TS(10,2)
       ENDIF
      DO 59 I=1,IS
        IF(TS(I,J).LT.0.) TS(I,J)=TSAV
 59   CONTINUE
 60   CONTINUE
      CALL PRXY(' TS ',ZS(K),TS,IS,1,JS,1,0.01)
C
c     WRITE(6,61) k,ZS(K),TSAV
c61   format(2X,'LEVEL K= ',I3,'  DEPTH Z= ',F8.0,' TAV= ',F8.2)
C
c ********* find data points and interpolate horizontally******
c
c     TPS=1. 
      do 30 i=1,im
      do 30 j=1,jm
        IF(H(I,J).LE.1.)GO TO 30
      x=alon(i,j)
      y=alat(i,j)
      ii=(x-wlon)+1
      jj=(y-slat)+1
c
      x1=wlon+(II-1)  
      y1=slat+(JJ-1)  
      f1=ts(ii,jj)
c
      x2=x1
      y2=y1+1.  
      f2=ts(ii,jj+1)
c
      x3=x1+1.  
      y3=y1+1.  
      f3=ts(ii+1,jj+1)
c
      x4=x1+1.  
      y4=y1
      f4=ts(ii+1,jj)
C
      CALL BLINT(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2,f3,f4,x,y,f)
C
      IF(ITS.EQ.1)TB(I,J,K)=f*FSM(I,J)
      IF(ITS.EQ.2)SB(I,J,K)=f*FSM(I,J)
c
  30  CONTINUE
  10  CONTINUE
C
C     WRITE(12) BLON,BLAT,TTB,SSB
C     IF(2.GT.1) STOP
C
      RETURN
      END
C
C
C ------------------------------------------------------------------
      SUBROUTINE WIND(UW,VW,IMM,JMM)
C
C read wind velocity climatology (1x1 deg. 79W-70W,31N-40N) 
C and interpolate to model grid
C
       INCLUDE 'gridcom'
C --- data grid
      PARAMETER(IS=10,JS=10)
      COMMON/OUT40/Z(KB),ZZ(KB),DZ(KB),DZZ(KB),ALON(IM,JM),ALAT(IM,JM),
     1    DX(IM,JM),DY(IM,JM),H(IM,JM),FSM(IM,JM),DUM(IM,JM),
     2    DVM(IM,JM),ART(IM,JM),ARU(IM,JM),ARV(IM,JM),COR(IM,JM),
     3    T(IM,JM,KB),S(IM,JM,KB),TMEAN(IM,JM,KB),SMEAN(IM,JM,KB),
     4    RMEAN(IM,JM,KB),TBE(JM,KB),TBW(JM,KB),TBS(IM,KB),
     5    SBE(JM,KB),SBW(JM,KB),SBS(IM,KB)
      DIMENSION BLON(IS,JS),BLAT(IS,JS),TS(IS,JS),FTS(IS,JS)
      DIMENSION UUB(IS,JS),VVB(IS,JS),UW(IMM,JMM),VW(IMM,JMM)
      CHARACTER*6 TIT
C ---- south-western most data point of input data
      wlon=-79.
      slat=31.
C
C
      WRITE(6,'(//)')
C
C read LON/LAT of data 
C
      READ(12,'(A6)') TIT
       DO 101 J=JS,1,-1
 101  READ(12,'(10F7.2)') (BLON(I,J),I=1,IS)
      READ(12,'(A6)') TIT
       DO 102 J=JS,1,-1
 102  READ(12,'(10F7.2)') (BLAT(I,J),I=1,IS)
C
C read U,V data (land=-100) 
C
      READ(12,'(A6)') TIT
       DO 103 J=JS,1,-1
 103  READ(12,'(10F7.2)') (UUB(I,J),I=1,IS)
      READ(12,'(A6)') TIT
       DO 104 J=JS,1,-1
 104  READ(12,'(10F7.2)') (VVB(I,J),I=1,IS)
C
      DO 10 ITS=1,2
        DO I=1,IS
        DO J=1,JS
      IF(ITS.EQ.1)TS(I,J)=UUB(I,J)
      IF(ITS.EQ.2)TS(I,J)=VVB(I,J)
        END DO
        END DO
C
C put latitudal mean value on land points, 
C to eliminate bound. prob. when interp. to model grid.
C 
      DO 60 J=1,JS
       TSAV=0.
       ISAV=0
      DO 58 I=1,IS
        IF(TS(I,J).GT.-99.) THEN
       ISAV=ISAV+1
       TSAV=TSAV+TS(I,J)
        ENDIF
 58   CONTINUE
       IF(ISAV.NE.0) THEN 
        TSAV=TSAV/ISAV
       ELSE
c --- if no data available at this latitude
        TSAV=0.
       ENDIF
      DO 59 I=1,IS
        IF(TS(I,J).LT.-99.) TS(I,J)=TSAV
 59   CONTINUE
 60   CONTINUE
      CALL PRXY(' UV ',1.*ITS,TS,IS,1,JS,1,0.01)
C
c ********* find data points and interpolate horizontally******
c
c     TPS=1. 
      do 30 i=1,im
      do 30 j=1,jm
        IF(H(I,J).LE.1.)GO TO 30
      x=alon(i,j)
      y=alat(i,j)
      ii=(x-wlon)+1
      jj=(y-slat)+1
c
      x1=wlon+(II-1)  
      y1=slat+(JJ-1)  
      f1=ts(ii,jj)
c
      x2=x1
      y2=y1+1.  
      f2=ts(ii,jj+1)
c
      x3=x1+1.  
      y3=y1+1.  
      f3=ts(ii+1,jj+1)
c
      x4=x1+1.  
      y4=y1
      f4=ts(ii+1,jj)
C
      CALL BLINT(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2,f3,f4,x,y,f)
C
      IF(ITS.EQ.1)UW(I,J)=f*FSM(I,J)
      IF(ITS.EQ.2)VW(I,J)=f*FSM(I,J)
c
  30  CONTINUE
  10  CONTINUE
C
      RETURN
      END
C ------------------------------------------------------------------
      SUBROUTINE ZTOSIG(ZS,TB,ZZ,H,T,IM,JM,KS,KB)
C
C   VERTICAL INTERPOLATION
C
      DIMENSION ZS(KS),TB(IM,JM,KS),ZZ(KB),H(IM,JM),T(IM,JM,KB),
     1          TIN(KS),TOUT(KB),ZZH(KB)
      DO 40 I=1,IM
      DO 40 J=1,JM
        IF(H(I,J).LE.1.0)GO TO 40
C       IF(H(I,J).LT.1.00001)GO TO 40
C-- special interp on z-lev for cases of no data because H smoothing
      DO 45 K=1,KS
      TIN(K)=TB(I,J,K)
      IF(ZS(K).LE.H(I,J).AND.TIN(K).LT.0.01)THEN   
      TMAX=AMAX1(TB(I-1,J,K),TB(I+1,J,K),TB(I,J-1,K),TB(I,J+1,K))
      TIN(K)=TMAX
         ENDIF
      IF(TIN(K).LT.0.01.AND.K.NE.1)TIN(K)=TIN(K-1)
   45 CONTINUE
C
      DO 50 K=1,KB
   50 ZZH(K)=-ZZ(K)*H(I,J)
C
C        VERTICAL SPLINE INTERP. 
      CALL SPLINC(ZS,TIN,KS  ,2.e30,2.e30,ZZH,TOUT,KB)
C
      IPR=IM/2
      JPR=JM/2
      IPR=11 
      JPR=18     
      IF(I.EQ.IPR.AND.J.EQ.JPR) THEN    
      WRITE(6,'(//,'' Data interpolated from z to sigma grid at I,J =''
     1    ,I5,'','',I5)') IPR,JPR
      WRITE(6,'(''    H ='',F10.1)') H(IPR,JPR)
      WRITE(6,'(1X,I5,4F10.4)') (K,ZS(K),TIN(K),ZZH(K),TOUT(K),K=1,KB)
      WRITE(6,'(1X,I5,2F10.4)') (K,ZS(K),TIN(K),K=KB+1,KS)
      ENDIF    
C
       DO K=1,KB
      T(I,J,K)=TOUT(K)
       END DO
C
   40 CONTINUE
C
      RETURN 
      END
C
C ------------------------------------------------------------------
C
      SUBROUTINE DEPTH(Z,ZZ,DZ,DZZ,KB)
CC
      DIMENSION Z(KB),ZZ(KB),DZ(KB),DZZ(KB)
      KBM1=KB-1
      KL1=4
      KL2=KB-2
C***********************************************************************
C   THIS SUBROUTINE ESTABLISHES THE VERTICAL RESOLUTION WITH LOG
C   DISTRIBUTIONS  AT THE TOP AND BOTTOM AND A LINEAR DISTRIBUTION
C   BETWEEN. CHOOSE KL1 AND KL2. THE DEFAULT KL1 = .3*KB AND KL2 = KB-2
C   YIELDS A LOG DISTRIBUTION AT THE TOP AND NONE AT THE BOTTOM.
C***********************************************************************
      BB=FLOAT(KL2-KL1)+4.
      CC=FLOAT(KL1)-2.
      DEL1=2.0/BB/EXP(.693147*FLOAT(KL1-2))
      DEL2=2.0/BB/EXP(.693147*FLOAT(KB-KL2-1))
      Z(1)=0.0
      ZZ(1)=-DEL1/2.0
      DO 3 K=2,KL1
      Z(K)=-DEL1*EXP(.693147*FLOAT(K-2))
    3 ZZ(K)=-DEL1*EXP(.693147*(FLOAT(K)-1.50))
      DO 4 K=KL1,KL2
      Z(K)=-(FLOAT(K)-CC)/BB
    4 ZZ(K)=-(FLOAT(K)-CC+0.50)/BB
      DO 5 K=KL2,KBM1
      Z(K)=(1.0-DEL2*EXP(.693147*FLOAT(KB-K-1)))*(-1.0)
    5 ZZ(K)=(1.0-DEL2*EXP(.693147*(FLOAT(KB-K)-1.50)))*(-1.0)
      Z(KB)=-1.0
      ZZ(KBM1)=-1.0*(1.0-DEL2/2.0)
      ZZ(KB)=-1.0*(1.0+DEL2/2.0)
      DO 6 K=1,KBM1
      DZ(K)=Z(K)-Z(K+1)
    6 DZZ(K)=ZZ(K)-ZZ(K+1)
      RETURN
      END
C
      FUNCTION POTEM(T,S,P)
C***********************************************************************
C
C     POTENTIAL TEMPERATURE FUNCTION
C     BASED ON FOFONOFF AND FROESE (1958) AS SHOWN IN "THE SEA" VOL. 1,
C     PAGE 17, TABLE IV
C     INPUT IS TEMPERATURE, SALINITY, PRESSURE (OR DEPTH)
C     UNITS ARE DEG.C., PPT, DBARS (OR METERS)
C
C***********************************************************************
CC
      B1=-1.60E-5*P
      B2=1.014E-5*P*T
      T2=T*T
      T3=T2*T
      B3=-1.27E-7*P*T2
      B4=2.7E-9*P*T3
      B5=1.322E-6*P*S
      B6=-2.62E-8*P*S*T
      S2=S*S
      P2=P*P
      B7=4.1E-9*P*S2
      B8=9.14E-9*P2
      B9=-2.77E-10*P2*T
      B10=9.5E-13*P2*T2
      B11=-1.557E-13*P2*P
      POTEM=B1+B2+B3+B4+B5+B6+B7+B8+B9+B10+B11
      POTEM=T-POTEM
      RETURN
      END
C
C ------------------------------------------------------------------
      SUBROUTINE DENS(S,T,RHO,FSM,ZZ,H,IM,JM,KB)
      DIMENSION S(IM,JM,KB),T(IM,JM,KB),RHO(IM,JM,KB),FSM(IM,JM),
     1   ZZ(KB),H(IM,JM)
C
C
C         THIS SUBROUTINE COMPUTES DENSITY-1.0  
C         T = POTENTIAL TEMPERATURE
C
        GRAV=9.807
      DO 1 K=1,KB-1
      DO 1 J=1,JM
      DO 1 I=1,IM
      TR=T(I,J,K)
      SR=S(I,J,K)
C         Here, the (approximate) pressure is in units of bars.
      P=-GRAV*1.025*ZZ(K)*H(I,J)*0.01
C
      RHOR = 999.842594 + 6.793952E-2*TR
     1        - 9.095290E-3*TR**2 + 1.001685E-4*TR**3
     2        - 1.120083E-6*TR**4 + 6.536332E-9*TR**5
C
      RHOR = RHOR + (0.824493 - 4.0899E-3*TR
     1       + 7.6438E-5*TR**2 - 8.2467E-7*TR**3
     2       + 5.3875E-9*TR**4) * SR
     3       + (-5.72466E-3 + 1.0227E-4*TR
     4       - 1.6546E-6*TR**2) *(ABS(SR))**1.5
     5       + 4.8314E-4 * SR**2
C
      CR =1449.1+.0821*P+4.55*TR-.045*TR**2+1.34*(SR-35.)
      RHOR=RHOR + 1.E5*P/CR**2*(1.-2.0*P/CR**2)
C    
      RHO(I,J,K)=(RHOR-1000.)*1.E-3*FSM(I,J)
    1 CONTINUE
C
      DO 3 J=1,JM
      DO 3 I=1,IM
    3 RHO(I,J,KB)=0.E0
      RETURN
      END
C
C ------------------------------------------------------------------
      SUBROUTINE BORDER(IMm,JMm,X,Y)
C
C------------------------------------------------------------------------
C      A FEW POINTS OF THE BORDER ARE SPECIFIED IN DATA STATEMENTS.
C      THE REMAINING POINTS ARE INTERPOLATED
C------------------------------------------------------------------------
      INCLUDE 'gridcom'
C     PARAMETER (NT=6,NB=2,NL=3,NR=3)
      PARAMETER (Nn=10)                
      REAL LFTX,LFTY
      DIMENSION ITOP(Nn),       IBOT(Nn),
     &          TOPX(Nn),TOPY(Nn),BOTX(Nn),BOTY(Nn),
     &          JRHT(Nn),        JLFT(Nn),
     &          RHTX(Nn),RHTY(Nn),LFTX(Nn),LFTY(Nn)
      DIMENSION RITOP(Nn),RIBOT(Nn),RJLFT(Nn),RJRHT(Nn),
     &          RI(IM),RJ(JM)
C --- Fix an old elusive bug (original code gave error if IM<JM)
C     DIMENSION XXX(IM),YYY(IM)
      DIMENSION XXI(IM),YYI(IM)
      DIMENSION XXJ(JM),YYJ(JM)
C
      REAL X(IMm,JMm),Y(IMm,JMm)
C
C
      DATA ZER0/1.E-38/
      DATA PI/3.141593/
      rad=pi/180.
C
C ---- read boundary parameters from file -------------------
C
       INCLUDE 'gridborder'
C
C
      DO 1000 J=1,JM
      DO 1000 I=1,IM
      X(I,J)=0.
 1000 Y(I,J)=0.
      DO 1002 I=1,IM
      XXI(I)=0.
 1002 YYI(I)=0.
C
C------ FILL IN BORDER POINTS ------------------------------
      DO 10 I=1,IM
   10 RI(I)=I
      DO 20 I=1,NT
c---- write location of grid points for MATLAB grid plotting
      WRITE(45,'(2F10.3)') TOPX(I),TOPY(I)
   20 RITOP(I)=1.+FLOAT(ITOP(I)-1)*FLOAT(IM-1)/FLOAT(100-1)
      WRITE(6,'(''TOP '',I5,2F10.3)') (J,RITOP(J),TOPX(J),J=1,NT)
      CALL SPLINC (RITOP,TOPX,NT,2.E30,2.E30,RI,XXI,IM)
      CALL SPLINC (RITOP,TOPY,NT,2.E30,2.E30,RI,YYI,IM)
      DO 30 I=1,IM
         X(I,JM)=XXI(I)
         Y(I,JM)=YYI(I)
  30  CONTINUE
C
      DO 40 I=1,NB
c---- write location of grid points for MATLAB grid plotting
      WRITE(45,'(2F10.3)') BOTX(I),BOTY(I)
   40 RIBOT(I)=1.+FLOAT(IBOT(I)-1)*FLOAT(IM-1)/FLOAT(100-1)
      WRITE(6,'(''BOT '',I5,2F10.3)') (J,RIBOT(J),BOTX(J),J=1,NB)
      CALL SPLINC (RIBOT,BOTX,NB,2.E30,2.E30,RI,XXI,IM)
      CALL SPLINC (RIBOT,BOTY,NB,2.E30,2.E30,RI,YYI,IM)
      DO 50 I=1,IM
         X(I,1)=XXI(I)
         Y(I,1)=YYI(I)
  50  CONTINUE
C
      DO 60 J=1,JM
   60 RJ(J)=J
      DO 70 J=1,NL
c---- write location of grid points for MATLAB grid plotting
      WRITE(45,'(2F10.3)') LFTX(J),LFTY(J)
   70 RJLFT(J)=1.+FLOAT(JLFT(J)-1)*FLOAT(JM-1)/FLOAT(100-1)
      WRITE(6,'(''LFT '',I5,2F10.3)') (J,RJLFT(J),LFTX(J),J=1,NL)
      CALL SPLINC (RJLFT,LFTX,NL,2.E30,2.E30,RJ,XXJ,JM)
      CALL SPLINC (RJLFT,LFTY,NL,2.E30,2.E30,RJ,YYJ,JM)
      DO 80 I=1,JM
         X(1,I)=XXJ(I)
         Y(1,I)=YYJ(I)
  80  CONTINUE
C
      XTIM=X(IM,JM)-X(IM-1,JM) 
      YTIM=Y(IM,JM)-Y(IM-1,JM) 
      DO 90 J=1,NR
c---- write location of grid points for MATLAB grid plotting
      WRITE(45,'(2F10.3)') RHTX(J),RHTY(J)
   90 RJRHT(J)=1.+FLOAT(JRHT(J)-1)*FLOAT(JM-1)/FLOAT(100-1)
      WRITE(6,'(''RIT '',I5,2F10.3)') (J,RJRHT(J),RHTX(J),J=1,NR)
      CALL SPLINC (RJRHT,RHTX,NR,2.E30,-YTIM,RJ,XXJ,JM)
      CALL SPLINC (RJRHT,RHTY,NR,2.E30,XTIM,RJ,YYJ,JM)
      DO 100 I=1,JM
         X(IM,I)=XXJ(I)
         Y(IM,I)=YYJ(I)
  100 CONTINUE 
      RETURN
      END
C
C ------------------------------------------------------------------
      SUBROUTINE SPLINC(X,Y,N,YP1,YPN,XNEW,YNEW,M)
      PARAMETER (NMAX=300)
      DIMENSION X(N),Y(N),Y2(NMAX),U(NMAX),XNEW(M),YNEW(M)
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
C        
      DO 20 I =1,M
      CALL SPLINT(X,Y,Y2,N,XNEW(I),YNEW(I))
  20  CONTINUE
      RETURN
      END
C---------------------------------------------------------------------
C
      SUBROUTINE CRASH(ICRASH,IERR)
      IF(ICRASH.EQ.1)WRITE(6,10)IERR
 10   FORMAT(' SEPELI FAILED WHILE COMPUTING X SOLUTION: IERR =',I3)
      IF(ICRASH.EQ.2)WRITE(6,20)IERR
 20   FORMAT(' SEPELI FAILED WHILE COMPUTING Y SOLUTION: IERR =',I3)
      IF(ICRASH.EQ.3)WRITE(6,30)
 30   FORMAT(' WRITE ERROR WHILE OUTPUTING SOLUTION')
      IF(ICRASH.EQ.4)WRITE(6,40)
 40   FORMAT(' TRACKING ERROR IN RECT')
      IF(ICRASH.EQ.5)WRITE(6,50)
 50   FORMAT(' BAD XA INPUT IN SPLINT')
      CALL EXIT
      RETURN
      END
C
C ------------------------------------------------------------------
      SUBROUTINE COFX(X,AFUN,BFUN,CFUN)
C      SUBROUTINE TO COMPUTE THE COEFFICIENTS OF THE ELLIPTIC EQUATION
C      SOLVED TO FILL IN THE GRID.  IT IS PASSED (ALONG WITH COFY) TO
C      THE ELLIPTIC SOLVER SEPELI.
       INCLUDE 'gridcom'
      PARAMETER (     L2=IM-1   , M2=JM-1    )
      PARAMETER (     MM=M2/2   , LM=L2/2    )
      COMMON / XIEJ / SXI(0:L2),SETA(0:M2)
      I = NINT(X)
      DXDI = 0.5*(SXI(I+1)-SXI(I-1))
      AFUN = 1./DXDI**2
      BFUN = (-SXI(I+1)+2.*SXI(I)-SXI(I-1))/DXDI**3
      CFUN = 0.
      RETURN
      END
C ------------------------------------------------------------------
      SUBROUTINE COFY(Y,DFUN,EFUN,FFUN)
       INCLUDE 'gridcom'
      PARAMETER (     L2=IM-1   , M2=JM-1    )
      PARAMETER (     MM=M2/2   , LM=L2/2    )
      COMMON / XIEJ / SXI(0:L2),SETA(0:M2)
      J = NINT(Y)
      DEDJ = 0.5*(SETA(J+1)-SETA(J-1))
      DFUN = 1./DEDJ**2
      EFUN = (-SETA(J+1)+2.*SETA(J)-SETA(J-1))/DEDJ**3
      FFUN = 0.
      RETURN
      END
C     
C ------------------------------------------------------------------
      SUBROUTINE PRXY (LABEL,TIME,A,IM,ISKP,JM,JSKP,SCALA)
C
C     THIS WRITES A 2-D FIELD
C     TIME=TIME IN DAYS
C     A = ARRAY(IM,JM) TO BE PRINTED
C     ISKP=PRINT SKIP FOR I
C     JSKP=PRINT SKIP FOR J
C     SCALE=DIVISOR FOR VALUES OF A
      DIMENSION A(IM,JM),NUM(IM),LINE(IM)
      CHARACTER LABEL*(*)
C
      SCALE=SCALA
      IF (SCALE.GT.1.E-10) GO TO 160
      AMX=1.E-10
      DO 150 J=1,JM,JSKP
      DO 150 I=1,IM,ISKP
      AMX=MAX1(ABS(A(I,J)),AMX)
  150 CONTINUE
      SCALE=10.**(INT(LOG10(AMX)+1.E2)-103)
  160 CONTINUE
      SCALEI=1./SCALE
      WRITE(6,170) LABEL
  170 FORMAT(1X,A40)
      WRITE(6,180) TIME,SCALE
  180 FORMAT(' TIME =',F9.4,' DAYS     MULTIPLY ALL VALUES BY',1PE10.3)
      DO 190 I=1,IM
  190 NUM(I)=I
      IB=1
C
  200 CONTINUE
      IE=IB+23*ISKP
      IF(IE.GT.IM) IE=IM
      WRITE(6,210) (NUM(I),I=IB,IE,ISKP)
  210 FORMAT(/,2X,24I5,/)
      DO 260 J=1,JM,JSKP
      JWR=JM+1-J
      DO 220 I=IB,IE,ISKP
  220 LINE(I)=INT(SCALEI*A(I,JWR))
        IF(IE.le.IM) THEN
        WRITE(6,240) JWR,(LINE(I),I=IB,IE,ISKP)
        ELSE
        LINE(IM)=INT(SCALEI*A(IM,JWR))
        WRITE(6,240) JWR,(LINE(I),I=IB,IE,ISKP),LINE(IM)
        ENDIF
  240 FORMAT(1X,I3,24I5)
  260 CONTINUE
      IF(JWR.NE.1) THEN
      JWR=1
      DO 270 I=IB,IE,ISKP
  270 LINE(I)=INT(SCALEI*A(I,JWR))
        IF(IE.le.IM) THEN
        WRITE(6,240) JWR,(LINE(I),I=IB,IE,ISKP)
        ELSE
        LINE(IM)=INT(SCALEI*A(IM,JWR))
        WRITE(6,240) JWR,(LINE(I),I=IB,IE,ISKP),LINE(IM)
        ENDIF
      ENDIF
      WRITE(6,280)
  280 FORMAT(//)
      IF(IE.GE.IM) RETURN
      IB=IB+24*ISKP
      GO TO 200
      END
C
C ------------------------------------------------------------------
      SUBROUTINE PRXYZ (LABEL,TIME,A,IUO,IM,ISKP,JM,JSKP,KB,SCALA)
CC
C     WRITE HORIZONTAL LAYERS OF A 3-D FIELD (9/24/87)
C       TIME=TIME IN DAYS
C       A= ARRAY(IM,JM,KB) TO BE PRINTED
C       ISPL=PRINT SKIP FOR I
C       JSPL=PRINT SKIP FOR J
C       SCALE=DIVISOR FOR VALUES OF A
      DIMENSION A(IM,JM,KB),NUM(IM),LINE(IM),KP(4)
      CHARACTER*(*) LABEL
      DATA KE,KP/4,1,3,5,9/
      SCALE=SCALA
      IF (SCALE.GT.1.E-10) GO TO 160
      AMX=1.E-10
      DO 150 KM=1, KE
      K=KP(KM)
      DO 150 J=1,JM,JSKP
      DO 150 I=1,IM,ISKP
      AMX=MAX1(ABS(A(I,J,K)),AMX)
  150 CONTINUE
      SCALE=10.**(INT(LOG10(AMX)+1.E2)-103)
  160 CONTINUE
      SCALEI=1.0/SCALE
      WRITE(IUO,51) LABEL
   51 FORMAT(1X,A40)
      WRITE(IUO,52) TIME,SCALE
   52 FORMAT(' TIME = ',F8.2,'    MULTIPLY ALL VALUES BY',1PE10.3)
      DO 190 I=1,IM
  190 NUM(I)=I
      DO 499 KM=1,KE
      K=KP(KM)
      WRITE(IUO,53) K
   53 FORMAT(3X,/7H LAYER ,I2)
      IB=1
      IE=IB+23*ISKP
  200 CONTINUE
      IF(IE.GT.IM) IE=IM
      WRITE(IUO,54) (NUM(I),I=IB,IE,ISKP)
   54 FORMAT(/,2X,24I5,/)
      DO 230 J=JM,4,-JSKP
      DO 220 I=IB,IE,ISKP
  220 LINE(I)=AINT(SCALEI*A(I,J,K))
      WRITE(IUO,'(1X,I3,24(I5))') J,(LINE(I),I=IB,IE,ISKP)
  230 CONTINUE
      DO 250 J=3,1,-1
      DO 240 I=IB,IE,ISKP
  240 LINE(I)=AINT(SCALEI*A(I,J,K))
      WRITE(IUO,'(1X,I3,24(I5))') J,(LINE(I),I=IB,IE,ISKP)
  250 CONTINUE
      WRITE(IUO,56)
   56 FORMAT(//)
      IF(IE.EQ.IM) GO TO 499
      IB=IB+24*ISKP
      IE=IB+23*ISKP
      GO TO 200
  499 CONTINUE
      RETURN
      END
C
C ------------------------------------------------------------------
      SUBROUTINE RECT(Z,N,N1,N2,N3,N4)
C
C     THIS SUBROUTINE IS TAKEN DIRECTLY FROM IVES,D.C. AND
C     R.M.ZACHARIAS "CONFORMAL MAPPING AND ORTHOGONAL GRID GENERATION"
C     AIAA-87-2057.
C
      IMPLICIT REAL (A-H,O-Z)
c     COMPLEX Z(1),Z0,ZD
c     REAL R(1000),T(1000)
      COMPLEX Z(N),Z0,ZD
      REAL R(N),T(N)
      PI = 4.D0*DATAN(1.D0)
      DO 7 I=1,N
      IM = N-MOD(N-I+1,N)
      IP = 1+MOD(I,N)
      ZD = (Z(IM)-Z(I))/(Z(IP)-Z(I))
      ALPHA = ATAN2(AIMAG(ZD),REAL(ZD))
C     ALPHA = ATAN2(AIMAG(ZD),REAL(ZD))
      IF(ALPHA.LT.0)ALPHA = ALPHA+PI+PI
      PWR = PI/ALPHA
      IF(I.NE.N1.AND.I.NE.N2.AND.I.NE.N3.AND.I.NE.N4)GOTO 2
      ZD = (Z(IM)-Z(I))/CABS(Z(IM)-Z(I))
      DO 1 J=1,N
c1       Z(J) = CMPLX(0.E0,1.E0)*Z(J)/ZD
 1       Z(J) = CMPLX(0.E0,1.E0)*Z(J)/(ZD+1.E-10)
      PWR = PWR/2.
 2    PMIN = 100.
      PMAX =-100.
      TP   = 0.
      DO 5 J=2,N
         ZD = Z(MOD(J+I-2,N)+1)-Z(I)
         R(J) = CABS(ZD)
         T(J) = ATAN2(AIMAG(ZD),REAL(ZD))-PI-PI-PI-PI-PI-PI
         DO 3 K=1,7
            IF(ABS(T(J)-TP).LT.PI)GOTO 4
 3          T(J) = T(J)+PI+PI
C        PAUSE ' WARNING - TRACKING ERROR '
         CALL CRASH(4,0)
 4       TP = T(J)
         PMAX = AMAX1(PMAX,T(J)*PWR)
 5       PMIN = AMIN1(PMIN,T(J)*PWR)
      PWR = AMIN1(PWR,1.98E0*PI*PWR/(PMAX-PMIN))
      Z(I) = CMPLX(0.E0,0.E0)
      DO 6 J=2,N
 6       Z(MOD(J+I-2,N)+1) = R(J)**PWR*CEXP(CMPLX(0.E0,T(J)*PWR))
      ZD = 1./(Z(N2)-Z(N1))
      Z0 = Z(N1)
      DO 7 J=1,N
 7       Z(J) = (Z(J)-Z0)*ZD
      RETURN
      END       
C ------------------------------------------------------------------
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
C
C     THE FOLLOWING TWO SUBROUTINES ARE USED TO PERFORM THE CUBIC SPLINE
C     INTERPOLATION REQUIRED TO MATCH UP THE DISTRIBUTION OF POINTS ON
C     OPPOSITE SIDES OF THE TRANSFORMED PLANE RECTANGLE.  THE ROUTINES
C     ARE TAKEN FROM  PRESS,W.H., B.P.FLANNERY, S.A.TEUKOLSKY AND
C     W.T.VETTERLING: "NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING"
C     CAMBRIDGE UNIVERSITY PRESS, 1986.
      INCLUDE 'gridcom'
      PARAMETER (     L2=IM-1   , M2=JM-1    )
      PARAMETER (     MM=M2/2   , LM=L2/2    )
      PARAMETER (     NMAX=M2+L2             )
C UNIX 6 indicated 3 errors in the next line
C     DIMENSION X(N),Y(N),Y2(N),U(NMAX)
C
      DIMENSION X(NMAX),Y(NMAX),Y2(NMAX),U(NMAX)
C
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END          
C ------------------------------------------------------------------
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) CALL CRASH(5,0)
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END
C ------------------------------------------------------------------
      SUBROUTINE Z1(S1,X,Y,Z)
C
C     SUBROUTINES WHICH SPECIFY THE BOUNDARIES OF THE PHYSICAL
C     DOMAIN.  THEY ARE DEFINED AS FUNCTIONS OF THE VARIABLE S WHICH
C     RANGES FROM 0 TO 1 ON EACH BOUNDARY, PROCEEDING ANTI-CLOCKWISE.
C     COMPLEX*16 Z
      COMPLEX   Z
c     X = 0.
      pi=2.*4.*atan(1.)
      Y = 400.E+03 - 300.E+03*S1
      x=-50.E+03*sin((y-250.E+03)*pi/300.E+03)
      Z = CMPLX(X,Y)
      RETURN
C
      ENTRY Z2(S2,X,Y,Z)
      PI = 4.*ATAN(1.)
      X = 600.E+03*S2
      IF(X.LT.150.E+03)THEN
         Y = 100.E+03
      ELSE IF(X.GT.400.E+03)THEN
         Y = 0.
      ELSE
         Y = 50.E+03*(1.+COS((X-150.E+03)*PI/250.E+03))
      ENDIF
      Z = CMPLX(X,Y)
      RETURN
C
      ENTRY Z3(S3,X,Y,Z)
      PI = 4.*ATAN(1.)
      X = 500.E+03+100.E+03*COS(PI*S3)
      Y = 500.E+03*S3
      Z = CMPLX(X,Y)
      RETURN
C
      ENTRY Z4(S4,X,Y,Z)
      PI = 4.*ATAN(1.)
      X = 400.E+03*(1.-S4)
      Y = 450.E+03+50.E+03*COS(PI*S4)
      Z = CMPLX(X,Y)
      RETURN
      END
C
C ------------------------------------------------------------------
      SUBROUTINE BATH
C **********************************************************
C Subroutine to read bathymetry data on 1/12 deg. intervals
C and interpolate into the model grid.
C This example: ETOPO5 data from NGDC for 70-79W,31-40N
c ***********************************************************
      INCLUDE 'gridcom'
C --- grid of data
      PARAMETER (IB=109,JB=109)
      COMMON/OUT40/AZ(KB),ZZ(KB),DZ(KB),DZZ(KB),ALON(IM,JM),ALAT(IM,JM),
     1    DX(IM,JM),DY(IM,JM),H(IM,JM),FSM(IM,JM),DUM(IM,JM),
     2    DVM(IM,JM),ART(IM,JM),ARU(IM,JM),ARV(IM,JM),COR(IM,JM),
     3    T(IM,JM,KB),S(IM,JM,KB),TMEAN(IM,JM,KB),SMEAN(IM,JM,KB),
     4    RMEAN(IM,JM,KB),TBE(JM,KB),TBW(JM,KB),TBS(IM,KB),
     5    SBE(JM,KB),SBW(JM,KB),SBS(IM,KB)
      DIMENSION HB(IB,JB),NHB(IB)
C
c ************** read bathymetry ******************************
c
C --- south-west point of data
      wlon=-79.                                               
      slat=31.
c --- this data file starts at north-west corner (+/-=land/ocean)
        do 10 j=jb,1,-1
         read (75,*)NHB
        do 10 i=1,ib
c        HB(I,J)=AMIN1(-1.*NHB(I),0.)
         HB(I,J)=-1.*NHB(I)
  10   continue
C --- print original data
c     CALL PRXY(' TOPO DATA ',TIME,HB,IM,1,JM,1,1.)
C
c ************** find data points for interpolation ***********
c
      do 30 i=1,im
      do 30 j=1,jm
      x=alon(i,j)
      y=alat(i,j)
      ii=12.*(x-wlon)+1
      jj=12.*(y-slat)+1
cc                          if no interpolation
cc    h(i,j)=hb(ii,jj)
cc    go to 30
c
      x1=wlon+(II-1)/12.
      y1=slat+(JJ-1)/12.
      f1=hb(ii,jj)
c
      x2=x1
      y2=y1+1./12.
      f2=hb(ii,jj+1)
c
      x3=x1+1./12.
      y3=y1+1./12.
      f3=hb(ii+1,jj+1)
c
      x4=x1+1./12.
      y4=y1
      f4=hb(ii+1,jj)
c
c  Defind the coast if distance < 1/12 deg.
c
      if(f1.eq.0..or.f2.eq.0..or.f3.eq.0..or.f4.eq.0.)then
      h(i,j)=0.        
      go to 30
      endif
c
      CALL BLINT(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2,f3,f4,x,y,f)
      H(i,j)=f
c
  30  CONTINUE
c
      return
      end
C
C ------------------------------------------------------------------
      SUBROUTINE BLINT(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2,f3,f4,x,y,f)
C
C Bilinear interpolation subroutine.
C (Xi,Yi,fi) = data grid & values surounding model point (x,y)
C f = interpolated value at the model grid point.
C
      a1=x1-x2+x3-x4
      a2=-x1+x4
      a3=-x1+x2
      a4=x1-x
      b1=y1-y2+y3-y4
      b2=-y1+y4
      b3=-y1+y2
      b4=y1-y 
      A=a3*b1-a1*b3
      B=b2*a3+b1*a4-a1*b4-a2*b3
      C=-a2*b4+a4*b2
      if(ABS(A*C).gt.0.002*B**2) then
         t=(-B-sqrt(B*B-4.*A*C))/(2.*A)
         else
	 t=C/ABS(B)
      endif
  10  CONTINUE
      A=a2*b1-a1*b2
      B=b3*a2+b1*a4-a1*b4-a3*b2
      C=-a3*b4+a4*b3
      if(ABS(A*C).gt.0.002*B**2) then
         s=(-B+sqrt(B*B-4.*A*C))/(2.*A)
         else
	 s=-C/ABS(B)
      endif
  20  CONTINUE
      f=f1*(1.-t)*(1.-s)+f2*t*(1.-s)+f3*s*t+f4*(1.-t)*s
      return
      end
C
C ------------------------------------------------------------------
      SUBROUTINE SMOOTH(A,FSM,IM,JM,NITS)
C-------------------------------------------------------------------
C     THIS ROUTINE SMOOTHS DATA WITH A FIVE POINT LAPLACIAN FILTER.
C-------------------------------------------------------------------
      DIMENSION A(IM,JM),FSM(IM,JM)
      DO 100 N=1,NITS
      DO 200 J=2,JM-1
      DO 200 I=2,IM-1
C
      SMFAC=FSM(I+1,J)+FSM(I,J-1)+FSM(I-1,J)+FSM(I,J+1)+1.E-10
      A(I,J)=A(I,J)+(.5/SMFAC)
     1               *(A(I+1,J)*FSM(I+1,J)+A(I,J-1)*FSM(I,J-1)
     2                +A(I-1,J)*FSM(I-1,J)+A(I,J+1)*FSM(I,J+1)
     3                -SMFAC*A(I,J))
      A(I,J)=A(I,J)*FSM(I,J)
  200 CONTINUE
  100 CONTINUE
      RETURN
      END
C ------------------------------------------------------------------
C
      SUBROUTINE SLPMIN(H,IM,JM,FSM,SLMIN)
C
C     TOPOGRAPHIC SMOOTHER - DH/H < 2*SLMIN
C     The topography is smoothed such that the difference in 
C     depth of adjacent grid points divided by the mean depth
C     is less than than 2*SLMIN where SLMIN is a constant
C     in the range, 0.0 < SLMIN < 2.0. The adjustment is made
C     so as to conserve the combined volume of the adjacent 
C     grid points.
C
      DIMENSION H(IM,JM),FSM(IM,JM)
      DIMENSION SL(IM,JM)
      LOOPMAX=10
C
      DO 100 LOOP=1,LOOPMAX
C    sweep right
      DO 3 J=2,JM-1
      DO 1 I=2,IM-1
      IF(FSM(I,J).EQ.0..OR.FSM(I+1,J).EQ.0.) GOTO 1
      SL(I,J)=ABS(H(I+1,J)-H(I,J))/(H(I,J)+H(I+1,J))
      IF(SL(I,J).LT.SLMIN) GOTO 1
      DH=0.5*(SL(I,J)-SLMIN)*(H(I,J)+H(I+1,J))
      SN=-1.
      IF(H(I+1,J).GT.H(I,J)) SN=1.
      H(I+1,J)=H(I+1,J)-SN*DH
      H(I,J)=H(I,J)+SN*DH
   1  CONTINUE
C    sweep left
      DO 2 I=IM-1,2,-1
      IF(FSM(I,J).EQ.0..OR.FSM(I+1,J).EQ.0.) GOTO 2
      SL(I,J)=ABS(H(I+1,J)-H(I,J))/(H(I,J)+H(I+1,J))
      IF(SL(I,J).LT.SLMIN) GOTO 2
      DH=0.5*(SL(I,J)-SLMIN)*(H(I,J)+H(I+1,J))
      SN=-1.
      IF(H(I+1,J).GT.H(I,J)) SN=1.
      H(I+1,J)=H(I+1,J)-SN*DH
      H(I,J)=H(I,J)+SN*DH
   2  CONTINUE
   3  CONTINUE
C     CALL PRXY(' H AFTER left  ',0.  ,H,IM,3,JM,2,1)
C   sweep up     
      DO 13 I=2,IM-1
      DO 11 J=2,JM-1
      IF(FSM(I,J).EQ.0..OR.FSM(I,J+1).EQ.0.) GOTO 11
      SL(I,J)=ABS(H(I,J+1)-H(I,J))/(H(I,J)+H(I,J+1))
      IF(SL(I,J).LT.SLMIN) GOTO 11
      DH=0.5*(SL(I,J)-SLMIN)*(H(I,J)+H(I,J+1))
      SN=-1.
      IF(H(I,J+1).GT.H(I,J)) SN=1.
      H(I,J+1)=H(I,J+1)-SN*DH
      H(I,J)=H(I,J)+SN*DH
   11 CONTINUE
C   sweep down
      DO 12 J=JM-1,2,-1
      IF(FSM(I,J).EQ.0..OR.FSM(I,J+1).EQ.0.) GOTO 12
      SL(I,J)=ABS(H(I,J+1)-H(I,J))/(H(I,J)+H(I,J+1))
      IF(SL(I,J).LT.SLMIN) GOTO 12
      DH=0.5*(SL(I,J)-SLMIN)*(H(I,J)+H(I,J+1))
      SN=-1.
      IF(H(I,J+1).GT.H(I,J)) SN=1.
      H(I,J+1)=H(I,J+1)-SN*DH
      H(I,J)=H(I,J)+SN*DH
  12  CONTINUE
  13  CONTINUE
C     CALL PRXY(' H AFTER down ',0.  ,H,IM,3,JM,2,1.)
C
  100 CONTINUE 
      RETURN
      END
C ----------------------------------------------------------
      SUBROUTINE SEPELI (INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,
     1                   D,N,NBDCND,BDC,GAMA,BDD,XNU,COFX,COFY,GRHS,
     2                   USOL,IDMN,W,PERTRB,IERROR)
C
C DIMENSION OF           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
C ARGUMENTS              USOL(IDMN,N+1),GRHS(IDMN,N+1),
C                        W (SEE ARGUMENT LIST)
C
C LATEST REVISION        MARCH 1985
C
C PURPOSE                SEPELI SOLVES FOR EITHER THE SECOND-ORDER
C                        FINITE DIFFERENCE APPROXIMATION OR A
C                        FOURTH-ORDER APPROXIMATION TO A SEPARABLE
C                        ELLIPTIC EQUATION
C
C                          2    2
C                          AF(X)*D U/DX + BF(X)*DU/DX  + CF(X)*U +
C                          2    2
C                          DF(Y)*D U/DY  + EF(Y)*DU/DY + FF(Y)*U
C
C                          = G(X,Y)
C
C                        ON A RECTANGLE (X GREATER THAN OR EQUAL TO A
C                        AND LESS THAN OR EQUAL TO B; Y GREATER THAN
C                        OR EQUAL TO C AND LESS THAN OR EQUAL TO D).
C                        ANY COMBINATION OF PERIODIC OR MIXED BOUNDARY
C                        CONDITIONS IS ALLOWED.
C
C                        THE POSSIBLE BOUNDARY CONDITIONS ARE:
C                        IN THE X-DIRECTION:
C                        (0) PERIODIC, U(X+B-A,Y)=U(X,Y) FOR ALL
C                            Y,X (1) U(A,Y), U(B,Y) ARE SPECIFIED FOR
C                            ALL Y
C                        (2) U(A,Y), DU(B,Y)/DX+BETA*U(B,Y) ARE
C                            SPECIFIED FOR ALL Y
C                        (3) DU(A,Y)/DX+ALPHA*U(A,Y),DU(B,Y)/DX+
C                            BETA*U(B,Y) ARE SPECIFIED FOR ALL Y
C                        (4) DU(A,Y)/DX+ALPHA*U(A,Y),U(B,Y) ARE
C                            SPECIFIED FOR ALL Y
C
C                        IN THE Y-DIRECTION:
C                        (0) PERIODIC, U(X,Y+D-C)=U(X,Y) FOR ALL X,Y
C                        (1) U(X,C),U(X,D) ARE SPECIFIED FOR ALL X
C                        (2) U(X,C),DU(X,D)/DY+XNU*U(X,D) ARE
C                            SPECIFIED FOR ALL X
C                        (3) DU(X,C)/DY+GAMA*U(X,C),DU(X,D)/DY+
C                            XNU*U(X,D) ARE SPECIFIED FOR ALL X
C                        (4) DU(X,C)/DY+GAMA*U(X,C),U(X,D) ARE
C                            SPECIFIED FOR ALL X
C
C USAGE                  CALL SEPELI (INTL,IORDER,A,B,M,MBDCND,BDA,
C                                     ALPHA,BDB,BETA,C,D,N,NBDCND,BDC,
C                                     GAMA,BDD,XNU,COFX,COFY,GRHS,USOL,
C                                     IDMN,W,PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               INTL
C                          = 0 ON INITIAL ENTRY TO SEPELI OR IF ANY
C                              OF THE ARGUMENTS C,D, N, NBDCND, COFY
C                              ARE CHANGED FROM A PREVIOUS CALL
C                          = 1 IF C, D, N, NBDCND, COFY ARE UNCHANGED
C                              FROM THE PREVIOUS CALL.
C
C                        IORDER
C                          = 2 IF A SECOND-ORDER APPROXIMATION
C                              IS SOUGHT
C                          = 4 IF A FOURTH-ORDER APPROXIMATION
C                              IS SOUGHT
C
C                        A,B
C                          THE RANGE OF THE X-INDEPENDENT VARIABLE,
C                          I.E., X IS GREATER THAN OR EQUAL TO A
C                          AND LESS THAN OR EQUAL TO B.  A MUST BE
C                          LESS THAN B.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL [A,B] IS SUBDIVIDED. HENCE,
C                          THERE WILL BE M+1 GRID POINTS IN THE X-
C                          DIRECTION GIVEN BY XI=A+(I-1)*DLX
C                          FOR I=1,2,...,M+1 WHERE DLX=(B-A)/M IS
C                          THE PANEL WIDTH.  M MUST BE LESS THAN
C                          IDMN AND GREATER THAN 5.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT X=A AND X=B
C
C                          = 0 IF THE SOLUTION IS PERIODIC IN X, I.E.,
C                              U(X+B-A,Y)=U(X,Y)  FOR ALL Y,X
C                          = 1 IF THE SOLUTION IS SPECIFIED AT X=A
C                              AND X=B, I.E., U(A,Y) AND U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C                          = 2 IF THE SOLUTION IS SPECIFIED AT X=A AND
C                              THE BOUNDARY CONDITION IS MIXED AT X=B,
C                              I.E., U(A,Y) AND DU(B,Y)/DX+BETA*U(B,Y)
C                              ARE SPECIFIED FOR ALL Y
C                          = 3 IF THE BOUNDARY CONDITIONS AT X=A AND
C                              X=B ARE MIXED, I.E.,
C                              DU(A,Y)/DX+ALPHA*U(A,Y) AND
C                              DU(B,Y)/DX+BETA*U(B,Y) ARE SPECIFIED
C                              FOR ALL Y
C                          = 4 IF THE BOUNDARY CONDITION AT X=A IS
C                              MIXED AND THE SOLUTION IS SPECIFIED
C                              AT X=B, I.E., DU(A,Y)/DX+ALPHA*U(A,Y)
C                              AND U(B,Y) ARE SPECIFIED FOR ALL Y
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
C                          THAT SPECIFIES THE VALUES OF
C                          DU(A,Y)/DX+ ALPHA*U(A,Y) AT X=A, WHEN
C                          MBDCND=3 OR 4.
C                          BDA(J) = DU(A,YJ)/DX+ALPHA*U(A,YJ),
C                          J=1,2,...,N+1. WHEN MBDCND HAS ANY OTHER
C                          OTHER VALUE, BDA IS A DUMMY PARAMETER.
C
C                        ALPHA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT X=A
C                          (SEE ARGUMENT BDA).  IF MBDCND IS NOT
C                          EQUAL TO 3 OR 4 THEN ALPHA IS A DUMMY
C                          PARAMETER.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
C                          THAT SPECIFIES THE VALUES OF
C                          DU(B,Y)/DX+ BETA*U(B,Y) AT X=B.
C                          WHEN MBDCND=2 OR 3
C                          BDB(J) = DU(B,YJ)/DX+BETA*U(B,YJ),
C                          J=1,2,...,N+1. WHEN MBDCND HAS ANY OTHER
C                          OTHER VALUE, BDB IS A DUMMY PARAMETER.
C
C                        BETA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT
C                          X=B (SEE ARGUMENT BDB).  IF MBDCND IS
C                          NOT EQUAL TO 2 OR 3 THEN BETA IS A DUMMY
C                          PARAMETER.
C
C                        C,D
C                          THE RANGE OF THE Y-INDEPENDENT VARIABLE,
C                          I.E., Y IS GREATER THAN OR EQUAL TO C
C                          AND LESS THAN OR EQUAL TO D.  C MUST BE
C                          LESS THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL [C,D] IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS
C                          IN THE Y-DIRECTION GIVEN BY
C                          YJ=C+(J-1)*DLY FOR J=1,2,...,N+1 WHERE
C                          DLY=(D-C)/N IS THE PANEL WIDTH.
C                          IN ADDITION, N MUST BE GREATER THAN 4.
C
C                        NBDCND
C                          INDICATES THE TYPES OF BOUNDARY CONDITIONS
C                          AT Y=C AND Y=D
C
C                          = 0 IF THE SOLUTION IS PERIODIC IN Y,
C                              I.E., U(X,Y+D-C)=U(X,Y)  FOR ALL X,Y
C                          = 1 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND Y = D, I.E., U(X,C) AND U(X,D)
C                              ARE SPECIFIED FOR ALL X
C                          = 2 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND THE BOUNDARY CONDITION IS MIXED
C                              AT Y=D, I.E., U(X,C) AND
C                              DU(X,D)/DY+XNU*U(X,D) ARE SPECIFIED
C                              FOR ALL X
C                          = 3 IF THE BOUNDARY CONDITIONS ARE MIXED
C                              AT Y=C AND Y=D, I.E.,
C                              DU(X,D)/DY+GAMA*U(X,C) AND
C                              DU(X,D)/DY+XNU*U(X,D) ARE SPECIFIED
C                              FOR ALL X
C                          = 4 IF THE BOUNDARY CONDITION IS MIXED
C                              AT Y=C AND THE SOLUTION IS SPECIFIED
C                              AT Y=D, I.E. DU(X,C)/DY+GAMA*U(X,C)
C                              AND U(X,D) ARE SPECIFIED FOR ALL X
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
C                          THAT SPECIFIES THE VALUE OF
C                          DU(X,C)/DY+GAMA*U(X,C) AT Y=C.
C                          WHEN NBDCND=3 OR 4 BDC(I) = DU(XI,C)/DY +
C                          GAMA*U(XI,C), I=1,2,...,M+1.
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC
C                          IS A DUMMY PARAMETER.
C
C                        GAMA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT
C                          Y=C (SEE ARGUMENT BDC).  IF NBDCND IS
C                          NOT EQUAL TO 3 OR 4 THEN GAMA IS A DUMMY
C                          PARAMETER.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
C                          THAT SPECIFIES THE VALUE OF
C                          DU(X,D)/DY + XNU*U(X,D) AT Y=C.
C                          WHEN NBDCND=2 OR 3 BDD(I) = DU(XI,D)/DY +
C                          XNU*U(XI,D), I=1,2,...,M+1.
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD
C                          IS A DUMMY PARAMETER.
C
C                        XNU
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT
C                          Y=D (SEE ARGUMENT BDD).  IF NBDCND IS
C                          NOT EQUAL TO 2 OR 3 THEN XNU IS A
C                          DUMMY PARAMETER.
C
C                        COFX
C                          A USER-SUPPLIED SUBPROGRAM WITH
C                          PARAMETERS X, AFUN, BFUN, CFUN WHICH
C                          RETURNS THE VALUES OF THE X-DEPENDENT
C                          COEFFICIENTS AF(X), BF(X), CF(X) IN THE
C                          ELLIPTIC EQUATION AT X.
C
C                        COFY
C                          A USER-SUPPLIED SUBPROGRAM WITH PARAMETERS
C                          Y, DFUN, EFUN, FFUN WHICH RETURNS THE
C                          VALUES OF THE Y-DEPENDENT COEFFICIENTS
C                          DF(Y), EF(Y), FF(Y) IN THE ELLIPTIC
C                          EQUATION AT Y.
C
C                          NOTE:  COFX AND COFY MUST BE DECLARED
C                          EXTERNAL IN THE CALLING ROUTINE.
C                          THE VALUES RETURNED IN AFUN AND DFUN
C                          MUST SATISFY AFUN*DFUN GREATER THAN 0
C                          FOR A LESS THAN X LESS THAN B, C LESS
C                          THAN Y LESS THAN D (SEE IERROR=10).
C                          THE COEFFICIENTS PROVIDED MAY LEAD TO A
C                          MATRIX EQUATION WHICH IS NOT DIAGONALLY
C                          DOMINANT IN WHICH CASE SOLUTION MAY FAIL
C                          (SEE IERROR=4).
C
C                        GRHS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT-HAND SIDE OF THE
C                          ELLIPTIC EQUATION, I.E.,
C                          GRHS(I,J)=G(XI,YI), FOR I=2,...,M,
C                          J=2,...,N.  AT THE BOUNDARIES, GRHS IS
C                          DEFINED BY
C
C                          MBDCND   GRHS(1,J)   GRHS(M+1,J)
C                          ------   ---------   -----------
C                            0      G(A,YJ)     G(B,YJ)
C                            1         *           *
C                            2         *        G(B,YJ)  J=1,2,...,N+1
C                            3      G(A,YJ)     G(B,YJ)
C                            4      G(A,YJ)        *
C
C                          NBDCND   GRHS(I,1)   GRHS(I,N+1)
C                          ------   ---------   -----------
C                            0      G(XI,C)     G(XI,D)
C                            1         *           *
C                            2         *        G(XI,D)  I=1,2,...,M+1
C                            3      G(XI,C)     G(XI,D)
C                            4      G(XI,C)        *
C
C                          WHERE * MEANS THESE QUANTITIES ARE NOT USED.
C                          GRHS SHOULD BE DIMENSIONED IDMN BY AT LEAST
C                          N+1 IN THE CALLING ROUTINE.
C
C                        USOL
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE SOLUTION ALONG THE BOUNDARIES.
C                          AT THE BOUNDARIES, USOL IS DEFINED BY
C
C                          MBDCND   USOL(1,J)   USOL(M+1,J)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(A,YJ)     U(B,YJ)
C                            2      U(A,YJ)        *     J=1,2,...,N+1
C                            3         *           *
C                            4         *        U(B,YJ)
C
C                          NBDCND   USOL(I,1)   USOL(I,N+1)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(XI,C)     U(XI,D)
C                            2      U(XI,C)        *     I=1,2,...,M+1
C                            3         *           *
C                            4         *        U(XI,D)
C
C                          WHERE * MEANS THE QUANTITIES ARE NOT USED
C                          IN THE SOLUTION.
C
C                          IF IORDER=2, THE USER MAY EQUIVALENCE GRHS
C                          AND USOL TO SAVE SPACE.  NOTE THAT IN THIS
C                          CASE THE TABLES SPECIFYING THE BOUNDARIES
C                          OF THE GRHS AND USOL ARRAYS DETERMINE THE
C                          BOUNDARIES UNIQUELY EXCEPT AT THE CORNERS.
C                          IF THE TABLES CALL FOR BOTH G(X,Y) AND
C                          U(X,Y) AT A CORNER THEN THE SOLUTION MUST
C                          BE CHOSEN.  FOR EXAMPLE, IF MBDCND=2 AND
C                          NBDCND=4, THEN U(A,C), U(A,D), U(B,D) MUST
C                          BE CHOSEN AT THE CORNERS IN ADDITION
C                          TO G(B,C).
C
C                          IF IORDER=4, THEN THE TWO ARRAYS, USOL AND
C                          GRHS, MUST BE DISTINCT.
C
C                          USOL SHOULD BE DIMENSIONED IDMN BY AT LEAST
C                          N+1 IN THE CALLING ROUTINE.
C
C                        IDMN
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAYS
C                          GRHS AND USOL AS IT APPEARS IN THE PROGRAM
C                          CALLING SEPELI.  THIS PARAMETER IS USED
C                          TO SPECIFY THE VARIABLE DIMENSION OF GRHS
C                          AND USOL.  IDMN MUST BE AT LEAST 7 AND
C                          GREATER THAN OR EQUAL TO M+1.
C
C                        W
C                          A ONE-DIMENSIONAL ARRAY THAT MUST BE
C                          PROVIDED BY THE USER FOR WORK SPACE.
C                          LET K=INT(LOG2(N+1))+1 AND SET L=2**(K+1).
C                          THEN (K-2)*L+K+10*N+12*M+27 WILL SUFFICE
C                          AS A LENGTH OF W.  THE ACTUAL LENGTH OF W
C                          IN THE CALLING ROUTINE MUST BE SET IN W(1)
C                          (SEE IERROR=11).
C
C ON OUTPUT              USOL
C                          CONTAINS THE APPROXIMATE SOLUTION TO THE
C                          ELLIPTIC EQUATION.
C                          USOL(I,J) IS THE APPROXIMATION TO U(XI,YJ)
C                          FOR I=1,2...,M+1 AND J=1,2,...,N+1.
C                          THE APPROXIMATION HAS ERROR
C                          O(DLX**2+DLY**2) IF CALLED WITH IORDER=2
C                          AND O(DLX**4+DLY**4) IF CALLED WITH
C                          IORDER=4.
C
C                        W
C                          CONTAINS INTERMEDIATE VALUES THAT MUST NOT
C                          BE DESTROYED IF SEPELI IS CALLED AGAIN WITH
C                          INTL=1.  IN ADDITION W(1) CONTAINS THE
C                          EXACT MINIMAL LENGTH (IN FLOATING POINT)
C                          REQUIRED FOR THE WORK SPACE (SEE IERROR=11).
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS
C                          (I.E., ALPHA=BETA=0 IF MBDCND=3;
C                          GAMA=XNU=0 IF NBDCND=3) IS SPECIFIED
C                          AND IF THE COEFFICIENTS OF U(X,Y) IN THE
C                          SEPARABLE ELLIPTIC EQUATION ARE ZERO
C                          (I.E., CF(X)=0 FOR X GREATER THAN OR EQUAL
C                          TO A AND LESS THAN OR EQUAL TO B;
C                          FF(Y)=0 FOR Y GREATER THAN OR EQUAL TO C
C                          AND LESS THAN OR EQUAL TO D) THEN A
C                          SOLUTION MAY NOT EXIST.  PERTRB IS A
C                          CONSTANT CALCULATED AND SUBTRACTED FROM
C                          THE RIGHT-HAND SIDE OF THE MATRIX EQUATIONS
C                          GENERATED BY SEPELI WHICH INSURES THAT A
C                          SOLUTION EXISTS. SEPELI THEN COMPUTES THIS
C                          SOLUTION WHICH IS A WEIGHTED MINIMAL LEAST
C                          SQUARES SOLUTION TO THE ORIGINAL PROBLEM.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS OR FAILURE TO FIND A SOLUTION
C                          = 0 NO ERROR
C                          = 1 IF A GREATER THAN B OR C GREATER THAN D
C                          = 2 IF MBDCND LESS THAN 0 OR MBDCND GREATER
C                              THAN 4
C                          = 3 IF NBDCND LESS THAN 0 OR NBDCND GREATER
C                              THAN 4
C                          = 4 IF ATTEMPT TO FIND A SOLUTION FAILS.
C                              (THE LINEAR SYSTEM GENERATED IS NOT
C                              DIAGONALLY DOMINANT.)
C                          = 5 IF IDMN IS TOO SMALL
C                              (SEE DISCUSSION OF IDMN)
C                          = 6 IF M IS TOO SMALL OR TOO LARGE
C                              (SEE DISCUSSION OF M)
C                          = 7 IF N IS TOO SMALL (SEE DISCUSSION OF N)
C                          = 8 IF IORDER IS NOT 2 OR 4
C                          = 9 IF INTL IS NOT 0 OR 1
C                          = 10 IF AFUN*DFUN LESS THAN OR EQUAL TO 0
C                               FOR SOME INTERIOR MESH POINT (XI,YJ)
C                          = 11 IF THE WORK SPACE LENGTH INPUT IN W(1)
C                               IS LESS THAN THE EXACT MINIMAL WORK
C                               SPACE LENGTH REQUIRED OUTPUT IN W(1).
C
C                          NOTE (CONCERNING IERROR=4):  FOR THE
C                          COEFFICIENTS INPUT THROUGH COFX, COFY,
C                          THE DISCRETIZATION MAY LEAD TO A BLOCK
C                          TRIDIAGONAL LINEAR SYSTEM WHICH IS NOT
C                          DIAGONALLY DOMINANT (FOR EXAMPLE, THIS
C                          HAPPENS IF CFUN=0 AND BFUN/(2.*DLX) GREATER
C                          THAN AFUN/DLX**2).  IN THIS CASE SOLUTION
C                          MAY FAIL.  THIS CANNOT HAPPEN IN THE LIMIT
C                          AS DLX, DLY APPROACH ZERO.  HENCE, THE
C                          CONDITION MAY BE REMEDIED BY TAKING LARGER
C                          VALUES FOR M OR N.
C
C SPECIAL CONDITIONS     SEE COFX, COFY ARGUMENT DESCRIPTIONS ABOVE.
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       BLKTRI, COMF, Q8QST4, AND SEPAUX, WHICH ARE
C FILES                  LOADED BY DEFAULT ON NCAR'S CRAY MACHINES.
C
C LANGUAGE               FORTRAN
C
C HISTORY                DEVELOPED AT NCAR DURING 1975-76 BY
C                        JOHN C. ADAMS OF THE SCIENTIFIC COMPUTING
C                        DIVISION.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980.
C
C PORTABILITY            FORTRAN 66
C
C ALGORITHM              SEPELI AUTOMATICALLY DISCRETIZES THE
C                        SEPARABLE ELLIPTIC EQUATION WHICH IS THEN
C                        SOLVED BY A GENERALIZED CYCLIC REDUCTION
C                        ALGORITHM IN THE SUBROUTINE, BLKTRI.  THE
C                        FOURTH-ORDER SOLUTION IS OBTAINED USING
C                        'DEFERRED CORRECTIONS' WHICH IS DESCRIBED
C                        AND REFERENCED IN SECTIONS, REFERENCES AND
C                        METHOD.
C
C TIMING                 THE OPERATIONAL COUNT IS PROPORTIONAL TO
C                        M*N*LOG2(N).
C
C ACCURACY               THE FOLLOWING ACCURACY RESULTS WERE OBTAINED
C                        ON A CDC 7600.  NOTE THAT THE FOURTH-ORDER
C                        ACCURACY IS NOT REALIZED UNTIL THE MESH IS
C                        SUFFICIENTLY REFINED.
C
C                                     SECOND-ORDER  FOURTH-ORDER
C                            M    N     ERROR         ERROR
C
C                             6    6    6.8E-1        1.2E0
C                            14   14    1.4E-1        1.8E-1
C                            30   30    3.2E-2        9.7E-3
C                            62   62    7.5E-3        3.0E-4
C                           126  126    1.8E-3        3.5E-6
C
C PORTABILITY            FORTRAN 66
C
C REFERENCES             KELLER, H.B., NUMERICAL METHODS FOR TWO-POINT
C                        BOUNDARY-VALUE PROBLEMS, BLAISDEL (1968),
C                        WALTHAM, MASS.
C
C                        SWARZTRAUBER, P., AND R. SWEET (1975):
C                        EFFICIENT FORTRAN SUBPROGRAMS FOR THE
C                        SOLUTION OF ELLIPTIC PARTIAL DIFFERENTIAL
C                        EQUATIONS.  NCAR TECHNICAL NOTE
C                        NCAR-TN/IA-109, PP. 135-137.
C***********************************************************************
      DIMENSION       GRHS(IDMN,1)           ,USOL(IDMN,1)
      DIMENSION       BDA(1)     ,BDB(1)     ,BDC(1)     ,BDD(1)     ,
     1                W(1)
      EXTERNAL        COFX       ,COFY
      LOGICAL Q8Q4
      SAVE Q8Q4
      DATA Q8Q4 /.TRUE./
C
      IF (Q8Q4) THEN
C         CALL Q8QST4('LOCLIB','SEPELI','SEPELI','VERSION 01')
          Q8Q4 = .FALSE.
      ENDIF
C
C     CHECK INPUT PARAMETERS
C
      CALL CHKPRM (INTL,IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,COFY,
     1             IDMN,IERROR)
      IF (IERROR .NE. 0) RETURN
C
C     COMPUTE MINIMUM WORK SPACE AND CHECK WORK SPACE LENGTH INPUT
C
      L = N+1
      IF (NBDCND .EQ. 0) L = N
      LOGB2N = INT(ALOG(FLOAT(L)+0.5)/ALOG(2.0))+1
      LL = 2**(LOGB2N+1)
      K = M+1
      L = N+1
      LENGTH = (LOGB2N-2)*LL+LOGB2N+MAX0(2*L,6*K)+5
      IF (NBDCND .EQ. 0) LENGTH = LENGTH+2*L
      IERROR = 11
      LINPUT = INT(W(1)+0.5)
      LOUTPT = LENGTH+6*(K+L)+1
      W(1) = FLOAT(LOUTPT)
      IF (LOUTPT .GT. LINPUT) RETURN
      IERROR = 0
C
C     SET WORK SPACE INDICES
C
      I1 = LENGTH+2
      I2 = I1+L
      I3 = I2+L
      I4 = I3+L
      I5 = I4+L
      I6 = I5+L
      I7 = I6+L
      I8 = I7+K
      I9 = I8+K
      I10 = I9+K
      I11 = I10+K
      I12 = I11+K
      I13 = 2
      CALL SPELIP (INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     1             NBDCND,BDC,GAMA,BDD,XNU,COFX,COFY,W(I1),W(I2),W(I3),
     2             W(I4),W(I5),W(I6),W(I7),W(I8),W(I9),W(I10),W(I11),
     3             W(I12),GRHS,USOL,IDMN,W(I13),PERTRB,IERROR)
      RETURN
      END
      SUBROUTINE SPELIP (INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,
     1                   D,N,NBDCND,BDC,GAMA,BDD,XNU,COFX,COFY,AN,BN,
     2                   CN,DN,UN,ZN,AM,BM,CM,DM,UM,ZM,GRHS,USOL,IDMN,
     3                   W,PERTRB,IERROR)
C
C     SPELIP SETS UP VECTORS AND ARRAYS FOR INPUT TO BLKTRI
C     AND COMPUTES A SECOND ORDER SOLUTION IN USOL.  A RETURN JUMP TO
C     SEPELI OCCURRS IF IORDER=2.  IF IORDER=4 A FOURTH ORDER
C     SOLUTION IS GENERATED IN USOL.
C
      DIMENSION       BDA(1)     ,BDB(1)     ,BDC(1)     ,BDD(1)     ,
     1                W(1)
      DIMENSION       GRHS(IDMN,1)           ,USOL(IDMN,1)
      DIMENSION       AN(1)      ,BN(1)      ,CN(1)      ,DN(1)      ,
     1                UN(1)      ,ZN(1)
      DIMENSION       AM(1)      ,BM(1)      ,CM(1)      ,DM(1)      ,
     1                UM(1)      ,ZM(1)
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      LOGICAL         SINGLR
      EXTERNAL        COFX       ,COFY
C
C     SET PARAMETERS INTERNALLY
C
      KSWX = MBDCND+1
      KSWY = NBDCND+1
      K = M+1
      L = N+1
      AIT = A
      BIT = B
      CIT = C
      DIT = D
C
C     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
C     AND NON-SPECIFIED BOUNDARIES.
C
      DO  20 I=2,M
         DO  10 J=2,N
            USOL(I,J) = GRHS(I,J)
   10    CONTINUE
   20 CONTINUE
      IF (KSWX.EQ.2 .OR. KSWX.EQ.3) GO TO  40
      DO  30 J=2,N
         USOL(1,J) = GRHS(1,J)
   30 CONTINUE
   40 CONTINUE
      IF (KSWX.EQ.2 .OR. KSWX.EQ.5) GO TO  60
      DO  50 J=2,N
         USOL(K,J) = GRHS(K,J)
   50 CONTINUE
   60 CONTINUE
      IF (KSWY.EQ.2 .OR. KSWY.EQ.3) GO TO  80
      DO  70 I=2,M
         USOL(I,1) = GRHS(I,1)
   70 CONTINUE
   80 CONTINUE
      IF (KSWY.EQ.2 .OR. KSWY.EQ.5) GO TO 100
      DO  90 I=2,M
         USOL(I,L) = GRHS(I,L)
   90 CONTINUE
  100 CONTINUE
      IF (KSWX.NE.2 .AND. KSWX.NE.3 .AND. KSWY.NE.2 .AND. KSWY.NE.3)
     1    USOL(1,1) = GRHS(1,1)
      IF (KSWX.NE.2 .AND. KSWX.NE.5 .AND. KSWY.NE.2 .AND. KSWY.NE.3)
     1    USOL(K,1) = GRHS(K,1)
      IF (KSWX.NE.2 .AND. KSWX.NE.3 .AND. KSWY.NE.2 .AND. KSWY.NE.5)
     1    USOL(1,L) = GRHS(1,L)
      IF (KSWX.NE.2 .AND. KSWX.NE.5 .AND. KSWY.NE.2 .AND. KSWY.NE.5)
     1    USOL(K,L) = GRHS(K,L)
      I1 = 1
C
C     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
C
      MP = 1
      NP = 1
      IF (KSWX .EQ. 1) MP = 0
      IF (KSWY .EQ. 1) NP = 0
C
C     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
C     IN NINT,MINT
C
      DLX = (BIT-AIT)/FLOAT(M)
      MIT = K-1
      IF (KSWX .EQ. 2) MIT = K-2
      IF (KSWX .EQ. 4) MIT = K
      DLY = (DIT-CIT)/FLOAT(N)
      NIT = L-1
      IF (KSWY .EQ. 2) NIT = L-2
      IF (KSWY .EQ. 4) NIT = L
      TDLX3 = 2.0*DLX**3
      DLX4 = DLX**4
      TDLY3 = 2.0*DLY**3
      DLY4 = DLY**4
C
C     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI
C
      IS = 1
      JS = 1
      IF (KSWX.EQ.2 .OR. KSWX.EQ.3) IS = 2
      IF (KSWY.EQ.2 .OR. KSWY.EQ.3) JS = 2
      NS = NIT+JS-1
      MS = MIT+IS-1
C
C     SET X - DIRECTION
C
      DO 110 I=1,MIT
         XI = AIT+FLOAT(IS+I-2)*DLX
         CALL COFX (XI,AI,BI,CI)
         AXI = (AI/DLX-0.5*BI)/DLX
         BXI = -2.*AI/DLX**2+CI
         CXI = (AI/DLX+0.5*BI)/DLX
         AM(I) = AXI
         BM(I) = BXI
         CM(I) = CXI
  110 CONTINUE
C
C     SET Y DIRECTION
C
      DO 120 J=1,NIT
         YJ = CIT+FLOAT(JS+J-2)*DLY
         CALL COFY (YJ,DJ,EJ,FJ)
         DYJ = (DJ/DLY-0.5*EJ)/DLY
         EYJ = (-2.*DJ/DLY**2+FJ)
         FYJ = (DJ/DLY+0.5*EJ)/DLY
         AN(J) = DYJ
         BN(J) = EYJ
         CN(J) = FYJ
  120 CONTINUE
C
C     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
C
      AX1 = AM(1)
      CXM = CM(MIT)
      GO TO (170,130,150,160,140),KSWX
C
C     DIRICHLET-DIRICHLET IN X DIRECTION
C
  130 AM(1) = 0.0
      CM(MIT) = 0.0
      GO TO 170
C
C     MIXED-DIRICHLET IN X DIRECTION
C
  140 AM(1) = 0.0
      BM(1) = BM(1)+2.*ALPHA*DLX*AX1
      CM(1) = CM(1)+AX1
      CM(MIT) = 0.0
      GO TO 170
C
C     DIRICHLET-MIXED IN X DIRECTION
C
  150 AM(1) = 0.0
      AM(MIT) = AM(MIT)+CXM
      BM(MIT) = BM(MIT)-2.*BETA*DLX*CXM
      CM(MIT) = 0.0
      GO TO 170
C
C     MIXED - MIXED IN X DIRECTION
C
  160 CONTINUE
      AM(1) = 0.0
      BM(1) = BM(1)+2.*DLX*ALPHA*AX1
      CM(1) = CM(1)+AX1
      AM(MIT) = AM(MIT)+CXM
      BM(MIT) = BM(MIT)-2.*DLX*BETA*CXM
      CM(MIT) = 0.0
  170 CONTINUE
C
C     ADJUST IN Y DIRECTION UNLESS PERIODIC
C
      DY1 = AN(1)
      FYN = CN(NIT)
      GO TO (220,180,200,210,190),KSWY
C
C     DIRICHLET-DIRICHLET IN Y DIRECTION
C
  180 CONTINUE
      AN(1) = 0.0
      CN(NIT) = 0.0
      GO TO 220
C
C     MIXED-DIRICHLET IN Y DIRECTION
C
  190 CONTINUE
      AN(1) = 0.0
      BN(1) = BN(1)+2.*DLY*GAMA*DY1
      CN(1) = CN(1)+DY1
      CN(NIT) = 0.0
      GO TO 220
C
C     DIRICHLET-MIXED IN Y DIRECTION
C
  200 AN(1) = 0.0
      AN(NIT) = AN(NIT)+FYN
      BN(NIT) = BN(NIT)-2.*DLY*XNU*FYN
      CN(NIT) = 0.0
      GO TO 220
C
C     MIXED - MIXED DIRECTION IN Y DIRECTION
C
  210 CONTINUE
      AN(1) = 0.0
      BN(1) = BN(1)+2.*DLY*GAMA*DY1
      CN(1) = CN(1)+DY1
      AN(NIT) = AN(NIT)+FYN
      BN(NIT) = BN(NIT)-2.0*DLY*XNU*FYN
      CN(NIT) = 0.0
  220 IF (KSWX .EQ. 1) GO TO 270
C
C     ADJUST USOL ALONG X EDGE
C
      DO 260 J=JS,NS
         IF (KSWX.NE.2 .AND. KSWX.NE.3) GO TO 230
         USOL(IS,J) = USOL(IS,J)-AX1*USOL(1,J)
         GO TO 240
  230    USOL(IS,J) = USOL(IS,J)+2.0*DLX*AX1*BDA(J)
  240    IF (KSWX.NE.2 .AND. KSWX.NE.5) GO TO 250
         USOL(MS,J) = USOL(MS,J)-CXM*USOL(K,J)
         GO TO 260
  250    USOL(MS,J) = USOL(MS,J)-2.0*DLX*CXM*BDB(J)
  260 CONTINUE
  270 IF (KSWY .EQ. 1) GO TO 320
C
C     ADJUST USOL ALONG Y EDGE
C
      DO 310 I=IS,MS
         IF (KSWY.NE.2 .AND. KSWY.NE.3) GO TO 280
         USOL(I,JS) = USOL(I,JS)-DY1*USOL(I,1)
         GO TO 290
  280    USOL(I,JS) = USOL(I,JS)+2.0*DLY*DY1*BDC(I)
  290    IF (KSWY.NE.2 .AND. KSWY.NE.5) GO TO 300
         USOL(I,NS) = USOL(I,NS)-FYN*USOL(I,L)
         GO TO 310
  300    USOL(I,NS) = USOL(I,NS)-2.0*DLY*FYN*BDD(I)
  310 CONTINUE
  320 CONTINUE
C
C     SAVE ADJUSTED EDGES IN GRHS IF IORDER=4
C
      IF (IORDER .NE. 4) GO TO 350
      DO 330 J=JS,NS
         GRHS(IS,J) = USOL(IS,J)
         GRHS(MS,J) = USOL(MS,J)
  330 CONTINUE
      DO 340 I=IS,MS
         GRHS(I,JS) = USOL(I,JS)
         GRHS(I,NS) = USOL(I,NS)
  340 CONTINUE
  350 CONTINUE
      IORD = IORDER
      PERTRB = 0.0
C
C     CHECK IF OPERATOR IS SINGULAR
C
      CALL CHKSNG (MBDCND,NBDCND,ALPHA,BETA,GAMA,XNU,COFX,COFY,SINGLR)
C
C     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
C     IF SINGULAR
C
      IF (SINGLR) CALL SEPTRI (MIT,AM,BM,CM,DM,UM,ZM)
      IF (SINGLR) CALL SEPTRI (NIT,AN,BN,CN,DN,UN,ZN)
C
C     MAKE INITIALIZATION CALL TO BLKTRI
C
      IF (INTL .EQ. 0)
     1    CALL BLKTRI (INTL,NP,NIT,AN,BN,CN,MP,MIT,AM,BM,CM,IDMN,
     2                 USOL(IS,JS),IERROR,W)
      IF (IERROR .NE. 0) RETURN
C
C     ADJUST RIGHT HAND SIDE IF NECESSARY
C
  360 CONTINUE
      IF (SINGLR) CALL SEPORT (USOL,IDMN,ZN,ZM,PERTRB)
C
C     COMPUTE SOLUTION
C
      CALL BLKTRI (I1,NP,NIT,AN,BN,CN,MP,MIT,AM,BM,CM,IDMN,USOL(IS,JS),
     1             IERROR,W)
      IF (IERROR .NE. 0) RETURN
C
C     SET PERIODIC BOUNDARIES IF NECESSARY
C
      IF (KSWX .NE. 1) GO TO 380
      DO 370 J=1,L
         USOL(K,J) = USOL(1,J)
  370 CONTINUE
  380 IF (KSWY .NE. 1) GO TO 400
      DO 390 I=1,K
         USOL(I,L) = USOL(I,1)
  390 CONTINUE
  400 CONTINUE
C
C     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
C     NORM IF OPERATOR IS SINGULAR
C
      IF (SINGLR) CALL SEPMIN (USOL,IDMN,ZN,ZM,PRTRB)
C
C     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
C     NOT FLAGGED
C
      IF (IORD .EQ. 2) RETURN
      IORD = 2
C
C     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
C
      CALL DEFER (COFX,COFY,IDMN,USOL,GRHS)
      GO TO 360
      END
      SUBROUTINE CHKPRM (INTL,IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,
     1                   COFY,IDMN,IERROR)
C
C     THIS PROGRAM CHECKS THE INPUT PARAMETERS FOR ERRORS
C
      EXTERNAL        COFX       ,COFY
C
C     CHECK DEFINITION OF SOLUTION REGION
C
      IERROR = 1
      IF (A.GE.B .OR. C.GE.D) RETURN
C
C     CHECK BOUNDARY SWITCHES
C
      IERROR = 2
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) RETURN
      IERROR = 3
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) RETURN
C
C     CHECK FIRST DIMENSION IN CALLING ROUTINE
C
      IERROR = 5
      IF (IDMN .LT. 7) RETURN
C
C     CHECK M
C
      IERROR = 6
      IF (M.GT.(IDMN-1) .OR. M.LT.6) RETURN
C
C     CHECK N
C
      IERROR = 7
      IF (N .LT. 5) RETURN
C
C     CHECK IORDER
C
      IERROR = 8
      IF (IORDER.NE.2 .AND. IORDER.NE.4) RETURN
C
C     CHECK INTL
C
      IERROR = 9
      IF (INTL.NE.0 .AND. INTL.NE.1) RETURN
C
C     CHECK THAT EQUATION IS ELLIPTIC
C
      DLX = (B-A)/FLOAT(M)
      DLY = (D-C)/FLOAT(N)
      DO  30 I=2,M
         XI = A+FLOAT(I-1)*DLX
         CALL COFX (XI,AI,BI,CI)
         DO  20 J=2,N
            YJ = C+FLOAT(J-1)*DLY
            CALL COFY (YJ,DJ,EJ,FJ)
            IF (AI*DJ .GT. 0.0) GO TO  10
            IERROR = 10
            RETURN
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
C
C     NO ERROR FOUND
C
      IERROR = 0
      RETURN
      END
      SUBROUTINE CHKSNG (MBDCND,NBDCND,ALPHA,BETA,GAMA,XNU,COFX,COFY,
     1                   SINGLR)
C
C     THIS SUBROUTINE CHECKS IF THE PDE   SEPELI
C     MUST SOLVE IS A SINGULAR OPERATOR
C
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      LOGICAL         SINGLR
      SINGLR = .FALSE.
C
C     CHECK IF THE BOUNDARY CONDITIONS ARE
C     ENTIRELY PERIODIC AND/OR MIXED
C
      IF ((MBDCND.NE.0 .AND. MBDCND.NE.3) .OR.
     1    (NBDCND.NE.0 .AND. NBDCND.NE.3)) RETURN
C
C     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
C
      IF (MBDCND .NE. 3) GO TO  10
      IF (ALPHA.NE.0.0 .OR. BETA.NE.0.0) RETURN
   10 IF (NBDCND .NE. 3) GO TO  20
      IF (GAMA.NE.0.0 .OR. XNU.NE.0.0) RETURN
   20 CONTINUE
C
C     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
C     ARE ZERO
C
      DO  30 I=IS,MS
         XI = AIT+FLOAT(I-1)*DLX
         CALL COFX (XI,AI,BI,CI)
         IF (CI .NE. 0.0) RETURN
   30 CONTINUE
      DO  40 J=JS,NS
         YJ = CIT+FLOAT(J-1)*DLY
         CALL COFY (YJ,DJ,EJ,FJ)
         IF (FJ .NE. 0.0) RETURN
   40 CONTINUE
C
C     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
C
      SINGLR = .TRUE.
      RETURN
      END
      SUBROUTINE DEFER (COFX,COFY,IDMN,USOL,GRHS)
C
C     THIS SUBROUTINE FIRST APPROXIMATES THE TRUNCATION ERROR GIVEN BY
C     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY WHERE
C     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 ON THE INTERIOR AND
C     AT THE BOUNDARIES IF PERIODIC(HERE UXXX,UXXXX ARE THE THIRD
C     AND FOURTH PARTIAL DERIVATIVES OF U WITH RESPECT TO X).
C     TX IS OF THE FORM AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
C     AT X=A OR X=B IF THE BOUNDARY CONDITION THERE IS MIXED.
C     TX=0.0 ALONG SPECIFIED BOUNDARIES.  TY HAS SYMMETRIC FORM
C     IN Y WITH X,AFUN(X),BFUN(X) REPLACED BY Y,DFUN(Y),EFUN(Y).
C     THE SECOND ORDER SOLUTION IN USOL IS USED TO APPROXIMATE
C     (VIA SECOND ORDER FINITE DIFFERENCING) THE TRUNCATION ERROR
C     AND THE RESULT IS ADDED TO THE RIGHT HAND SIDE IN GRHS
C     AND THEN TRANSFERRED TO USOL TO BE USED AS A NEW RIGHT
C     HAND SIDE WHEN CALLING BLKTRI FOR A FOURTH ORDER SOLUTION.
C
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       GRHS(IDMN,1)           ,USOL(IDMN,1)
      EXTERNAL        COFX       ,COFY
C
C     COMPUTE TRUNCATION ERROR APPROXIMATION OVER THE ENTIRE MESH
C
      DO  40 J=JS,NS
         YJ = CIT+FLOAT(J-1)*DLY
         CALL COFY (YJ,DJ,EJ,FJ)
         DO  30 I=IS,MS
            XI = AIT+FLOAT(I-1)*DLX
            CALL COFX (XI,AI,BI,CI)
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
C
            CALL SEPDX (USOL,IDMN,I,J,UXXX,UXXXX)
            CALL SEPDY (USOL,IDMN,I,J,UYYY,UYYYY)
            TX = AI*UXXXX/12.0+BI*UXXX/6.0
            TY = DJ*UYYYY/12.0+EJ*UYYY/6.0
C
C     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
C
            IF (KSWX.EQ.1 .OR. (I.GT.1 .AND. I.LT.K)) GO TO  10
            TX = AI/3.0*(UXXXX/4.0+UXXX/DLX)
   10       IF (KSWY.EQ.1 .OR. (J.GT.1 .AND. J.LT.L)) GO TO  20
            TY = DJ/3.0*(UYYYY/4.0+UYYY/DLY)
   20       GRHS(I,J) = GRHS(I,J)+DLX**2*TX+DLY**2*TY
   30    CONTINUE
   40 CONTINUE
C
C     RESET THE RIGHT HAND SIDE IN USOL
C
      DO  60 I=IS,MS
         DO  50 J=JS,NS
            USOL(I,J) = GRHS(I,J)
   50    CONTINUE
   60 CONTINUE
      RETURN
C
C REVISION HISTORY---
C
C DECEMBER 1979    FIRST ADDED TO NSSL
C-----------------------------------------------------------------------
      END
      SUBROUTINE BLKTRI (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,
     1                   IERROR,W)
C
C DIMENSION OF           AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N),
C ARGUMENTS              W(SEE ARGUMENT LIST)
C
C LATEST REVISION        JANUARY 1985
C
C USAGE                  CALL BLKTRI (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,
C                                     CM,IDIMY,Y,IERROR,W)
C
C PURPOSE                BLKTRI SOLVES A SYSTEM OF LINEAR EQUATIONS
C                        OF THE FORM
C
C                        AN(J)*X(I,J-1) + AM(I)*X(I-1,J) +
C                        (BN(J)+BM(I))*X(I,J) + CN(J)*X(I,J+1) +
C                        CM(I)*X(I+1,J) = Y(I,J)
C
C                        FOR I = 1,2,...,M  AND  J = 1,2,...,N.
C
C                        I+1 AND I-1 ARE EVALUATED MODULO M AND
C                        J+1 AND J-1 MODULO N, I.E.,
C
C                        X(I,0) = X(I,N),  X(I,N+1) = X(I,1),
C                        X(0,J) = X(M,J),  X(M+1,J) = X(1,J).
C
C                        THESE EQUATIONS USUALLY RESULT FROM THE
C                        DISCRETIZATION OF SEPARABLE ELLIPTIC
C                        EQUATIONS.  BOUNDARY CONDITIONS MAY BE
C                        DIRICHLET, NEUMANN, OR PERIODIC.
C
C ARGUMENTS
C
C ON INPUT               IFLG
C
C                          = 0  INITIALIZATION ONLY.
C                               CERTAIN QUANTITIES THAT DEPEND ON NP,
C                               N, AN, BN, AND CN ARE COMPUTED AND
C                               STORED IN THE WORK ARRAY W.
C
C                          = 1  THE QUANTITIES THAT WERE COMPUTED
C                               IN THE INITIALIZATION ARE USED
C                               TO OBTAIN THE SOLUTION X(I,J).
C
C                               NOTE:
C                               A CALL WITH IFLG=0 TAKES
C                               APPROXIMATELY ONE HALF THE TIME
C                               AS A CALL WITH IFLG = 1.
C                               HOWEVER, THE INITIALIZATION DOES
C                               NOT HAVE TO BE REPEATED UNLESS NP,
C                               N, AN, BN, OR CN CHANGE.
C
C                        NP
C                          = 0  IF AN(1) AND CN(N) ARE NOT ZERO,
C                               WHICH CORRESPONDS TO PERIODIC
C                               BOUNARY CONDITIONS.
C
C                          = 1  IF AN(1) AND CN(N) ARE ZERO.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
C                          N MUST BE GREATER THAN 4.
C                          THE OPERATION COUNT IS PROPORTIONAL TO
C                          MNLOG2(N), HENCE N SHOULD BE SELECTED
C                          LESS THAN OR EQUAL TO M.
C
C                        AN,BN,CN
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH N
C                          THAT SPECIFY THE COEFFICIENTS IN THE
C                          LINEAR EQUATIONS GIVEN ABOVE.
C
C                        MP
C                          = 0  IF AM(1) AND CM(M) ARE NOT ZERO,
C                               WHICH CORRESPONDS TO PERIODIC
C                               BOUNDARY CONDITIONS.
C
C                          = 1  IF AM(1) = CM(M) = 0  .
C
C                        M
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
C                           M MUST BE GREATER THAN 4.
C
C                        AM,BM,CM
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
C                          SPECIFY THE COEFFICIENTS IN THE LINEAR
C                          EQUATIONS GIVEN ABOVE.
C
C                        IDIMY
C                          THE ROW (OR FIRST) DIMENSION OF THE
C                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS
C                          IN THE PROGRAM CALLING BLKTRI.
C                          THIS PARAMETER IS USED TO SPECIFY THE
C                          VARIABLE DIMENSION OF Y.
C                          IDIMY MUST BE AT LEAST M.
C
C                        Y
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE RIGHT SIDE OF THE LINEAR
C                          SYSTEM OF EQUATIONS GIVEN ABOVE.
C                          Y MUST BE DIMENSIONED AT LEAST M*N.
C
C                        W
C                          A ONE-DIMENSIONAL ARRAY THAT MUST BE
C                          PROVIDED BY THE USER FOR WORK SPACE.
C                          IF NP=1 DEFINE K=INT(LOG2(N))+1 AND
C                          SET L=2**(K+1) THEN W MUST HAVE DIMENSION
C                          (K-2)*L+K+5+MAX(2N,6M)
C
C                          IF NP=0 DEFINE K=INT(LOG2(N-1))+1 AND
C                          SET L=2**(K+1) THEN W MUST HAVE DIMENSION
C                          (K-2)*L+K+5+2N+MAX(2N,6M)
C
C                          **IMPORTANT**
C                          FOR PURPOSES OF CHECKING, THE REQUIRED
C                          DIMENSION OF W IS COMPUTED BY BLKTRI AND
C                          STORED IN W(1) IN FLOATING POINT FORMAT.
C
C ARGUMENTS
C
C ON OUTPUT              Y
C                          CONTAINS THE SOLUTION X.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID
C                          INPUT PARAMETERS.  EXCEPT FOR NUMBER ZER0,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                        = 0  NO ERROR.
C                        = 1  M IS LESS THAN 5
C                        = 2  N IS LESS THAN 5
C                        = 3  IDIMY IS LESS THAN M.
C                        = 4  BLKTRI FAILED WHILE COMPUTING RESULTS
C                             THAT DEPEND ON THE COEFFICIENT ARRAYS
C                             AN, BN, CN.  CHECK THESE ARRAYS.
C                        = 5  AN(J)*CN(J-1) IS LESS THAN 0 FOR SOME J.
C
C                             POSSIBLE REASONS FOR THIS CONDITION ARE
C                             1. THE ARRAYS AN AND CN ARE NOT CORRECT
C                             2. TOO LARGE A GRID SPACING WAS USED
C                                IN THE DISCRETIZATION OF THE ELLIPTIC
C                                EQUATION.
C                             3. THE LINEAR EQUATIONS RESULTED FROM A
C                                PARTIAL DIFFERENTIAL EQUATION WHICH
C                                WAS NOT ELLIPTIC.
C
C                        W
C                           CONTAINS INTERMEDIATE VALUES THAT MUST
C                           NOT BE DESTROYED IF BLKTRI WILL BE CALLED
C                           AGAIN WITH IFLG=1. W(1) CONTAINS THE
C                           NUMBER OF LOCATIONS REQUIRED BY W IN
C                           FLOATING POINT FORMAT.
C
C
C SPECIAL CONDITIONS     THE ALGORITHM MAY FAIL IF ABS(BM(I)+BN(J))
C                        IS LESS THAN ABS(AM(I))+ABS(AN(J))+
C                        ABS(CM(I))+ABS(CN(J))
C                        FOR SOME I AND J. THE ALGORITHM WILL ALSO
C                        FAIL IF AN(J)*CN(J-1) IS LESS THAN ZERO FOR
C                        SOME J.
C                        SEE THE DESCRIPTION OF THE OUTPUT PARAMETER
C                        IERROR.
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       COMF AND Q8QST4, WHICH ARE AUTOMATICALLY LOADED
C FILES                  NCAR'S CRAY MACHINES.
C
C LANGUAGE               FORTRAN
C
C HISTORY                WRITTEN BY PAUL SWARZTRAUBER AT NCAR IN THE
C                        EARLY 1970'S.  REWRITTEN AND RELEASED IN
C                        JANUARY, 1980.
C
C ALGORITHM              GENERALIZED CYCLIC REDUCTION
C
C PORTABILITY            FORTRAN 66.  APPROXIMATE MACHINE ACCURACY
C                        IS COMPUTED IN FUNCTION EPMACH.
C
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, 'EFFICIENT
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
C                        ELLIPTIC EQUATIONS'
C                        NCAR TN/IA-109, JULY, 1975, 138 PP.
C
C                        SWARZTRAUBER P. N.,A DIRECT METHOD FOR
C                        THE DISCRETE SOLUTION OF SEPARABLE
C                        ELLIPTIC EQUATIONS, S.I.A.M.
C                        J. NUMER. ANAL.,11(1974) PP. 1136-1150.
C***********************************************************************
      DIMENSION       AN(1)      ,BN(1)      ,CN(1)      ,AM(1)      ,
     1                BM(1)      ,CM(1)      ,Y(IDIMY,1) ,W(2)
      EXTERNAL        PROD       ,PRODP      ,CPROD      ,CPRODP
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
C
C THE FOLLOWING CALL IS FOR MONITORING LIBRARY USE AT NCAR
C
      LOGICAL Q8Q4
      SAVE Q8Q4
      DATA Q8Q4 /.TRUE./
      IF (Q8Q4) THEN
C         CALL Q8QST4('LOCLIB','BLKTRI','BLKTRI','VERSION 01')
          Q8Q4 = .FALSE.
      ENDIF
C
C TEST M AND N FOR THE PROPER FORM
C
      NM = N
      IERROR = 0
      IF (M-5) 101,102,102
  101 IERROR = 1
      GO TO 119
  102 IF (NM-3) 103,104,104
  103 IERROR = 2
      GO TO 119
  104 IF (IDIMY-M) 105,106,106
  105 IERROR = 3
      GO TO 119
  106 NH = N
      NPP = NP
      IF (NPP) 107,108,107
  107 NH = NH+1
  108 IK = 2
      K = 1
  109 IK = IK+IK
      K = K+1
      IF (NH-IK) 110,110,109
  110 NL = IK
      IK = IK+IK
      NL = NL-1
      IWAH = (K-2)*IK+K+6
      IF (NPP) 111,112,111
C
C     DIVIDE W INTO WORKING SUB ARRAYS
C
  111 IW1 = IWAH
      IWBH = IW1+NM
      W(1) = FLOAT(IW1-1+MAX0(2*NM,6*M))
      GO TO 113
  112 IWBH = IWAH+NM+NM
      IW1 = IWBH
      W(1) = FLOAT(IW1-1+MAX0(2*NM,6*M))
      NM = NM-1
C
C SUBROUTINE COMP B COMPUTES THE ROOTS OF THE B POLYNOMIALS
C
  113 IF (IERROR) 119,114,119
  114 IW2 = IW1+M
      IW3 = IW2+M
      IWD = IW3+M
      IWW = IWD+M
      IWU = IWW+M
      IF (IFLG) 116,115,116
  115 CALL COMPB (NL,IERROR,AN,BN,CN,W(2),W(IWAH),W(IWBH))
      GO TO 119
  116 IF (MP) 117,118,117
C
C SUBROUTINE BLKTR1 SOLVES THE LINEAR SYSTEM
C
  117 CALL BLKTR1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),W(IW1),W(IW2),
     1             W(IW3),W(IWD),W(IWW),W(IWU),PROD,CPROD)
      GO TO 119
  118 CALL BLKTR1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),W(IW1),W(IW2),
     1             W(IW3),W(IWD),W(IWW),W(IWU),PRODP,CPRODP)
  119 CONTINUE
      RETURN
      END
      SUBROUTINE BLKTR1 (N,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,B,W1,W2,W3,WD,
     1                   WW,WU,PRDCT,CPRDCT)
C
C BLKTR1 SOLVES THE LINEAR SYSTEM
C
C B  CONTAINS THE ROOTS OF ALL THE B POLYNOMIALS
C W1,W2,W3,WD,WW,WU  ARE ALL WORKING ARRAYS
C PRDCT  IS EITHER PRODP OR PROD DEPENDING ON WHETHER THE BOUNDARY
C CONDITIONS IN THE M DIRECTION ARE PERIODIC OR NOT
C CPRDCT IS EITHER CPRODP OR CPROD WHICH ARE THE COMPLEX VERSIONS
C OF PRODP AND PROD. THESE ARE CALLED IN THE EVENT THAT SOME
C OF THE ROOTS OF THE B SUB P POLYNOMIAL ARE COMPLEX
C
C
      DIMENSION       AN(1)      ,BN(1)      ,CN(1)      ,AM(1)      ,
     1                BM(1)      ,CM(1)      ,B(1)       ,W1(1)      ,
     2                W2(1)      ,W3(1)      ,WD(1)      ,WW(1)      ,
     3                WU(1)      ,Y(IDIMY,1)
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
C
C BEGIN REDUCTION PHASE
C
      KDO = K-1
      DO 109 L=1,KDO
         IR = L-1
         I2 = 2**IR
         I1 = I2/2
         I3 = I2+I1
         I4 = I2+I2
         IRM1 = IR-1
         CALL INDXB (I2,IR,IM2,NM2)
         CALL INDXB (I1,IRM1,IM3,NM3)
         CALL INDXB (I3,IRM1,IM1,NM1)
         CALL PRDCT (NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,Y(1,I2),W3,
     1               M,AM,BM,CM,WD,WW,WU)
         IF = 2**K
         DO 108 I=I4,IF,I4
            IF (I-NM) 101,101,108
  101       IPI1 = I+I1
            IPI2 = I+I2
            IPI3 = I+I3
            CALL INDXC (I,IR,IDXC,NC)
            IF (I-IF) 102,108,108
  102       CALL INDXA (I,IR,IDXA,NA)
            CALL INDXB (I-I1,IRM1,IM1,NM1)
            CALL INDXB (IPI2,IR,IP2,NP2)
            CALL INDXB (IPI1,IRM1,IP1,NP1)
            CALL INDXB (IPI3,IRM1,IP3,NP3)
            CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W3,W1,M,AM,
     1                  BM,CM,WD,WW,WU)
            IF (IPI2-NM) 105,105,103
  103       DO 104 J=1,M
               W3(J) = 0.
               W2(J) = 0.
  104       CONTINUE
            GO TO 106
  105       CALL PRDCT (NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,
     1                  Y(1,IPI2),W3,M,AM,BM,CM,WD,WW,WU)
            CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W3,W2,M,AM,
     1                  BM,CM,WD,WW,WU)
  106       DO 107 J=1,M
               Y(J,I) = W1(J)+W2(J)+Y(J,I)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      IF (NPP) 132,110,132
C
C     THE PERIODIC CASE IS TREATED USING THE CAPACITANCE MATRIX METHOD
C
  110 IF = 2**K
      I = IF/2
      I1 = I/2
      CALL INDXB (I-I1,K-2,IM1,NM1)
      CALL INDXB (I+I1,K-2,IP1,NP1)
      CALL INDXB (I,K-1,IZ,NZ)
      CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,Y(1,I),W1,M,AM,
     1            BM,CM,WD,WW,WU)
      IZR = I
      DO 111 J=1,M
         W2(J) = W1(J)
  111 CONTINUE
      DO 113 LL=2,K
         L = K-LL+1
         IR = L-1
         I2 = 2**IR
         I1 = I2/2
         I = I2
         CALL INDXC (I,IR,IDXC,NC)
         CALL INDXB (I,IR,IZ,NZ)
         CALL INDXB (I-I1,IR-1,IM1,NM1)
         CALL INDXB (I+I1,IR-1,IP1,NP1)
         CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W1,W1,M,AM,BM,
     1               CM,WD,WW,WU)
         DO 112 J=1,M
            W1(J) = Y(J,I)+W1(J)
  112    CONTINUE
         CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,W1,M,AM,
     1               BM,CM,WD,WW,WU)
  113 CONTINUE
      DO 118 LL=2,K
         L = K-LL+1
         IR = L-1
         I2 = 2**IR
         I1 = I2/2
         I4 = I2+I2
         IFD = IF-I2
         DO 117 I=I2,IFD,I4
            IF (I-I2-IZR) 117,114,117
  114       IF (I-NM) 115,115,118
  115       CALL INDXA (I,IR,IDXA,NA)
            CALL INDXB (I,IR,IZ,NZ)
            CALL INDXB (I-I1,IR-1,IM1,NM1)
            CALL INDXB (I+I1,IR-1,IP1,NP1)
            CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W2,W2,M,AM,
     1                  BM,CM,WD,WW,WU)
            DO 116 J=1,M
               W2(J) = Y(J,I)+W2(J)
  116       CONTINUE
            CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W2,W2,M,
     1                  AM,BM,CM,WD,WW,WU)
            IZR = I
            IF (I-NM) 117,119,117
  117    CONTINUE
  118 CONTINUE
  119 DO 120 J=1,M
         Y(J,NM+1) = Y(J,NM+1)-CN(NM+1)*W1(J)-AN(NM+1)*W2(J)
  120 CONTINUE
      CALL INDXB (IF/2,K-1,IM1,NM1)
      CALL INDXB (IF,K-1,IP,NP)
      IF (NCMPLX) 121,122,121
  121 CALL CPRDCT (NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1),
     1             Y(1,NM+1),M,AM,BM,CM,W1,W3,WW)
      GO TO 123
  122 CALL PRDCT (NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1),
     1            Y(1,NM+1),M,AM,BM,CM,WD,WW,WU)
  123 DO 124 J=1,M
         W1(J) = AN(1)*Y(J,NM+1)
         W2(J) = CN(NM)*Y(J,NM+1)
         Y(J,1) = Y(J,1)-W1(J)
         Y(J,NM) = Y(J,NM)-W2(J)
  124 CONTINUE
      DO 126 L=1,KDO
         IR = L-1
         I2 = 2**IR
         I4 = I2+I2
         I1 = I2/2
         I = I4
         CALL INDXA (I,IR,IDXA,NA)
         CALL INDXB (I-I2,IR,IM2,NM2)
         CALL INDXB (I-I2-I1,IR-1,IM3,NM3)
         CALL INDXB (I-I1,IR-1,IM1,NM1)
         CALL PRDCT (NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,W1,W1,M,AM,
     1               BM,CM,WD,WW,WU)
         CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W1,W1,M,AM,BM,
     1               CM,WD,WW,WU)
         DO 125 J=1,M
            Y(J,I) = Y(J,I)-W1(J)
  125    CONTINUE
  126 CONTINUE
C
      IZR = NM
      DO 131 L=1,KDO
         IR = L-1
         I2 = 2**IR
         I1 = I2/2
         I3 = I2+I1
         I4 = I2+I2
         IRM1 = IR-1
         DO 130 I=I4,IF,I4
            IPI1 = I+I1
            IPI2 = I+I2
            IPI3 = I+I3
            IF (IPI2-IZR) 127,128,127
  127       IF (I-IZR) 130,131,130
  128       CALL INDXC (I,IR,IDXC,NC)
            CALL INDXB (IPI2,IR,IP2,NP2)
            CALL INDXB (IPI1,IRM1,IP1,NP1)
            CALL INDXB (IPI3,IRM1,IP3,NP3)
            CALL PRDCT (NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,W2,W2,M,
     1                  AM,BM,CM,WD,WW,WU)
            CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W2,W2,M,AM,
     1                  BM,CM,WD,WW,WU)
            DO 129 J=1,M
               Y(J,I) = Y(J,I)-W2(J)
  129       CONTINUE
            IZR = I
            GO TO 131
  130    CONTINUE
  131 CONTINUE
C
C BEGIN BACK SUBSTITUTION PHASE
C
  132 DO 144 LL=1,K
         L = K-LL+1
         IR = L-1
         IRM1 = IR-1
         I2 = 2**IR
         I1 = I2/2
         I4 = I2+I2
         IFD = IF-I2
         DO 143 I=I2,IFD,I4
            IF (I-NM) 133,133,143
  133       IMI1 = I-I1
            IMI2 = I-I2
            IPI1 = I+I1
            IPI2 = I+I2
            CALL INDXA (I,IR,IDXA,NA)
            CALL INDXC (I,IR,IDXC,NC)
            CALL INDXB (I,IR,IZ,NZ)
            CALL INDXB (IMI1,IRM1,IM1,NM1)
            CALL INDXB (IPI1,IRM1,IP1,NP1)
            IF (I-I2) 134,134,136
  134       DO 135 J=1,M
               W1(J) = 0.
  135       CONTINUE
            GO TO 137
  136       CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),Y(1,IMI2),
     1                  W1,M,AM,BM,CM,WD,WW,WU)
  137       IF (IPI2-NM) 140,140,138
  138       DO 139 J=1,M
               W2(J) = 0.
  139       CONTINUE
            GO TO 141
  140       CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),Y(1,IPI2),
     1                  W2,M,AM,BM,CM,WD,WW,WU)
  141       DO 142 J=1,M
               W1(J) = Y(J,I)+W1(J)+W2(J)
  142       CONTINUE
            CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,Y(1,I),
     1                  M,AM,BM,CM,WD,WW,WU)
  143    CONTINUE
  144 CONTINUE
      RETURN
      END
      SUBROUTINE INDXB (I,IR,IDX,IDP)
C
C B(IDX) IS THE LOCATION OF THE FIRST ROOT OF THE B(I,IR) POLYNOMIAL
C
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
      IDP = 0
      IF (IR) 107,101,103
  101 IF (I-NM) 102,102,107
  102 IDX = I
      IDP = 1
      RETURN
  103 IZH = 2**IR
      ID = I-IZH-IZH
      IDX = ID+ID+(IR-1)*IK+IR+(IK-I)/IZH+4
      IPL = IZH-1
      IDP = IZH+IZH-1
      IF (I-IPL-NM) 105,105,104
  104 IDP = 0
      RETURN
  105 IF (I+IPL-NM) 107,107,106
  106 IDP = NM+IPL-I+1
  107 RETURN
      END
      SUBROUTINE INDXA (I,IR,IDXA,NA)
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
      NA = 2**IR
      IDXA = I-NA+1
      IF (I-NM) 102,102,101
  101 NA = 0
  102 RETURN
      END
      SUBROUTINE INDXC (I,IR,IDXC,NC)
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
      NC = 2**IR
      IDXC = I
      IF (IDXC+NC-1-NM) 102,102,101
  101 NC = 0
  102 RETURN
      END
      SUBROUTINE PROD (ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,W,U)
C
C PROD APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,W,U ARE WORKING ARRAYS
C IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
c     DIMENSION       A(1)       ,B(1)       ,C(1)       ,X(1)       ,
c    1                Y(1)       ,D(1)       ,W(1)       ,BD(1)      ,
c    2                BM1(1)     ,BM2(1)     ,AA(1)      ,U(1)
      DIMENSION       A(M)       ,B(M)       ,C(M)       ,X(M)       ,
     1                Y(M)       ,D(M)       ,W(M)       ,BD(M)      ,
     2                BM1(M)     ,BM2(M)     ,AA(M)      ,U(M)
      DO 101 J=1,M
         W(J) = X(J)
         Y(J) = W(J)
  101 CONTINUE
      MM = M-1
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 IF (IA) 105,105,103
  103 RT = AA(IA)
      IF (ND .EQ. 0) RT = -RT
      IA = IA-1
C
C SCALAR MULTIPLICATION
C
      DO 104 J=1,M
         Y(J) = RT*W(J)
  104 CONTINUE
  105 IF (ID) 125,125,106
  106 RT = BD(ID)
      ID = ID-1
      IF (ID .EQ. 0) IBR = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      D(M) = A(M)/(B(M)-RT)
      W(M) = Y(M)/(B(M)-RT)
      DO 107 J=2,MM
         K = M-J
         DEN = B(K+1)-RT-C(K+1)*D(K+2)
         D(K+1) = A(K+1)/DEN
         W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
  107 CONTINUE
      DEN = B(1)-RT-C(1)*D(2)
      W(1) = 1.
      IF (DEN) 108,109,108
  108 W(1) = (Y(1)-C(1)*W(2))/DEN
  109 DO 110 J=2,M
         W(J) = W(J)-D(J)*W(J-1)
  110 CONTINUE
      IF (NA) 113,113,102
  111 DO 112 J=1,M
         Y(J) = W(J)
  112 CONTINUE
      IBR = 1
      GO TO 102
  113 IF (M1) 114,114,115
  114 IF (M2) 111,111,120
  115 IF (M2) 117,117,116
  116 IF (ABS(BM1(M1))-ABS(BM2(M2))) 120,120,117
  117 IF (IBR) 118,118,119
  118 IF (ABS(BM1(M1)-BD(ID))-ABS(BM1(M1)-RT)) 111,119,119
  119 RT = RT-BM1(M1)
      M1 = M1-1
      GO TO 123
  120 IF (IBR) 121,121,122
  121 IF (ABS(BM2(M2)-BD(ID))-ABS(BM2(M2)-RT)) 111,122,122
  122 RT = RT-BM2(M2)
      M2 = M2-1
  123 DO 124 J=1,M
         Y(J) = Y(J)+RT*W(J)
  124 CONTINUE
      GO TO 102
  125 RETURN
      END
      SUBROUTINE PRODP (ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,U,W)
C
C PRODP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y        PERIODIC BOUNDARY CONDITIONS
C
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,U,W ARE WORKING ARRAYS
C IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      DIMENSION       A(1)       ,B(1)       ,C(1)       ,X(1)       ,
     1                Y(1)       ,D(1)       ,U(1)       ,BD(1)      ,
     2                BM1(1)     ,BM2(1)     ,AA(1)      ,W(1)
      DO 101 J=1,M
         Y(J) = X(J)
         W(J) = Y(J)
  101 CONTINUE
      MM = M-1
      MM2 = M-2
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 IF (IA) 105,105,103
  103 RT = AA(IA)
      IF (ND .EQ. 0) RT = -RT
      IA = IA-1
      DO 104 J=1,M
         Y(J) = RT*W(J)
  104 CONTINUE
  105 IF (ID) 128,128,106
  106 RT = BD(ID)
      ID = ID-1
      IF (ID .EQ. 0) IBR = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      BH = B(M)-RT
      YM = Y(M)
      DEN = B(1)-RT
      D(1) = C(1)/DEN
      U(1) = A(1)/DEN
      W(1) = Y(1)/DEN
      V = C(M)
      IF (MM2-2) 109,107,107
  107 DO 108 J=2,MM2
         DEN = B(J)-RT-A(J)*D(J-1)
         D(J) = C(J)/DEN
         U(J) = -A(J)*U(J-1)/DEN
         W(J) = (Y(J)-A(J)*W(J-1))/DEN
         BH = BH-V*U(J-1)
         YM = YM-V*W(J-1)
         V = -V*D(J-1)
  108 CONTINUE
  109 DEN = B(M-1)-RT-A(M-1)*D(M-2)
      D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
      W(M-1) = (Y(M-1)-A(M-1)*W(M-2))/DEN
      AM = A(M)-V*D(M-2)
      BH = BH-V*U(M-2)
      YM = YM-V*W(M-2)
      DEN = BH-AM*D(M-1)
      IF (DEN) 110,111,110
  110 W(M) = (YM-AM*W(M-1))/DEN
      GO TO 112
  111 W(M) = 1.
  112 W(M-1) = W(M-1)-D(M-1)*W(M)
      DO 113 J=2,MM
         K = M-J
         W(K) = W(K)-D(K)*W(K+1)-U(K)*W(M)
  113 CONTINUE
      IF (NA) 116,116,102
  114 DO 115 J=1,M
         Y(J) = W(J)
  115 CONTINUE
      IBR = 1
      GO TO 102
  116 IF (M1) 117,117,118
  117 IF (M2) 114,114,123
  118 IF (M2) 120,120,119
  119 IF (ABS(BM1(M1))-ABS(BM2(M2))) 123,123,120
  120 IF (IBR) 121,121,122
  121 IF (ABS(BM1(M1)-BD(ID))-ABS(BM1(M1)-RT)) 114,122,122
  122 RT = RT-BM1(M1)
      M1 = M1-1
      GO TO 126
  123 IF (IBR) 124,124,125
  124 IF (ABS(BM2(M2)-BD(ID))-ABS(BM2(M2)-RT)) 114,125,125
  125 RT = RT-BM2(M2)
      M2 = M2-1
  126 DO 127 J=1,M
         Y(J) = Y(J)+RT*W(J)
  127 CONTINUE
      GO TO 102
  128 RETURN
      END
      SUBROUTINE CPROD (ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,YY,M,A,B,C,D,W,Y)
C
C PROD APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN YY           (COMPLEX CASE)
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C NA IS THE LENGTH OF THE ARRAY AA
C X,YY THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS YY
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,W,Y ARE WORKING ARRAYS
C ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      COMPLEX         Y          ,D          ,W          ,BD         ,
     1                CRT        ,DEN        ,Y1         ,Y2
c     DIMENSION       A(1)       ,B(1)       ,C(1)       ,X(1)       ,
c    1                Y(1)       ,D(1)       ,W(1)       ,BD(1)      ,
c    2                BM1(1)     ,BM2(1)     ,AA(1)      ,YY(1)
      DIMENSION       A(M)       ,B(M)       ,C(M)       ,X(M)       ,
     1                Y(M)       ,D(M)       ,W(M)       ,BD(M)      ,
     2                BM1(M)     ,BM2(M)     ,AA(M)      ,YY(M)
      DO 101 J=1,M
         Y(J) = CMPLX(X(J),0.)
  101 CONTINUE
      MM = M-1
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
  102 IFLG = 0
      IF (ID) 109,109,103
  103 CRT = BD(ID)
      ID = ID-1
C
C BEGIN SOLUTION TO SYSTEM
C
      D(M) = A(M)/(B(M)-CRT)
      W(M) = Y(M)/(B(M)-CRT)
      DO 104 J=2,MM
         K = M-J
         DEN = B(K+1)-CRT-C(K+1)*D(K+2)
         D(K+1) = A(K+1)/DEN
         W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
  104 CONTINUE
      DEN = B(1)-CRT-C(1)*D(2)
      IF (CABS(DEN)) 105,106,105
  105 Y(1) = (Y(1)-C(1)*W(2))/DEN
      GO TO 107
  106 Y(1) = (1.,0.)
  107 DO 108 J=2,M
         Y(J) = W(J)-D(J)*Y(J-1)
  108 CONTINUE
  109 IF (M1) 110,110,112
  110 IF (M2) 121,121,111
  111 RT = BM2(M2)
      M2 = M2-1
      GO TO 117
  112 IF (M2) 113,113,114
  113 RT = BM1(M1)
      M1 = M1-1
      GO TO 117
  114 IF (ABS(BM1(M1))-ABS(BM2(M2))) 116,116,115
  115 RT = BM1(M1)
      M1 = M1-1
      GO TO 117
  116 RT = BM2(M2)
      M2 = M2-1
  117 Y1 = (B(1)-RT)*Y(1)+C(1)*Y(2)
      IF (MM-2) 120,118,118
C
C MATRIX MULTIPLICATION
C
  118 DO 119 J=2,MM
         Y2 = A(J)*Y(J-1)+(B(J)-RT)*Y(J)+C(J)*Y(J+1)
         Y(J-1) = Y1
         Y1 = Y2
  119 CONTINUE
  120 Y(M) = A(M)*Y(M-1)+(B(M)-RT)*Y(M)
      Y(M-1) = Y1
      IFLG = 1
      GO TO 102
  121 IF (IA) 124,124,122
  122 RT = AA(IA)
      IA = IA-1
      IFLG = 1
C
C SCALAR MULTIPLICATION
C
      DO 123 J=1,M
         Y(J) = RT*Y(J)
  123 CONTINUE
  124 IF (IFLG) 125,125,102
  125 DO 126 J=1,M
         YY(J) = REAL(Y(J))
  126 CONTINUE
      RETURN
      END
      SUBROUTINE CPRODP (ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,YY,M,A,B,C,D,U,Y)
C
C PRODP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN YY       PERIODIC BOUNDARY CONDITIONS
C AND  COMPLEX  CASE
C
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,YY THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS YY
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,U,Y ARE WORKING ARRAYS
C ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      COMPLEX         Y          ,D          ,U          ,V          ,
     1                DEN        ,BH         ,YM         ,AM         ,
     2                Y1         ,Y2         ,YH         ,BD         ,
     3                CRT
c     DIMENSION       A(1)       ,B(1)       ,C(1)       ,X(1)       ,
c    1                Y(1)       ,D(1)       ,U(1)       ,BD(1)      ,
c    2                BM1(1)     ,BM2(1)     ,AA(1)      ,YY(1)
      DIMENSION       A(M)       ,B(M)       ,C(M)       ,X(M)       ,
     1                Y(M)       ,D(M)       ,U(M)       ,BD(M)      ,
     2                BM1(M)     ,BM2(M)     ,AA(M)      ,YY(M)
      DO 101 J=1,M
         Y(J) = CMPLX(X(J),0.)
  101 CONTINUE
      MM = M-1
      MM2 = M-2
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
  102 IFLG = 0
      IF (ID) 111,111,103
  103 CRT = BD(ID)
      ID = ID-1
      IFLG = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      BH = B(M)-CRT
      YM = Y(M)
      DEN = B(1)-CRT
      D(1) = C(1)/DEN
      U(1) = A(1)/DEN
      Y(1) = Y(1)/DEN
      V = CMPLX(C(M),0.)
      IF (MM2-2) 106,104,104
  104 DO 105 J=2,MM2
         DEN = B(J)-CRT-A(J)*D(J-1)
         D(J) = C(J)/DEN
         U(J) = -A(J)*U(J-1)/DEN
         Y(J) = (Y(J)-A(J)*Y(J-1))/DEN
         BH = BH-V*U(J-1)
         YM = YM-V*Y(J-1)
         V = -V*D(J-1)
  105 CONTINUE
  106 DEN = B(M-1)-CRT-A(M-1)*D(M-2)
      D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
      Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2))/DEN
      AM = A(M)-V*D(M-2)
      BH = BH-V*U(M-2)
      YM = YM-V*Y(M-2)
      DEN = BH-AM*D(M-1)
      IF (CABS(DEN)) 107,108,107
  107 Y(M) = (YM-AM*Y(M-1))/DEN
      GO TO 109
  108 Y(M) = (1.,0.)
  109 Y(M-1) = Y(M-1)-D(M-1)*Y(M)
      DO 110 J=2,MM
         K = M-J
         Y(K) = Y(K)-D(K)*Y(K+1)-U(K)*Y(M)
  110 CONTINUE
  111 IF (M1) 112,112,114
  112 IF (M2) 123,123,113
  113 RT = BM2(M2)
      M2 = M2-1
      GO TO 119
  114 IF (M2) 115,115,116
  115 RT = BM1(M1)
      M1 = M1-1
      GO TO 119
  116 IF (ABS(BM1(M1))-ABS(BM2(M2))) 118,118,117
  117 RT = BM1(M1)
      M1 = M1-1
      GO TO 119
  118 RT = BM2(M2)
      M2 = M2-1
C
C MATRIX MULTIPLICATION
C
  119 YH = Y(1)
      Y1 = (B(1)-RT)*Y(1)+C(1)*Y(2)+A(1)*Y(M)
      IF (MM-2) 122,120,120
  120 DO 121 J=2,MM
         Y2 = A(J)*Y(J-1)+(B(J)-RT)*Y(J)+C(J)*Y(J+1)
         Y(J-1) = Y1
         Y1 = Y2
  121 CONTINUE
  122 Y(M) = A(M)*Y(M-1)+(B(M)-RT)*Y(M)+C(M)*YH
      Y(M-1) = Y1
      IFLG = 1
      GO TO 102
  123 IF (IA) 126,126,124
  124 RT = AA(IA)
      IA = IA-1
      IFLG = 1
C
C SCALAR MULTIPLICATION
C
      DO 125 J=1,M
         Y(J) = RT*Y(J)
  125 CONTINUE
  126 IF (IFLG) 127,127,102
  127 DO 128 J=1,M
         YY(J) = REAL(Y(J))
  128 CONTINUE
      RETURN
      END
      SUBROUTINE PPADD (N,IERROR,A,C,    BP,BH)
c     SUBROUTINE PPADD (N,IERROR,A,C,CBP,BP,BH)
C
C     PPADD COMPUTES THE EIGENVALUES OF THE PERIODIC TRIDIAGONAL MATRIX
C     WITH COEFFICIENTS AN,BN,CN
C
C N IS THE ORDER OF THE BH AND BP POLYNOMIALS
C     ON OUTPUT BP CONTIANS THE EIGENVALUES
C CBP IS THE SAME AS BP EXCEPT TYPE COMPLEX
C BH IS USED TO TEMPORARILY STORE THE ROOTS OF THE B HAT POLYNOMIAL
C WHICH ENTERS THROUGH BP
C
      COMPLEX         CF         ,CX         ,FSG        ,HSG        ,
     1                DD         ,F          ,FP         ,FPP        ,
     2                CDIS       ,R1         ,R2         ,R3         ,
     3                CBP
c     DIMENSION       A(1)       ,C(1)       ,BP(1)      ,BH(1)      ,
c    1                CBP(1)
      DIMENSION       A(N)       ,C(N)       ,BP(N)      ,BH(N)      ,
     1                CBP(N)
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
      EXTERNAL        PSGF       ,PPSPF      ,PPSGF
      SCNV = SQRT(CNV)
      IZ = N
      IZM = IZ-1
      IZM2 = IZ-2
      IF (BP(N)-BP(1)) 101,142,103
  101 DO 102 J=1,N
         NT = N-J
         BH(J) = BP(NT+1)
  102 CONTINUE
      GO TO 105
  103 DO 104 J=1,N
         BH(J) = BP(J)
  104 CONTINUE
  105 NCMPLX = 0
      MODIZ = MOD(IZ,2)
      IS = 1
      IF (MODIZ) 106,107,106
  106 IF (A(1)) 110,142,107
  107 XL = BH(1)
      DB = BH(3)-BH(1)
  108 XL = XL-DB
      IF (PSGF(XL,IZ,C,A,BH)) 108,108,109
  109 SGN = -1.
      CBP(1) = CMPLX(BSRH(XL,BH(1),IZ,C,A,BH,PSGF,SGN),0.)
      IS = 2
  110 IF = IZ-1
      IF (MODIZ) 111,112,111
  111 IF (A(1)) 112,142,115
  112 XR = BH(IZ)
      DB = BH(IZ)-BH(IZ-2)
  113 XR = XR+DB
      IF (PSGF(XR,IZ,C,A,BH)) 113,114,114
  114 SGN = 1.
      CBP(IZ) = CMPLX(BSRH(BH(IZ),XR,IZ,C,A,BH,PSGF,SGN),0.)
      IF = IZ-2
  115 DO 136 IG=IS,IF,2
         XL = BH(IG)
         XR = BH(IG+1)
         SGN = -1.
         XM = BSRH(XL,XR,IZ,C,A,BH,PPSPF,SGN)
         PSG = PSGF(XM,IZ,C,A,BH)
         IF (ABS(PSG)-EPS) 118,118,116
  116    IF (PSG*PPSGF(XM,IZ,C,A,BH)) 117,118,119
C
C     CASE OF A REAL ZERO
C
  117    SGN = 1.
         CBP(IG) = CMPLX(BSRH(BH(IG),XM,IZ,C,A,BH,PSGF,SGN),0.)
         SGN = -1.
         CBP(IG+1) = CMPLX(BSRH(XM,BH(IG+1),IZ,C,A,BH,PSGF,SGN),0.)
         GO TO 136
C
C     CASE OF A MULTIPLE ZERO
C
  118    CBP(IG) = CMPLX(XM,0.)
         CBP(IG+1) = CMPLX(XM,0.)
         GO TO 136
C
C     CASE OF A COMPLEX ZERO
C
  119    IT = 0
         ICV = 0
         CX = CMPLX(XM,0.)
  120    FSG = (1.,0.)
         HSG = (1.,0.)
         FP = (0.,0.)
         FPP = (0.,0.)
         DO 121 J=1,IZ
            DD = 1./(CX-BH(J))
            FSG = FSG*A(J)*DD
            HSG = HSG*C(J)*DD
            FP = FP+DD
            FPP = FPP-DD*DD
  121    CONTINUE
         IF (MODIZ) 123,122,123
  122    F = (1.,0.)-FSG-HSG
         GO TO 124
  123    F = (1.,0.)+FSG+HSG
  124    I3 = 0
         IF (CABS(FP)) 126,126,125
  125    I3 = 1
         R3 = -F/FP
  126    I2 = 0
         IF (CABS(FPP)) 132,132,127
  127    I2 = 1
         CDIS = CSQRT(FP**2-2.*F*FPP)
         R1 = CDIS-FP
         R2 = -FP-CDIS
         IF (CABS(R1)-CABS(R2)) 129,129,128
  128    R1 = R1/FPP
         GO TO 130
  129    R1 = R2/FPP
  130    R2 = 2.*F/FPP/R1
         IF (CABS(R2) .LT. CABS(R1)) R1 = R2
         IF (I3) 133,133,131
  131    IF (CABS(R3) .LT. CABS(R1)) R1 = R3
         GO TO 133
  132    R1 = R3
  133    CX = CX+R1
         IT = IT+1
         IF (IT .GT. 50) GO TO 142
         IF (CABS(R1) .GT. SCNV) GO TO 120
         IF (ICV) 134,134,135
  134    ICV = 1
         GO TO 120
  135    CBP(IG) = CX
         CBP(IG+1) = CONJG(CX)
  136 CONTINUE
      IF (CABS(CBP(N))-CABS(CBP(1))) 137,142,139
  137 NHALF = N/2
      DO 138 J=1,NHALF
         NT = N-J
         CX = CBP(J)
         CBP(J) = CBP(NT+1)
         CBP(NT+1) = CX
  138 CONTINUE
  139 NCMPLX = 1
      DO 140 J=2,IZ
         IF (AIMAG(CBP(J))) 143,140,143
  140 CONTINUE
      NCMPLX = 0
      DO 141 J=2,IZ
         BP(J) = REAL(CBP(J))
  141 CONTINUE
      GO TO 143
  142 IERROR = 4
  143 CONTINUE
      RETURN
      END
      FUNCTION PSGF (X,IZ,C,A,BH)
      DIMENSION       A(1)       ,C(1)       ,BH(1)
      FSG = 1.
      HSG = 1.
      DO 101 J=1,IZ
         DD = 1./(X-BH(J))
         FSG = FSG*A(J)*DD
         HSG = HSG*C(J)*DD
  101 CONTINUE
      IF (MOD(IZ,2)) 103,102,103
  102 PSGF = 1.-FSG-HSG
      RETURN
  103 PSGF = 1.+FSG+HSG
      RETURN
      END
      FUNCTION BSRH (XLL,XRR,IZ,C,A,BH,F,SGN)
      DIMENSION       A(1)       ,C(1)       ,BH(1)
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
      XL = XLL
      XR = XRR
      DX = .5*ABS(XR-XL)
  101 X = .5*(XL+XR)
      IF (SGN*F(X,IZ,C,A,BH)) 103,105,102
  102 XR = X
      GO TO 104
  103 XL = X
  104 DX = .5*DX
      IF (DX-CNV) 105,105,101
  105 BSRH = .5*(XL+XR)
      RETURN
      END
      FUNCTION PPSGF (X,IZ,C,A,BH)
      DIMENSION       A(1)       ,C(1)       ,BH(1)
      SUM = 0.
      DO 101 J=1,IZ
         SUM = SUM-1./(X-BH(J))**2
  101 CONTINUE
      PPSGF = SUM
      RETURN
      END
      FUNCTION PPSPF (X,IZ,C,A,BH)
      DIMENSION       A(1)       ,C(1)       ,BH(1)
      SUM = 0.
      DO 101 J=1,IZ
         SUM = SUM+1./(X-BH(J))
  101 CONTINUE
      PPSPF = SUM
      RETURN
      END
      SUBROUTINE COMPB (N,IERROR,AN,BN,CN,B,AH,BH)
C
C     COMPB COMPUTES THE ROOTS OF THE B POLYNOMIALS USING SUBROUTINE
C     TEVLS WHICH IS A MODIFICATION THE EISPACK PROGRAM TQLRAT.
C     IERROR IS SET TO 4 IF EITHER TEVLS FAILS OR IF A(J+1)*C(J) IS
C     LESS THAN ZERO FOR SOME J.  AH,BH ARE TEMPORARY WORK ARRAYS.
C
      DIMENSION       AN(1)      ,BN(1)      ,CN(1)      ,B(1)       ,
     1                AH(1)      ,BH(1)
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
      EPS = EPMACH(DUM)
      BNORM = ABS(BN(1))
      DO 102 J=2,NM
         BNORM = AMAX1(BNORM,ABS(BN(J)))
         ARG = AN(J)*CN(J-1)
         IF (ARG) 119,101,101
  101    B(J) = SIGN(SQRT(ARG),AN(J))
  102 CONTINUE
      CNV = EPS*BNORM
      IF = 2**K
      KDO = K-1
      DO 108 L=1,KDO
         IR = L-1
         I2 = 2**IR
         I4 = I2+I2
         IPL = I4-1
         IFD = IF-I4
         DO 107 I=I4,IFD,I4
            CALL INDXB (I,L,IB,NB)
            IF (NB) 108,108,103
  103       JS = I-IPL
            JF = JS+NB-1
            LS = 0
            DO 104 J=JS,JF
               LS = LS+1
               BH(LS) = BN(J)
               AH(LS) = B(J)
  104       CONTINUE
            CALL TEVLS (NB,BH,AH,IERROR)
            IF (IERROR) 118,105,118
  105       LH = IB-1
            DO 106 J=1,NB
               LH = LH+1
               B(LH) = -BH(J)
  106       CONTINUE
  107    CONTINUE
  108 CONTINUE
      DO 109 J=1,NM
         B(J) = -BN(J)
  109 CONTINUE
      IF (NPP) 117,110,117
  110 NMP = NM+1
      NB = NM+NMP
      DO 112 J=1,NB
         L1 = MOD(J-1,NMP)+1
         L2 = MOD(J+NM-1,NMP)+1
         ARG = AN(L1)*CN(L2)
         IF (ARG) 119,111,111
  111    BH(J) = SIGN(SQRT(ARG),-AN(L1))
         AH(J) = -BN(L1)
  112 CONTINUE
      CALL TEVLS (NB,AH,BH,IERROR)
      IF (IERROR) 118,113,118
  113 CALL INDXB (IF,K-1,J2,LH)
      CALL INDXB (IF/2,K-1,J1,LH)
      J2 = J2+1
      LH = J2
      N2M2 = J2+NM+NM-2
  114 D1 = ABS(B(J1)-B(J2-1))
      D2 = ABS(B(J1)-B(J2))
      D3 = ABS(B(J1)-B(J2+1))
      IF ((D2 .LT. D1) .AND. (D2 .LT. D3)) GO TO 115
      B(LH) = B(J2)
      J2 = J2+1
      LH = LH+1
      IF (J2-N2M2) 114,114,116
  115 J2 = J2+1
      J1 = J1+1
      IF (J2-N2M2) 114,114,116
  116 B(LH) = B(N2M2+1)
      CALL INDXB (IF,K-1,J1,J2)
      J2 = J1+NMP+NMP
      CALL PPADD (NM+1,IERROR,AN,CN,      B(J1),B(J2))
c     CALL PPADD (NM+1,IERROR,AN,CN,B(J1),B(J1),B(J2))
  117 RETURN
  118 IERROR = 4
      RETURN
  119 IERROR = 5
      RETURN
      END
      SUBROUTINE TEVLS (N,D,E2,IERR)
C
      INTEGER         I          ,J          ,L          ,M          ,
     1                N          ,II         ,L1         ,MML        ,
     2                IERR
      REAL            D(N)       ,E2(N)
      REAL            B          ,C          ,F          ,G          ,
     1                H          ,P          ,R          ,S          ,
     2                MACHEP
C
C     REAL SQRT,ABS,SIGN
C
      COMMON /CBLKT/  NPP        ,K          ,MACHEP     ,CNV        ,
     1                NM         ,NCMPLX     ,IK
C
C     THIS SUBROUTINE IS A MODIFICATION OF THE EISPACK SUBROUTINE TQLRAT
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E2 CONTAINS THE                SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES,
C
C        E2 HAS BEEN DESTROYED,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
C
      IERR = 0
      IF (N .EQ. 1) GO TO 115
C
      DO 101 I=2,N
         E2(I-1) = E2(I)*E2(I)
  101 CONTINUE
C
      F = 0.0
      B = 0.0
      E2(N) = 0.0
C
      DO 112 L=1,N
         J = 0
         H = MACHEP*(ABS(D(L))+SQRT(E2(L)))
         IF (B .GT. H) GO TO 102
         B = H
         C = B*B
C
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
C
  102    DO 103 M=L,N
            IF (E2(M) .LE. C) GO TO 104
C
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
C
  103    CONTINUE
C
  104    IF (M .EQ. L) GO TO 108
  105    IF (J .EQ. 30) GO TO 114
         J = J+1
C
C     ********** FORM SHIFT **********
C
         L1 = L+1
         S = SQRT(E2(L))
         G = D(L)
         P = (D(L1)-G)/(2.0*S)
         R = SQRT(P*P+1.0)
         D(L) = S/(P+SIGN(R,P))
         H = G-D(L)
C
         DO 106 I=L1,N
            D(I) = D(I)-H
  106    CONTINUE
C
         F = F+H
C
C     ********** RATIONAL QL TRANSFORMATION **********
C
         G = D(M)
         IF (G .EQ. 0.0) G = B
         H = G
         S = 0.0
         MML = M-L
C
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
C
         DO 107 II=1,MML
            I = M-II
            P = G*H
            R = P+E2(I)
            E2(I+1) = S*R
            S = E2(I)/R
            D(I+1) = H+S*(H+D(I))
            G = D(I)-E2(I)/G
            IF (G .EQ. 0.0) G = B
            H = G*P/R
  107    CONTINUE
C
         E2(L) = S*G
         D(L) = H
C
C     ********** GUARD AGAINST UNDERFLOWED H **********
C
         IF (H .EQ. 0.0) GO TO 108
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 108
         E2(L) = H*E2(L)
         IF (E2(L) .NE. 0.0) GO TO 105
  108    P = D(L)+F
C
C     ********** ORDER EIGENVALUES **********
C
         IF (L .EQ. 1) GO TO 110
C
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
C
         DO 109 II=2,L
            I = L+2-II
            IF (P .GE. D(I-1)) GO TO 111
            D(I) = D(I-1)
  109    CONTINUE
C
  110    I = 1
  111    D(I) = P
  112 CONTINUE
C
      IF (ABS(D(N)) .GE. ABS(D(1))) GO TO 115
      NHALF = N/2
      DO 113 I=1,NHALF
         NTOP = N-I
         DHOLD = D(I)
         D(I) = D(NTOP+1)
         D(NTOP+1) = DHOLD
  113 CONTINUE
      GO TO 115
C
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
C
  114 IERR = L
  115 RETURN
C
C     ********** LAST CARD OF TQLRAT **********
C
C
C REVISION HISTORY---
C
C DECEMBER 1979    FIRST ADDED TO NSSL
C-----------------------------------------------------------------------
      END
C PACKAGE COMF           THE ENTRIES IN THIS PACKAGE ARE LOWLEVEL
C                        ENTRIES, SUPPORTING ULIB PACKAGES BLKTRI
c                        AND CBLKTRI. THAT IS, THESE ROUTINES ARE
C                        NOT CALLED DIRECTLY BY USERS, BUT RATHER
C                        BY ENTRIES WITHIN BLKTRI AND CBLKTRI.
C                        DESCRIPTION OF ENTRIES EPMACH AND PIMACH
C                        FOLLOW BELOW.
C
C LATEST REVISION        JANUARY 1985
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       NONE
C FILES
C
C LANGUAGE               FORTRAN
C ********************************************************************
C
C FUNCTION EPMACH (DUM)
C
C PURPOSE                TO COMPUTE AN APPROXIMATE MACHINE ACCURACY
C                        EPSILON ACCORDING TO THE FOLLOWING DEFINITION:
C                        EPSILON IS THE SMALLEST NUMBER SUCH THAT
C                        (1.+EPSILON).GT.1.)
C
C USAGE                  EPS = EPMACH (DUM)
C
C ARGUMENTS
C ON INPUT               DUM
C                          DUMMY VALUE
C
C ARGUMENTS
C ON OUTPUT              NONE
C
C HISTORY                THE ORIGINAL VERSION, WRITTEN WHEN THE
C                        BLKTRI PACKAGE WAS CONVERTED FROM THE
C                        CDC 7600 TO RUN ON THE CRAY-1, CALCULATED
C                        MACHINE ACCURACY BY SUCCESSIVE DIVISIONS
C                        BY 10.  USE OF THIS CONSTANT CAUSED BLKTRI
C                        TO COMPUTE SOLUTIONS ON THE CRAY-1 WITH FOUR
C                        FEWER PLACES OF ACCURACY THAN THE VERSION
C                        ON THE 7600.  IT WAS FOUND THAT COMPUTING
C                        MACHINE ACCURACY BY SUCCESSIVE DIVISIONS
C                        OF 2 PRODUCED A MACHINE ACCURACY 29% LESS
C                        THAN THE VALUE GENERATED BY SUCCESSIVE
C                        DIVISIONS BY 10, AND THAT USE OF THIS
C                        MACHINE CONSTANT IN THE BLKTRI PACKAGE
C                        RECOVERED THE ACCURACY THAT APPEARED TO
C                        BE LOST ON CONVERSION.
C
C ALGORITHM              COMPUTES MACHINE ACCURACY BY SUCCESSIVE
C                        DIVISIONS OF TWO.
C
C PORTABILITY            THIS CODE WILL EXECUTE ON MACHINES OTHER
C                        THAN THE CRAY1, BUT THE RETURNED VALUE MAY
C                        BE UNSATISFACTORY.  SEE HISTORY ABOVE.
C ********************************************************************
C
C FUNCTION PIMACH (DUM)
C
C PURPOSE                TO SUPPLY THE VALUE OF THE CONSTANT PI
C                        CORRECT TO MACHINE PRECISION WHERE
C                        PI=3.141592653589793238462643383279502884197
C                             1693993751058209749446
C
C USAGE                  PI = PIMACH (DUM)
C
C ARGUMENTS
C ON INPUT               DUM
C                          DUMMY VALUE
C
C ARGUMENTS
C ON OUTPUT              NONE
C
C ALGORITHM              THE VALUE OF PI IS SET IN A CONSTANT.
C
C PORTABILITY            THIS ENTRY IS PORTABLE, BUT USERS SHOULD
C                        CHECK TO SEE WHETHER GREATER ACCURACY IS
C                        REQUIRED.
C
C***********************************************************************
      FUNCTION EPMACH (DUM)
      COMMON /VALUE/  V
      EPS = 1.
  101 EPS = EPS/2.
      CALL STORE (EPS+1.)
      IF (V-1.) 102,102,101
  102 EPMACH = 100.*EPS
      RETURN
      END
      SUBROUTINE STORE (X)
      COMMON /VALUE/  V
      V = X
      RETURN
      END
      FUNCTION PIMACH (DUM)
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C
      PIMACH = 3.14159265358979
      RETURN
      END
C PACKAGE SEPAUX         CONTAINS NO USER ENTRY POINTS.
C
C LATEST REVISION        MARCH 1985
C
C PURPOSE                THIS PACKAGE CONTAINS AUXILIARY ROUTINES FOR
C                        NCAR PUBLIC SOFTWARE PACKAGES SUCH AS SEPELI
C                        AND SEPX4.
C
C USAGE                  SINCE THIS PACKAGE CONTAINS NO USER ENTRIES,
C                        NO USAGE INSTRUCTIONS OR ARGUMENT DESCRIPTIONS
C                        ARE GIVEN HERE.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       NONE
C FILES
C
C LANGUAGE               FORTRAN
C
C HISTORY                DEVELOPED IN THE LATE 1970'S BY JOHN C. ADAMS
C                        OF NCAR'S SCIENTTIFIC COMPUTING DIVISION.
C
C PORTABILITY            FORTRAN 66
C **********************************************************************
      SUBROUTINE SEPORT (USOL,IDMN,ZN,ZM,PERTRB)
C
C     THIS SUBROUTINE ORTHOGANALIZES THE ARRAY USOL WITH RESPECT TO
C     THE CONSTANT ARRAY IN A WEIGHTED LEAST SQUARES NORM
C
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       USOL(IDMN,1)           ,ZN(1)      ,ZM(1)
      ISTR = IS
      IFNL = MS
      JSTR = JS
      JFNL = NS
C
C     COMPUTE WEIGHTED INNER PRODUCTS
C
      UTE = 0.0
      ETE = 0.0
      DO  20 I=IS,MS
         II = I-IS+1
         DO  10 J=JS,NS
            JJ = J-JS+1
            ETE = ETE+ZM(II)*ZN(JJ)
            UTE = UTE+USOL(I,J)*ZM(II)*ZN(JJ)
   10    CONTINUE
   20 CONTINUE
C
C     SET PERTURBATION PARAMETER
C
      PERTRB = UTE/ETE
C
C     SUBTRACT OFF CONSTANT PERTRB
C
      DO  40 I=ISTR,IFNL
         DO  30 J=JSTR,JFNL
            USOL(I,J) = USOL(I,J)-PERTRB
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SEPMIN (USOL,IDMN,ZN,ZM,PERTB)
C
C     THIS SUBROUTINE ORHTOGONALIZES THE ARRAY USOL WITH RESPECT TO
C     THE CONSTANT ARRAY IN A WEIGHTED LEAST SQUARES NORM
C
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       USOL(IDMN,1)           ,ZN(1)      ,ZM(1)
C
C     ENTRY AT SEPMIN OCCURRS WHEN THE FINAL SOLUTION IS
C     TO BE MINIMIZED WITH RESPECT TO THE WEIGHTED
C     LEAST SQUARES NORM
C
      ISTR = 1
      IFNL = K
      JSTR = 1
      JFNL = L
C
C     COMPUTE WEIGHTED INNER PRODUCTS
C
      UTE = 0.0
      ETE = 0.0
      DO  20 I=IS,MS
         II = I-IS+1
         DO  10 J=JS,NS
            JJ = J-JS+1
            ETE = ETE+ZM(II)*ZN(JJ)
            UTE = UTE+USOL(I,J)*ZM(II)*ZN(JJ)
   10    CONTINUE
   20 CONTINUE
C
C     SET PERTURBATION PARAMETER
C
      PERTRB = UTE/ETE
C
C     SUBTRACT OFF CONSTANT PERTRB
C
      DO  40 I=ISTR,IFNL
         DO  30 J=JSTR,JFNL
            USOL(I,J) = USOL(I,J)-PERTRB
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SEPTRI (N,A,B,C,D,U,Z)
C
C     THIS SUBROUTINE SOLVES FOR A NON-ZERO EIGENVECTOR CORRESPONDING
C     TO THE ZERO EIGENVALUE OF THE TRANSPOSE OF THE RANK
C     DEFICIENT ONE MATRIX WITH SUBDIAGONAL A, DIAGONAL B, AND
C     SUPERDIAGONAL C , WITH A(1) IN THE (1,N) POSITION, WITH
C     C(N) IN THE (N,1) POSITION, AND ALL OTHER ELEMENTS ZERO.
C
      DIMENSION       A(N)       ,B(N)       ,C(N)       ,D(N)       ,
     1                U(N)       ,Z(N)
      BN = B(N)
      D(1) = A(2)/B(1)
      V = A(1)
      U(1) = C(N)/B(1)
      NM2 = N-2
      DO  10 J=2,NM2
         DEN = B(J)-C(J-1)*D(J-1)
         D(J) = A(J+1)/DEN
         U(J) = -C(J-1)*U(J-1)/DEN
         BN = BN-V*U(J-1)
         V = -V*D(J-1)
   10 CONTINUE
      DEN = B(N-1)-C(N-2)*D(N-2)
      D(N-1) = (A(N)-C(N-2)*U(N-2))/DEN
      AN = C(N-1)-V*D(N-2)
      BN = BN-V*U(N-2)
      DEN = BN-AN*D(N-1)
C
C     SET LAST COMPONENT EQUAL TO ONE
C
      Z(N) = 1.0
      Z(N-1) = -D(N-1)
      NM1 = N-1
      DO  20 J=2,NM1
         K = N-J
         Z(K) = -D(K)*Z(K+1)-U(K)*Z(N)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE SEPDX (U,IDMN,I,J,UXXX,UXXXX)
C
C     THIS PROGRAM COMPUTES SECOND ORDER FINITE DIFFERENCE
C     APPROXIMATIONS TO THE THIRD AND FOURTH X
C     PARTIAL DERIVATIVES OF U AT THE (I,J) MESH POINT
C
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
c     DIMENSION       U(IDMN,1)
      DIMENSION       U(IDMN,NS)
      IF (I.GT.2 .AND. I.LT.(K-1)) GO TO  50
      IF (I .EQ. 1) GO TO  10
      IF (I .EQ. 2) GO TO  30
      IF (I .EQ. K-1) GO TO  60
      IF (I .EQ. K) GO TO  80
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A
C
   10 IF (KSWX .EQ. 1) GO TO  20
      UXXX = (-5.0*U(1,J)+18.0*U(2,J)-24.0*U(3,J)+14.0*U(4,J)-
     1                                               3.0*U(5,J))/(TDLX3)
      UXXXX = (3.0*U(1,J)-14.0*U(2,J)+26.0*U(3,J)-24.0*U(4,J)+
     1                                      11.0*U(5,J)-2.0*U(6,J))/DLX4
      RETURN
C
C     PERIODIC AT X=A
C
   20 UXXX = (-U(K-2,J)+2.0*U(K-1,J)-2.0*U(2,J)+U(3,J))/(TDLX3)
      UXXXX = (U(K-2,J)-4.0*U(K-1,J)+6.0*U(1,J)-4.0*U(2,J)+U(3,J))/DLX4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A+DLX
C
   30 IF (KSWX .EQ. 1) GO TO  40
      UXXX = (-3.0*U(1,J)+10.0*U(2,J)-12.0*U(3,J)+6.0*U(4,J)-U(5,J))/
     1       TDLX3
      UXXXX = (2.0*U(1,J)-9.0*U(2,J)+16.0*U(3,J)-14.0*U(4,J)+6.0*U(5,J)-
     1                                                      U(6,J))/DLX4
      RETURN
C
C     PERIODIC AT X=A+DLX
C
   40 UXXX = (-U(K-1,J)+2.0*U(1,J)-2.0*U(3,J)+U(4,J))/(TDLX3)
      UXXXX = (U(K-1,J)-4.0*U(1,J)+6.0*U(2,J)-4.0*U(3,J)+U(4,J))/DLX4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
C
   50 CONTINUE
      UXXX = (-U(I-2,J)+2.0*U(I-1,J)-2.0*U(I+1,J)+U(I+2,J))/TDLX3
      UXXXX = (U(I-2,J)-4.0*U(I-1,J)+6.0*U(I,J)-4.0*U(I+1,J)+U(I+2,J))/
     1        DLX4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B-DLX
C
   60 IF (KSWX .EQ. 1) GO TO  70
      UXXX = (U(K-4,J)-6.0*U(K-3,J)+12.0*U(K-2,J)-10.0*U(K-1,J)+
     1                                                 3.0*U(K,J))/TDLX3
      UXXXX = (-U(K-5,J)+6.0*U(K-4,J)-14.0*U(K-3,J)+16.0*U(K-2,J)-
     1                                     9.0*U(K-1,J)+2.0*U(K,J))/DLX4
      RETURN
C
C     PERIODIC AT X=B-DLX
C
   70 UXXX = (-U(K-3,J)+2.0*U(K-2,J)-2.0*U(1,J)+U(2,J))/TDLX3
      UXXXX = (U(K-3,J)-4.0*U(K-2,J)+6.0*U(K-1,J)-4.0*U(1,J)+U(2,J))/
     1        DLX4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B
C
   80 UXXX = -(3.0*U(K-4,J)-14.0*U(K-3,J)+24.0*U(K-2,J)-18.0*U(K-1,J)+
     1                                                 5.0*U(K,J))/TDLX3
      UXXXX = (-2.0*U(K-5,J)+11.0*U(K-4,J)-24.0*U(K-3,J)+26.0*U(K-2,J)-
     1                                    14.0*U(K-1,J)+3.0*U(K,J))/DLX4
      RETURN
      END
      SUBROUTINE SEPDY (U,IDMN,I,J,UYYY,UYYYY)
C
C     THIS PROGRAM COMPUTES SECOND ORDER FINITE DIFFERENCE
C     APPROXIMATIONS TO THE THIRD AND FOURTH Y
C     PARTIAL DERIVATIVES OF U AT THE (I,J) MESH POINT
C
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
c     DIMENSION       U(IDMN,1)
      DIMENSION       U(IDMN,NS)
      IF (J.GT.2 .AND. J.LT.(L-1)) GO TO  50
      IF (J .EQ. 1) GO TO  10
      IF (J .EQ. 2) GO TO  30
      IF (J .EQ. L-1) GO TO  60
      IF (J .EQ. L) GO TO  80
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C
C
   10 IF (KSWY .EQ. 1) GO TO  20
      UYYY = (-5.0*U(I,1)+18.0*U(I,2)-24.0*U(I,3)+14.0*U(I,4)-
     1                                                 3.0*U(I,5))/TDLY3
      UYYYY = (3.0*U(I,1)-14.0*U(I,2)+26.0*U(I,3)-24.0*U(I,4)+
     1                                      11.0*U(I,5)-2.0*U(I,6))/DLY4
      RETURN
C
C     PERIODIC AT X=A
C
   20 UYYY = (-U(I,L-2)+2.0*U(I,L-1)-2.0*U(I,2)+U(I,3))/TDLY3
      UYYYY = (U(I,L-2)-4.0*U(I,L-1)+6.0*U(I,1)-4.0*U(I,2)+U(I,3))/DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C+DLY
C
   30 IF (KSWY .EQ. 1) GO TO  40
      UYYY = (-3.0*U(I,1)+10.0*U(I,2)-12.0*U(I,3)+6.0*U(I,4)-U(I,5))/
     1       TDLY3
      UYYYY = (2.0*U(I,1)-9.0*U(I,2)+16.0*U(I,3)-14.0*U(I,4)+6.0*U(I,5)-
     1                                                      U(I,6))/DLY4
      RETURN
C
C     PERIODIC AT Y=C+DLY
C
   40 UYYY = (-U(I,L-1)+2.0*U(I,1)-2.0*U(I,3)+U(I,4))/TDLY3
      UYYYY = (U(I,L-1)-4.0*U(I,1)+6.0*U(I,2)-4.0*U(I,3)+U(I,4))/DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
C
   50 CONTINUE
      UYYY = (-U(I,J-2)+2.0*U(I,J-1)-2.0*U(I,J+1)+U(I,J+2))/TDLY3
      UYYYY = (U(I,J-2)-4.0*U(I,J-1)+6.0*U(I,J)-4.0*U(I,J+1)+U(I,J+2))/
     1        DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D-DLY
C
   60 IF (KSWY .EQ. 1) GO TO  70
      UYYY = (U(I,L-4)-6.0*U(I,L-3)+12.0*U(I,L-2)-10.0*U(I,L-1)+
     1                                                 3.0*U(I,L))/TDLY3
      UYYYY = (-U(I,L-5)+6.0*U(I,L-4)-14.0*U(I,L-3)+16.0*U(I,L-2)-
     1                                     9.0*U(I,L-1)+2.0*U(I,L))/DLY4
      RETURN
C
C     PERIODIC AT Y=D-DLY
C
   70 CONTINUE
      UYYY = (-U(I,L-3)+2.0*U(I,L-2)-2.0*U(I,1)+U(I,2))/TDLY3
      UYYYY = (U(I,L-3)-4.0*U(I,L-2)+6.0*U(I,L-1)-4.0*U(I,1)+U(I,2))/
     1        DLY4
      RETURN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D
C
   80 UYYY = -(3.0*U(I,L-4)-14.0*U(I,L-3)+24.0*U(I,L-2)-18.0*U(I,L-1)+
     1                                                 5.0*U(I,L))/TDLY3
      UYYYY = (-2.0*U(I,L-5)+11.0*U(I,L-4)-24.0*U(I,L-3)+26.0*U(I,L-2)-
     1                                    14.0*U(I,L-1)+3.0*U(I,L))/DLY4
      RETURN
C
C REVISION HISTORY---
C
C DECEMBER 1979    FIRST ADDED TO NSSL
C
      END

