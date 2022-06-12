C********************** GETGRID ***********************
      PROGRAM POST_ONE_POINT
      IMPLICIT NONE
      INCLUDE 'PARAM.H'
      COMMON/MESH/Y(0:NY),DX,DZ
      COMMON/DIRECTORY/DIR_INS,DIR_SAVE

      REAL*8 Y, DX, DZ
      CHARACTER*100 DIR_INS
      CHARACTER*100 DIR_SAVE

      ! SET DIRECTORY TO READ INS FIELD
      DIR_INS='/home/kimsh8037/vscode/post_sample/INSFIELD6/'
      DIR_SAVE='meanresult/OUTPUT_'

      CALL GET_GRID
      CALL INDICES
      CALL CALC_AVG
    
      STOP
      END

C********************** GETGRID ***********************     
      SUBROUTINE GET_GRID
      IMPLICIT NONE

      INCLUDE 'PARAM.H'
      COMMON/MESH/Y(0:NY),DX,DZ
      COMMON/MESH2/DY(0:NY),H(NY)

      REAL*8 Y, DX, DZ,DY,H
      REAL*8 PI
      INTEGER J

      PI = ACOS(-1.)
      ! CALCULATE DX, DZ
      DX = LX*PI/REAL(NX-1)
      DZ = LZ*PI/REAL(NZ-1)
      
      ! OPEN GRID FILE
      OPEN(2,FILE='grid.dat',ACTION='READ')
      ! READ DATA
      READ(2,*) (Y(J),J=0,NY)
      ! CLOSE FILE
      CLOSE(2)

      ! CALCULATE DY, H(SEE CODE NOTE)
      DY(1)=Y(2)
      DO 20 J=2,(N2-1)
      DY(J)=Y(J+1)-Y(J)
      H(J)=0.5*(DY(J)+DY(J-1))
 20   CONTINUE
      H(1)=0.5*DY(1)
      H(NY)=0.5*DY(N2)

      RETURN
      END

C********************** INDICES ***********************     
      SUBROUTINE INDICES
      IMPLICIT NONE
      INCLUDE 'PARAM.H'
      COMMON/INDX/IPA(NX),IMA(NX),KPA(NZ),KMA(NZ)
      COMMON/INDX2/JPA(NY),JMU(NY),JMV(NY)
      COMMON/FINDX2/FJPA(NY),FJMU(NY),FJMV(NY)

      INTEGER IC,KC,JC

      INTEGER IPA,IMA,KPA,KMA,JPA,JMU,JMV,FJPA,FJMU,FJMV

      DO 10 IC=1,N1
      IPA(IC)=IC+1
 10   IMA(IC)=IC-1
      IPA(N1)=1
      IMA(1)=N1

      DO 20 KC=1,N3
      KPA(KC)=KC+1
 20   KMA(KC)=KC-1
      KPA(N3)=1
      KMA(1)=N3

      DO 30 JC=1,N2
      JPA(JC)=JC+1
      JMU(JC)=JC-1
 30   JMV(JC)=JC-1
      JPA(N2)=N2
      JMU(1)=1
      JMV(2)=2

      DO 40 JC=1,N2
      FJPA(JC)=JPA(JC)-JC   ! AT BOUNDARY, THESE VALUES ARE "0".
      FJMU(JC)=JC-JMU(JC)   ! OTHERWISE, THESE ARE "1".
 40   FJMV(JC)=JC-JMV(JC)

      RETURN
      END

C********************** INDICES ***********************     
      SUBROUTINE GET_FILENAME_DISK(FILENAME,NTIME,NV)
C     NV=1 : U
C        2 : V
C        3 : W
      IMPLICIT NONE

      COMMON/DIRECTORY/DIR_INS,DIR_SAVE
      INTEGER NV, NN, NTIME

      CHARACTER*100 FILENAME
      CHARACTER*100 DIR_INS
      CHARACTER*100 DIR_SAVE

      IF (NV.EQ.1) THEN

      FILENAME=DIR_INS
      NN=INDEX(FILENAME,'6')
      WRITE(UNIT=FILENAME(NN+2:),FMT='(BN,A6)') 'INS_U.'
      WRITE(UNIT=FILENAME(NN+8:),FMT='(BN,I6.6)') NTIME

      ENDIF

      IF (NV.EQ.2) THEN


      FILENAME=DIR_INS
      NN=INDEX(FILENAME,'6')
      WRITE(UNIT=FILENAME(NN+2:),FMT='(BN,A6)') 'INS_V.'
      WRITE(UNIT=FILENAME(NN+8:),FMT='(BN,I6.6)') NTIME

      ENDIF

      IF (NV.EQ.3) THEN


      FILENAME=DIR_INS
      NN=INDEX(FILENAME,'6')
      WRITE(UNIT=FILENAME(NN+2:),FMT='(BN,A6)') 'INS_W.'
      WRITE(UNIT=FILENAME(NN+8:),FMT='(BN,I6.6)') NTIME

      ENDIF

      RETURN
      END

C***************** GET VELOCITY ***********************
      SUBROUTINE GET_VEL(Q,NTIME)
      IMPLICIT NONE

      INCLUDE 'PARAM.H'

      REAL*8 Q(3,0:NX,0:NY,0:NZ)
      CHARACTER*100 FILENAME
      INTEGER NTIME,I,J,K,NV

C---  READ U

      CALL GET_FILENAME_DISK(FILENAME,NTIME,1)
      WRITE(*,*) FILENAME
      OPEN(11,FILE=FILENAME,STATUS='OLD',FORM='UNFORMATTED'
     &     ,ACTION='READ')
      READ(11) (((Q(1,I,J,K),I=0,NX),J=0,NY),K=0,NZ)
      CLOSE(11)

C---  READ V

      CALL GET_FILENAME_DISK(FILENAME,NTIME,2)
      WRITE(*,*) FILENAME
      OPEN(11,FILE=FILENAME,STATUS='OLD',FORM='UNFORMATTED'
     &     ,ACTION='READ')
      READ(11) (((Q(2,I,J,K),I=0,NX),J=0,NY),K=0,NZ)
      CLOSE(11)

C---  READ W

      CALL GET_FILENAME_DISK(FILENAME,NTIME,3)
      WRITE(*,*) FILENAME
      OPEN(11,FILE=FILENAME,STATUS='OLD',FORM='UNFORMATTED'
     &     ,ACTION='READ')
      READ(11) (((Q(3,I,J,K),I=0,NX),J=0,NY),K=0,NZ)
      CLOSE(11)

      RETURN
      END

c******************** CALC AVERAGE ********************
C     CODE BUILT BY : SEONGHWAN KIM
C     UPDEATED DATE : 2022.06.02.
c******************************************************
      SUBROUTINE CALC_AVG
      IMPLICIT NONE
      INCLUDE 'PARAM.H'
      COMMON/MESH/Y(0:NY),DX,DZ
      COMMON/MESH2/DY(0:NY),H(NY)
      COMMON/INDX/IPA(NX),IMA(NX),KPA(NZ),KMA(NZ)
      COMMON/INDX2/JPA(NY),JMU(NY),JMV(NY)
      COMMON/FINDX2/FJPA(NY),FJMU(NY),FJMV(NY)
      COMMON/DIRECTORY/DIR_INS,DIR_SAVE

      REAL*8 Y,DX, DZ,DY,H

      INTEGER IPA,IMA,KPA,KMA,JPA,JMU,JMV,FJPA,FJMU,FJMV

      INTEGER NUM_INS,NTIME,I,J,K,L,NV,NTIME2
      REAL*8 YC,YP

      REAL*8 U(3,0:NX,0:NY,0:NZ)

      REAL*8 U1(3,0:N2) ! MEAN VELOCITY
      REAL*8 U2(6,0:N2) ! TMP_VALUE AT INSTANTANEOUS TIME
      REAL*8 U3(3,0:N2) ! RMS FLUCTUATION
      REAL*8 U4(6,0:N2) ! MINUS FLUCTUATION PRODUCT(STRESS)

      REAL*8 URMS(3,0:N2)

      INTEGER NN

      CHARACTER*100 FILENAME

      REAL*8 DUDY_WALL
      REAL*8 TAU_WALL
      REAL*8 UTAU
      REAL*8 RETAU
      REAL*8 DEL_V
      REAL*8 N2M
      REAL*8 DUDY(N2)
C     CUSTOM VARIABLE DEFINITON

      REAL T_START
      REAL T_INS
      REAL T_END
      REAL T_INS_START
      REAL T_INS_END
      REAL INS_PERCENT

      INTEGER N_INS

      CHARACTER*100 DIR_INS, DIR_SAVE

C     CALCULATION START
      T_START = SECNDS(0.0)

      N_INS = (N_END-N_INI)/N_SKIP + 1

      U1(:,:)=0.0
      U2(:,:)=0.0

C     INSTANTANEOUS TIME CALCULATION
      NUM_INS = 0
      DO 5 NTIME=N_INI,N_END,N_SKIP
      
      T_INS_START = SECNDS(0.0)

      NTIME2=NTIME
      NUM_INS = NUM_INS + 1
      ! GET INSTANTANEOUS VELOCITY
      CALL GET_VEL(U,NTIME2)
      ! OPERATE INSTANTANEOUS SPATIAL AVERAGING
      DO J=0,N2
      DO K=1,N3
      DO I=1,N1

      DO L = 1,3
      U1(L,J)=U1(L,J)+U(L,I,J,K)/DBLE(N1*N3*N_INS)
      U2(L,J)=U2(L,J)+U(L,I,J,K)**2/DBLE(N1*N3*N_INS)
      ENDDO
      U2(4,J)=U2(4,J)+U(1,I,J,K)*U(2,I,J,K)/DBLE(N1*N3*N_INS)
      U2(5,J)=U2(5,J)+U(1,I,J,K)*U(3,I,J,K)/DBLE(N1*N3*N_INS)
      U2(6,J)=U2(6,J)+U(2,I,J,K)*U(3,I,J,K)/DBLE(N1*N3*N_INS)

      ENDDO  ! K,I,J
      ENDDO
      ENDDO

C     TIME CHECK
      INS_PERCENT = 100.0*DBLE(NUM_INS)/(N_INS)
      T_INS_END = SECNDS(T_INS_START)
      T_INS = SECNDS(T_START)

      WRITE(*,*) ''
      WRITE(*,*) '************************************************'
      WRITE(*,*) 'INS FILE ST/EN : ',N_INI,'/',N_END
      WRITE(*,*) 'INS FILE INDEX : ',NTIME
      WRITE(*,*) 'TOTAL STEP     : ',N_INS
      WRITE(*,*) 'STEPS FINISHED : ',NUM_INS,'(',INS_PERCENT,'%)'
      WRITE(*,*) 'TIME TAKEN INS : ',T_INS_END,'SEC'
      WRITE(*,*) 'TIME TAKEN TTL : ',T_INS,'SEC'
      WRITE(*,*) '************************************************'
      WRITE(*,*) ''
   5  CONTINUE  ! NTIME

C----+-------------------------------------------------------------
C     INITIALIZE AVERAGED VALUE
      U3(:,:)=0.0
      U4(:,:)=0.0

C     AVERAGE CALCULATION
      ! U3 : FOR FLUCTUATION PRODUCTION
      ! 1-UU, 2-VV, 3-WW, 4-UV, 5-UW, 6-VW
      DO 200 J=0,N2
      DO L = 1,3
      U3(L,J)=U2(L,J)-U1(L,J)*U1(L,J)
      U4(L,J)=-U3(L,J)
      U3(L,J)=SQRT(U3(L,J))
      ENDDO
      U4(4,J) = U2(4,J)-U1(1,J)*U1(2,J)
      U4(5,J) = U2(5,J)-U1(1,J)*U1(3,J)
      U4(6,J) = U2(6,J)-U1(2,J)*U1(3,J)

      U4(4,J) = -U4(4,J)
      U4(5,J) = -U4(5,J)
      U4(6,J) = -U4(6,J)
 200  CONTINUE

C     CALCULATE COMPUTE AVERAGE PROPERTIES
      DUDY_WALL = (U1(1,1)-U1(1,0))/(Y(2)*0.5)
      N2M=N2-1
      DUDY(:)=0.0
      DUDY(1)=DUDY_WALL
      DO 300 J=2,N2
      DUDY(J)=(U1(1,J)-U1(1,J-1))/(0.5*(Y(J+1)-Y(J-1)))
 300  CONTINUE

      DUDY(N2) = U1(1,N2)/((2.0-Y(N2))*0.5)
      TAU_WALL=DUDY_WALL/RE
      UTAU=SQRT(TAU_WALL)
      RETAU=RE*UTAU
      DEL_V=1.0D0/RETAU


C     EXPORT TO CONSOLE
      WRITE(*,*) '##### FLOW INFORMATION #####'
      WRITE(*,*) 'TAU_WALL = ', TAU_WALL
      WRITE(*,*) 'UTAU = ', UTAU
      WRITE(*,*) 'RETAU = ', RETAU
      WRITE(*,*) 'DELV = ', DEL_V
      WRITE(*,*) ''
      WRITE(*,*) '#####  GRID INFORMATION  #####'
      WRITE(*,*) 'NX = ', NX
      WRITE(*,*) 'NY = ', NY
      WRITE(*,*) 'NZ = ', NZ
      WRITE(*,*) ''

C     EXPORT TO FILE(INFO ABOUT FLOW FIELD)
      FILENAME='meanresult/POST_INFO.txt'
      FILENAME=DIR_SAVE
      NN=INDEX(FILENAME,'_')
      WRITE(UNIT=FILENAME(NN+1:),FMT='(BN,A13)') 'POST_INFO.txt'
      WRITE(*,*) FILENAME
      OPEN(1997, FILE=FILENAME, ACTION='WRITE')
      WRITE(1997,*) '##### FLOW INFORMATION #####'
      WRITE(1997,*) 'TAU_SHEAR = ', TAU_WALL
      WRITE(1997,*) 'UTAU = ', UTAU
      WRITE(1997,*) 'RETAU = ', RETAU
      WRITE(1997,*) 'DELV = ', DEL_V

      WRITE(1997,*) ''
      WRITE(1997,*) '##### GRID INFORMATION #####'
      WRITE(1997,*) 'NX = ', NX
      WRITE(1997,*) 'NY = ', NY
      WRITE(1997,*) 'NZ = ', NZ

C     WRITE TIME-AVERAGED INFORMATION
      FILENAME=DIR_SAVE
      NN=INDEX(FILENAME,'_')
      WRITE(UNIT=FILENAME(NN+1:),FMT='(BN,I6.6)') N_INI
      WRITE(UNIT=FILENAME(NN+7:),FMT='(BN,A1)') '-'
      WRITE(UNIT=FILENAME(NN+8:),FMT='(BN,I6.6)') N_END
      WRITE(UNIT=FILENAME(NN+14:),FMT='(BN,A5)') '_MEAN'     
      WRITE(*,*) FILENAME

      OPEN(11,FILE=FILENAME,ACTION='WRITE')
      WRITE(11,*)'VARIABLES=YC,YP,U,V,W,URMS,VRMS,WRMS,
     >T11,T22,T33,T12,T13,T23'
      
      DO 400 J=0,N2
      YC = 0.5*(Y(J)+Y(J+1))
      YP = YC/DEL_V
      WRITE(11,120)
     &YC,YP,(U1(NV,J),NV=1,3),(U3(NV,J),NV=1,3),(U4(NV,J),NV=1,6)
 400  CONTINUE

      CLOSE(11)

C     NORMALIZE VALUE CALCULATION
      DO 305 J=0,N2
      DO L=1,3
      U1(L,J)  = U1(L,J)/UTAU
      U3(L,J)  = U3(L,J)/UTAU
      U4(L,J)  = U4(L,J)/(UTAU**2)
      U4(L+3,J)= U4(L+3,J)/(UTAU**2)
      ENDDO
 305  CONTINUE

C     WRITE NORMALIZED TIME-AVERAGED INFORMATION
      FILENAME=DIR_SAVE
      NN=INDEX(FILENAME,'_')
      WRITE(UNIT=FILENAME(NN+1:),FMT='(BN,I6.6)') N_INI
      WRITE(UNIT=FILENAME(NN+7:),FMT='(BN,A1)') '-'
      WRITE(UNIT=FILENAME(NN+8:),FMT='(BN,I6.6)') N_END
      WRITE(UNIT=FILENAME(NN+14:),FMT='(BN,A12)') '_MEAN_NORMAL'     
      WRITE(*,*) FILENAME

      OPEN(12,FILE=FILENAME,ACTION='WRITE')
      WRITE(12,*)'VARIABLES=YC,YP,UN,VN,WN,URMSN,VRMSN,WRMSN,
     &T11N,T22N,T33N,T12N,T13N,T23N'
      
      DO 405 J=0,N2
      YC = 0.5*(Y(J)+Y(J+1))
      YP = YC/DEL_V
      WRITE(12,1200)
     &YC,YP,(U1(NV,J),NV=1,3),(U3(NV,J),NV=1,3),(U4(NV,J),NV=1,6)
 405  ENDDO

      CLOSE(12)

      T_END = SECNDS(T_START)

      WRITE(*,*) ''
      WRITE(*,*) '************************************************'
      WRITE(*,*) 'TIME TAKEN : ',T_END,'SEC'
      WRITE(*,*) '************************************************'

120   FORMAT(9(E12.5,X))
1200  FORMAT(13(E12.5,X))
      RETURN
      END
