C****************** DOMAIN SIZE ***********************
      REAL,PARAMETER :: LX = 2.0
      REAL,PARAMETER :: LZ = 2.0

C********************* GRID INFO **********************

      INTEGER,PARAMETER :: NX = 65
      INTEGER,PARAMETER :: NY = 65
      INTEGER,PARAMETER :: NZ = 49
      REAL*8, PARAMETER :: RE = 10500

      INTEGER,PARAMETER :: N1 = NX-1
      INTEGER,PARAMETER :: N2 = NY-1
      INTEGER,PARAMETER :: N3 = NZ-1
     
      INTEGER,PARAMETER :: N_INI=022100
      INTEGER,PARAMETER :: N_END=025000
      INTEGER,PARAMETER :: N_SKIP=100