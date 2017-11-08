      PROGRAM LATLON

      REAL LAT(64,32)
      REAL LON(64,32)
      REAL LTT

      INTEGER I,J

      OPEN(10,FILE="lat.txt",FORM="FORMATTED")

      DO 100 J=1,32
      READ(10,*) LTT
      DO 100 I=1,64
      LAT(I,J)=LTT   
      LON(I,J)=(REAL(I)-1.0)*5.625
100   CONTINUE
      CLOSE(10)
 
      OPEN(20,FILE="latlon.dat",FORM="UNFORMATTED")
      WRITE(20) LAT
      WRITE(20) LON
      CLOSE(20)

      END 
