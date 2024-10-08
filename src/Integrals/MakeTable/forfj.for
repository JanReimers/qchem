C**********************************************************************
C
C  CALCULATE THE ERROR FUNCTION OR ITS Jth DERIVATIVE
C
	PROGRAM TESTFJ

	INTEGER JMAX
	PARAMETER (JMAX=9)
	REAL*8  FJ,FJ1(0:JMAX),FJ2(0:JMAX),T,DF
	INTEGER J,IT

	EXTERNAL FJ

C	CALL MKFJTB(JMAX)
	J=0
	DO IT=1,1000
	  T=RND()*60.0D0
	  DO J=0,JMAX
	    CALL QROMB1(FJ,0.0D0,1.0D0,FJ1(J),T,J)
	  ENDDO
	  CALL CALCFJ(T,FJ2)
	  DF=0.0D0
	  DO J=0,JMAX
	    DF=DF+(FJ1(J)-FJ2(J))**2
	  ENDDO
	  DF=SQRT(DF)
	  PRINT *,'T,DF',T,DF
	ENDDO

	STOP
	END
