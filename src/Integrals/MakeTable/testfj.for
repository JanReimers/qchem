C**********************************************************************
C
C  CALCULATE THE ERROR FUNCTION OR ITS Jth DERIVATIVE
C
      FUNCTION FORTFJ(T,J)

      INTEGER J
      REAL*8 T,FORTFJ
      EXTERNAL FJ

      CALL QROMB(FJ,0.0D0,1.0D0,FORTFJ,T,J)
      RETURN
      END
C**********************************************************************
C
C  INTEGRAND OF THE ERROR FUNCTION INTEGRAL
C
	FUNCTION FJ(U,T,J)

        REAL*8  FJ,T,U,U2
	INTEGER J

C        print *,J,T,U

        U2=U**2
	IF (J .EQ. 0) THEN
	  FJ=EXP(-T*U2)
        ELSE
          IF(U2 .EQ. 0) THEN
            U2=0
           ELSE
            FJ=EXP(J*LOG(U2)-T*U2)
          ENDIF
	ENDIF

	RETURN
	END


C**********************************************************************
C
C  ROMBERG INTEGRATION
C
      SUBROUTINE QROMB(FUNC,A,B,SS,T,JJ)

      REAL*8 EPS,A,B,FUNC,SS,DSS,T
      INTEGER JMAX,JMAXP,K,KM,J,JJ
      EXTERNAL FUNC

      PARAMETER (EPS=1.E-13, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)

      REAL*8 S(JMAXP),H(JMAXP)
      H(1)=1.
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,S(J),J,T,JJ)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.0D0,SS,DSS)
          IF (ABS(DSS) .LT. EPS) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
11    CONTINUE
      PRINT *,'Too many steps in QROMB.'
      RETURN
      END

C**********************************************************************
C
C  TRAPAZOIDAL RULE INTEGRATION
C
      SUBROUTINE TRAPZD(FUNC,A,B,S,N,T,JJ)

      REAL*8 FUNC,A,B,S,TNM,DEL,X,SUM,T
      INTEGER N,IT,J,JJ
      EXTERNAL FUNC
      SAVE IT

      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A,T,JJ)+FUNC(B,T,JJ))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X,T,JJ)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
C**********************************************************************
C
C  POLYNOMIAL INTERPOLATION FUNCTION
C
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)

      INTEGER NMAX
      PARAMETER (NMAX=10) 

      REAL*8 X,Y,DY,DIF,DIFT,HO,HP,W,DEN
      INTEGER N,NS,I,M
      REAL*8 XA(N),YA(N),C(NMAX),D(NMAX)

      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
