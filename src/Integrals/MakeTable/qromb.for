C**********************************************************************
C
C  ROMBERG INTEGRATION
C
      SUBROUTINE QROMB(FUNC,A,B,SS,I,IA,ALPHA)

      REAL*8 EPS,A,B,FUNC,SS,DSS,ALPHA
      INTEGER JMAX,JMAXP,K,KM,J,I,IA
      EXTERNAL FUNC

      PARAMETER (EPS=1.E-7, JMAX=5, JMAXP=JMAX+1, K=5, KM=K-1)

      REAL*8 S(JMAXP),H(JMAXP)
      H(1)=1.
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,S(J),J,I,IA,ALPHA)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.0D0,SS,DSS)
C          IF (ABS(DSS) .LT. EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
11    CONTINUE
C      PRINT *,'Too many steps in QROMB.'
C      PRINT *,'ALPHA',ALPHA
C      PRINT *,'RMAX',B
C       PRINT *,'QROMB',ABS(DSS)
      RETURN
      END
C**********************************************************************
C
C  TRAPAZOIDAL RULE INTEGRATION
C
      SUBROUTINE TRAPZD(FUNC,A,B,S,N,I,IA,ALPHA)

      REAL*8 FUNC,A,B,S,ALPHA,TNM,DEL,X,SUM
      INTEGER N,IT,J,I,IA
      EXTERNAL FUNC

      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A,I,IA,ALPHA)+FUNC(B,I,IA,ALPHA))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X,I,IA,ALPHA)
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
