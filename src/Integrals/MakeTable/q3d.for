      SUBROUTINE QROMBX(FUNC,A,B,SS)

      REAL*8 EPS,DSS,FUNC,A,B,SS
      INTEGER JMAX,JMAXP,K,KM,J
      PARAMETER (EPS=1.E-5, JMAX=6, JMAXP=JMAX+1, K=5, KM=K-1)
      REAL*8 S(JMAXP),H(JMAXP)
      EXTERNAL FUNC

      IF (A .EQ. B) THEN
	SS=0.0D0
        RETURN
      ENDIF
      H(1)=1.
      DO 11 J=1,JMAX
        CALL TRAPZX(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.0D0,SS,DSS)
          IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
C	PRINT *,'QROMBX',J
11    CONTINUE
C      PAUSE 'Too many steps.'
      END

      SUBROUTINE TRAPZX(FUNC,A,B,S,N)

      REAL*8 FUNC,A,B,S,X,SUM,DEL,TNM
      INTEGER N,IT,J
      EXTERNAL FUNC

      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A)+FUNC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END

      SUBROUTINE QROMBY(FUNC,A,B,SS)

      REAL*8 EPS,DSS,FUNC,A,B,SS
      INTEGER JMAX,JMAXP,K,KM,J
      PARAMETER (EPS=1.E-5, JMAX=6, JMAXP=JMAX+1, K=5, KM=K-1)
      REAL*8 S(JMAXP),H(JMAXP)
      EXTERNAL FUNC

      IF (A .EQ. B) THEN
	SS=0.0D0
        RETURN
      ENDIF
      H(1)=1.
      DO 11 J=1,JMAX
        CALL TRAPZY(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.0D0,SS,DSS)
          IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
C	PRINT *,'   QROMBY',J
11    CONTINUE
C      PAUSE 'Too many steps.'
      END

      SUBROUTINE TRAPZY(FUNC,A,B,S,N)

      REAL*8 FUNC,A,B,S,X,SUM,DEL,TNM
      INTEGER N,IT,J
      EXTERNAL FUNC

      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A)+FUNC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END

      SUBROUTINE QROMBZ(FUNC,A,B,SS)

      REAL*8 EPS,DSS,FUNC,A,B,SS
      INTEGER JMAX,JMAXP,K,KM,J
      PARAMETER (EPS=1.E-5, JMAX=6, JMAXP=JMAX+1, K=5, KM=K-1)
      REAL*8 S(JMAXP),H(JMAXP)
      EXTERNAL FUNC

      IF (A .EQ. B) THEN
	SS=0.0D0
        RETURN
      ENDIF
      H(1)=1.
      DO 11 J=1,JMAX
        CALL TRAPZZ(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.0D0,SS,DSS)
          IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
11    CONTINUE
C      PAUSE 'Too many steps.'
      END

      SUBROUTINE TRAPZZ(FUNC,A,B,S,N)

      REAL*8 FUNC,A,B,S,X,SUM,DEL,TNM
      INTEGER N,IT,J
      EXTERNAL FUNC

      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A)+FUNC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END

C      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
C
C     INTEGER NMAX,N,NS,I,M
C      PARAMETER (NMAX=10) 
C      REAL*8 XA(N),YA(N),C(NMAX),D(NMAX)
C      REAL*8 X,Y,DY,DIFT,DIF,W,HO,HP,DEN
C
C      NS=1
C      DIF=ABS(X-XA(1))
C      DO 11 I=1,N 
C        DIFT=ABS(X-XA(I))
C        IF (DIFT.LT.DIF) THEN
C          NS=I
C          DIF=DIFT
C        ENDIF
C        C(I)=YA(I)
C        D(I)=YA(I)
C11    CONTINUE
C      Y=YA(NS)
C      NS=NS-1
C      DO 13 M=1,N-1
C        DO 12 I=1,N-M
C          HO=XA(I)-X
C          HP=XA(I+M)-X
C          W=C(I+1)-D(I)
C          DEN=HO-HP
C          IF(DEN.EQ.0.)PAUSE
C          DEN=W/DEN
C          D(I)=HP*DEN
C          C(I)=HO*DEN
C12      CONTINUE
C        IF (2*NS.LT.N-M)THEN
C          DY=C(NS+1)
C        ELSE
C          DY=D(NS)
C          NS=NS-1
C        ENDIF
C        Y=Y+DY
C13    CONTINUE
C      RETURN
C      END
