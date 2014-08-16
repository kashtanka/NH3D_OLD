      SUBROUTINE MATRIXSUM(a,b,c,k)
      implicit none

      integer(4), parameter :: vector_length  = 350    
 
!     MATRIXES: C=A+B
      real(8), dimension (vector_length,2,2)::a,b,c
      integer(4) k,j,i

      do i=1,2
       do j=1,2
        c(k,i,j)=a(k,i,j)+b(k,i,j)
       enddo
      enddo
      
      END SUBROUTINE MATRIXSUM


      SUBROUTINE MATRIXMULT(a,b,c,k)
      implicit none

      integer(4), parameter :: vector_length = 350

!     MATRIXES: C=A*B      

      real(8), dimension (vector_length,2,2)::a,b,c
      integer(4) k

      c(k,1,1)=a(k,1,1)*b(k,1,1)+a(k,1,2)*b(k,2,1)
      c(k,1,2)=a(k,1,1)*b(k,1,2)+a(k,1,2)*b(k,2,2)
      c(k,2,1)=a(k,2,1)*b(k,1,1)+a(k,2,2)*b(k,2,1)
      c(k,2,2)=a(k,2,1)*b(k,1,2)+a(k,2,2)*b(k,2,2)

      END SUBROUTINE MATRIXMULT


      SUBROUTINE MATRIXMULTVECTOR(a,g,f,k)
      implicit none

      integer(4), parameter :: vector_length = 350      

!     MATRIX A, VECTORS g, f: Ag=f

      real(8) a(vector_length,2,2),f(vector_length,2),
     & g(vector_length,2)
      integer(4) k

      f(k,1)=a(k,1,1)*g(k,1)+a(k,1,2)*g(k,2)
      f(k,2)=a(k,2,1)*g(k,1)+a(k,2,2)*g(k,2)

      return
      END SUBROUTINE MATRIXMULTVECTOR


      SUBROUTINE VECTORSUM(a,b,c,k)
      implicit none

      integer(4), parameter :: vector_length = 350     

!     VECTORS: C=A+B

      real(8), dimension(vector_length,2)::a,b,c
      integer(4) k
      c(k,1)=a(k,1)+b(k,1)
      c(k,2)=a(k,2)+b(k,2)
      END SUBROUTINE VECTORSUM


      SUBROUTINE INVERSMATRIX(a,a1,k)
      implicit none 

      integer(4), parameter :: vector_length = 350
      
!     MATRIXES: A1=A*(-1)

      real(8), dimension(vector_length,2,2)::a,a1
      integer(4) k
      
      a1(k,1,1)=a(k,2,2)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1)) 
      a1(k,1,2)=-a(k,1,2)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1))
      a1(k,2,1)=-a(k,2,1)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1))
      a1(k,2,2)=a(k,1,1)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1))  

      END SUBROUTINE INVERSMATRIX

      
      SUBROUTINE MATRIXPROGONKA(a,b,c,d,y,N)
      
!     MATRIXPROGONKA solves the set of MATRIX three-point diference equations 
      
      implicit none

      integer(4), parameter :: vector_length = 350

      real(8), dimension(vector_length,2,2):: a,b,c,x3,x32,x31,alpha
      real(8), dimension(vector_length,2):: y,d,x2,beta,x21
      integer(4) N,i,j,k

      call INVERSMATRIX(c,x32,1)
      call MATRIXMULT(x32,b,x3,1)
      do i=1,2
       do j=1,2
        alpha(2,i,j)=x3(1,i,j)
       enddo
      enddo
      call MATRIXMULTVECTOR(x32,d,x2,1)
      do i=1,2
       beta(2,i)=x2(1,i)
      enddo

      do k=3,N
       CALL MATRIXMULT(A,ALPHA,X3,k-1)
       CALL MATRIXSUM(C,-X3,X31,k-1)
       CALL INVERSMATRIX(X31,X32,k-1)
       CALL MATRIXMULT(X32,B,X3,k-1)
       do i=1,2
        do j=1,2
         alpha(k,i,j)=X3(k-1,i,j)
        enddo
       enddo
       !call matrixmult(x3,x31,x33,k-1)
       !call matrixsum(-b,x33,x3,k-1)
       CALL MATRIXMULTVECTOR(A,BETA,X2,K-1)
       CALL VECTORSUM(D,X2,X21,K-1)
       CALL MATRIXMULTVECTOR(X32,X21,X2,K-1)
       do i=1,2
        beta(k,i)=X2(k-1,i)
       enddo
      enddo

      CALL MATRIXMULT(A,ALPHA,X3,N)
      CALL MATRIXSUM(C,-X3,X31,N)
      CALL INVERSMATRIX(X31,X32,N)
      CALL MATRIXMULTVECTOR(A,BETA,X2,N)
      CALL VECTORSUM(D,X2,X21,N)
      CALL MATRIXMULTVECTOR(X32,X21,X2,N)
      do i=1,2
       Y(N,i)=X2(N,i)
      enddo

      do k=N-1,1,-1
       CALL MATRIXMULTVECTOR(ALPHA,Y,X2,K+1)
       CALL VECTORSUM(X2,BETA,X21,K+1)
       Y(K,1)=X21(K+1,1)
       Y(K,2)=X21(K+1,2)
      enddo

      return

      END SUBROUTINE MATRIXPROGONKA


      REAL(8) FUNCTION KRON(i,j)
      implicit none
      integer(4) i,j
      kron=0.
      if(i==j) kron=1.
      END FUNCTION 


      SUBROUTINE IND_STAB_FACT_DB (a,b,c,N,M,ind_stab,ind_bound)

      implicit none

      integer(4), parameter :: vector_length = 350

      real(8), dimension(1:vector_length):: a,b,c
      integer(4) M,i,N
      logical ind_stab, ind_bound 

      SAVE

      ind_stab=.true.
      if (ind_bound .eqv. .true.) then 
       if (dabs(b(N))>=dabs(c(N)).or.dabs(a(M))>=dabs(c(M))) then
        ind_stab=.false.
        RETURN
       endif
      endif
      do i=N+1,M-1
       if (dabs(a(i))+dabs(b(i))>=dabs(c(i))) then
        ind_stab=.false.
        RETURN
       endif
      enddo
      
      END SUBROUTINE IND_STAB_FACT_DB



      LOGICAL FUNCTION CHECK_PROGONKA(N,a,b,c,d,y)

!     Function CHECK_PROGONKA checks the accuracy
!     of tridiagonal matrix system solution

      implicit none

      integer(4), intent(in):: N
      real(8), intent(in):: a(1:N)
      real(8), intent(in):: b(1:N)
      real(8), intent(in):: c(1:N)
      real(8), intent(in):: d(1:N)
      real(8), intent(in):: y(1:N)

      real(8), parameter:: del0 = 1.0d-13
      real(8) del

      integer(4) i

      del = 0.d0
      del = dmax1(c(1)*y(1)-b(1)*y(2)-d(1),del)
      del = dmax1(c(N)*y(N)-a(N)*y(N-1)-d(N),del)
      do i = 2, N-1
        del = dmax1(-a(i)*y(i-1)+c(i)*y(i)-b(i)*y(i+1)-d(i),del)
      enddo

      CHECK_PROGONKA = del < del0

      END FUNCTION CHECK_PROGONKA



      SUBROUTINE progonka(a, b, c, f, y, K, N)
      implicit none
!FACTORIZATION METHOD FOR THE FOLLOWING SYSTEM OF LINEAR EQUATIONS:
!     -a(i)*y(i-1)+c(i)*y(i)-b(i)*y(i+1)=f(i) i=K+1,N-1
!      c(K)*y(K)-b(K)*y(K+1)=f(K)
!     -a(N)*y(N-1)+c(N)*y(N)=f(N)
!
      integer(4), parameter :: M=350

      integer(4) K, N, i
      real(8) a(M), b(M), c(M), f(M), y(M)
      real(8) alpha(M+2), beta(M+2) 
      SAVE
                  
      alpha(K+1) = b(K)/c(K)
      beta(K+1) = f(K)/c(K)
      do i = K+2, N+1 
       alpha(i) = b(i-1)/(c(i-1)-a(i-1)*alpha(i-1))
       beta(i) = (f(i-1)+a(i-1)*beta(i-1))/ 
     & (c(i-1)-a(i-1)*alpha(i-1))
      end do
      y(N) = beta(N+1)
      do i = N-1, K, -1
       y(i) = alpha(i+1)*y(i+1)+beta(i+1)
      end do
       
      END 


!     The subroutine AEH1D is from SRCC MSU numerical analysis library, file AEH1D.FOR

      SUBROUTINE AEH1D(N,A,EV,V,RAB1,IERR)
      DIMENSION A(N,N),EV(N),V(N,N),RAB1(N)
      INTEGER N,IERR,B,C,D,E,F,G,M,O,P,Q,R,S,T,U
      DOUBLE PRECISION A,EV,V,RAB1,H,I,J,K,L,DSQRT,DABS,
     1DSIGN,W,X,Y,Z,BA,BB,BC,BD,SYS051
      DATA  SYS051/1.1107652D-16/
      DO 1 B=1,N
      DO 1 C=1,B
      V(B,C)=A(B,C)
    1 CONTINUE
      IF(N.EQ.1) GO TO 13
      DO 12 F=2,N
      B=N+2-F
      E=B-1
      J=0.0D0
      L=0.0D0
      IF(E.LT.2) GO TO 3
      DO 2 D=1,E
    2 L=L+DABS(V(B,D))
      IF(L.NE.0.0D0) GO TO 4
    3 RAB1(B)=V(B,E)
      GO TO 11
    4 DO 5 D=1,E
      V(B,D)=V(B,D)/L
      J=J+V(B,D)*V(B,D)
    5 CONTINUE
      H=V(B,E)
      I=-DSIGN(DSQRT(J),H)
      RAB1(B)=L*I
      J=J-H*I
      V(B,E)=H-I
      H=0.0D0
      DO 9 C=1,E
      V(C,B)=V(B,C)/J
      I=0.0D0
      DO 6 D=1,C
    6 I=I+V(C,D)*V(B,D)
      G=C+1
      IF(E.LT.G) GO TO 8
      DO 7 D=G,E
    7 I=I+V(D,C)*V(B,D)
    8 RAB1(C)=I/J
      H=H+RAB1(C)*V(B,C)
    9 CONTINUE
      K=H/(J+J)
      DO 10 C=1,E
      H=V(B,C)
      I=RAB1(C)-K*H
      RAB1(C)=I
      DO 10 D=1,C
      V(C,D)=V(C,D)-H*RAB1(D)-I*V(B,D)
   10 CONTINUE
   11 EV(B)=J
   12 CONTINUE
   13 EV(1)=0.0D0
      RAB1(1)=0.D0
      DO 18 B=1,N
      E=B-1
      IF(EV(B).EQ.0.0D0) GO TO 16
      DO 15 C=1,E
      I=0.0D0
      DO 14 D=1,E
   14 I=I+V(B,D)*V(D,C)
      DO 15 D=1,E
      V(D,C)=V(D,C)-I*V(D,B)
   15 CONTINUE
   16 EV(B)=V(B,B)
      V(B,B)=1.0D0
      IF(E.LT.1) GO TO 18
      DO 17 C=1,E
      V(B,C)=0.0D0
      V(C,B)=0.0D0
   17 CONTINUE
   18 CONTINUE
      IERR=0
      IF(N.EQ.1) GO TO 34
      DO 19 M=2,N
   19 RAB1(M-1)=RAB1(M)
      Y=0.0D0
      W=0.0D0
      RAB1(N)=0.0D0
      DO 29 Q=1,N
      O=0
      BA=SYS051*(DABS(EV(Q))+DABS(RAB1(Q)))
      IF(W.LT.BA) W=BA
      DO 20 R=Q,N
      IF(DABS(RAB1(R)).LE.W) GO TO 21
   20 CONTINUE
   21 IF(R.EQ.Q) GO TO 28
   22 IF(O.EQ.30) GO TO 33
      O=O+1
      T=Q+1
      Z=EV(Q)
      BB=(EV(T)-Z)/(2.0D0*RAB1(Q))
      IF(DABS(BB).GT.1.D0) BC=DABS(BB)*DSQRT(1.D0+(1/BB)**2)
      IF(DABS(BB).LE.1.D0) BC=DSQRT(BB*BB+1.D0)
      EV(Q)=RAB1(Q)/(BB+DSIGN(BC,BB))
      BA=Z-EV(Q)
      DO 23 M=T,N
   23 EV(M)=EV(M)-BA
      Y=Y+BA
      BB=EV(R)
      X=1.0D0
      BD=0.0D0
      U=R-Q
      DO 27 S=1,U
      M=R-S
      Z=X*RAB1(M)
      BA=X*BB
      IF(DABS(BB).LT.DABS(RAB1(M))) GO TO 24
      X=RAB1(M)/BB
      BC=DSQRT(X*X+1.0D0)
      RAB1(M+1)=BD*BB*BC
      BD=X/BC
      X=1.0D0/BC
      GO TO 25
   24 X=BB/RAB1(M)
      BC=DSQRT(X*X+1.0D0)
      RAB1(M+1)=BD*RAB1(M)*BC
      BD=1.0D0/BC
      X=X*BD
   25 BB=X*EV(M)-BD*Z
      EV(M+1)=BA+BD*(X*Z+BD*EV(M))
      DO 26 P=1,N
      BA=V(P,M+1)
      V(P,M+1)=BD*V(P,M)+X*BA
      V(P,M)=X*V(P,M)-BD*BA
   26 CONTINUE
   27 CONTINUE
      RAB1(Q)=BD*BB
      EV(Q)=X*BB
      IF(DABS(RAB1(Q)).GT.W) GO TO 22
   28 EV(Q)=EV(Q)+Y
   29 CONTINUE
      DO 32 S=2,N
      M=S-1
      P=M
      BB=EV(M)
      DO 30 O=S,N
      IF(EV(O).GE.BB) GO TO 30
      P=O
      BB=EV(O)
   30 CONTINUE
      IF(P.EQ.M) GO TO 32
      EV(P)=EV(M)
      EV(M)=BB
      DO 31 O=1,N
      BB=V(O,M)
      V(O,M)=V(O,P)
      V(O,P)=BB
   31 CONTINUE
   32 CONTINUE
      GO TO 34
   33 IERR=Q
   34 IF(IERR.NE.0) CALL UTAE10(IERR,N,13)
      RETURN
      END

!     The subroutine ASG0D is from SRCC MSU numerical analysis library, file ASG0D.FOR

      SUBROUTINE ASG0D(A,B,X,N,P)
      DIMENSION A(N,N),B(N),X(N)
      INTEGER N,P,F,L,M,I,J,K
      DOUBLE PRECISION A,B,X,S,R,W
      IF(N.GT.1) GO TO 15
      IF(N.LE.0) GO TO 10
      X(N)=B(N)/A(N,N)
      GO TO 10
   15 IF(P.NE.1) GO TO 4
C
      F=0
    1 F=F+1
      S=1.D0/A(F,F)
      L=F+1
      DO 2 I=L,N,1
    2 A(I,F)=-(S*A(I,F))
      DO 3 J=L,N,1
      DO 3 I=L,N,1
    3 A(I,J)=A(I,J)+A(F,J)*A(I,F)
      IF(L-N) 1,4,4
C
    4 DO 5 I=1,N
    5 X(I)=B(I)
      M=N-1
      DO 6 J=1,M
      W=X(J)
      K=J+1
      DO 6 I=K,N
    6 X(I)=X(I)+A(I,J)*W
C
      X(N)=X(N)/A(N,N)
      IF(N.EQ.1) GO TO 10
      K=N
      DO 9 I=1,M
      R=0.0D0
      L=K-1
      DO 8 J=K,N
    8 R=R+A(L,J)*X(J)
      K=N-I
    9 X(K)=(X(K)-R)/A(K,K)
   10 RETURN
      END

!     The subroutine UTAE10 is from SRCC MSU numerical analysis library, file UTAE10.FOR

      SUBROUTINE UTAE10(IERR,N,M)
      IF(M.EQ.10) GO TO 100
      IF(M.EQ.11) GO TO 110
      IF(M.EQ.12) GO TO 120
      IF(M.EQ.13) GO TO 130
      IF(M.EQ.14) GO TO 140
      IF(M.EQ.15) GO TO 150
      IF(M.EQ.18) GO TO 180
      IF(M.EQ.19) GO TO 190
      IF(M.EQ.20) GO TO 200
      IF(M.EQ.21) GO TO 210
      IF(M.EQ.22) GO TO 220
      IF(M.EQ.23) GO TO 230
      IF(M.EQ.25) GO TO 250
      IF(M.EQ.30) GO TO 300
      IF(M.EQ.31) GO TO 310
      IF(M.EQ.32) GO TO 320
      IF(M.EQ.33) GO TO 330
      IF(M.EQ.34) GO TO 340
      IF(M.EQ.35) GO TO 350
      RETURN
  100 PRINT 101
  101 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEG1C(AEG1P):',
     *'�ATA��HA� O���KA')
      GO TO 1
  110 PRINT 111
  111 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEH1C(AEH1P):',
     *'�ATA��HA� O���KA')
      GO TO 2
  120 PRINT 121
  121 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEG1R(AEG1D):',
     *'�ATA��HA� O���KA')
      GO TO 1
  130 PRINT 131
  131 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEH1R(AEH1D):',
     *'�ATA��HA� O���KA')
      GO TO 2
  140 PRINT 141
  141 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEE1R(AEE1D):',
     *'�ATA��HA� O���KA')
      GO TO 2
  150 PRINT 151
  151 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEJ1R(AEJ1D):',
     *'�ATA��HA� O���KA')
      GO TO 3
  300 PRINT 301
  301 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEG2C(AEG2P):',
     *'�ATA��HA� O���KA')
      GO TO 6
  310 PRINT 311
  311 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEH2C(AEH2P):',
     *'�ATA��HA� O���KA')
      GO TO 7
  320 PRINT 321
  321 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEG2R(AEG2D):',
     *'�ATA��HA� O���KA')
      GO TO 6
  330 PRINT 331
  331 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEH2R(AEH2D):',
     *'�ATA��HA� O���KA')
      GO TO 7
  340 PRINT 341
  341 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEE2R(AEE2D):',
     *'�ATA��HA� O���KA')
      GO TO 7
    1 PRINT 11
   11 FORMAT(' N1 - ��� B���C�EH�� CO�CTBEHHO�O �HA�EH�� C �H�EKCOM ',
     *'IERR TPE�YETC� �O�EE'/'     30 �TEPA���. �PAB���HO B���C�EH� ',
     *'CO�CTBEHH�E �HA�EH�� C �H�EKCAM�'/'     IERR+1,IERR+2,...,N,',
     *' HO CO�CTBEHH�E BEKTOP� HE B���C���TC�.')
      RETURN
    2 PRINT 12
   12 FORMAT(' N1 - ��� B���C�EH�� CO�CTBEHHO�O �HA�EH�� C �H�EKCOM ',
     *'IERR TPE�YETC� �O�EE'/'     30 �TEPA���. �PAB���HO B���C�EH� ',
     *'CO�CTBEHH�E �HA�EH�� � CO�CTBEHH�E'/'     BEKTOP� C �H�EKCAM� ',
     *'1,2,...,IERR-1, HO CO�CTBEHH�E �HA�EH��'/'     HEY�OP��O�EH�.')
      RETURN
    3 IF(IERR.GT.N.AND.IERR.LE.2*N) GO TO 4
      IF(IERR.GT.2*N) GO TO 5
      IF(IERR.LE.N) GO TO 2
    4 PRINT 14
   14 FORMAT(' N1 - HE BCE �O�APH�E �PO��BE�EH�� COOTBETCBY���X ',
     *'��EMEHTOB �O�O�H�X'/'     ��A�OHA�E� HEOTP��ATE��H�.')
      RETURN
    5 PRINT 15
   15 FORMAT(' N1 - �MEETC� PABHOE HY�� �PO��BE�EH�E COOTBETCTBY���X ',
     *'��EMEHTOB �O�O�H�X'/'     ��A�OHA�E�, �P��EM COMHO��TE�� ',
     *'PABH� HY�� HEO�HOBPEMEHHO.B �TOM'/'     C�Y�AE HE CY�ECTBYET ',
     *'C�MMETP��A���, HEO�XO��MO� ��� �PAB���HO�O'/'     B���C�EH�� ',
     *'CO�CTBEHH�X BEKTOPOB.')
      RETURN
    6 PRINT 16
   16 FORMAT(' N1 - ��� B���C�EH�� CO�CTBEHHO�O �HA�EH�� C �H�EKCOM ',
     *'IERR TPE�YETC� �O�EE'/'     30 �TEPA���. �PAB���HO B���C�EH� ',
     *'CO�CTBEHH�E �HA�EH�� C �H�EKCAM�'/'     IERR+1,IERR+2,...,N.')
      RETURN
    7 PRINT 17
   17 FORMAT(' N1 - ��� B���C�EH�� CO�CTBEHHO�O �HA�EH�� C �H�EKCOM ',
     *'IERR TPE�YETC� �O�EE'/'     30 �TEPA���. �PAB���HO B���C�EH� � ',
     *'Y�OP��O�EH� CO�CTBEHH�E �HA�EH��'/'    C �H�EKCAM� 1,2,...,'
     *'IERR-1. OH� HE O���ATE��HO �B���TC� CAM�M�'/'     MEH���M� �� ',
     *'���� N CO�CTBEHH�X �HA�EH��.')
      RETURN
  350 IF(IERR.LT.0) GO TO 351
      GO TO 353
  351 PRINT 352
  352 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEJ2R(AEJ2D):',
     *'HE�ATA��HA� O���KA.'/' N1 - �MEETC� O�HO ��� HECKO��KO PABH�X ',
     *'HY�� �PO��BE�EH�� COOTBETCTBY���X'/'     ��EMEHTOB �O�O�H�X ',
     *'��A�OHA�E�, �P��EM COMHO��TE�� PABH� HY��',
     */'     ��������������. ',
     *'B �TOM C�Y�AE BCE CO�CTBEHH�E �HA�EH�� ',
     *'B���C�EH� �PAB���HO.')
      RETURN
  353 PRINT 354
  354 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEJ2R(AEJ2D):',
     *'�ATA��HA� O���KA')
      IF(IERR.GT.N)GOTO 4
      IF(IERR.GT.0.AND.IERR.LE.N) GO TO 7
  180 PRINT 181
  181 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEH3C(AEH3P):',
     *'�ATA��HA� O���KA')
      IF(IERR.GT.3*N) GO TO 8
    9 PRINT 19
   19 FORMAT(' N1 - ��� B���C�EH�� CO�CTBEHHO�O BEKTOPA ',
     *'C �H�EKCOM IERR TPE�YETC� �O�EE'/'     5 �TEPA���,',
     *' �P� �TOM KOM�OHEHT� �TO�O BEKTOPA'/'     �P�PABH�BA�TC� HY��.'
     *,
     *'EC�� TAK�X CO�CTBEHH�X BEKTOPOB'/'      HECKO��KO, TO ',
     *'�HA�EH�E IERR �O�A�AETC� PABH�M �H�EKCY'/'     �OC�E�HE�O ',
     *'�� ���.')
      RETURN
    8 PRINT 18
   18 FORMAT(/'N1 - �HA�EH�E MM MEH��E �CT�HHO�O ��C�A CO�CTBEHH�X ',
     *'�HA�EH�� M HA �HTEPBA�E.')
      RETURN
  190 PRINT 191
  191 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEH4C(AEH4P):',
     *'�ATA��HA� O���KA')
      GO TO 8
  200 PRINT 201
  201 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEH3R(AEH3D):',
     *'�ATA��HA� O���KA')
      IF(IERR.GT.3*N) GO TO 8
      GO TO 9
  210 PRINT 211
  211 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEH4R(AEH4D):',
     *'�ATA��HA� O���KA')
      GO TO 8
  220 PRINT 221
  221 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEE9R(AEE9D):',
     *'�ATA��HA� O���KA')
      IF(IERR.GT.3*N) GO TO 8
      GO TO 9
  230 PRINT 231
  231 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEJ3R(AEJ3D):',
     *'�ATA��HA� O���KA')
      IF(IERR.LT.0) GO TO 10
      IF(IERR.GT.2*N.AND.IERR.LE.3*N) GO TO 10
      IF(IERR.GT.3*N) GO TO 8
      IF(IERR.LE.N) GO TO 9
      GO TO 1111
   10 PRINT 20
   20 FORMAT(' N1 - �A�AHHA� MATP��A HE �B��ETC� MATP��E� �KO��, ',
     *'T.K. �MEETC� PABHOE'/'     HY�� �PO��BE�EH�E ',
     *'COOTBETCTBY���X ��EMEHTOB �O�O�H�X'/'     ��A�OHA�E�, �P��EM ',
     *'COMHO��TE�� PABH� HY�� HEO�HOBPEMEHHO.'/'     B �TOM C�Y�AE ',
     *'HET C�MMETP��A���, HEO�XO��MO� ��� �PAB���HO�O'/
     *'     ���������� CO�CTBEHH�X BEKTOPOB.')
      RETURN
 1111 PRINT 21
   21 FORMAT(' N1 - �A�AHHA� MATP��A HE  �B��ETC� MATP��E� �KO��, ',
     *'TAK KAK HE BCE'/'     �O�APH�E �PO��BE�EH�� COOTBETCTBY���X ',
     *'��EMEHTOB �O�O�H�X ��A�OHA�E�'/'     HEOTP��ATE��H�.')
      RETURN
  250 IF(IERR.LT.0) GO TO 251
      GO TO 254
  251 PRINT 252
  252 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEJ4R(AEJ4D):',
     *'HE�ATA��HA� O���KA')
      PRINT 253
  253 FORMAT(' N1 - �A�AHHA� MATP��A HE �B��ETC� MATP��E� �KO��, ',
     *'TAK KAK �MEETC� PABHOE'/'     HY�� �PO��BE�EH�E ',
     *'COOTBETCTBY���X ��EMEHTOB �O�O�H�X'/'    ��A�OHA�E�, �P��EM ',
     *'COMHO��TE�� PABH� HY�� HEO�HOBPEMEHHO.'/'     B �TOM ',
     *'C�Y�AE CO�CTBEHH�E �HA�EH�� B���C�EH� �PAB���HO.')
      RETURN
  254 PRINT 255
  255 FORMAT(' �����OTEKA H�B� M�Y, �O��PO�PAMMA AEJ4R(AEJ4D):',
     *'�ATA��HA� O���KA')
      IF(IERR.GT.N.AND.IERR.LE.2*N) GO TO 1111
      IF(IERR.GT.3*N) GO TO 8
      RETURN
      END


      subroutine BICGSTAB(n, matrix, rhs, solution, aquracy)
      
c     Bicojugate Gradient Stabilized solver
c     Solves dense linear system : matrix * solution = rhs
c
c     where :
c     n         - the matrix & vectors dimensions
c     matrix    - the dense n x n matrix
c     rhs       - the right hand side vector
c     solution  - the vector linear system solution to place in
c                 (may contain an initial guess)
c     aquracy   - the solver aquracy
      
      implicit none
      
      integer :: n
      real(8) :: matrix(1:n, 1:n), rhs(1:n), solution(1:n), aquracy
      
c     Internal variables      
      real(8) :: residual(1:n)
      real(8) :: rtilde(1:n), p(1:n), v(1:n), s(1:n), t(1:n)
      real(8) :: alpha, beta, omega, rho_1, rho_2
      integer :: iterationsNumber
      
c     BODY

c     Setup the residual vector.
      residual = rhs - matmul(matrix, solution)

c     Setup rtilde, p and v.
      rtilde = residual
      p = residual
      v = residual

c     Initialize parameters.
      alpha = 0; omega = 0; rho_2 = 0; 

c     Repeat until converged with the specified accuracy.
      iterationsNumber = 0
      do while (maxval(rhs - matmul(matrix, solution)) .GE. aquracy)
        print*, 'residual', maxval(rhs - matmul(matrix, solution))
        iterationsNumber = iterationsNumber + 1

        rho_1 = dot_product(rtilde, residual)

        if (iterationsNumber .NE. 1) then

            beta = (rho_1 / rho_2) * (alpha / omega)
            
            p = residual + beta * (p - omega * v)
            
        endif
            
        v = matmul(matrix, p)

          alpha = rho_1 / dot_product(rtilde, v)

          s = residual - alpha * v

        t = matmul(matrix, s)

          omega = dot_product(t, s) / dot_product(t, t)

          solution = solution + alpha * p + omega * s
          residual = s - omega * t

          rho_2 = rho_1
          
      enddo

      end subroutine BICGSTAB
