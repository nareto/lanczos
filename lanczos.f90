PROGRAM lanczos
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp=KIND(0.d0)
  INTEGER :: k,j, n, left_index, right_index, iterations, uscito = 0
  REAL(dp) :: epsilon =1.e0
  REAL(dp), DIMENSION(:), allocatable :: rnd, t, alfa, beta, r, &
       v, w, z, eig, alfa_tmp, beta_tmp

  WRITE(*,*) "Assegna la dimensione n della matrice A:"
  READ(*,*) n
  WRITE(*,*) "Assegna il numero di iterazioni K da effettuare"
  READ(*,*) iterations
  iterations = min(iterations,n)

  ALLOCATE(rnd(n),t(n),alfa(iterations),beta(iterations),r(n),&
       v(n),w(n),z(n),eig(iterations),alfa_tmp(n),beta_tmp(n))

  OPEN(unit=4, file="errors.txt")
  OPEN(unit=5, file="eigv.txt")
  OPEN(unit=6, file="correcteigv.txt")
  WRITE(6,*) "#Correct eigenvalues:"
  DO j=1,n
     WRITE(6,*) correctegv(j,n)
  END DO
  WRITE(5,*) "#Approximated eigenvalues"

  CALL random_seed
  CALL random_number(rnd)
  v = rnd/NORM2(rnd)
  w = 0
  DO k=1, iterations
     CALL prodotto(v, n, t)
     alfa(k) = dot_product(v, t)
     r = t - alfa(k)*v - beta(k-1)*w
     IF (k < n) then
        beta(k) = norm2(r)
        IF (abs(beta(k)) <= epsilon) then
           WRITE(0,*) "|beta(k)| <= epsilon =", epsilon, " per k = ", k
           uscito = 1
           exit
        END IF
        z = r/beta(k)
     END IF
     w = v
     v = z

     !la seguente subroutine scrive in eig i k autovalori di T_k
     !siccome DSTEVR, la routine di LAPACK utilizzata, potrebbe moltiplicare
     !i vettori diagonali e sopradiagonali per evitare instabilità numerica,
     !gliene passo delle copie per lasciare i vettori alfa e beta inalterati
     alfa_tmp = alfa
     beta_tmp = beta
     CALL eigenvalues(alfa_tmp(1:k), beta_tmp(1:k), eig, k)

     WRITE(5,"(I4 A)", advance = "no") k, " " 
     WRITE(4,"(I4 A)", advance = "no") k, " " 
     DO j=1, k
        IF (mod(j,2)==1) then
           right_index = k - ceiling(REAL(j)/2.0) + 1
           WRITE(4, "(F12.9 A)", advance="no") abs(correctegv(right_index + n -k,n)&
                - eig(right_index)), " "
           WRITE(5,"(F12.9 A)", advance="no") eig(right_index), " "
        else
           left_index = j/2
           WRITE(4, "(F12.9 A)", advance="no") abs(correctegv(left_index,n)&
                - eig(left_index)), " "
           WRITE(5, "(F12.9 A)", advance="no") eig(left_index), " "
        END IF
     END DO
     WRITE(5,*)
     WRITE(4,*)
  END DO
  
  IF (uscito == 1) then
     k = k + 1
  END IF

  WRITE(5,*) "#iterations: ", k - 1
  WRITE(4,*) "#iterations: ", k - 1
  WRITE(0,"(A I6 A)") "ho effettuato", k-1, " iterazioni, esco"

CONTAINS
  
  REAL FUNCTION correctegv(i,n)
    IMPLICIT NONE
    INTEGER :: i,n
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    REAL(dp) :: pi = acos(-1.0)
    correctegv = (2+2* cos(pi *(REAL(n -i +1)/REAL((n+1)))))**2
    RETURN
  END FUNCTION correctegv
  
  SUBROUTINE prodotto(invec, n, outvec) !calcola outvec=A invec
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    INTEGER :: n, i
    REAL(dp), DIMENSION(n) :: invec,outvec
    
    outvec(1) = 5*invec(1) + 4*invec(2) + invec(3)
    outvec(2) = 4*invec(1) + 6*invec(2) + 4*invec(3) + invec(4)
    DO i = 3, n-2
       outvec(i) = invec(i-2) + 4*invec(i-1) + 6*invec(i) + 4*invec(i+1) + invec(i+2)
    END DO
    outvec(n-1) = invec(n-3) + 4*invec(n-2) + 6*invec(n-1) + 4*invec(n)
    outvec(n) = invec(n-2) + 4*invec(n-1) + 5*invec(n)
    return
  END SUBROUTINE prodotto

  SUBROUTINE eigenvalues(d, u, eig, n)
    implicit none
    INTEGER :: n
    CHARACTER :: JOBZ = 'N', RANGE = 'A'
    !lower and upper bounds of wanted eigenvalues 
    !used only IF RANGE = 'V':
    REAL(dp) :: VL, VU 
    !lower and upper index (ascending order) of wanted eigenvalues
    !used only IF RANGE = 'I':
    INTEGER :: IL, IU 
    REAL(dp) :: ABSTOL = 1.e-5 !tolerance on approximation error of eigenvalues
    REAL(dp), dimension(n):: eig, d !array with eigenvalues, diagonal
    REAL(dp), dimension(n - 1):: u !upper diagonal
    REAL(dp), dimension(:), allocatable :: Z !not referenced IF JOBZ='N'
    INTEGER, dimension(2*n) :: ISUPPZ
    INTEGER :: LWORK, LIWORK
    REAL(dp), dimension(:), allocatable :: WORK
    INTEGER, dimension(:), allocatable :: IWORK
    INTEGER :: INFO, M
    
    ALLOCATE(WORK(1), IWORK(1))
    CALL DSTEVR(JOBZ, RANGE, n, d, u, VL, VU, IL, IU, ABSTOL, n, eig,Z, n,&
         ISUPPZ,WORK, -1, IWORK, -1,INFO)
    LWORK = WORK(1)
    LIWORK = IWORK(1)
    DEALLOCATE(WORK, IWORK)
    ALLOCATE(WORK(LWORK), IWORK(LIWORK))
    CALL DSTEVR(JOBZ, RANGE, n, d, u, VL, VU,IL,IU,ABSTOL,M,eig,Z,n,&
         ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    DEALLOCATE(WORK, IWORK)
  END SUBROUTINE eigenvalues
END PROGRAM lanczos
