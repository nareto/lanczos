PROGRAM lanczos
  IMPLICIT NONE
  INTEGER, PARAMETER :: dim = 500
  !INTEGER, DIMENSION(N,N) :: A
  INTEGER, DIMENSION(dim) :: v
  INTEGER :: i

  
  !v = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
  !do i = 1, dim
  !   v(i) = i
  !end do
  !CALL prodotto(v, dim)
  
  call lanczos_naive(dim)

  CONTAINS

  SUBROUTINE prodotto(invec, dim, outvec)
    !calcola outvec=A invec
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    INTEGER :: dim, i
    REAL(dp), DIMENSION(dim) :: invec,outvec

    outvec(1) = 5*invec(1) + 4*invec(2) + 1*invec(3)
    outvec(2) = 4*invec(1) + 6*invec(2) + 4*invec(3) + invec(4)
    do i = 3, dim-2
       outvec(i) = invec(i-2) + 4*invec(i-1) + 6*invec(i) + 4*invec(i+1) + invec(i+2)
    end do
    outvec(dim-1) = invec(dim-3) + 4*invec(dim-2) + 6*invec(dim-1) + 4*invec(dim)
    outvec(dim) = invec(dim-2) + 4*invec(dim-1) + 5*invec(dim)

    !stampa
    !do i = 1, dim
    !   write(*,*) y(i)
    !end do
    return
  END SUBROUTINE prodotto

  SUBROUTINE lanczos_naive(dim)
    INTEGER :: i, dim
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    REAL(dp), DIMENSION(dim,dim) :: Q
    REAL(dp), DIMENSION(dim) :: rnd, tmp, alfa, beta, r
    character :: JOBZ = 'N', RANGE = 'A', UPLO = 'U'

    call random_seed
    call random_number(rnd)
    Q(:,1) = rnd/NORM2(rnd)
    
       !write(*,*), Q(:,1)
    !write(*,*), tmp
    
    call prodotto(Q(:,1), dim, tmp)
    alfa(1) = dot_product(Q(:,1), tmp)
    r = tmp - alfa(1)*Q(:,1)
    beta(1) = norm2(r)
    Q(:,1+1) = r/beta(1+1)
    
    do i=2, dim
       call prodotto(Q(:,i), dim, tmp)
       alfa(i) = dot_product(Q(:,i), tmp)
       r = tmp - alfa(i)*Q(:,i) - beta(i-1)*Q(:,i-1)
       beta(i) = norm2(r)
       if (i < dim) then
          Q(:,i+1) = r/beta(i)
       end if
    end do

    

    !call dsyevr( JOBZ, RANGE, UPLO, dim, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )
       
  END SUBROUTINE lanczos_naive

END PROGRAM lanczos
