PROGRAM lanczos
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp=KIND(0.d0), ndiprova = 5
  INTEGER :: i
  REAL(dp), DIMENSION(ndiprova, ndiprova) ::  T
  INTEGER, PARAMETER :: dim = 500
  !INTEGER, DIMENSION(N,N) :: A
  REAL(dp), DIMENSION(ndiprova) :: v
  REAL(dp), DIMENSION(ndiprova - 1) :: w

  !do i = 1, ndiprova - 1
  !   v(i) = real(i)
  !   w(i) = real(5)
  !end do
  !v(ndiprova) = ndiprova
  !
  !call tridiag(T, v, w, ndiprova)
  !call m_print(T)

  call lanczos_naive(dim)

  CONTAINS

    subroutine m_print(A)
      INTEGER, PARAMETER :: dp=KIND(0.d0)
      real(dp), dimension(:,:), intent(in) :: A
      integer :: m, n, i, j
      m=size(A,dim=1)
      n=size(A,dim=2)
      do i=1,m
         print "(9f10.4)", (A(i,j),j=1,n)
      end do
      print *
    end subroutine m_print

  SUBROUTINE tridiag(T, D, U, dim) !creates symmetric tridiagonal matrix (stores in T) from vectors D (diagonal) and U (upper diagonal) of dimension dim and dim-1
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    INTEGER :: i, j, dim
    REAL(dp), DIMENSION(dim,dim) :: T
    REAL(dp), DIMENSION(dim) :: D
    REAL(dp), DIMENSION(dim - 1) :: U

    T(1,1) = 1.0
    T(1,2) = 5.0
    do j=3, dim
       T(1,j) = 0
    end do
    do i =2, dim-1 !righe
       do j=1, i-2
          T(i,j) = 0
       end do
       T(i,i-1) = U(i-1) 
       T(i,i) = D(i)
       T(i,i+1) = U(i)
       do j=i+2, dim
          T(i,j) = 0
       end do
    end do
    do j=1, dim-2
       T(dim,j) = 0
    end do
    T(dim,dim-1) = U(dim-1)
    T(dim,dim) = D(dim)
  END SUBROUTINE tridiag

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
    INTEGER :: i, j, dim
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    REAL(dp), DIMENSION(dim,dim) :: Q, T
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
    

    !call dsyevr( JOBZ, RANGE, UPLO, dim, T, dim, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )
       
  END SUBROUTINE lanczos_naive


END PROGRAM lanczos
