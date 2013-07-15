PROGRAM lanczos
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp=KIND(0.d0), ndiprova = 5
  INTEGER :: i
  INTEGER, PARAMETER :: dim = 1000
  REAL(dp), DIMENSION(dim) ::  D
  REAL(dp), DIMENSION(dim - 1) ::  U
  !INTEGER, DIMENSION(N,N) :: A
  REAL(dp), DIMENSION(ndiprova) :: v
  REAL(dp), DIMENSION(ndiprova - 1) :: w
  !for LAPACK's dsyevr
  character :: JOBZ = 'N', RANGE = 'A'
  real(dp) :: VL, VU !lower and upper bounds of wanted eigenvalues - used if RANGE = 'V'
  integer :: IL, IU !lower and upper index (ascending order) of wanted eigenvalues - used if RANGE = 'I'
  real(dp) :: ABSTOL = 0.001 !tolerance on approximation error of eigenvalues
  integer :: TOT_EGV = dim !total number of eigenvalues found
  real(dp), dimension(dim) :: EGV !eigenvalues
  real(dp) :: Z !not referenced if JOBZ='N'
  integer, dimension(2*dim) :: ISUPPZ
  !integer :: LWORK = -1, LIWORK = -1
  integer :: LWORK, LIWORK
  real(dp), dimension(:), allocatable :: WORK
  integer, dimension(:), allocatable :: IWORK
  integer :: INFO

  !do i = 1, ndiprova - 1
  !   v(i) = real(i)
  !   w(i) = real(5)
  !end do
  !v(ndiprova) = ndiprova
  !
  !call tridiag(T, v, w, ndiprova)
  !call m_print(T)
  
  call lanczos_naive(dim, D, U)
  if (RANGE == 'I') then
     TOT_EGV = IU - IL +1
  end if
  allocate(WORK(1), IWORK(1))
  call dstevr(JOBZ, RANGE, dim, D, U, VL, VU, IL, IU, ABSTOL, TOT_EGV, EGV, Z, dim, ISUPPZ, WORK, -1, IWORK, -1,INFO)
  LWORK = WORK(1)
  LIWORK = IWORK(1)
  deallocate(WORK, IWORK)
  allocate(WORK(LWORK), IWORK(LIWORK))
  call dstevr(JOBZ, RANGE, dim, D, U, VL, VU, IL, IU, ABSTOL, TOT_EGV, EGV, Z, dim, ISUPPZ, WORK, LWORK, IWORK, LIWORK,INFO)

  call compare_egv(EGV, dim)

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

    subroutine compare_egv(egv, dim)
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    INTEGER :: i, dim
    REAL(dp), DIMENSION(dim) :: egv
    real(dp) :: pi = acos(-1.0), correctegv

    do i=1, dim
       !print "(9f10.4)", (2+2* cos(pi *(real(i)/real((dim+1)))))**2, "(9f10.4)", egv(i)
       !print "(9f10.6)", real(2**4)
       correctegv = (2+2* cos(pi *(real(i)/real((dim+1)))))**2
       write(*,*), correctegv, egv(dim - i +1), correctegv - egv(dim - i +1)
    end do

    end subroutine compare_egv

  SUBROUTINE tridiag(T, D, U, dim) !creates symmetric tridiagonal matrix (stores in T) from vectors D (diagonal) and U (upper diagonal) of dimension dim and dim-1
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    INTEGER :: i, j, dim
    REAL(dp), DIMENSION(dim,dim) :: T
    REAL(dp), DIMENSION(dim) :: D
    REAL(dp), DIMENSION(dim - 1) :: U

    T(1,1) = D(1)
    T(1,2) = U(1)
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

  SUBROUTINE lanczos_naive(dim, D, U) !dimension, diagonal and upper diagonal of resulting symmetric tridiagonal matrix
    INTEGER :: i, j, dim
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    REAL(dp), DIMENSION(dim,dim) :: Q, T
    REAL(dp), DIMENSION(dim) :: rnd, tmp, D, r
    real(dp), dimension(dim-1) :: U

    call random_seed
    call random_number(rnd)
    Q(:,1) = rnd/NORM2(rnd)
    
    call prodotto(Q(:,1), dim, tmp)
    D(1) = dot_product(Q(:,1), tmp)
    r = tmp - D(1)*Q(:,1)
    U(1) = norm2(r)
    Q(:,2) = r/U(1)
    do i=2, dim
       call prodotto(Q(:,i), dim, tmp)
       D(i) = dot_product(Q(:,i), tmp)
       r = tmp - D(i)*Q(:,i) - U(i-1)*Q(:,i-1)
       if (i < dim) then
          U(i) = norm2(r)
          Q(:,i+1) = r/U(i)
       end if
    end do

  END SUBROUTINE lanczos_naive


END PROGRAM lanczos
