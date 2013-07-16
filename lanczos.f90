PROGRAM lanczos
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp=KIND(0.d0), ndiprova = 5
  INTEGER :: i, extreme =3 !how many extreme eigenvalues to print
  INTEGER, PARAMETER :: dim = 1.e1
  REAL(dp), DIMENSION(dim) ::  D, eig
  REAL(dp), DIMENSION(dim - 1) ::  U

  !call lanczos_naive(dim, D, U)
  call compute_lanczos(dim, D, U)
  !call eigenvalues(D, U, eig, dim)
  !call compare_egv(eig, extreme, dim)

  CONTAINS
    
    real function correctegv(i)
      implicit none
      integer :: i
      INTEGER, PARAMETER :: dp=KIND(0.d0)
      real(dp) :: pi = acos(-1.0)
      correctegv = (2+2* cos(pi *(real(dim -i +1)/real((dim+1)))))**2
      return
    end function correctegv

    subroutine compare_egv(egv, extreme, dim)
      implicit none
      INTEGER, PARAMETER :: dp=KIND(0.d0)
      INTEGER :: i, dim, extreme
      REAL(dp), DIMENSION(dim) :: egv
      real(dp) :: pi = acos(-1.0)

      write(*,"(A16, A16, A16)"), "Correct", "Approximation", "Relative Error"
      do i=1, extreme
         write(*,"(E16.8E2, E16.8E2, E16.8E2)"), correctegv(i), egv(i), (correctegv(i) - egv(i))/correctegv(i)
      end do
      write(*,*), "...."
      do i=dim-extreme+1, dim
       write(*,"(E16.8E2, E16.8E2, E16.8E2)"), correctegv(i), egv(i), (correctegv(i) - egv(i))/correctegv(i)
    end do
  end subroutine compare_egv

    subroutine eigenvalues(d, u, eig, dim)
      !for LAPACK's dsyevr
      implicit none
      integer :: dim
      character :: JOBZ = 'N', RANGE = 'A'
      real(dp) :: VL, VU !lower and upper bounds of wanted eigenvalues - used if RANGE = 'V'
      integer :: IL, IU !lower and upper index (ascending order) of wanted eigenvalues - used if RANGE = 'I'
      real(dp) :: ABSTOL = 0.001 !tolerance on approximation error of eigenvalues
      real(dp), dimension(dim):: eig, d !array with eigenvalues, diagonal
      real(dp), dimension(dim - 1):: u !upper diagonal
      real(dp), dimension(:), allocatable :: Z !not referenced if JOBZ='N'
      integer, dimension(2*dim) :: ISUPPZ
      integer :: LWORK, LIWORK
      real(dp), dimension(:), allocatable :: WORK
      integer, dimension(:), allocatable :: IWORK
      integer :: INFO, M

      allocate(WORK(1), IWORK(1))
      call dstevr(JOBZ, RANGE, dim, d, u, VL, VU, IL, IU, ABSTOL, dim, eig,Z, dim,ISUPPZ,WORK, -1, IWORK, -1,INFO)
      LWORK = WORK(1)
      LIWORK = IWORK(1)
      deallocate(WORK, IWORK)
      allocate(WORK(LWORK), IWORK(LIWORK))
      call dstevr(JOBZ, RANGE, dim, d, u, VL, VU,IL,IU,ABSTOL,M,eig,Z,dim,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
      deallocate(WORK, IWORK)
    end subroutine eigenvalues

  SUBROUTINE prodotto(invec, dim, outvec) !calcola outvec=A invec
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
    return
  END SUBROUTINE prodotto

  SUBROUTINE lanczos_naive(dim, D, U) !dimension, diagonal and upper diagonal of resulting symmetric tridiagonal matrix
    INTEGER :: i, dim
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    REAL(dp), DIMENSION(dim,dim) :: Q
    REAL(dp), DIMENSION(dim) :: rnd, tmp, D, r, dtmp
    real(dp), dimension(dim-1) :: U, utmp
    open(unit=4, file="eigv.txt")

    call random_seed
    call random_number(rnd)
    Q(:,1) = rnd/NORM2(rnd)
    
    call prodotto(Q(:,1), dim, tmp)
    D(1) = dot_product(Q(:,1), tmp)
    write(4,*), 0, abs(correctegv(dim) - d(1)), abs(correctegv(dim-1) - d(1)),&
         abs(correctegv(2) - d(1)), abs(correctegv(1) - d(1))
    
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
       dtmp = d
       utmp = u
       call eigenvalues(dtmp(1:i), utmp(1:i), eig, i)
       write(4,*), i-1, abs(correctegv(dim) - eig(i)), abs(correctegv(dim-1) - eig(i-1)),&
            abs(correctegv(2) - eig(2)), abs(correctegv(1) - eig(1))
    end do
  END SUBROUTINE lanczos_naive


  SUBROUTINE compute_lanczos(dim, D, U) !dimension, diagonal and upper diagonal of resulting symmetric tridiagonal matrix
    INTEGER :: i, j, dim
    INTEGER, PARAMETER :: dp=KIND(0.d0)
    REAL(dp), DIMENSION(dim) :: rnd, tmp, D, r, dtmp, v, w, z
    real(dp), dimension(dim-1) :: U, utmp
    open(unit=4, file="eigv.txt")

    call random_seed
    call random_number(rnd)
    v = rnd/NORM2(rnd)
    w = 0
    do i=1, dim
       call prodotto(v, dim, tmp)
       D(i) = dot_product(v, tmp)
       r = tmp - D(i)*v - U(i)*w
       if (i < dim) then
          U(i) = norm2(r)
          z = r/U(i)
       end if
       v = z
       w = v
       dtmp = d
       utmp = u
       call eigenvalues(dtmp(1:i), utmp(1:i), eig, i)
       write(4,*), i-1, abs(correctegv(dim) - eig(i)), abs(correctegv(dim-1) - eig(i-1)),&
            abs(correctegv(2) - eig(2)), abs(correctegv(1) - eig(1))
    end do
  END SUBROUTINE compute_lanczos

END PROGRAM lanczos
