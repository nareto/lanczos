PROGRAM lanczos
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp=KIND(0.d0), ndiprova = 5
  INTEGER :: i,j, iterations, extreme = 3 !how many extreme eigenvalues to print
  INTEGER, PARAMETER :: dim = 1.e2
  REAL(dp), DIMENSION(dim) :: rnd, tmp, D, r, dtmp, v, w, z, eig
  real(dp), dimension(dim-1) :: U, utmp
  open(unit=4, file="eigv.txt")

  call random_seed
  call random_number(rnd)
  v = rnd/NORM2(rnd)
  w = 0
  do i=1, dim 
     call prodotto(v, dim, tmp)
     D(i) = dot_product(v, tmp)
     r = tmp - D(i)*v - U(i-1)*w
     if (i < dim) then
        U(i) = norm2(r)
        if (U(i) == 0) then
           exit
        end if
        z = r/U(i)
     end if
     w = v
     v = z
     !write(*,*), d !d(i), eig(i)
     dtmp = d
     utmp = u
     call eigenvalues(dtmp(1:i), utmp(1:i), eig, i)
     !write(*,*), eig
     write(4,*), i-1, abs(correctegv(dim) - eig(i)), abs(correctegv(dim-1) - eig(i-1)),&
          abs(correctegv(2) - eig(2)), abs(correctegv(1) - eig(1))
     !write(4,*), i-1, correctegv(dim), eig(1)
  end do

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
END PROGRAM lanczos
