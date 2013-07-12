PROGRAM lanczos
  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 10
  !INTEGER, DIMENSION(N,N) :: A
  INTEGER, DIMENSION(N) :: q

  
  q = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
  CALL prodotto(q, N)

  CONTAINS

  SUBROUTINE prodotto(x, dim)
    !calcola y=Ax
    IMPLICIT NONE
    INTEGER :: dim
    INTEGER :: i
    INTEGER, DIMENSION(dim) :: x,y

    y(1) = 5*x(1) + 4*x(2) + 1*x(3)
    y(2) = 4*x(1) + 6*x(2) + 4*x(3) + x(4)
    do i = 3, dim-2
       y(i) = x(i-2) + 4*x(i-1) + 6*x(i) + 4*x(i+1) + x(i+2)
    end do
    y(n-1) = x(n-3) + 4*x(n-2) + 6*x(n-1) + 4*x(n)
    y(n) = x(n-2) + 4*x(n-1) + 5*x(n)

    !stampa
    do i = 1, dim
       write(*,*) y(i)
    end do

  END SUBROUTINE prodotto

END PROGRAM lanczos
