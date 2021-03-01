module util_math
  implicit none
  !> parameters (going to be added in const)
  integer,  parameter :: dp = 8
  real(dp), parameter :: pi = 3.1415926535_dp
  real(dp), parameter :: one = 1._dp
  complex(dp), parameter :: cone = (1._dp, 0._dp)
  complex(dp), parameter :: czero = (0._dp, 0._dp)
contains

  !> transform degree to rad
  subroutine util_degree2rad( degree, rad )
    implicit none
    !> input
    real(dp), intent(in) :: degree
    real(dp), intent(inout) :: rad
    rad = degree * pi / 180._dp
  end subroutine util_degree2rad

  !> transform rad to degree
  subroutine util_rad2degree( rad, degree )
    implicit none
    !> input
    real(dp), intent(in) :: rad
    real(dp), intent(inout) :: degree
    degree = rad * 180._dp / pi
  end subroutine util_rad2degree

  !> generate a rotation matrix in two dimension
  subroutine util_rot2d( rad, rMat )
    implicit none
    !> rotation angle, + for anti-clockwise
    !> input in rad
    real(dp), intent(in) :: rad
    !> rotation matrix acting on target vectors
    real(dp), intent(out):: rMat(2,2)

    rMat(1,:) = [ dcos(rad), -dsin(rad) ]
    rMat(2,:) = [ dsin(rad),  dcos(rad) ]
  end subroutine util_rot2d

  !> to be added
  !! subroutine util_rot3d

  !> obtain inversion of certain vectors / matrices
  subroutine util_real_matInv( Mat, MatInv, nDim )
    implicit none
    !> dimension of the matrix
    integer, intent(in) :: nDim
    !> matrix to process
    real(dp), intent(in) :: Mat(nDim, nDim)
    real(dp), intent(out):: MatInv(nDim, nDim)
    !> local variables
    integer :: i
    real(dp):: tmp1(nDim, nDim)
    real(dp):: tmp2(nDim, nDim)
    !> workspace for lapack
    integer :: info
    integer :: ipiv(nDim)
    !> initialization
    ipiv = 0._dp; tmp1 = 0._dp; tmp2 = 0._dp
    tmp1 = Mat
    do i = 1, nDim
      tmp2(i,i) = one
    enddo
    !> obtain inversion of matrix mat via calling lapack
    call dgesv(nDim, nDim, tmp1, nDim, ipiv, tmp2, nDim, info)
    !> if no inversion
    if (info /= 0) print*, 'something wrong with zgesv.'
    !> result
    MatInv = tmp2
  end subroutine util_real_matInv

  !> to be added
  !! subroutine util_complex_matInv

  !> diagonalize a matrix (Hermitian)
  subroutine util_hermitian_diag( Mat, rMat, Eig, nDim )
    implicit none
    !> input & output
    !> dimension of the matrix
    integer, intent(in) :: nDim
    !> matrix to be diagonalized
    complex(dp), intent(in) :: Mat(nDim, nDim)
    !> unitary rotation of the matrix
    !> to obtain the coefficients of basis
    complex(dp), intent(out) :: rMat(nDim, nDim)
    !> eigenvalues of the matrix
    real(dp), intent(out) :: Eig(nDim)
    !> local variables
    integer :: i, j
    complex(dp) :: tmp(nDim, nDim)
    !> workspace for lapack
    real(dp) :: rwork(nDim*7)
    complex(dp) :: mat_pack((nDim*(nDim+1))/2), cwork(nDim*2)
    !> zhpevx
    integer  :: info, nfound, iwork(nDim*5), ifail(nDim)
    integer  :: lwork, work(nDim*2)

    !> zheev
!   complex(dp) :: work(nDim*2)
!   real(dp) :: rwork(nDim*3)

    !> initialization
    do j = 1, nDim
    do i = 1, j
      mat_pack(i+((j-1)*j)/2) = Mat(i,j)
    enddo
    enddo
    !> zhpevx
    cwork = czero; rwork = 0._dp; iwork = 0
    !> zheev
!   lwork = nDim*2
    tmp = czero; Eig = 0._dp
    !> diagonalization
    call zhpevx('V', 'A', 'U', &
      nDim, mat_pack, 0._dp, 0._dp, 0, 0, -1._dp, &
      nfound, Eig(1), tmp, nDim, &
      cwork, rwork, iwork, ifail, info)

!   call zheev('V', 'U', &
!     nDim, Mat, nDim, eig, work, lwork, rwork, info)
    !> if fail
    if (info < 0) write(*,'(a,i3,a)') &
      'The ', -info, ' argument of zhpevx had an illigal value.'
    if (info > 0) write(*,'(i3,a)') &
      info, ' Eigenvectors fail to converge.'
    !> else
    rMat = tmp
  end subroutine util_hermitian_diag

  !> obtain trace of a matrix
  subroutine util_complex_trace( Mat, nDim, trace )
    implicit none
    !> input & output
    !> dimension
    integer, intent(in) :: nDim
    !> matrix
    complex(dp), intent(in) :: Mat(nDim, nDim)
    complex(dp), intent(out):: trace
    !> local variables
    integer :: i
    complex(dp) :: tmp
    !> initialization
    tmp = 0._dp
    do i = 1, nDim
      tmp = tmp + Mat(i,i)
    enddo
    !> result
    trace = tmp
  end subroutine util_complex_trace

end module util_math
