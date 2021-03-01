module util_band
  use util_math
  implicit none
contains

  !> generate k-path and k-list
  subroutine obtain_Klist_bnd_2d( Klist, nK, kpath, ab )
    implicit none
    !> input & output
    !> number of kpoint of each kpath segment
    integer, intent(in) :: nK
    !> coordinates of kpoint included
!   real(dp), allocatable, intent(out) :: Klist(:,:)
    real(dp), intent(out) :: Klist(2,201)
    !> kpath for band plotting
!   real(dp), allocatable, intent(out) :: kpath(:)
    real(dp), intent(out) :: kpath(201)
    !> basis vectors in reciprocal space in cartesian
    real(dp), intent(in) :: ab(2,2)
    !> local variables
    integer :: i, j
    !> number of kpath segment
    integer :: nk_seg
    real(dp) :: dk_frac(2,10), dk_cart(2,10), dkabs(10)
    real(dp) :: K_highsym(2,10), k_path(2,10)
    !> list of high symmetry kpoint
    !> written in fractional
    !> full path
    K_highsym = 0._dp
    K_highsym(:,1) = [ 1._dp/3._dp,-1._dp/3._dp ] !> A, bot K
    K_highsym(:,2) = [ 0._dp/3._dp, 0._dp/3._dp ] !> B, top K
    K_highsym(:,3) = [-1._dp/3._dp, 1._dp/3._dp ] !> C
    K_highsym(:,4) = [ 2._dp/3._dp, 1._dp/3._dp ] !> D
    K_highsym(:,5) = [ 1._dp/3._dp,-1._dp/3._dp ] !> A, bot K

    !> full path
    K_path = 0._dp
    K_path(:,1) = matmul([ 1._dp/3._dp,-1._dp/3._dp ],ab)!> A, bot K
    K_path(:,2) = matmul([ 0._dp/3._dp, 0._dp/3._dp ],ab)!> B, top K
    K_path(:,3) = matmul([-1._dp/3._dp, 1._dp/3._dp ],ab)!> C
    K_path(:,4) = matmul([ 2._dp/3._dp, 1._dp/3._dp ],ab)!> D
    K_path(:,5) = matmul([ 1._dp/3._dp,-1._dp/3._dp ],ab)!> A, bot K

    nk_seg = 4 !> = nk_highsym - 1
!   allocate(Klist(2,nk_seg*nK+1)) !> stat=status?
!   allocate(kpath(nk_seg*nK+1)) !> stat=status?

    !> generate Klist
    !> generate delta_k
    do i = 1, nk_seg !> variable
    dk_frac(:,i) = (K_highsym(:,i+1) - K_highsym(:,i)) / real(nK,dp)
    dk_cart(:,i) = (K_path(:,i+1) - K_path(:,i)) / real(nK,dp)
    dkabs(i) = dsqrt(dot_product(dk_cart(:,i),dk_cart(:,i)))
    enddo
    !>  
    !>
    Klist = 0._dp; kpath = 0._dp
!   Klist(:,1) = K_highsym(:,1)
    Klist(:,1) = K_path(:,1)
    do j = 1, nk_seg
    do i = (j-1)*nK+1, j*nK
!     Klist(:,i+1) = Klist(:,i) + dk_frac(:,j)
      Klist(:,i+1) = Klist(:,i) + dk_cart(:,j)
      kpath(i+1) = kpath(i) + dkabs(j)
    enddo
    enddo

    !> generate kpath for band plotting
!   deallocate(Klist,kpath)
  end subroutine obtain_Klist_bnd_2d

end module util_band
