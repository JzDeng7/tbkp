program test
  use util_math
  use util_tbg
  use util_band
  implicit none
  !> const
  !> rotation angle from bot to top & lattice constant
  real(dp),parameter :: theta = 0.5_dp
! real(dp),parameter :: theta = 1.05_dp
! real(dp),parameter :: theta = 5._dp
  real(dp),parameter :: latt_const = 1.42_dp
  real(dp),parameter :: unitt = 0.0054508_dp

  !> basis vectors of direct lattice & reciprocal lattice
  real(dp) :: oab(2,2), oabinv(2,2) !> origin layer
! real(dp) :: tab(2,2), tabinv(2,2) !> top layer
! real(dp) :: bab(2,2), babinv(2,2) !> bottom layer
  real(dp) :: mab(2,2), mabinv(2,2) !> moire layer
  real(dp) :: mab_(2,2) !> moire layer, dimensionless

  !> variables for Hk generating
  !>   
  real(dp) :: k_D, alpha
  integer :: nDim
! complex(dp), allocatable :: Hk(:,:), rot(:,:)
  complex(dp) :: Hk(148,148), rot(148,148)
! real(dp), allocatable :: eig(:)
  real(dp) :: eig(148)
  !>   
  !> off-diagonal terms generating
  integer :: Glist(2,100), nG

  !> variables for band plotting
! real(dp),allocatable :: Klist(:,:), kpath(:)
  real(dp) :: Klist(2,201), kpath(201)
  real(dp) :: kp(2,1)
  integer :: nK = 50

  !> local variable
  integer :: i, j
  real(dp):: rad
  !> initialization
  oab = 0._dp; mab = 0._dp; k_D = 0._dp

  oab(:,1) = [ dsqrt(3._dp)/2._dp, 3._dp/2._dp ] * latt_const
  oab(:,2) = [-dsqrt(3._dp)/2._dp, 3._dp/2._dp ] * latt_const
! call plot_BZ(theta, oab, k_D, mab, mab_)
  call util_real_matinv(oab, oabinv, 2)

  oabinv = oabinv * 2._dp * pi
  call obtain_kD( oabinv, theta, k_D )
  print*,k_D

  mab(1,:) = [ dsqrt(3._dp)/2._dp,-3._dp/2._dp ] * k_D
  mab(2,:) = [ dsqrt(3._dp)/2._dp, 3._dp/2._dp ] * k_D

  mab_(1,:) = [ dsqrt(3._dp)/2._dp,-3._dp/2._dp ]
  mab_(2,:) = [ dsqrt(3._dp)/2._dp, 3._dp/2._dp ]


  call util_degree2rad(theta, rad)
  alpha = unitt / dsin(rad/2._dp)
  print*, dsin(rad)/dsin(rad/2._dp)

  call obtain_Glist( mab_, 5.55_dp, nG, Glist )

  nDim = nG*4
  print*, nDim
! allocate(Hk(nDim,nDim),rot(nDim,nDim),eig(nDim))

  kpath = 0._dp; Klist = 0._dp
  call obtain_Klist_bnd_2d(Klist, nK, kpath, mab)
  print*, Kpath(1)
  print*, Kpath(51)
  print*, Kpath(101)
  print*, Kpath(151)
  print*, Kpath(201)

  open(321,file='bnd.dat',status='replace')
  do i = 1, nDim
  do j = 1, nK*4+1
    kp(:,1) = Klist(:,j)
    call obtain_hk(Hk, alpha, kp, Glist, mab, nDim)
    call util_hermitian_diag(Hk, rot, eig, nDim)
    write(321,"(2F12.6)") kpath(j), eig(i)
  enddo
  write(321,"(2F12.6)")
  write(321,"(2F12.6)")
  enddo
  close(321)

! deallocate(Hk,rot,eig)
end program test
