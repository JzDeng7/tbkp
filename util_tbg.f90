module util_tbg
  use util_math! only: dp, pi
  implicit none
contains

  !> obtain the absolute value of K_top - K_bot, i.e., |k_D|
  !> with rotation angle begin theta
  subroutine obtain_kD( abinv, theta, k_D )
    implicit none
    !> input & output
    !> basis vectors of original graphene BZ
    !> information of ab(2,2) & latt_const is included in abinv(2,2)
    real(dp), intent(in) :: abinv(2,2)
    !> rotation angle of tbg
    real(dp), intent(in) :: theta
    !> absolute value of \bt_{\bq}
    real(dp), intent(out):: k_D
    !> local variables
    real(dp) :: rad, topK(2), Kabs
    !> initialization
    rad = 0._dp; topK = 0._dp
    topK = matmul( [ 1._dp/3._dp, -1._dp/3._dp ], abinv )
    Kabs = dsqrt(topK(1)**2 + topK(2)**2)
    call util_degree2rad( theta/2._dp, rad )
    !> result
    k_D = 2._dp * dsin(rad) * Kabs
  end subroutine obtain_kD

  !> obtain sites on moire lattice included in kp Hamiltonian
  subroutine obtain_Glist( ab, qmax, nG, Glist )
    implicit none
    !> input & output
    !> moire superlattice basis in fractional
    real(dp), intent(in) :: ab(2,2)
    !> radius of the range of sites included on moire reciprocal lattice
    real(dp), intent(in) :: qmax
    !> number of sites (K-points) ~, i.e., dimension of kp Hamiltonian
    integer, intent(out) :: nG
    integer, intent(out) :: Glist(2,100)
    !> local variables
    integer :: i, j, k
    real(dp):: gxy(2), gxyabs
    !> initialization
    k = 0
    !> output file
!   open(101, file='topKlist.dat', status='replace')
!   open(102, file='botKlist.dat', status='replace')
    do i = -10, 10
    do j = -10, 10
      !> fractional distance to top K
      gxy(:) = matmul( [i, j], ab )
      gxyabs = dsqrt(gxy(1)**2 + gxy(2)**2)
      !>
      if (gxyabs > qmax) cycle
      k = k + 1
      Glist(:,k) = [i, j]
      !> one can check coordinates of all Kpoints included
!     write(101,"(3I5,2F8.4)") k, Glist(:,k), gxy(:)
!     write(102,"(3I5,2F8.4)") k, Glist(:,k), &
!       gxy(:) + matmul( [ 1._dp/3._dp, -1._dp/3._dp ], ab )
    enddo
    enddo
!   close(101);close(102)
    !> result: dimension of kp Hamiltonian
    nG = k
  end subroutine obtain_Glist

  !> obtain kp Hamiltonian at a given kpoint
  subroutine obtain_hk( Hk, alpha, Kpt, Glist, mab, nDim )
    implicit none
    !> input & output
    !> dimension of Hk
    integer, intent(in) :: nDim
    !> kpoint
    real(dp), intent(in) :: Kpt(2,1)
    !> dimensionless variable
    real(dp), intent(in) :: alpha
    !>  
    integer, intent(in) :: Glist(2,100)
    !> basis
    real(dp), intent(in) :: mab(2,2)
    !> Hk
    complex(dp), intent(out) :: Hk(nDim, nDim)
    !> local variables
    !> local indices
    integer :: i, j, k
    real(dp):: Gi(2), Gj(2), gxy(2), gxyabs
    !> local parameters
    integer :: nG !> nKpts
    real(dp):: q(2,3), qabs
    real(dp):: topK(2,1), botK(2,1)

    !> to be moved into const.f90   
    !> Pauli Sigma matrix
    complex(dp) :: s_x(2,2), s_y(2,2), s_z(2,2), s_0(2,2), tmat(2,2,3)
    !> initialization
    s_x = 0._dp; s_y = 0._dp; s_z = 0._dp; s_0 = 0._dp
    tmat= 0._dp; Hk = 0._dp;   q = 0._dp
    !> set value
    s_0(1,1) = 1._dp; s_0(2,2) = 1._dp
    s_x(1,2) = 1._dp; s_x(2,1) = 1._dp
    s_y(1,2) = (0._dp,-1._dp); s_y(2,1) = (0._dp,1._dp)
    !> to be moved part ends
    !> initialization
    nG = nDim / 4
    !> inter-layer hopping term
    !> t_{\bq}
    tmat(:,:,1) = s_0 + s_x
    tmat(:,:,2) = s_0 + dcos(2._dp*pi/3._dp)*s_x + dsin(2._dp*pi/3._dp)*s_y
    tmat(:,:,3) = s_0 + dcos(2._dp*pi/3._dp)*s_x - dsin(2._dp*pi/3._dp)*s_y
!   write(*,"(4F12.6)") tmat
    !> \bq_i
    q(:,1) = matmul( [ 1._dp/3._dp, -1._dp/3._dp ], mab )
    q(:,2) = matmul( [ 1._dp/3._dp,  2._dp/3._dp ], mab )
    q(:,3) = matmul( [-2._dp/3._dp, -1._dp/3._dp ], mab )
    qabs = dsqrt(q(1,1)**2 + q(2,1)**2)

    !> generate Hk
    !> diagonal elements, Dirac fermion?
    do i = 1, nG
      topK(:,1) =-(Kpt(:,1)-matmul(real(Glist(:,i),dp),mab)) / qabs
      botK(:,1) =-(Kpt(:,1)-matmul(real(Glist(:,i),dp),mab)-q(:,1)) / qabs
      hk([-1:0]+2*i, [-1:0]+2*i) = topK(1,1)*s_x + topK(2,1)*s_y
      hk([-1:0]+2*(i+nG), [-1:0]+2*(i+nG)) = botK(1,1)*s_x + botK(2,1)*s_y
    enddo

    !> off-diagonal elements
    do k = 1, 3
    do i = 1, nG
    Gi(:) = matmul(real(Glist(:,i),dp),mab)
    do j = 1, nG
      Gj(:) = matmul(real(Glist(:,j),dp),mab) + q(:,1)
      gxy(:)= Gi(:) + q(:,k) - Gj(:)
      gxyabs= dsqrt(gxy(1)**2 + gxy(2)**2)
      if (gxyabs > qabs/100._dp) cycle
      hk([-1:0]+2*(i), [-1:0]+2*(j+nG)) = tmat(:,:,k) * alpha
      hk([-1:0]+2*(j+nG), [-1:0]+2*(i)) = tmat(:,:,k) * alpha
    enddo
    enddo
    enddo
  end subroutine obtain_hk

  !> to be added and adjusted into allocatable
  !>   

  !> obtain moire BZ, top BZ & bot BZ
  subroutine plot_BZ( theta, oab, k_D, mab, mab_ )
    implicit none
    !> input & output
    !> rotation angle in degree (to be transformed to rad)
    real(dp), intent(in) :: theta
    !> original basis vectors of graphene
    real(dp), intent(in) :: oab(2,2)
    !> |k_D| = |t_\bq|
    real(dp), intent(out):: k_D
    !> moire superlattice basis
    !> mab (with unit), mab_ (dimensionless)
    !> mab = mab_ * k_D
    real(dp), intent(out):: mab(2,2), mab_(2,2)

    !> local variables
    !> inversion of oab, basis of original reciprocal space
    real(dp) :: oabinv(2,2)
    !> basis of rotated lattice
    real(dp) :: tab(2,2), bab(2,2)
    real(dp) :: tabinv(2,2), babinv(2,2) !> t: top, b: bot
    !> rotation angle of each layer (in rad)
    real(dp) :: rad, trMat(2,2), brMat(2,2)

    !> initialization
    rad = 0._dp; trMat = 0._dp; brMat = 0._dp
    tab = 0._dp; bab = 0._dp
    k_D = 0._dp

    !> pre
    !>> rotation
    call util_degree2rad( theta/2._dp, rad )
    call util_rot2d( rad, trMat)
    call util_rot2d(-rad, brMat)
    tab = matmul(trMat,oab)
    bab = matmul(brMat,oab)
    call util_real_matinv(tab, tabinv, 2)
    call util_real_matinv(bab, babinv, 2)
    tabinv = tabinv * 2._dp * pi
    babinv = babinv * 2._dp * pi

    !> generate moire basis
    call util_real_matinv(oab, oabinv, 2)
    oabinv = oabinv * 2._dp * pi

    call obtain_kD( oabinv, theta, k_D )

    mab_(1,:) = [ dsqrt(3._dp)/2._dp,-3._dp/2._dp ]
    mab_(2,:) = [ dsqrt(3._dp)/2._dp, 3._dp/2._dp ]
    mab = mab_ * k_D
  end subroutine plot_BZ

end module util_tbg
