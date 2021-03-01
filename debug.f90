program plot
  use util_math
  use util_tbg
  implicit none
! integer,  parameter :: dp = 8
! real(dp), parameter :: pi = 3.141592653589793238_dp
! real(dp), parameter :: th = 5.00_dp
! real(dp), parameter :: th = 1.05_dp
  real(dp), parameter :: th = 0.50_dp
  real(dp) :: rad

  real(dp) :: ab(2,2),  invab(2,2)
  real(dp) :: abt(2,2), invabt(2,2)
  real(dp) :: abb(2,2), invabb(2,2)
  real(dp) :: Kuc(2), Kabs, qabs
  real(dp) :: mab(2,2), mbz(2,7)
  real(dp) :: fbz(2,7)
  real(dp) :: ng(2), ngabs
  real(dp) :: ng_1(2), ng_2(2), ng_3(2)
  integer  :: ii, mm, nn, i, j, kk

  integer :: ndim
  real(dp):: alpha

! complex(dp):: hk(37*4,37*4)
  complex(dp), allocatable :: hk(:,:), rot(:,:)
  real(dp),    allocatable :: eigenv(:)
  real(dp)   :: kp(2,1), kpt(2,81), kpath(81)
  complex(dp):: sum_hk
  integer    :: nGtot, iind
  real(dp)   :: mab_(2,2), gxy(2), absgxy, listbG(2,100)
  integer :: Klist(2,100)

  real(dp) :: G1(2,1)
  real(dp) :: MX(2,1)
  real(dp) :: G2(2,1)
  real(dp) :: XM(2,1)

  !> lattice constant a = 1 = |\bt_{\alpha} - \bt_{\beta}|
  !> original lattice without rotation
  ab(:,1) = (/ sqrt(3._dp)*0.5_dp, 3._dp*0.5_dp /) ! * a
  ab(:,2) = (/-sqrt(3._dp)*0.5_dp, 3._dp*0.5_dp /) ! * a

  !> add rotation
  rad = th*pi / 360._dp
  alpha = 0.0054508_dp/dsin(rad) ! v_F = 0.9*10**6
! alpha = 0.0049058_dp/dsin(rad) ! v_F = 10**6
! alpha = 0.1_dp
  print"('alpha =',F8.4,'; ','theta =',F8.4)", alpha, th

  abt(1,1) = cos(rad)*ab(1,1)-sin(rad)*ab(2,1) ! top layer \theta/2
  abt(2,1) = sin(rad)*ab(1,1)+cos(rad)*ab(2,1) ! a_1

  abt(1,2) = cos(rad)*ab(1,2)-sin(rad)*ab(2,2)
  abt(2,2) = sin(rad)*ab(1,2)+cos(rad)*ab(2,2) ! a_2

  abb(1,1) = cos(rad)*ab(1,1)+sin(rad)*ab(2,1) ! bottom layer -\theta/2
  abb(2,1) =-sin(rad)*ab(1,1)+cos(rad)*ab(2,1) ! a_1

  abb(1,2) = cos(rad)*ab(1,2)+sin(rad)*ab(2,2)
  abb(2,2) =-sin(rad)*ab(1,2)+cos(rad)*ab(2,2) ! a_2

  invab  = 0._dp
  invabt = 0._dp
  invabb = 0._dp
  !> reciprocal lattice bases
  call invreal22(ab,  invab)
  call invreal22(abt, invabt)
  call invreal22(abb, invabb)
  invab  = invab  * 2._dp * pi
  invabt = invabt * 2._dp * pi
  invabb = invabb * 2._dp * pi

  !> \bbK (1/3,-1/3) of top layer FBZ
  Kuc = 0._dp
  Kuc  = matmul((/ 1._dp/3.,-1._dp/3. /), invabt(:,:))
  !> module of \bbK
  Kabs = 0._dp
  Kabs = sqrt(Kuc(1)*Kuc(1)+Kuc(2)*Kuc(2))
  ! \bq_\da = \bbKpr - \bbK
  qabs = 0._dp
  qabs = 2. * dsin(rad) * Kabs
! write(*,*) qabs
! write(*,"(2F12.6)") Kuc

  ! moire super lattice reciprocal bases
  ! mab is mabinv
  mab(:,:) = 0._dp
  mab(1,:) = (/ sqrt(3._dp)/2._dp,-3._dp/2._dp /) * qabs
  mab(2,:) = (/ sqrt(3._dp)/2._dp, 3._dp/2._dp /) * qabs

  !-- generate listbG, i.e., all reciprocal vectors included
  mbz(:,:) = 0._dp
  !> plot start form \bbK
  mbz(1,:) = Kuc(1)
  mbz(2,:) = Kuc(2)
  do i = 1, 7
  mbz(:,i) = mbz(:,i)-matmul((/ 1._dp/3  ,-1._dp/3+1/),mab(:,:))
  enddo
! mbz(:,1) = mbz(:,1)+matmul((/ 1._dp/3  ,-1._dp/3  /),mab(:,:))
! mbz(:,2) = mbz(:,2)+matmul((/-1._dp/3+1, 1._dp/3  /),mab(:,:))
  mbz(:,3) = mbz(:,3)+matmul((/ 1._dp/3  ,-1._dp/3+1/),mab(:,:))
! mbz(:,4) = mbz(:,4)+matmul((/-1._dp/3  , 1._dp/3  /),mab(:,:))
! mbz(:,5) = mbz(:,5)+matmul((/ 1._dp/3-1,-1._dp/3  /),mab(:,:))
! mbz(:,6) = mbz(:,6)+matmul((/-1._dp/3  , 1._dp/3-1/),mab(:,:))
! mbz(:,7) = mbz(:,7)+matmul((/ 1._dp/3  ,-1._dp/3  /),mab(:,:))

  mab_(1,:) = (/ 0.5d0,-sqrt(3.d0)*0.5d0 /)
  mab_(2,:) = (/ 0.5d0, sqrt(3.d0)*0.5d0 /)
  !> determine the dimension of Hk
  iind=0
  do i=-7,7
  do j=-7,7
    gxy(:)=matmul((/i,j/),mab_)
!   write(2,"(2F8.4)") gxy
    absgxy=dsqrt(gxy(1)**2+gxy(2)**2)
    if (absgxy > 3.d0) cycle
    iind=iind+1
  enddo
  enddo
  !> dimension / 2
  nGtot=iind
  write(*,*) nGtot
  allocate( hk(nGtot*4, nGtot*4) )
  allocate( rot(nGtot*4, nGtot*4) )
  allocate( eigenv(nGtot*4) )
  !-- generate listbG ends

  mbz(:,:) = 0._dp
  !> plot start form \bbK
  mbz(1,:) = Kuc(1)
  mbz(2,:) = Kuc(2)
  do i = 1, 7
  mbz(:,i) = mbz(:,i)-matmul((/ 1._dp/3  ,-1._dp/3+1/),mab(:,:))
  enddo
  mbz(:,1) = mbz(:,1)+matmul((/ 1._dp/3  ,-1._dp/3  /),mab(:,:))
  mbz(:,2) = mbz(:,2)+matmul((/-1._dp/3+1, 1._dp/3  /),mab(:,:))
  mbz(:,3) = mbz(:,3)+matmul((/ 1._dp/3  ,-1._dp/3+1/),mab(:,:))
  mbz(:,4) = mbz(:,4)+matmul((/-1._dp/3  , 1._dp/3  /),mab(:,:))
  mbz(:,5) = mbz(:,5)+matmul((/ 1._dp/3-1,-1._dp/3  /),mab(:,:))
  mbz(:,6) = mbz(:,6)+matmul((/-1._dp/3  , 1._dp/3-1/),mab(:,:))
  mbz(:,7) = mbz(:,7)+matmul((/ 1._dp/3  ,-1._dp/3  /),mab(:,:))
  !> plot MBZ
  open(123,file='mbz.dat',status='replace')
  do i = -10,10
  do j = -10,10
  ng   = 0._dp
  ng_1 = 0._dp
  ng_2 = 0._dp
  ng_3 = 0._dp
  ng(:)   = matmul((/ 1._dp*i  , 1._dp*j  /),mab(:,:)) ! (m,n) in Cartesian
  ng_1(:) = matmul((/ 0._dp    , 0._dp    /),mab(:,:)) ! coordinates in \bq_\da representation
  ng_2(:) = matmul((/ 0._dp    , 1._dp    /),mab(:,:))
  ng_3(:) = matmul((/-1._dp    , 0._dp    /),mab(:,:))
  ngabs = sqrt(ng(1)*ng(1)+ng(2)*ng(2))
  gxy(:)=matmul((/i,j/),mab_)
  absgxy=dsqrt(gxy(1)**2+gxy(2)**2)
  if (absgxy > 3.d0) cycle
!   listbG(:,iind)=(/real(i,dp),real(j,dp)/)
! if (ngabs < 1.8_dp) then !> th = 6.008
! if (ngabs < 0.35_dp) then !> th = 1.05
! if (sqrt(i*1._dp*i+j*1._dp*j) < 5) then
    write(123,"(2F12.6,3I5)") mbz(:,3)+ng(:),         i, j, 0 !\bbK
    write(123,"(2F12.6,3I5)") mbz(:,2)+ng_1(:)+ng(:), i, j, 1 !\bbKpr
    write(123,"(2F12.6,3I5)")
    write(123,"(2F12.6,3I5)")
    write(123,"(2F12.6,3I5)") mbz(:,3)+ng(:),         i, j,   0
    write(123,"(2F12.6,3I5)") mbz(:,2)+ng_2(:)+ng(:), i, j+1, 2
    write(123,"(2F12.6,3I5)")
    write(123,"(2F12.6,3I5)")
    write(123,"(2F12.6,3I5)") mbz(:,3)+ng(:),         i,   j, 0
    write(123,"(2F12.6,3I5)") mbz(:,2)+ng_3(:)+ng(:), i-1, j, 3
    write(123,"(2F12.6,3I5)")
    write(123,"(2F12.6,3I5)")
! endif
  enddo
  enddo
  close(123)

  !> bottom layer FBZ
  open(124,file='bbz.dat',status='replace')
  fbz(:,:)=0._dp
  fbz(:,1) = matmul((/ 1._dp/3  ,-1._dp/3  /),invabb(:,:))
  fbz(:,2) = matmul((/-1._dp/3+1, 1._dp/3  /),invabb(:,:))
  fbz(:,3) = matmul((/ 1._dp/3  ,-1._dp/3+1/),invabb(:,:))
  fbz(:,4) = matmul((/-1._dp/3  , 1._dp/3  /),invabb(:,:))
  fbz(:,5) = matmul((/ 1._dp/3-1,-1._dp/3  /),invabb(:,:))
  fbz(:,6) = matmul((/-1._dp/3  , 1._dp/3-1/),invabb(:,:))
  fbz(:,7) = matmul((/ 1._dp/3  ,-1._dp/3  /),invabb(:,:))
  do i = -5, 5
  do j = -5, 5
  ng=0._dp
  ng(:) = matmul((/ 1._dp*i  , 1._dp*j  /),invabb(:,:))
  ngabs = sqrt(ng(1)*ng(1)+ng(2)*ng(2))
  if (ngabs <  9._dp) then
    write(124,"(2F12.6)") (fbz(:,ii)+ng(:),ii=1,7)
    write(124,"(2F12.6)")
    write(124,"(2F12.6)")
  endif
  enddo
  enddo
  close(124)
  !> top layer FBZ
  open(125,file='tbz.dat',status='replace')
  fbz(:,:)=0._dp
  fbz(:,1) = matmul((/ 1._dp/3  ,-1._dp/3  /),invabt(:,:))
  fbz(:,2) = matmul((/-1._dp/3+1, 1._dp/3  /),invabt(:,:))
  fbz(:,3) = matmul((/ 1._dp/3  ,-1._dp/3+1/),invabt(:,:))
  fbz(:,4) = matmul((/-1._dp/3  , 1._dp/3  /),invabt(:,:))
  fbz(:,5) = matmul((/ 1._dp/3-1,-1._dp/3  /),invabt(:,:))
  fbz(:,6) = matmul((/-1._dp/3  , 1._dp/3-1/),invabt(:,:))
  fbz(:,7) = matmul((/ 1._dp/3  ,-1._dp/3  /),invabt(:,:))
  do i = -5, 5
  do j = -5, 5
  ng=0._dp
  ng(:) = matmul((/ 1._dp*i  , 1._dp*j  /),invabt(:,:))
  ngabs = sqrt(ng(1)*ng(1)+ng(2)*ng(2))
  if (ngabs <  9._dp) then
    write(125,"(2F12.6)") (fbz(:,ii)+ng(:),ii=1,7)
    write(125,"(2F12.6)")
    write(125,"(2F12.6)")
  endif
  enddo
  enddo
  close(125)

! !> rearrange part
! !> rearrange part ends, output: mbz.dat and tbz.dat, bbz.dat

  !> generate k-path
  !> H(k) at any k-point

  call obtain_Klist(mab_, 3.05_dp, nGtot, Klist)
  ndim =nGtot*2*2  ! nGtot * 2 * 2
  write(*,*) ndim
  write(*,"(2I5)") Klist(:,:)
  write(1238,"(2F8.5)") real(Klist(:,:),dp)

  !> check Hermicity
  MX(:,1) = matmul( [ 1.d0/7.d0, 3.d0/7.d0 ],mab(:,:)) ! top K , B
  kp(:,1) = MX(:,1)
  call obtain_hk(hk, alpha, kp, Klist, mab, ndim)
! call hmnk( hk, ndim, kp, mab, alpha )
  do i = 1, ndim
  do j = 1, ndim
  sum_hk = hk(i,j) + hk(j,i)
  if (aimag(sum_hk) /= 0.0d0) then
    print*,"Hk is not Hermitian."
    stop
  endif
  enddo
  enddo
  print*,"Hk is Hermitian."

  !> diagonalization
! !> output: bnd.dat (band structure)
  G1(:,1) = matmul( [ 1.d0/3.d0,-1.d0/3.d0 ],mab(:,:)) ! A
  MX(:,1) = matmul( [ 0.d0/3.d0, 0.d0/3.d0 ],mab(:,:)) ! top K , B
  G2(:,1) = matmul( [-1.d0/3.d0, 1.d0/3.d0 ],mab(:,:)) ! C
  XM(:,1) = matmul( [ 2.d0/3.d0, 1.d0/3.d0 ],mab(:,:)) ! D

  open(321,file='bnd.dat',status='replace')
  kpath(:)=0._dp
  call bnd_plt(kpath, G1, G2, MX, XM, mab)
! print*, kpath(1)
! print*, kpath(21)
! print*, kpath(41)
! print*, kpath(61)
! print*, kpath(81)
  call k_path(kpt, G1, G2, MX, XM)
  do i = 1, ndim
  do j = 1, 81
    kp(:,1) = kpt(:,j)
    call obtain_hk(hk, alpha, kp, Klist, mab, ndim)
! call hmnk( hk, ndim, kp, mab, alpha )
    call util_diag( hk, ndim, eigenv, rot)
    write(321,"(2F12.6)") kpath(j), eigenv(i)
  enddo
  write(321,"(2F12.6)")
  write(321,"(2F12.6)")
  enddo
  close(321)

  !> clean up
! call execute_command_line("rm fort.*")
end program plot


!> inversion of a real matrix
subroutine  invreal22(mat2,invmat2)
implicit none
real(8) ,intent(in):: Mat2(2,2)
real(8) ,intent(out):: Invmat2(2,2)

real(8) :: Mat(3,3)
real(8) :: Invmat(3,3)

real(8) :: pa(3),pb(3),pc(3),volumn
real(8) :: pra(3),prb(3),prc(3),t
Mat(:,:)=0._8
Mat(1:2,1:2)=Mat2(:,:)
Mat(3,3)=1._8

pa(:)=Mat(:,1); pb(:)=Mat(:,2); pc(:)=Mat(:,3)

volumn =  pa(1)*(pb(2)*pc(3)-pb(3)*pc(2)) &
        + pa(2)*(pb(3)*pc(1)-pb(1)*pc(3)) &
        + pa(3)*(pb(1)*pc(2)-pb(2)*pc(1))

t=1._8/volumn

pra(1)=t*(pb(2)*pc(3)-pb(3)*pc(2))
pra(2)=t*(pb(3)*pc(1)-pb(1)*pc(3))
pra(3)=t*(pb(1)*pc(2)-pb(2)*pc(1))

prb(1)=t*(pc(2)*pa(3)-pc(3)*pa(2))
prb(2)=t*(pc(3)*pa(1)-pc(1)*pa(3))
prb(3)=t*(pc(1)*pa(2)-pc(2)*pa(1))

prc(1)=t*(pa(2)*pb(3)-pa(3)*pb(2))
prc(2)=t*(pa(3)*pb(1)-pa(1)*pb(3))
prc(3)=t*(pa(1)*pb(2)-pa(2)*pb(1))

invMat(1,:) = pra(:); invMat(2,:) = prb(:); invMat(3,:) = prc(:)

invMat2(:,:)=invMat(1:2,1:2)
end subroutine  invreal22
