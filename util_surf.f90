! reference:    M P Lopez Sancho et al 1985 J. Phys. F: Met. Phys. 15 851
! date     : 2011-1-20
! author   : zjwang@iphy.ac.cn
! modified : 2012-9-19
SUBROUTINE green(lda,H00,H01,rho,omega,info)
use control, only :dp,eta,maxloop,CI,ifsurf,sfE
implicit none

integer    , intent(in)  :: lda
real(dp)   , intent(in)  :: omega
complex(dp), intent(in)  :: H00(lda,lda),H01(lda,lda)
real(dp)   , intent(out) :: rho 
integer    , intent(out) :: info 

!##########################################
complex(dp)  :: G00(lda,lda),ttmmpp
complex(dp)  :: ti(lda,lda),ti_1(lda,lda)
complex(dp)  :: tj(lda,lda),tj_1(lda,lda)
complex(dp)  :: omega_h00(lda,lda),Cone(lda,lda)
complex(dp)  :: omega_hs(lda,lda),Sigma(lda,lda)
!##########################################

integer  ::i,j,iter
complex(dp)  :: temp(lda,lda),inver_omega_h(lda,lda),tij(lda,lda),tji(lda,lda),ti2(lda,lda),tj2(lda,lda)
complex(dp)  :: epti(lda,lda),eptj(lda,lda),T(lda,lda),T_tilde(lda,lda)
real(dp)     :: diff
info=0

!###define################
Cone=0.d0
do i=1,lda
   Cone(i,i)=(1.d0,0.d0)
enddo
!###define################

!print*,(omega-eta*CI)
omega_h00=(omega+eta*CI)*Cone-H00



inver_omega_h=omega_h00
call inv(lda,inver_omega_h)

! add for bulk 2012.11.28
IF(.not.ifsurf) THEN
G00=inver_omega_h
rho=0.d0
  ttmmpp=0.d0
    do j=1,lda
      ttmmpp=ttmmpp+G00(j,j)
    enddo
 rho=-aimag(ttmmpp)
return
ENDIF
! add for bulk 2012.11.28

temp=conjg(transpose(H01))
call mat_mul(lda,inver_omega_h,temp,ti_1)
call mat_mul(lda,inver_omega_h,H01,tj_1)

epti=Cone
eptj=Cone
T=ti_1
T_tilde=tj_1
iter_loop:do iter=1,maxloop
          !DO WHILE(abs(diff)>0.1d-8)
          !--------ti_1,tj_1 >>>>> ti, tj--------
          call mat_mul(lda,ti_1,tj_1,tij)
          call mat_mul(lda,tj_1,ti_1,tji)
          call mat_mul(lda,ti_1,ti_1,ti2)
          call mat_mul(lda,tj_1,tj_1,tj2)
          
          temp=Cone-tij-tji
          call inv(lda,temp)
          call mat_mul(lda,temp,ti2,ti)
          call mat_mul(lda,temp,tj2,tj)
          
          
          temp=epti
          call mat_mul(lda,temp,tj_1,epti)
          call mat_mul(lda,epti,ti,temp)
          T=T+temp

          temp=eptj
          call mat_mul(lda,temp,ti_1,eptj)
          call mat_mul(lda,eptj,tj,temp)
          T_tilde=T_tilde+temp

          diff=0.d0
          do i=1,lda
             do j=1,lda
              diff= diff+abs(ti(j,i))+abs(tj(j,i)) 
             enddo
          enddo
          if(abs(diff)<0.1d-8) EXIT iter_loop
          
          ti_1=ti
          tj_1=tj
          !----------------------------------------
          !print*,diff
          END DO iter_loop

          !------error--------------
          if(iter>maxloop) then
           info=1
           write(*,*) "ERROR: T_matrix is not convergent"
           return
          endif
          !------error-------------

call mat_mul(lda,H01,T,temp)
!call mat_mul(lda,conjg(transpose(H01)),T_tilde,temp)
temp=omega_h00-temp
call inv(lda,temp)

   !for surface modification
    omega_hs=(omega+eta*CI)*Cone-H00
    omega_hs=omega_hs-cmplx(sfE,0.d0,dp)*Cone

    Sigma=(0.d0,0.d0)
    call mat_mul(lda,H01,temp,Sigma)
    call mat_mul(lda,Sigma,conjg(transpose(H01)),temp)
    temp=omega_hs-temp
    call inv(lda,temp)
   !for surface modification


G00=temp
 ttmmpp=0.d0
   do j=1,lda
     ttmmpp=ttmmpp+G00(j,j)
   enddo

rho = -aimag(ttmmpp)

!print*, rho, ttmmpp
!rho(1)=-aimag(ttmmpp)
!  temp=0.d0
!  call mat_mul(lda,G00,sigmax,temp)
!   ttmmpp=0.d0
!     do j=1,lda
!       ttmmpp=ttmmpp+temp(j,j)
!     enddo
!rho(2)=-aimag(ttmmpp)
!  call mat_mul(lda,G00,sigmay,temp)
!   ttmmpp=0.d0
!     do j=1,lda
!       ttmmpp=ttmmpp+temp(j,j)
!     enddo
!rho(3)=-aimag(ttmmpp)
!  call mat_mul(lda,G00,sigmaz,temp)
!   ttmmpp=0.d0
!     do j=1,lda
!       ttmmpp=ttmmpp+temp(j,j)
!     enddo
!rho(4)=-aimag(ttmmpp)
End SUBROUTINE

