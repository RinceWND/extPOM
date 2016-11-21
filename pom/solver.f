! solver.f

! main subroutines for POM

!_______________________________________________________________________
      subroutine advave
! calculate horizontal advection and diffusion
      implicit none
      include 'pom.h'
      double precision curv2d(im,jm)
      integer i,j

! u-advection and diffusion

! advective fluxes
      advua  = 0.
      fluxua = 0.
      fluxva = 0.

      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=.125d0*((d(i+1,j)+d(i,j))*ua(i+1,j)
     $                       +(d(i,j)+d(i-1,j))*ua(i,j))
     $                      *(ua(i+1,j)+ua(i,j))
        end do
      end do

      do j=2,jm
        do i=2,im
          fluxva(i,j)=.125d0*((d(i,j)+d(i,j-1))*va(i,j)
     $                       +(d(i-1,j)+d(i-1,j-1))*va(i-1,j))
     $                      *(ua(i,j)+ua(i,j-1))
        end do
      end do

! add viscous fluxes
      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=fluxua(i,j)
     $                 -d(i,j)*2.d0*aam2d(i,j)*(uab(i+1,j)-uab(i,j))
     $                   /dx(i,j)
        end do
      end do

      do j=2,jm
        do i=2,im
          tps(i,j)=.25d0*(d(i,j)+d(i-1,j)+d(i,j-1)+d(i-1,j-1))
     $              *(aam2d(i,j)+aam2d(i,j-1)
     $                +aam2d(i-1,j)+aam2d(i-1,j-1))
     $              *((uab(i,j)-uab(i,j-1))
     $                 /(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
     $               +(vab(i,j)-vab(i-1,j))
     $                 /(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1)))
          fluxua(i,j)=fluxua(i,j)*dy(i,j)
          fluxva(i,j)=(fluxva(i,j)-tps(i,j))*.25d0
     $                 *(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1))
        end do
      end do

      call exchange2d_mpi(fluxua,im_local,jm_local)
      call exchange2d_mpi(fluxva,im_local,jm_local)

      do j=2,jmm1
        do i=2,imm1
          advua(i,j)=fluxua(i,j)-fluxua(i-1,j)
     $                +fluxva(i,j+1)-fluxva(i,j)
        end do
      end do

      call exchange2d_mpi(advua,im_local,jm_local)

! v-advection and diffusion
      advva  = 0.
      fluxua = 0.
      fluxva = 0.

! advective fluxes
      do j=2,jm
        do i=2,im
          fluxua(i,j)=.125d0*((d(i,j)+d(i-1,j))*ua(i,j)
     $                       +(d(i,j-1)+d(i-1,j-1))*ua(i,j-1))
     $                      *(va(i-1,j)+va(i,j))
        end do
      end do

      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=.125d0*((d(i,j+1)+d(i,j))*va(i,j+1)
     $                       +(d(i,j)+d(i,j-1))*va(i,j))
     $                      *(va(i,j+1)+va(i,j))
        end do
      end do

! add viscous fluxes
      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=fluxva(i,j)
     $                 -d(i,j)*2.d0*aam2d(i,j)*(vab(i,j+1)-vab(i,j))
     $                   /dy(i,j)
        end do
      end do

      do j=2,jm
        do i=2,im
          fluxva(i,j)=fluxva(i,j)*dx(i,j)
          fluxua(i,j)=(fluxua(i,j)-tps(i,j))*.25d0
     $                 *(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
        end do
      end do

      call exchange2d_mpi(fluxua,im_local,jm_local)
      call exchange2d_mpi(fluxva,im_local,jm_local)

      do j=2,jmm1
        do i=2,imm1
          advva(i,j)=fluxua(i+1,j)-fluxua(i,j)
     $                +fluxva(i,j)-fluxva(i,j-1)
        end do
      end do

      call exchange2d_mpi(advva,im_local,jm_local)

      if(mode.eq.2) then

        do j=2,jmm1
          do i=2,imm1
            wubot(i,j)=-0.5d0*(cbc(i,j)+cbc(i-1,j))
     $                  *sqrt(uab(i,j)**2
     $                        +(.25d0*(vab(i,j)+vab(i,j+1)
     $                                 +vab(i-1,j)+vab(i-1,j+1)))**2)
     $                  *uab(i,j)
          end do
        end do

        do j=2,jmm1
          do i=2,imm1
            wvbot(i,j)=-0.5d0*(cbc(i,j)+cbc(i,j-1))
     $                  *sqrt(vab(i,j)**2
     $                        +(.25d0*(uab(i,j)+uab(i+1,j)
     $                                +uab(i,j-1)+uab(i+1,j-1)))**2)
     $                  *vab(i,j)
          end do
        end do

        do j=2,jmm1
          do i=2,imm1
            curv2d(i,j)=.25d0
     $                   *((va(i,j+1)+va(i,j))*(dy(i+1,j)-dy(i-1,j))
     $                    -(ua(i+1,j)+ua(i,j))*(dx(i,j+1)-dx(i,j-1)))
     $                   /(dx(i,j)*dy(i,j))
          end do
        end do
        call exchange2d_mpi(curv2d,im,jm)

        do j=2,jmm1
          if(n_west.eq.-1) then
          do i=3,imm1
            advua(i,j)=advua(i,j)-aru(i,j)*.25d0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(va(i,j+1)+va(i,j))
     $                    +curv2d(i-1,j)*d(i-1,j)
     $                    *(va(i-1,j+1)+va(i-1,j)))
          end do
          else
          do i=2,imm1
            advua(i,j)=advua(i,j)-aru(i,j)*.25d0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(va(i,j+1)+va(i,j))
     $                    +curv2d(i-1,j)*d(i-1,j)
     $                    *(va(i-1,j+1)+va(i-1,j)))
          end do
          end if
        end do

        do i=2,imm1
          if(n_south.eq.-1) then
          do j=3,jmm1
            advva(i,j)=advva(i,j)+arv(i,j)*.25d0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(ua(i+1,j)+ua(i,j))
     $                    +curv2d(i,j-1)*d(i,j-1)
     $                    *(ua(i+1,j-1)+ua(i,j-1)))
          end do
          else
          do j=2,jmm1
            advva(i,j)=advva(i,j)+arv(i,j)*.25d0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(ua(i+1,j)+ua(i,j))
     $                    +curv2d(i,j-1)*d(i,j-1)
     $                    *(ua(i+1,j-1)+ua(i,j-1)))
          end do
          end if
        end do

      end if

      return
      end

!_______________________________________________________________________
      subroutine advct
! calculate the horizontal portions of momentum advection well in
! advance of their use in advu and advv so that their vertical integrals
! (created in the main program) may be used in the external (2-D) mode
! calculation
      implicit none
      include 'pom.h'
      double precision xflux(im,jm,kb),yflux(im,jm,kb)
      double precision curv(im,jm,kb)
      double precision dtaam
      integer i,j,k

      curv  = 0.
      advx  = 0.
      xflux = 0.
      yflux = 0.

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            curv(i,j,k)=.25d0*((v(i,j+1,k)+v(i,j,k))
     $                         *(dy(i+1,j)-dy(i-1,j))
     $                         -(u(i+1,j,k)+u(i,j,k))
     $                         *(dx(i,j+1)-dx(i,j-1)))
     $                       /(dx(i,j)*dy(i,j))
          end do
        end do
      end do
      call exchange3d_mpi(curv(:,:,1:kbm1),im,jm,kbm1)

! calculate x-component of velocity advection

! calculate horizontal advective fluxes
      do k=1,kbm1
        do j=1,jm
          do i=2,imm1
            xflux(i,j,k)=.125d0*((dt(i+1,j)+dt(i,j))*u(i+1,j,k)
     $                           +(dt(i,j)+dt(i-1,j))*u(i,j,k))
     $                         *(u(i+1,j,k)+u(i,j,k))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.125d0*((dt(i,j)+dt(i,j-1))*v(i,j,k)
     $                           +(dt(i-1,j)+dt(i-1,j-1))*v(i-1,j,k))
     $                         *(u(i,j,k)+u(i,j-1,k))
          end do
        end do
      end do

! add horizontal diffusive fluxes
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            xflux(i,j,k)=xflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.d0
     $                    *(ub(i+1,j,k)-ub(i,j,k))/dx(i,j)
            dtaam=.25d0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $             *(aam(i,j,k)+aam(i-1,j,k)
     $               +aam(i,j-1,k)+aam(i-1,j-1,k))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))
     $                            /(dy(i,j)+dy(i-1,j)
     $                              +dy(i,j-1)+dy(i-1,j-1))
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                            /(dx(i,j)+dx(i-1,j)
     $                              +dx(i,j-1)+dx(i-1,j-1)))

            xflux(i,j,k)=dy(i,j)*xflux(i,j,k)
            yflux(i,j,k)=.25d0*(dx(i,j)+dx(i-1,j)
     $                          +dx(i,j-1)+dx(i-1,j-1))*yflux(i,j,k)
          end do
        end do
      end do

      call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1)

! do horizontal advection
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advx(i,j,k)=xflux(i,j,k)-xflux(i-1,j,k)
     $                   +yflux(i,j+1,k)-yflux(i,j,k)
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          if(n_west.eq.-1) then
          do i=3,imm1
            advx(i,j,k)=advx(i,j,k)
     $                   -aru(i,j)*.25d0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(v(i,j+1,k)+v(i,j,k))
     $                       +curv(i-1,j,k)*dt(i-1,j)
     $                        *(v(i-1,j+1,k)+v(i-1,j,k)))
          end do
          else
          do i=2,imm1
            advx(i,j,k)=advx(i,j,k)
     $                   -aru(i,j)*.25d0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(v(i,j+1,k)+v(i,j,k))
     $                       +curv(i-1,j,k)*dt(i-1,j)
     $                        *(v(i-1,j+1,k)+v(i-1,j,k)))
          end do
          end if
        end do
      end do

      call exchange3d_mpi(advx,im_local,jm_local,kb)

! calculate y-component of velocity advection

      advy  = 0.
      xflux = 0.
      yflux = 0.

! calculate horizontal advective fluxes
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.125d0*((dt(i,j)+dt(i-1,j))*u(i,j,k)
     $                           +(dt(i,j-1)+dt(i-1,j-1))*u(i,j-1,k))
     $                         *(v(i,j,k)+v(i-1,j,k))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=1,im
            yflux(i,j,k)=.125d0*((dt(i,j+1)+dt(i,j))*v(i,j+1,k)
     $                           +(dt(i,j)+dt(i,j-1))*v(i,j,k))
     $                         *(v(i,j+1,k)+v(i,j,k))
          end do
        end do
      end do

! add horizontal diffusive fluxes
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            dtaam=.25d0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $             *(aam(i,j,k)+aam(i-1,j,k)
     $               +aam(i,j-1,k)+aam(i-1,j-1,k))
            xflux(i,j,k)=xflux(i,j,k)
     $                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))
     $                            /(dy(i,j)+dy(i-1,j)
     $                              +dy(i,j-1)+dy(i-1,j-1))
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                            /(dx(i,j)+dx(i-1,j)
     $                              +dx(i,j-1)+dx(i-1,j-1)))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.d0
     $                    *(vb(i,j+1,k)-vb(i,j,k))/dy(i,j)

            xflux(i,j,k)=.25d0*(dy(i,j)+dy(i-1,j)
     $                          +dy(i,j-1)+dy(i-1,j-1))*xflux(i,j,k)
            yflux(i,j,k)=dx(i,j)*yflux(i,j,k)
          end do
        end do
      end do

      call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1)

! do horizontal advection
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advy(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                   +yflux(i,j,k)-yflux(i,j-1,k)
          end do
        end do
      end do

      do k=1,kbm1
        do i=2,imm1
          if(n_south.eq.-1) then
          do j=3,jmm1
            advy(i,j,k)=advy(i,j,k)
     $                   +arv(i,j)*.25d0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(u(i+1,j,k)+u(i,j,k))
     $                       +curv(i,j-1,k)*dt(i,j-1)
     $                        *(u(i+1,j-1,k)+u(i,j-1,k)))
          end do
          else
          do j=2,jmm1
            advy(i,j,k)=advy(i,j,k)
     $                   +arv(i,j)*.25d0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(u(i+1,j,k)+u(i,j,k))
     $                       +curv(i,j-1,k)*dt(i,j-1)
     $                        *(u(i+1,j-1,k)+u(i,j-1,k)))
          end do
          end if
        end do
      end do

      call exchange3d_mpi(advy,im_local,jm_local,kb)

      return
      end

!_______________________________________________________________________
      subroutine advq(qb,q,qf)
! calculates horizontal advection and diffusion, and vertical advection
! for turbulent quantities
      implicit none
      include 'pom.h'
      double precision qb(im_local,jm_local,kb),q(im_local,jm_local,kb)
      double precision qf(im_local,jm_local,kb)
      double precision xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k

      xflux = 0.
      yflux = 0.

! do horizontal advection
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.125d0*(q(i,j,k)+q(i-1,j,k))
     $                    *(dt(i,j)+dt(i-1,j))*(u(i,j,k)+u(i,j,k-1))
            yflux(i,j,k)=.125d0*(q(i,j,k)+q(i,j-1,k))
     $                    *(dt(i,j)+dt(i,j-1))*(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do

! do horizontal diffusion
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.25d0*(aam(i,j,k)+aam(i-1,j,k)
     $                            +aam(i,j,k-1)+aam(i-1,j,k-1))
     $                          *(h(i,j)+h(i-1,j))
     $                          *(qb(i,j,k)-qb(i-1,j,k))*dum(i,j)
     $                          /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.25d0*(aam(i,j,k)+aam(i,j-1,k)
     $                            +aam(i,j,k-1)+aam(i,j-1,k-1))
     $                          *(h(i,j)+h(i,j-1))
     $                          *(qb(i,j,k)-qb(i,j-1,k))*dvm(i,j)
     $                          /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5d0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5d0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do

      call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1)
      call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1)

! do vertical advection, add flux terms, then step forward in time
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            qf(i,j,k)=(w(i,j,k-1)*q(i,j,k-1)-w(i,j,k+1)*q(i,j,k+1))
     $                 *art(i,j)/(dz(k)+dz(k-1))
     $                 +xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
            qf(i,j,k)=((h(i,j)+etb(i,j))*art(i,j)
     $                 *qb(i,j,k)-dti2*qf(i,j,k))
     $                /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine advt1(fb,f,fclim,ff)
! integrate conservative scalar equations
! this is centred scheme, as originally provide in POM (previously
! called advt)
      implicit none
      include 'pom.h'
      double precision fb(im_local,jm_local,kb),f(im_local,jm_local,kb)
      double precision fclim(im_local,jm_local,kb),
     $                 ff(im_local,jm_local,kb)
      double precision xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k

      xflux = 0.
      yflux = 0.

      f(:,:,kb)  = f(:,:,kbm1)
      fb(:,:,kb) = fb(:,:,kbm1)

! do advective fluxes
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25d0*((dt(i,j)+dt(i-1,j))
     $                          *(f(i,j,k)+f(i-1,j,k))*u(i,j,k))
            yflux(i,j,k)=.25d0*((dt(i,j)+dt(i,j-1))
     $                          *(f(i,j,k)+f(i,j-1,k))*v(i,j,k))
          end do
        end do
      end do

! add diffusive fluxes
      fb = fb-fclim

      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.5d0*(aam(i,j,k)+aam(i-1,j,k))
     $                         *(h(i,j)+h(i-1,j))*tprni
     $                         *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                         /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.5d0*(aam(i,j,k)+aam(i,j-1,k))
     $                         *(h(i,j)+h(i,j-1))*tprni
     $                         *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                         /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5d0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5d0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do

      fb = fb+fclim

! do vertical advection
      do j=2,jmm1
        do i=2,imm1
          zflux(i,j,1)=f(i,j,1)*w(i,j,1)*art(i,j)
          zflux(i,j,kb)=0.d0
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,k)=.5d0*(f(i,j,k-1)+f(i,j,k))*w(i,j,k)*art(i,j)
          end do
        end do
      end do

! add net horizontal fluxes and then step forward in time
!      do j=2,jmm1
!        do i=2,imm1
!          do k=1,kbm1
!            ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
!     $                 +yflux(i,j+1,k)-yflux(i,j,k)
!     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
!            ff(i,j,k)=(fb(i,j,k)*dble((h(i,j)+etb(i,j))*art(i,j))
!     $                 -dti2*ff(i,j,k))/dble((h(i,j)+etf(i,j))*art(i,j))
!          end do
!        end do
!      end do
      do k=1,kbm1
        ff(2:imm1,2:jmm1,k)=xflux(3:im,2:jmm1,k)-xflux(2:imm1,2:jmm1,k)
     $          +yflux(2:imm1,3:jm,k)-yflux(2:imm1,2:jmm1,k)
     $          +(zflux(2:imm1,2:jmm1,k)-zflux(2:imm1,2:jmm1,k+1))/dz(k)
        ff(2:imm1,2:jmm1,k)=(fb(2:imm1,2:jmm1,k)
     $         *(h(2:imm1,2:jmm1)+etb(2:imm1,2:jmm1))*art(2:imm1,2:jmm1)
     $         -dti2*ff(2:imm1,2:jmm1,k))
     $         /((h(2:imm1,2:jmm1)+etf(2:imm1,2:jmm1))
     $         *art(2:imm1,2:jmm1))
      end do

      return
      end

!_______________________________________________________________________
      subroutine advt2(fb,f,fclim,ff)
! integrate conservative scalar equations
! this is a first-order upstream scheme, which reduces implicit
! diffusion using the Smolarkiewicz iterative upstream scheme with an
! antidiffusive velocity
! it is based on the subroutines of Gianmaria Sannino (Inter-university
! Computing Consortium, Rome, Italy) and Vincenzo Artale (Italian
! National Agency for New Technology and Environment, Rome, Italy)
      implicit none
      include 'pom.h'
      double precision fb(im_local,jm_local,kb),f(im_local,jm_local,kb)
      double precision fclim(im_local,jm_local,kb),
     $                 ff(im_local,jm_local,kb)
      double precision xflux(im,jm,kb),yflux(im,jm,kb)
      double precision fbmem(im,jm,kb),eta(im,jm)
      double precision xmassflux(im,jm,kb),ymassflux(im,jm,kb),
     $                 zwflux(im,jm,kb)
      integer i,j,k,itera

! calculate horizontal mass fluxes
      xflux = 0.
      yflux = 0.
      xmassflux = 0.
      ymassflux = 0.

      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            xmassflux(i,j,k)=0.25d0*(dy(i-1,j)+dy(i,j))
     $                             *(dt(i-1,j)+dt(i,j))*u(i,j,k)
          end do
        end do

        do j=2,jm
          do i=2,imm1
            ymassflux(i,j,k)=0.25d0*(dx(i,j-1)+dx(i,j))
     $                             *(dt(i,j-1)+dt(i,j))*v(i,j,k)
          end do
        end do
      end do

      fb(:,:,kb) = fb(:,:,kbm1)
      eta = etb(1:im,1:jm)

      zwflux = w(1:im,1:jm,:)
      fbmem  = fb(1:im,1:jm,:)

! start Smolarkiewicz scheme
      do itera=1,nitera

! upwind advection scheme
        do k=1,kbm1
          do j=2,jm
            do i=2,im
              xflux(i,j,k)=0.5d0
     $                      *((xmassflux(i,j,k)+abs(xmassflux(i,j,k)))
     $                        *fbmem(i-1,j,k)+
     $                        (xmassflux(i,j,k)-abs(xmassflux(i,j,k)))
     $                        *fbmem(i,j,k))

              yflux(i,j,k)=0.5d0
     $                      *((ymassflux(i,j,k)+abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j-1,k)+
     $                        (ymassflux(i,j,k)-abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j,k))
            end do
          end do
        end do

        zflux(2:imm1,2:jmm1,1) = 0.
        if(itera.eq.1) then
          zflux(2:imm1,2:jmm1,1) = w(2:imm1,2:jmm1,1)*f(2:imm1,2:jmm1,1)
     $                            *art(2:imm1,2:jmm1)
        end if
        zflux(2:imm1,2:jmm1,kb) = 0.

        do k=2,kbm1
          do j=2,jmm1
            do i=2,imm1
              zflux(i,j,k)=0.5d0
     $                      *((zwflux(i,j,k)+abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k)+
     $                        (zwflux(i,j,k)-abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k-1))
              zflux(i,j,k)=zflux(i,j,k)*art(i,j)
            end do
          end do
        end do

! add net advective fluxes and step forward in time
      do j=2,jmm1
        do i=2,imm1
          do k=1,kbm1
              ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
              ff(i,j,k)=(fbmem(i,j,k)*dble((h(i,j)+eta(i,j))*art(i,j))
     $                 -dti2*ff(i,j,k))/dble((h(i,j)+etf(i,j))*art(i,j))
            end do
          end do
        end do
        ! next line added on 22-Jul-2009 by Raffaele Bernardello
        call exchange3d_mpi(ff(:,:,1:kbm1),im_local,jm_local,kbm1)

! calculate antidiffusion velocity
        call smol_adif(xmassflux,ymassflux,zwflux,ff)

        eta = etf(1:im,1:jm)
        fbmem = ff(1:im,1:jm,:)

! end of Smolarkiewicz scheme
      end do

! add horizontal diffusive fluxes
      fb = fb-fclim

      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xmassflux(i,j,k)=0.5d0*(aam(i,j,k)+aam(i-1,j,k))
            ymassflux(i,j,k)=0.5d0*(aam(i,j,k)+aam(i,j-1,k))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,im
           xflux(i,j,k)=-xmassflux(i,j,k)*(h(i,j)+h(i-1,j))*tprni
     $                   *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                   *(dy(i,j)+dy(i-1,j))*0.5d0/(dx(i,j)+dx(i-1,j))
           yflux(i,j,k)=-ymassflux(i,j,k)*(h(i,j)+h(i,j-1))*tprni
     $                   *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                   *(dx(i,j)+dx(i,j-1))*0.5d0/(dy(i,j)+dy(i,j-1))
          end do
        end do
      end do

      fb = fb+fclim

! add net horizontal fluxes and step forward in time
      do j=2,jmm1
        do i=2,imm1
          do k=1,kbm1
            ff(i,j,k)=ff(i,j,k)-dti2*(xflux(i+1,j,k)-xflux(i,j,k)
     $                               +yflux(i,j+1,k)-yflux(i,j,k))
     $                           /dble((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do

      call exchange3d_mpi(ff(:,:,1:kbm1),im_local,jm_local,kbm1)

      return
      end

!_______________________________________________________________________
      subroutine advu
! do horizontal and vertical advection of u-momentum, and includes
! coriolis, surface slope and baroclinic terms
      implicit none
      include 'pom.h'
      integer i,j,k

! do vertical advection
      uf = 0.

      do k=2,kbm1
        do j=1,jm
          do i=2,im
            uf(i,j,k)=.25d0*(w(i,j,k)+w(i-1,j,k))
     $                     *(u(i,j,k)+u(i,j,k-1))
          end do
        end do
      end do

! combine horizontal and vertical advection with coriolis, surface
! slope and baroclinic terms
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=advx(i,j,k)
     $                 +(uf(i,j,k)-uf(i,j,k+1))*aru(i,j)/dz(k)
     $                 -aru(i,j)*.25d0
     $                   *(cor(i,j)*dt(i,j)
     $                      *(v(i,j+1,k)+v(i,j,k))
     $                     +cor(i-1,j)*dt(i-1,j)
     $                       *(v(i-1,j+1,k)+v(i-1,j,k)))
     $                 +grav*.125d0*(dt(i,j)+dt(i-1,j))
     $                   *(egf(i,j)-egf(i-1,j)+egb(i,j)-egb(i-1,j)
     $                     +(e_atmos(i,j)-e_atmos(i-1,j))*2.d0)
     $                   *(dy(i,j)+dy(i-1,j))
     $                 +drhox(i,j,k)
          end do
        end do
      end do

!  step forward in time
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=((h(i,j)+etb(i,j)+h(i-1,j)+etb(i-1,j))
     $                 *aru(i,j)*ub(i,j,k)
     $                 -2.d0*dti2*uf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))
     $                  *aru(i,j))
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine advv
! do horizontal and vertical advection of v-momentum, and includes
! coriolis, surface slope and baroclinic terms
      implicit none
      include 'pom.h'
      integer i,j,k

! do vertical advection
      vf = 0.

      do k=2,kbm1
        do j=2,jm
          do i=1,im
            vf(i,j,k)=.25d0*(w(i,j,k)+w(i,j-1,k))
     $                     *(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do

! combine horizontal and vertical advection with coriolis, surface
! slope and baroclinic terms
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=advy(i,j,k)
     $                 +(vf(i,j,k)-vf(i,j,k+1))*arv(i,j)/dz(k)
     $                 +arv(i,j)*.25d0
     $                   *(cor(i,j)*dt(i,j)
     $                      *(u(i+1,j,k)+u(i,j,k))
     $                     +cor(i,j-1)*dt(i,j-1)
     $                       *(u(i+1,j-1,k)+u(i,j-1,k)))
     $                 +grav*.125d0*(dt(i,j)+dt(i,j-1))
     $                   *(egf(i,j)-egf(i,j-1)+egb(i,j)-egb(i,j-1)
     $                     +(e_atmos(i,j)-e_atmos(i,j-1))*2.d0)
     $                   *(dx(i,j)+dx(i,j-1))
     $                 +drhoy(i,j,k)
          end do
        end do
      end do

! step forward in time
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=((h(i,j)+etb(i,j)+h(i,j-1)+etb(i,j-1))
     $                 *arv(i,j)*vb(i,j,k)
     $                 -2.d0*dti2*vf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
     $                  *arv(i,j))
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine baropg
! calculate  baroclinic pressure gradient
      implicit none
      include 'pom.h'
      integer i,j,k

      rho = rho-rmean
      
! calculate x-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhox(i,j,1)=.5d0*grav*(-zz(1))*(dt(i,j)+dt(i-1,j))
     $                  *(rho(i,j,1)-rho(i-1,j,1))
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=drhox(i,j,k-1)
     $                    +grav*.25d0*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i-1,j))
     $                      *(rho(i,j,k)-rho(i-1,j,k)
     $                        +rho(i,j,k-1)-rho(i-1,j,k-1))
     $                    +grav*.25d0*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i-1,j))
     $                      *(rho(i,j,k)+rho(i-1,j,k)
     $                        -rho(i,j,k-1)-rho(i-1,j,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=.25d0*(dt(i,j)+dt(i-1,j))
     $                        *drhox(i,j,k)*dum(i,j)
     $                        *(dy(i,j)+dy(i-1,j))
          end do
        end do
      end do

!      call exchange3d_mpi(drhox,im_local,jm_local,kb)

! calculate y-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhoy(i,j,1)=.5d0*grav*(-zz(1))*(dt(i,j)+dt(i,j-1))
     $                  *(rho(i,j,1)-rho(i,j-1,1))
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)
     $                    +grav*.25d0*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i,j-1))
     $                      *(rho(i,j,k)-rho(i,j-1,k)
     $                        +rho(i,j,k-1)-rho(i,j-1,k-1))
     $                    +grav*.25d0*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i,j-1))
     $                      *(rho(i,j,k)+rho(i,j-1,k)
     $                        -rho(i,j,k-1)-rho(i,j-1,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=.25d0*(dt(i,j)+dt(i,j-1))
     $                        *drhoy(i,j,k)*dvm(i,j)
     $                        *(dx(i,j)+dx(i,j-1))
          end do
        end do
      end do

!      call exchange3d_mpi(drhoy,im_local,jm_local,kb)

      do k=1,kb
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=ramp*drhox(i,j,k)
            drhoy(i,j,k)=ramp*drhoy(i,j,k)
          end do
        end do
      end do

      rho = rho+rmean

      return
      end

!_______________________________________________________________________
      subroutine baropg_mcc
! calculate  baroclinic pressure gradient
! 4th order correction terms, following McCalpin
      implicit none
      include 'pom.h'
      integer i,j,k
      double precision d4(im,jm),ddx(im,jm),drho(im,jm,kb),
     $                 rhou(im,jm,kb)
      double precision rho4th(0:im_local,0:jm_local,kb)
     $                ,d4th(0:im_local,0:jm_local)

      rho = rho-rmean

! convert a 2nd order matrices to special 4th order
! special 4th order case
      call order2d_mpi(d,d4th,im_local,jm_local)
      call order3d_mpi(rho,rho4th,im_local,jm_local,kb)

! compute terms correct to 4th order
      ddx = 0.
      d4  = 0.
      rhou = 0.
      drho = 0.

! compute DRHO, RHOU, DDX and D4
      do j=1,jm
        do i=2,im
          do k=1,kbm1
            drho(i,j,k)=(rho(i,j,k)-rho(i-1,j,k))*dum(i,j)
            rhou(i,j,k)=0.5*(rho(i,j,k)+rho(i-1,j,k))*dum(i,j)
          end do
          ddx(i,j)=(d(i,j)-d(i-1,j))*dum(i,j)
          d4(i,j)=.5*(d(i,j)+d(i-1,j))*dum(i,j)
        end do
      end do

      if(n_west.eq.-1) then
        do j=1,jm
          do i=3,imm1
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k) - (1./24.)*
     $                    (dum(i+1,j)*(rho(i+1,j,k)-rho(i,j,k))-
     $                    2*(rho(i,j,k)-rho(i-1,j,k))+
     $                    dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
              rhou(i,j,k)=rhou(i,j,k) + (1./16.)*
     $                    (dum(i+1,j)*(rho(i,j,k)-rho(i+1,j,k))+
     $                    dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24.)*
     $               (dum(i+1,j)*(d(i+1,j)-d(i,j))-
     $               2*(d(i,j)-d(i-1,j))+
     $               dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dum(i+1,j)*(d(i,j)-d(i+1,j))+
     $              dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
          end do
        end do
      else
        do j=1,jm
          do i=2,imm1
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k) - (1./24.)*
     $                   (dum(i+1,j)*(rho(i+1,j,k)-rho(i,j,k))-
     $                    2*(rho(i,j,k)-rho(i-1,j,k))+
     $                    dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)))
              rhou(i,j,k)=rhou(i,j,k) + (1./16.)*
     $                    (dum(i+1,j)*(rho(i,j,k)-rho(i+1,j,k))+
     $                    dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24.)*
     $               (dum(i+1,j)*(d(i+1,j)-d(i,j))-
     $               2*(d(i,j)-d(i-1,j))+
     $               dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dum(i+1,j)*(d(i,j)-d(i+1,j))+
     $              dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)))
          end do
        end do
      end if

! calculate x-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhox(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1)
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=drhox(i,j,k-1)
     $                   +grav*0.5d0*dzz(k-1)*d4(i,j)
     $                   *(drho(i,j,k-1)+drho(i,j,k))
     $                   +grav*0.5d0*(zz(k-1)+zz(k))*ddx(i,j)
     $                   *(rhou(i,j,k)-rhou(i,j,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=.25d0*(dt(i,j)+dt(i-1,j))
     $                        *drhox(i,j,k)*dum(i,j)
     $                        *(dy(i,j)+dy(i-1,j))
          end do
        end do
      end do

!      call exchange3d_mpi(drhox,im_local,jm_local,kb)

! compute terms correct to 4th order
      ddx = 0.
      d4  = 0.
      rhou = 0.
      drho = 0.

! compute DRHO, RHOU, DDX and D4
      do j=2,jm
        do i=1,im
          do k=1,kbm1
            drho(i,j,k)=(rho(i,j,k)-rho(i,j-1,k))*dvm(i,j)
            rhou(i,j,k)=.5*(rho(i,j,k)+rho(i,j-1,k))*dvm(i,j)
          end do
          ddx(i,j)=(d(i,j)-d(i,j-1))*dvm(i,j)
          d4(i,j)=.5*(d(i,j)+d(i,j-1))*dvm(i,j)
        end do
      end do

      if(n_south.eq.-1) then
        do j=3,jmm1
          do i=1,im
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k)-(1./24.)*
     $                    (dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j,k))-
     $                    2*(rho(i,j,k)-rho(i,j-1,k))+
     $                    dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
              rhou(i,j,k)=rhou(i,j,k)+(1./16.)*
     $                    (dvm(i,j+1)*(rho(i,j,k)-rho(i,j+1,k))+
     $                    dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24)*
     $               (dvm(i,j+1)*(d(i,j+1)-d(i,j))-
     $               2*(d(i,j)-d(i,j-1))+
     $               dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dvm(i,j+1)*(d(i,j)-d(i,j+1))+
     $              dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
          end do
        end do
      else
        do j=2,jmm1
          do i=1,im
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k)-(1./24.)*
     $                    (dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j,k))-
     $                    2*(rho(i,j,k)-rho(i,j-1,k))+
     $                    dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)))
              rhou(i,j,k)=rhou(i,j,k)+(1./16.)*
     $                    (dvm(i,j+1)*(rho(i,j,k)-rho(i,j+1,k))+
     $                    dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24)*
     $               (dvm(i,j+1)*(d(i,j+1)-d(i,j))-
     $               2*(d(i,j)-d(i,j-1))+
     $               dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dvm(i,j+1)*(d(i,j)-d(i,j+1))+
     $              dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)))
          end do
        end do
      end if

! calculate y-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhoy(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1)
        end do
      end do
      
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)
     $                   +grav*0.5d0*dzz(k-1)*d4(i,j)
     $                   *(drho(i,j,k-1)+drho(i,j,k))
     $                   +grav*0.5d0*(zz(k-1)+zz(k))*ddx(i,j)
     $                   *(rhou(i,j,k)-rhou(i,j,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=.25d0*(dt(i,j)+dt(i,j-1))
     $                        *drhoy(i,j,k)*dvm(i,j)
     $                        *(dx(i,j)+dx(i,j-1))
          end do
        end do
      end do

!      call exchange3d_mpi(drhoy,im_local,jm_local,kb)

      do k=1,kb
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=ramp*drhox(i,j,k)
            drhoy(i,j,k)=ramp*drhoy(i,j,k)
          end do
        end do
      end do

      rho = rho+rmean

      return
      end

!_______________________________________________________________________
      subroutine dens(si,ti,rhoo)
! calculate (density-1000.)/rhoref.
! see: Mellor, G.L., 1991, J. Atmos. Oceanic Tech., 609-611
! note: if pressure is not used in dens, buoyancy term (boygr) in
! subroutine profq must be changed (see note in subroutine profq)
      implicit none
      include 'pom.h'
      double precision si(im_local,jm_local,kb),ti(im_local,jm_local,kb)
      double precision rhoo(im_local,jm_local,kb)
      integer i,j,k
      double precision cr,p,rhor,sr,tr,tr2,tr3,tr4

      do k=1,kbm1
        do j=1,jm
          do i=1,im

            tr=ti(i,j,k)+tbias
            sr=si(i,j,k)+sbias
            tr2=tr*tr
            tr3=tr2*tr
            tr4=tr3*tr

! approximate pressure in units of bars
            p=grav*rhoref*(-zz(k)* h(i,j))*1.d-5

            rhor=-0.157406d0+6.793952d-2*tr
     $           -9.095290d-3*tr2+1.001685d-4*tr3
     $           -1.120083d-6*tr4+6.536332d-9*tr4*tr

            rhor=rhor+(0.824493d0-4.0899d-3*tr
     $               +7.6438d-5*tr2-8.2467d-7*tr3
     $               +5.3875d-9*tr4)*sr
     $               +(-5.72466d-3+1.0227d-4*tr
     $               -1.6546d-6*tr2)*abs(sr)**1.5
     $               +4.8314d-4*sr*sr

            cr=1449.1d0+.0821d0*p+4.55d0*tr-.045d0*tr2
     $                 +1.34d0*(sr-35.d0)
            rhor=rhor+1.d5*p/(cr*cr)*(1.d0-2.d0*p/(cr*cr))

            rhoo(i,j,k)=rhor/rhoref*fsm(i,j)

          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine profq
! solve for q2 (twice the turbulent kinetic energy), q2l (q2 x turbulent
! length scale), km (vertical kinematic viscosity) and kh (vertical
! kinematic diffusivity), using a simplified version of the level 2 1/2
! model of Mellor and Yamada (1982)
! in this version, the Craig-Banner sub-model whereby breaking wave tke
! is injected into the surface is included. However, we use an
! analytical solution to the near surface tke equation to solve for q2
! at the surface giving the same result as C-B diffusion. The new scheme
! is simpler and more robust than the latter scheme
      implicit none
      include 'pom.h'
      double precision a(im,jm,kb),c(im,jm,kb)
      double precision ee(im,jm,kb),gg(im,jm,kb)
      double precision sm(im,jm,kb),sh(im,jm,kb)
      double precision cc(im,jm,kb)
      double precision gh(im,jm,kb),boygr(im,jm,kb),
     $                 dh(im,jm),stf(im,jm,kb)
      double precision prod(im,jm,kb)
      double precision a1,a2,b1,b2,c1
      double precision coef1,coef2,coef3,coef4,coef5
      double precision const1,e1,e2,ghc
      double precision p,sef,sp,tp
      double precision l0(im,jm)
      double precision cbcnst,surfl,shiw
      double precision utau2(im,jm)
!      double precision df0,df1,df2
      integer i,j,k,ki

      data a1,b1,a2,b2,c1/0.92d0,16.6d0,0.74d0,10.1d0,0.08d0/
      data e1/1.8d0/,e2/1.33d0/
      data sef/1.d0/
      data cbcnst/100.d0/,surfl/2.d5/,shiw/0.d0/

      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do

      a  = 0.
      c  = 0.
      ee = 0.
      gg = 0.
      utau2 = 0.

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k)=-dti2*(kq(i,j,k+1)+kq(i,j,k)+2.d0*umol)*.5d0
     $                /(dzz(k-1)*dz(k)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kq(i,j,k-1)+kq(i,j,k)+2.d0*umol)*.5d0
     $                /(dzz(k-1)*dz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do

! the following section solves the equation:
!     dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b

! surface and bottom boundary conditions
      const1=(16.6d0**(2.d0/3.d0))*sef

! initialize fields that are not calculated on all boundaries
! but are later used there
      l0    = 0.
      boygr = 0.
      prod  = 0.

      do j=1,jmm1
        do i=1,imm1
          utau2(i,j)=sqrt((.5d0*(wusurf(i,j)+wusurf(i+1,j)))**2
     $                   +(.5d0*(wvsurf(i,j)+wvsurf(i,j+1)))**2)
          uf(i,j,kb)=sqrt((.5d0*(wubot(i,j)+wubot(i+1,j)))**2
     $                   +(.5d0*(wvbot(i,j)+wvbot(i,j+1)))**2)*const1
        end do
      end do
      call exchange2d_mpi(utau2,im,jm)
      call exchange2d_mpi(uf(:,:,kb),im_local,jm_local)

      do j=1,jm
        do i=1,im
! wave breaking energy- a variant of Craig & Banner (1994)
! see Mellor and Blumberg, 2003.
          ee(i,j,1)=0.d0
          gg(i,j,1)=(15.8*cbcnst)**(2./3.)*utau2(i,j)
! surface length scale following Stacey (1999).
          l0(i,j)=surfl*utau2(i,j)/grav
        end do
      end do

! calculate speed of sound squared
      cc = 0.
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tp=t(i,j,k)+tbias
            sp=s(i,j,k)+sbias
! calculate pressure in units of decibars
            p=grav*rhoref*(-zz(k)*h(i,j))*1.d-4
            cc(i,j,k)=1449.1d0+.00821d0*p+4.55d0*tp-.045d0*tp**2
     $                 +1.34d0*(sp-35.0d0)
            cc(i,j,k)=cc(i,j,k)
     $                 /sqrt((1.d0-.01642d0*p/cc(i,j,k))
     $                   *(1.d0-0.40d0*p/cc(i,j,k)**2))
          end do
        end do
      end do

! calculate buoyancy gradient
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            q2b(i,j,k)=abs(q2b(i,j,k))
            q2lb(i,j,k)=abs(q2lb(i,j,k))
            boygr(i,j,k)=grav*(rho(i,j,k-1)-rho(i,j,k))
     $                    /(dzz(k-1)*h(i,j))
! *** note: comment out next line if dens does not include pressure
     $      +(grav**2)*2.d0/(cc(i,j,k-1)**2+cc(i,j,k)**2)
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k))
            if(z(k).gt.-0.5) l(i,j,k)=max(l(i,j,k),kappa*l0(i,j))
!            if (my_task == 0) then
!              write(*,*) i,j,k,":",l(i,j,k), boygr(i,j,k), q2b(i,j,k)
!            end if
            gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)/q2b(i,j,k)
            gh(i,j,k)=min(gh(i,j,k),.028d0)
          end do
        end do
      end do

      do j=1,jm
        do i=1,im
          l(i,j,1)=kappa*l0(i,j)
          l(i,j,kb)=0.d0
          gh(i,j,1)=0.d0
          gh(i,j,kb)=0.d0
        end do
      end do

! calculate production of turbulent kinetic energy:
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            prod(i,j,k)=km(i,j,k)*.25d0*sef
     $                   *((u(i,j,k)-u(i,j,k-1)
     $                      +u(i+1,j,k)-u(i+1,j,k-1))**2
     $                     +(v(i,j,k)-v(i,j,k-1)
     $                      +v(i,j+1,k)-v(i,j+1,k-1))**2)
     $                   /(dzz(k-1)*dh(i,j))**2
! add shear due to internal wave field
     $             -shiw*km(i,j,k)*boygr(i,j,k)
            prod(i,j,k)=prod(i,j,k)+kh(i,j,k)*boygr(i,j,k)
          end do
        end do
      end do
      call exchange3d_mpi(prod(:,:,2:kbm1),im,jm,kbm2)

! note: Richardson # dep. dissipation correction (Mellor, 2001; Ezer,
! 2000), depends on ghc the critical number (empirical -6 to -2) to
! increase mixing
      ghc=-6.0d0
      do k=1,kb
        do j=1,jm
          do i=1,im
            stf(i,j,k)=1.d0
! It is unclear yet if diss. corr. is needed when surf. waves are included.
!           if(gh(i,j,k).lt.0.d0)
!    $        stf(i,j,k)=1.0d0-0.9d0*(gh(i,j,k)/ghc)**1.5d0
!           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1d0
            dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)
     $                   /(b1*l(i,j,k)+small)
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1))
     $                      -(2.d0*dti2*dtef(i,j,k)+1.d0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(-2.d0*dti2*prod(i,j,k)+c(i,j,k)*gg(i,j,k-1)
     $                 -uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do k=1,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
            uf(i,j,ki)=ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do

! the following section solves the equation:
!     dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
      do j=1,jm
        do i=1,im
          vf(i,j,1)=0.
          vf(i,j,kb)=0.
          ee(i,j,2)=0.d0
          gg(i,j,2)=-kappa*z(2)*dh(i,j)*q2(i,j,2)
          vf(i,j,kb-1)=kappa*(1+z(kbm1))*dh(i,j)*q2(i,j,kbm1)
        end do
      end do
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            dtef(i,j,k)=dtef(i,j,k)
     $                   *(1.d0+e2*((1.d0/abs(z(k)-z(1))
     $                               +1.d0/abs(z(k)-z(kb)))
     $                                *l(i,j,k)/(dh(i,j)*kappa))**2)
          end do
        end do
      end do
      do k=3,kbm1
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1))
     $                      -(dti2*dtef(i,j,k)+1.d0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(dti2*(-prod(i,j,k)*l(i,j,k)*e1)
     $                 +c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do k=1,kb-2
        ki=kb-k
        do j=1,jm
          do i=1,im
            vf(i,j,ki)=ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do
! the following is to counter the problem of the ratio of two small
! numbers (l = q2l/q2) or one number becoming negative. Two options are
! included below. In this application, the second option, l was less
! noisy when uf or vf is small
      do k=2,kbm1
        do j=1,jm
          do i=1,im
!           if(uf(i,j,k).le.small.or.vf(i,j,k).le.small) then
!             uf(i,j,k)=small
!             vf(i,j,k)=0.1*dt(i,j)*small
!           end if
          uf(i,j,k)=abs(uf(i,j,k))
          vf(i,j,k)=abs(vf(i,j,k))
          end do
        end do
      end do

! the following section solves for km and kh
      coef4=18.d0*a1*a1+9.d0*a1*a2
      coef5=9.d0*a1*a2

! note that sm and sh limit to infinity when gh approaches 0.0288
      do k=1,kb
        do j=1,jm
          do i=1,im
            coef1=a2*(1.d0-6.d0*a1/b1*stf(i,j,k))
            coef2=3.d0*a2*b2/stf(i,j,k)+18.d0*a1*a2
            coef3=a1*(1.d0-3.d0*c1-6.d0*a1/b1*stf(i,j,k))
            sh(i,j,k)=coef1/(1.d0-coef2*gh(i,j,k))
            sm(i,j,k)=coef3+sh(i,j,k)*coef4*gh(i,j,k)
            sm(i,j,k)=sm(i,j,k)/(1.d0-coef5*gh(i,j,k))
          end do
        end do
      end do

! there are 2 options for kq which, unlike km and kh, was not derived by
! Mellor and Yamada but was purely empirical based on neutral boundary
! layer data. The choice is whether or not it should be subject to the
! stability factor, sh. Generally, there is not a great difference in
! output
      do k=1,kb
        do j=1,jm
          do i=1,im
            prod(i,j,k)=l(i,j,k)*sqrt(abs(q2(i,j,k)))
            kq(i,j,k)=(prod(i,j,k)*.41d0*sh(i,j,k)+kq(i,j,k))*.5d0
!            kq(i,j,k)=(prod(i,j,k)*.20+kq(i,j,k))*.5d0
            km(i,j,k)=(prod(i,j,k)*sm(i,j,k)+km(i,j,k))*.5d0
            kh(i,j,k)=(prod(i,j,k)*sh(i,j,k)+kh(i,j,k))*.5d0
          end do
        end do
      end do

! cosmetics: make boundr. values as interior (even if not used, printout
! may show strange values)
      if(n_north.eq.-1) then
        km(:,jm,:) = km(:,jmm1,:)
        kh(:,jm,:) = kh(:,jmm1,:)
        kq(:,jm,:) = kq(:,jmm1,:)
      end if
      if(n_south.eq.-1) then
        km(:,1,:) = km(:,2,:)
        kh(:,1,:) = kh(:,2,:)
        kq(:,1,:) = kq(:,2,:)
      end if
      if(n_east.eq.-1) then
        km(im,:,:) = km(imm1,:,:)
        kh(im,:,:) = kh(imm1,:,:)
        kq(im,:,:) = kq(imm1,:,:)
      end if
      if(n_west.eq.-1) then
        km(1,:,:) = km(2,:,:)
        kh(1,:,:) = kh(2,:,:)
        kq(1,:,:) = kq(2,:,:)
      end if

      do k=1,kb
        km(:,:,k) = km(:,:,k)*fsm
        kh(:,:,k) = kh(:,:,k)*fsm
        kq(:,:,k) = kq(:,:,k)*fsm
      end do

      return
      end

!_______________________________________________________________________
      subroutine proft(f,wfsurf,fsurf,nbc)
! solves for vertical diffusion of temperature and salinity using method
! described by Richmeyer and Morton (1967)
! note: wfsurf and swrad are negative values when water column is
! warming or salt is being added
      implicit none
      include 'pom.h'
      double precision f(im_local,jm_local,kb)
      double precision wfsurf(im_local,jm_local),
     $                 fsurf(im_local,jm_local)
      integer nbc
      double precision a(im,jm,kb),c(im,jm,kb)
      double precision ee(im,jm,kb),gg(im,jm,kb)
      double precision dh(im,jm),rad(im,jm,kb)
      double precision r(5),ad1(5),ad2(5)
      integer i,j,k,ki

! irradiance parameters after Paulson and Simpson (1977)
!       ntp               1      2       3       4       5
!   Jerlov type           i      ia      ib      ii     iii
      data r   /       .58d0,  .62d0,  .67d0,  .77d0,  .78d0 /
      data ad1 /       .35d0,  .60d0,  1.0d0,  1.5d0,  1.4d0 /
      data ad2 /       23.d0,  20.d0,  17.d0,  14.d0,  7.9d0 /

! surface boundary condition:
!       nbc   prescribed    prescribed   short wave
!             temperature      flux      penetration
!             or salinity               (temperature
!                                           only)
!        1        no           yes           no
!        2        no           yes           yes
!        3        yes          no            no
!        4        yes          no            yes
! note that only 1 and 3 are allowed for salinity

! the following section solves the equation
!     dti2*(kh*f')'-f=-fb
      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do

      a  = 0.
      c  = 0.
      ee = 0.
      gg = 0.

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do

! calculate penetrative radiation. At the bottom any unattenuated
! radiation is deposited in the bottom layer
      rad = 0.

      if(nbc.eq.2.or.nbc.eq.4) then
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              rad(i,j,k)=real(swrad(i,j)
     $                    *(r(ntp)*exp(real(z(k)*dh(i,j)/ad1(ntp),16))
     $              +(1.d0-r(ntp))*exp(real(z(k)*dh(i,j)/ad2(ntp),16)))
     $                        ,8)
            end do
          end do
        end do
      end if

      if(nbc.eq.1) then

        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.d0)
            gg(i,j,1)=dti2*wfsurf(i,j)/(dz(1)*dh(i,j))-f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.d0)
          end do
        end do

      else if(nbc.eq.2) then

        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.d0)
            gg(i,j,1)=dti2*(wfsurf(i,j)+rad(i,j,1)-rad(i,j,2))
     $                 /(dz(1)*dh(i,j))
     $                   -f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.d0)
          end do
        end do

      else if(nbc.eq.3.or.nbc.eq.4) then

        do j=1,jm
          do i=1,im
            ee(i,j,1)=0.d0
            gg(i,j,1)=fsurf(i,j)
          end do
        end do

      end if

      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1))-1.d0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-f(i,j,k)
     $                 +dti2*(rad(i,j,k)-rad(i,j,k+1))
     $                   /(dh(i,j)*dz(k)))
     $                 *gg(i,j,k)
          end do
        end do
      end do

! bottom adiabatic boundary condition
      do j=1,jm
        do i=1,im
          f(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-f(i,j,kbm1)
     $                 +dti2*(rad(i,j,kbm1)-rad(i,j,kb))
     $                   /(dh(i,j)*dz(kbm1)))
     $                 /(c(i,j,kbm1)*(1.d0-ee(i,j,kbm2))-1.d0)
        end do
      end do

      do k=2,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
            f(i,j,ki)=(ee(i,j,ki)*f(i,j,ki+1)+gg(i,j,ki))
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine profu
! solves for vertical diffusion of x-momentum using method described by
! Richmeyer and Morton (1967)
! note: wusurf has the opposite sign to the wind speed
      implicit none
      include 'pom.h'
      double precision a(im,jm,kb),c(im,jm,kb)
      double precision ee(im,jm,kb),gg(im,jm,kb)
      double precision dh(im,jm)
      integer i,j,k,ki

! the following section solves the equation
!   dti2*(km*u')'-u=-ub
      dh = 1.

      do j=2,jm
        do i=2,im
          dh(i,j)=(h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))*.5d0
        end do
      end do

      a  = 0.
      c  = 0.
      ee = 0.
      gg = 0.

      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i-1,j,k))*.5d0
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do

      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.d0)
          gg(i,j,1)=(-dti2*wusurf(i,j)/(-dz(1)*dh(i,j))
     $               -uf(i,j,1))
     $               /(a(i,j,1)-1.d0)
        end do
      end do

      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1))-1.d0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5d0*(cbc(i,j)+cbc(i-1,j))
     $              *sqrt(ub(i,j,kbm1)**2
     $                +(.25d0*(vb(i,j,kbm1)+vb(i,j+1,kbm1)
     $                         +vb(i-1,j,kbm1)+vb(i-1,j+1,kbm1)))**2)
          uf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-uf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.d0
     $                    -(ee(i,j,kbm2)-1.d0)*c(i,j,kbm1))
          uf(i,j,kbm1)=uf(i,j,kbm1)*dum(i,j)
        end do
      end do

      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,ki)=(ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki))*dum(i,j)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          wubot(i,j)=-tps(i,j)*uf(i,j,kbm1)
        end do
      end do
      call exchange2d_mpi(wubot,im_local,jm_local)

      return
      end

!_______________________________________________________________________
      subroutine profv
! solves for vertical diffusion of x-momentum using method described by
! Richmeyer and Morton (1967)
! note: wvsurf has the opposite sign to the wind speed
      implicit none
      include 'pom.h'
      double precision a(im,jm,kb),c(im,jm,kb)
      double precision ee(im,jm,kb),gg(im,jm,kb)
      double precision dh(im,jm)
      integer i,j,k,ki

! the following section solves the equation
!     dti2*(km*u')'-u=-ub

      dh = 1.

      do j=2,jm
        do i=2,im
          dh(i,j)=.5d0*(h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
        end do
      end do

      a  = 0.
      c  = 0.
      ee = 0.
      gg = 0.

      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i,j-1,k))*.5d0
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do

      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.d0)
          gg(i,j,1)=(-dti2*wvsurf(i,j)/(-dz(1)*dh(i,j))-vf(i,j,1))
     $               /(a(i,j,1)-1.d0)
        end do
      end do

      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1))-1.d0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5d0*(cbc(i,j)+cbc(i,j-1))
     $              *sqrt((.25d0*(ub(i,j,kbm1)+ub(i+1,j,kbm1)
     $                            +ub(i,j-1,kbm1)+ub(i+1,j-1,kbm1)))**2
     $                    +vb(i,j,kbm1)**2)
          vf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-vf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.d0
     $                    -(ee(i,j,kbm2)-1.d0)*c(i,j,kbm1))
          vf(i,j,kbm1)=vf(i,j,kbm1)*dvm(i,j)
        end do
      end do

      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,ki)=(ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki))*dvm(i,j)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          wvbot(i,j)=-tps(i,j)*vf(i,j,kbm1)
        end do
      end do
      call exchange2d_mpi(wvbot,im_local,jm_local)

      return
      end

!_______________________________________________________________________
      subroutine smol_adif(xmassflux,ymassflux,zwflux,ff)
! calculate the antidiffusive velocity used to reduce the numerical
! diffusion associated with the upstream differencing scheme
! this is based on a subroutine of Gianmaria Sannino (Inter-university
! Computing Consortium, Rome, Italy) and Vincenzo Artale (Italian
! National Agency for New Technology and Environment, Rome, Italy)
      implicit none
      include 'pom.h'
      double precision ff(im_local,jm_local,kb)
      double precision xmassflux(im,jm,kb),ymassflux(im,jm,kb),
     $                 zwflux(im,jm,kb)
      double precision mol,abs_1,abs_2
      double precision value_min,epsilon
      double precision udx,u2dt,vdy,v2dt,wdz,w2dt
      integer i,j,k
      parameter (value_min=1.d-9,epsilon=1.0d-14)

! apply temperature and salinity mask
      do k=1,kb
        ff(:,:,k) = ff(:,:,k)*fsm
      end do

! recalculate mass fluxes with antidiffusion velocity
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i-1,j,k).lt.value_min) then
              xmassflux(i,j,k)=0.d0
            else
              udx=abs(xmassflux(i,j,k))
              u2dt=dti2*xmassflux(i,j,k)*xmassflux(i,j,k)*2.d0
     $              /(aru(i,j)*(dt(i-1,j)+dt(i,j)))
              mol=(ff(i,j,k)-ff(i-1,j,k))
     $             /(ff(i-1,j,k)+ff(i,j,k)+epsilon)
              xmassflux(i,j,k)=(udx-u2dt)*mol*sw
              abs_1=abs(udx)
              abs_2=abs(u2dt)
              if(abs_1.lt.abs_2) xmassflux(i,j,k)=0.d0
            end if
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j-1,k).lt.value_min) then
              ymassflux(i,j,k)=0.d0
            else
             vdy=abs(ymassflux(i,j,k))
             v2dt=dti2*ymassflux(i,j,k)*ymassflux(i,j,k)*2.d0
     $             /(arv(i,j)*(dt(i,j-1)+dt(i,j)))
             mol=(ff(i,j,k)-ff(i,j-1,k))
     $            /(ff(i,j-1,k)+ff(i,j,k)+epsilon)
             ymassflux(i,j,k)=(vdy-v2dt)*mol*sw
             abs_1=abs(vdy)
             abs_2=abs(v2dt)
             if(abs_1.lt.abs_2) ymassflux(i,j,k)=0.d0
            end if
          end do
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j,k-1).lt.value_min) then
              zwflux(i,j,k)=0.d0
            else
              wdz=abs(zwflux(i,j,k))
              w2dt=dti2*zwflux(i,j,k)*zwflux(i,j,k)/
     $             (dzz(k-1)*dt(i,j))
              mol=(ff(i,j,k-1)-ff(i,j,k))
     $             /(ff(i,j,k)+ff(i,j,k-1)+epsilon)
              zwflux(i,j,k)=(wdz-w2dt)*mol*sw
              abs_1=abs(wdz)
              abs_2=abs(w2dt)
              if(abs_1.lt.abs_2)zwflux(i,j,k)=0.d0
            end if
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine vertvl
! calculates vertical velocity
      implicit none
      include 'pom.h'
      double precision xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k

      xflux = 0.
      yflux = 0.

! reestablish boundary conditions
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25d0*(dy(i,j)+dy(i-1,j))
     $                    *(dt(i,j)+dt(i-1,j))*u(i,j,k)
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.25d0*(dx(i,j)+dx(i,j-1))
     $                    *(dt(i,j)+dt(i,j-1))*v(i,j,k)
          end do
        end do
      end do

! note: if one wishes to include freshwater flux, the surface velocity
! should be set to vflux(i,j). See also change made to 2-D volume
! conservation equation which calculates elf
      do j=2,jmm1
        do i=2,imm1
          w(i,j,1)=0.5*(vfluxb(i,j)+vfluxf(i,j))
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            w(i,j,k+1)=w(i,j,k)
     $                +dz(k)*((xflux(i+1,j,k)-xflux(i,j,k)
     $                        +yflux(i,j+1,k)-yflux(i,j,k))
     $                        /(dx(i,j)*dy(i,j))
     $                        +(etf(i,j)-etb(i,j))/dti2)
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine realvertvl
! calculates real vertical velocity (wr)
      implicit none
      include 'pom.h'
      integer i,j,k
      double precision dxr,dxl,dyt,dyb

      wr = 0.

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tps(i,j)=zz(k)*dt(i,j) + et(i,j)
          end do
        end do
        do j=2,jmm1
          do i=2,imm1
            dxr=2.0/(dx(i+1,j)+dx(i,j))
            dxl=2.0/(dx(i,j)+dx(i-1,j))
            dyt=2.0/(dy(i,j+1)+dy(i,j))
            dyb=2.0/(dy(i,j)+dy(i,j-1))
            wr(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))+0.5*
     $                (u(i+1,j,k)*(tps(i+1,j)-tps(i,j))*dxr+
     $                 u(i,j,k)*(tps(i,j)-tps(i-1,j))*dxl+
     $                 v(i,j+1,k)*(tps(i,j+1)-tps(i,j))*dyt+
     $                 v(i,j,k)*(tps(i,j)-tps(i,j-1))*dyb)
     $                +(1.0+zz(k))*(etf(i,j)-etb(i,j))/dti2
          end do
        end do
      end do

      call exchange3d_mpi(wr(:,:,1:kbm1),im_local,jm_local,kbm1)

      if (n_south.eq.-1) wr(:,1,:) =wr(:,2,:)
      if (n_north.eq.-1) wr(:,jm,:)=wr(:,jmm1,:)
      if (n_west.eq.-1)  wr(1,:,:) =wr(2,:,:)
      if (n_east.eq.-1)  wr(im,:,:)=wr(imm1,:,:)

      do k=1,kbm1
        wr(:,:,k) = fsm(:,:)*wr(:,:,k)
      end do

      return
      end
