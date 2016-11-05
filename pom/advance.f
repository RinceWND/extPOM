! advance.f

! advance POM

!_______________________________________________________________________
      subroutine advance
! advance POM 1 step in time
      implicit none
      include 'pom.h'
      
! get time
      call get_time
      
! set time dependent surface boundary conditions
      call surface_forcing
      
! set time dependent lateral boundary conditions
      call lateral_bc
      
! set lateral viscosity
      call lateral_viscosity
      
! form vertical averages of 3-D fields for use in external (2-D) mode
      call mode_interaction
      
! external (2-D) mode calculation
      do iext=1,isplit
        call mode_external
      end do
      
! internal (3-D) mode calculation
      call mode_internal
      
! print section
      call print_section

! write output
      if(netcdf_file.ne.'nonetcdf' .and. mod(iint,iprint).eq.0)
     $                                         call write_output_pnetcdf

! write restart
      if(mod(iint,irestart).eq.0) call write_restart_pnetcdf

! check CFL condition
      call check_velocity

      return
      end

!_______________________________________________________________________
      subroutine get_time
! return the model time
      implicit none
      include 'pom.h'
      time=dti*float(iint)/86400.d0+time0
      if(iint.ge.iswtch) iprint=nint(prtd2*24.d0*3600.d0/dti)
      if(lramp) then
        ramp=time/period
        if(ramp.gt.1.d0) ramp=1.d0
      else
        ramp=1.d0
      endif
      return
      end

!_______________________________________________________________________
      subroutine surface_forcing
! set time dependent surface boundary conditions
      implicit none
      include 'pom.h'
      integer i,j
      double precision tatm,satm
! wind forcing is suplied in subroutine wind
      call wind
! heat fluxes are suplied in subroutine heat
!      call heat
! water fluxes are suplied in subroutine water
!      call water
      return
      end

!_______________________________________________________________________
      subroutine lateral_viscosity
! set the lateral viscosity
      implicit none
      include 'pom.h'
      integer i,j,k
! if mode=2 then initial values of aam2d are used. If one wishes
! to use Smagorinsky lateral viscosity and diffusion for an
! external (2-D) mode calculation, then appropiate code can be
! adapted from that below and installed just before the end of the
! "if(mode.eq.2)" loop in subroutine advave

! calculate Smagorinsky lateral viscosity:
! ( hor visc = horcon*dx*dy*sqrt((du/dx)**2+(dv/dy)**2
!                                +.5*(du/dy+dv/dx)**2) )
      if(mode.ne.2) then
        call advct

        if (npg.eq.1) then
          call baropg
        else if (npg.eq.2) then
          call baropg_mcc
        else
          error_status=1
          write(6,'(/''Error: invalid value for npg'')')
        end if

        do k=1,kbm1
          do j=2,jmm1
            do i=2,imm1
              aam(i,j,k)=horcon*dx(i,j)*dy(i,j)
     $                    *sqrt( ((u(i+1,j,k)-u(i,j,k))/dx(i,j))**2
     $                          +((v(i,j+1,k)-v(i,j,k))/dy(i,j))**2
     $                    +.5d0*(.25d0*(u(i,j+1,k)+u(i+1,j+1,k)
     $                                 -u(i,j-1,k)-u(i+1,j-1,k))
     $                    /dy(i,j)
     $                    +.25d0*(v(i+1,j,k)+v(i+1,j+1,k)
     $                           -v(i-1,j,k)-v(i-1,j+1,k))
     $                    /dx(i,j)) **2)
            end do
          end do
        end do
        call exchange3d_mpi(aam(:,:,1:kbm1),im,jm,kbm1)
      end if

      return
      end

!_______________________________________________________________________
      subroutine mode_interaction
! form vertical averages of 3-D fields for use in external (2-D) mode
      implicit none
      include 'pom.h'
      integer i,j,k

      if(mode.ne.2) then

        do j=1,jm
          do i=1,im
            adx2d(i,j)=0.d0
            ady2d(i,j)=0.d0
            drx2d(i,j)=0.d0
            dry2d(i,j)=0.d0
            aam2d(i,j)=0.d0
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              adx2d(i,j)=adx2d(i,j)+advx(i,j,k)*dz(k)
              ady2d(i,j)=ady2d(i,j)+advy(i,j,k)*dz(k)
              drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
              dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
              aam2d(i,j)=aam2d(i,j)+aam(i,j,k)*dz(k)
            end do
          end do
        end do

        call advave

        do j=1,jm
          do i=1,im
            adx2d(i,j)=adx2d(i,j)-advua(i,j)
            ady2d(i,j)=ady2d(i,j)-advva(i,j)
          end do
        end do

      end if

      do j=1,jm
        do i=1,im
          egf(i,j)=el(i,j)*ispi
        end do
      end do

      do j=1,jm
        do i=2,im
          utf(i,j)=ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
        end do
      end do
      do j=2,jm
        do i=1,im
          vtf(i,j)=va(i,j)*(d(i,j)+d(i,j-1))*isp2i
        end do
      end do
      
!      call exchange2d_mpi(utf, im, jm)
!      call exchange2d_mpi(vtf, im, jm)

      return
      end

!_______________________________________________________________________
      subroutine mode_external
! calculate the external (2-D) mode
      implicit none
      include 'pom.h'
      integer i,j

      do j=2,jm
        do i=2,im
          fluxua(i,j)=.25d0*(d(i,j)+d(i-1,j))
     $                 *(dy(i,j)+dy(i-1,j))*ua(i,j)
          fluxva(i,j)=.25d0*(d(i,j)+d(i,j-1))
     $                 *(dx(i,j)+dx(i,j-1))*va(i,j)
        end do
      end do
      
! NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
! with pom98.f. See also modifications to subroutine vertvl
      do j=2,jmm1
        do i=2,imm1
          elf(i,j)=elb(i,j)
     $              +dte2*(-(fluxua(i+1,j)-fluxua(i,j)
     $                      +fluxva(i,j+1)-fluxva(i,j))/art(i,j)
     $                      -vfluxf(i,j))
        end do
      end do
      
      call bcond(1)
      
      call exchange2d_mpi(elf,im,jm)
      
      if(mod(iext,ispadv).eq.0) call advave

      do j=2,jmm1
        do i=2,im
          uaf(i,j)=adx2d(i,j)+advua(i,j)
     $              -aru(i,j)*.25d0
     $                *(cor(i,j)*d(i,j)*(va(i,j+1)+va(i,j))
     $                 +cor(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j)))
     $              +.25d0*grav*(dy(i,j)+dy(i-1,j))
     $                *(d(i,j)+d(i-1,j))
     $                *((1.d0-2.d0*alpha)
     $                   *(el(i,j)-el(i-1,j))
     $                  +alpha*(elb(i,j)-elb(i-1,j)
     $                         +elf(i,j)-elf(i-1,j))
     $                  +e_atmos(i,j)-e_atmos(i-1,j))
     $              +drx2d(i,j)+aru(i,j)*(wusurf(i,j)-wubot(i,j))
        end do
      end do
      
      do j=2,jmm1
        do i=2,im
          uaf(i,j)=((h(i,j)+elb(i,j)+h(i-1,j)+elb(i-1,j))
     $                *aru(i,j)*uab(i,j)
     $              -4.d0*dte*uaf(i,j))
     $             /((h(i,j)+elf(i,j)+h(i-1,j)+elf(i-1,j))
     $                 *aru(i,j))
        end do
      end do
      
      do j=2,jm
        do i=2,imm1
          vaf(i,j)=ady2d(i,j)+advva(i,j)
     $              +arv(i,j)*.25d0
     $                *(cor(i,j)*d(i,j)*(ua(i+1,j)+ua(i,j))
     $               +cor(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1)))
     $              +.25d0*grav*(dx(i,j)+dx(i,j-1))
     $                *(d(i,j)+d(i,j-1))
     $                *((1.d0-2.d0*alpha)*(el(i,j)-el(i,j-1))
     $                  +alpha*(elb(i,j)-elb(i,j-1)
     $                         +elf(i,j)-elf(i,j-1))
     $                  +e_atmos(i,j)-e_atmos(i,j-1))
     $              +dry2d(i,j)+arv(i,j)*(wvsurf(i,j)-wvbot(i,j))
        end do
      end do
      
      do j=2,jm
        do i=2,imm1
          vaf(i,j)=((h(i,j)+elb(i,j)+h(i,j-1)+elb(i,j-1))
     $                *arv(i,j)*vab(i,j)
     $              -4.d0*dte*vaf(i,j))
     $             /((h(i,j)+elf(i,j)+h(i,j-1)+elf(i,j-1))
     $                 *arv(i,j))
        end do
      end do
      
      call bcond(2)
      
      call exchange2d_mpi(uaf,im,jm)
      call exchange2d_mpi(vaf,im,jm)
      
      if(iext.eq.(isplit-2))then
        do j=1,jm
          do i=1,im
            etf(i,j)=.25d0*smoth*elf(i,j)
          end do
        end do

      else if(iext.eq.(isplit-1)) then

        do j=1,jm
          do i=1,im
            etf(i,j)=etf(i,j)+.5d0*(1.-.5d0*smoth)*elf(i,j)
          end do
        end do

      else if(iext.eq.isplit) then

        do j=1,jm
          do i=1,im
            etf(i,j)=(etf(i,j)+.5d0*elf(i,j))*fsm(i,j)
          end do
        end do

      end if
      
! apply filter to remove time split
      do j=1,jm
        do i=1,im
          ua(i,j)=ua(i,j)+.5d0*smoth*(uab(i,j)-2.d0*ua(i,j)+uaf(i,j))
          va(i,j)=va(i,j)+.5d0*smoth*(vab(i,j)-2.d0*va(i,j)+vaf(i,j))
          el(i,j)=el(i,j)+.5d0*smoth*(elb(i,j)-2.d0*el(i,j)+elf(i,j))
          elb(i,j)=el(i,j)
          el(i,j)=elf(i,j)
          d(i,j)=h(i,j)+el(i,j)
          uab(i,j)=ua(i,j)
          ua(i,j)=uaf(i,j)
          vab(i,j)=va(i,j)
          va(i,j)=vaf(i,j)
        end do
      end do

      if(iext.ne.isplit) then
        do j=1,jm
          do i=1,im
            egf(i,j)=egf(i,j)+el(i,j)*ispi
          end do
        end do
        do j=1,jm
          do i=2,im
            utf(i,j)=utf(i,j)+ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
          end do
        end do
        do j=2,jm
          do i=1,im
            vtf(i,j)=vtf(i,j)+va(i,j)*(d(i,j)+d(i,j-1))*isp2i
          end do
        end do
       end if
       
       call exchange2d_mpi(utf, im, jm)
       call exchange2d_mpi(vtf, im, jm)

      return
      end

!_______________________________________________________________________
      subroutine mode_internal
! calculate the internal (3-D) mode
      implicit none
      include 'pom.h'
      integer i,j,k

      if((iint.ne.1.or.time0.ne.0.d0).and.mode.ne.2) then

! adjust u(z) and v(z) such that depth average of (u,v) = (ua,va)
        do j=1,jm
          do i=1,im
            tps(i,j)=0.d0
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)+u(i,j,k)*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=2,im
              u(i,j,k)=(u(i,j,k)-tps(i,j))+
     $                 (utb(i,j)+utf(i,j))/(dt(i,j)+dt(i-1,j))
            end do
          end do
        end do

        do j=1,jm
          do i=1,im
            tps(i,j)=0.d0
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)+v(i,j,k)*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=2,jm
            do i=1,im
              v(i,j,k)=(v(i,j,k)-tps(i,j))+
     $                 (vtb(i,j)+vtf(i,j))/(dt(i,j)+dt(i,j-1))
            end do
          end do
        end do

! calculate w from u, v, dt (h+et), etf and etb
        call vertvl

        call bcondorl(5)

        call exchange3d_mpi(w,im,jm,kb)

! set uf and vf to zero
        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=0.d0
              vf(i,j,k)=0.d0
            end do
          end do
        end do

! calculate q2f and q2lf using uf, vf, a and c as temporary variables
        call advq(q2b,q2,uf)
        call advq(q2lb,q2l,vf)
        call profq

        call bcond(6)

        call exchange3d_mpi(uf(:,:,2:kbm1),im,jm,kbm2)
        call exchange3d_mpi(vf(:,:,2:kbm1),im,jm,kbm2)

        do k=1,kb
          do j=1,jm
            do i=1,im
              q2(i,j,k)=q2(i,j,k)
     $                   +.5d0*smoth*(uf(i,j,k)+q2b(i,j,k)
     $                                -2.d0*q2(i,j,k))
              q2l(i,j,k)=q2l(i,j,k)
     $                   +.5d0*smoth*(vf(i,j,k)+q2lb(i,j,k)
     $                                -2.d0*q2l(i,j,k))
              q2b(i,j,k)=q2(i,j,k)
              q2(i,j,k)=uf(i,j,k)
              q2lb(i,j,k)=q2l(i,j,k)
              q2l(i,j,k)=vf(i,j,k)
            end do
          end do
        end do

! calculate tf and sf using uf, vf, a and c as temporary variables
        if(mode.ne.4) then
          if(nadv.eq.1) then
            call advt1(tb,t,tclim,uf)
            call advt1(sb,s,sclim,vf)
          else if(nadv.eq.2) then
            call advt2(tb,t,tclim,uf)
            call advt2(sb,s,sclim,vf)
          else
            error_status=1
            write(6,'(/''Error: invalid value for nadv'')')
          end if

          call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
          call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

          call proft(uf,wtsurf,tsurf,nbct)
          call proft(vf,wssurf,ssurf,nbcs)

          call bcond(4)

          do k=1,kb
            do j=1,jm
              do i=1,im
                t(i,j,k)=t(i,j,k)
     $                    +.5d0*smoth*(uf(i,j,k)+tb(i,j,k)
     $                                 -2.d0*t(i,j,k))
                s(i,j,k)=s(i,j,k)
     $                    +.5d0*smoth*(vf(i,j,k)+sb(i,j,k)
     $                                 -2.d0*s(i,j,k))
                tb(i,j,k)=t(i,j,k)
                t(i,j,k)=uf(i,j,k)
                sb(i,j,k)=s(i,j,k)
                s(i,j,k)=vf(i,j,k)
              end do
            end do
          end do

          ! restore temperature and salinity
          !call restore_interior ! TODO: Do not restore t and s yet

          call dens(s,t,rho)

        end if

! calculate uf and vf
        call advu
        call advv
        call profu
        call profv

        call bcondorl(3)

        call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
        call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

        do j=1,jm
          do i=1,im
            tps(i,j)=0.d0
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)
     $                  +(uf(i,j,k)+ub(i,j,k)-2.d0*u(i,j,k))*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              u(i,j,k)=u(i,j,k)
     $                  +.5d0*smoth*(uf(i,j,k)+ub(i,j,k)
     $                               -2.d0*u(i,j,k)-tps(i,j))
            end do
          end do
        end do

        do j=1,jm
          do i=1,im
            tps(i,j)=0.d0
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)
     $                  +(vf(i,j,k)+vb(i,j,k)-2.d0*v(i,j,k))*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              v(i,j,k)=v(i,j,k)
     $                  +.5d0*smoth*(vf(i,j,k)+vb(i,j,k)
     $                               -2.d0*v(i,j,k)-tps(i,j))
            end do
          end do
        end do
        
        call exchange3d_mpi(uf,im_local,jm_local,kb)
        call exchange3d_mpi(vf,im_local,jm_local,kb)
        call exchange3d_mpi(u ,im_local,jm_local,kb)
        call exchange3d_mpi(v ,im_local,jm_local,kb)
        
        do k=1,kb
          do j=1,jm
            do i=1,im
              ub(i,j,k)=u(i,j,k)
              u(i,j,k)=uf(i,j,k)
              vb(i,j,k)=v(i,j,k)
              v(i,j,k)=vf(i,j,k)
            end do
          end do
        end do

!        call exchange3d_mpi(ub,im_local,jm_local,kb)
!        call exchange3d_mpi(vb,im_local,jm_local,kb)

      end if

      do j=1,jm
        do i=1,im
          egb(i,j)=egf(i,j)
          etb(i,j)=et(i,j)
          et(i,j)=etf(i,j)
          dt(i,j)=h(i,j)+et(i,j)
          utb(i,j)=utf(i,j)
          vtb(i,j)=vtf(i,j)
          vfluxb(i,j)=vfluxf(i,j)
        end do
      end do

! calculate real w as wr
      call realvertvl

      return
      end

!_______________________________________________________________________
      subroutine print_section
! print output
      implicit none
      include 'pom.h'
      double precision atot,darea,dvol,eaver,saver,taver,vtot,tsalt
      integer i,j,k

      if(mod(iint,iprint).eq.0) then

! print time
        if(my_task.eq.master_task) write(6,'(/
     $    ''**********************************************************''
     $    /''time ='',f9.4,'', iint ='',i8,'', iext ='',i8,
     $    '', iprint ='',i8)') time,iint,iext,iprint

! check for errors
        call sum0d_mpi(error_status,master_task)
        call bcast0d_mpi(error_status,master_task)
        if(error_status.ne.0) then
          if(my_task.eq.master_task) write(*,'(/a)')
     $                                       'POM terminated with error'
          call finalize_mpi
          stop
        end if

! local averages
        vtot=0.d0
        atot=0.d0
        taver=0.d0
        saver=0.d0
        eaver=0.d0
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              darea=dx(i,j)*dy(i,j)*fsm(i,j)
              dvol=darea*dt(i,j)*dz(k)
              vtot=vtot+dvol
              taver=taver+tb(i,j,k)*dvol
              saver=saver+sb(i,j,k)*dvol
            end do
          end do
        end do

        do j=1,jm
          do i=1,im
            darea=dx(i,j)*dy(i,j)*fsm(i,j)
            atot=atot+darea
            eaver=eaver+et(i,j)*darea
          end do
        end do

        taver=taver/vtot
        saver=saver/vtot
        eaver=eaver/atot
        tsalt=(saver+sbias)*vtot

! print averages
! global averages requiere to transfer high amounts of data between
! processor - therefore, only local average for master_task is printed
        if(my_task.eq.master_task) write(6,'(/''vtot = '',e16.7,
     $    ''   atot = '',e16.7,''  eaver ='',e16.7/''taver ='',e16.7,
     $    ''   saver ='',e16.7,''  tsalt ='',e16.7)')
     $    vtot,atot,eaver,taver,saver,tsalt

      end if

      return
      end

!_______________________________________________________________________
      subroutine check_velocity
! check if velocity condition is violated
      implicit none
      include 'pom.h'
      double precision vamax,atot,darea,dvol,eaver,saver,taver,
     &                                                        vtot,tsalt
      integer i,j,k
      integer imax,jmax

      vamax=0.d0

      do j=1,jm
        do i=1,im
          if(abs(vaf(i,j)).ge.vamax) then
            vamax=abs(vaf(i,j))
            imax=i
            jmax=j
          end if
        end do
      end do

      if(vamax.gt.vmaxl) then
        if(my_task.eq.master_task.and.error_status.eq.0) write(6,'(/
     $    ''Error: velocity condition violated''/''time ='',f9.4,
     $    '', iint ='',i8,'', iext ='',i8,'', iprint ='',i8,/
     $    ''vamax ='',e12.3,''   imax,jmax ='',2i5)')
     $    time,iint,iext,iprint,vamax,imax,jmax
        error_status=1
      end if

      return
      end
