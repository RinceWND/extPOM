! bounds_forcing.f

! spcify variable boundary conditions, atmospheric forcing, restoring

!_______________________________________________________________________
      subroutine bcond(idx)
! apply open boundary conditions
! closed boundary conditions are automatically enabled through
! specification of the masks, dum, dvm and fsm, in which case the open
! boundary conditions, included below, will be overwritten
      implicit none
      include 'pom.h'
      integer idx
      integer i,j,k
      double precision ga,u1,wm
      double precision hmax

      if(idx.eq.1) then

! External (2-D) elevation boundary conditions (Clamped)
        if(n_west.eq.-1) then
!          elf(1,:) = elw
          elf(1,:) = elf(2,:)
        end if
        if(n_east.eq.-1) then
!          elf(im,:) = ele
          elf(im,:) = elf(imm1,:)
        end if

        if(n_south.eq.-1) then
!          elf(:,1) = els
          elf(:,1) = elf(:,2)
        end if
        if(n_north.eq.-1) then
!          elf(:,jm) = eln
          elf(:,jm) = elf(:,jmm1)
        end if

        elf = elf*fsm

        return

      else if(idx.eq.2) then

! external (2-D) velocity boundary conditions
        ! west
        if(n_west.eq.-1) then
          uaf(2,2:jmm1) = uabw(2:jmm1)
     $            -rfw*sqrt(grav/d(2,2:jmm1))*(el(2,2:jmm1)-elw(2:jmm1))
          uaf(2,2:jmm1) = ramp*uaf(2,2:jmm1)
          uaf(1,2:jmm1) = uaf(2,2:jmm1)
          vaf(1,2:jmm1) = vabw(2:jmm1)
        end if

        ! east
        if(n_east.eq.-1) then
          uaf(im,2:jmm1) = uabe(2:jmm1)
     $      +rfe*sqrt(grav/d(imm1,2:jmm1))*(el(imm1,2:jmm1)-ele(2:jmm1))
          uaf(im,2:jmm1) = ramp*uaf(im,2:jmm1)
          vaf(im,2:jmm1) = vabe(2:jmm1)
        end if

        ! south
        if(n_south.eq.-1) then
          vaf(2:imm1,2) = vabs(2:imm1)
     $            -rfs*sqrt(grav/d(2:imm1,2))*(el(2:imm1,2)-els(2:imm1))
          vaf(2:imm1,2) = ramp*vaf(2:imm1,2)
          vaf(2:imm1,1) = vaf(2:imm1,2)
          uaf(2:imm1,1) = uabs(2:imm1)
        end if

        ! north
        if(n_north.eq.-1) then
          vaf(2:imm1,jm) = vabn(2:imm1)
     $      +rfn*sqrt(grav/d(2:imm1,jmm1))*(el(2:imm1,jmm1)-eln(2:imm1))
          vaf(2:imm1,jm) = ramp*vaf(2:imm1,jm)
          uaf(2:imm1,jm) = uabn(2:imm1)
        end if

        uaf = uaf*dum
        vaf = vaf*dvm

        return

      else if(idx.eq.3) then

! internal (3-D) velocity boundary conditions
! radiation conditions
! smoothing is used in the direction tangential to the boundaries
        hmax=maxval(d) !4500.d0

        do k=1,kbm1
          do j=2,jmm1
            ! east
            if(n_east.eq.-1) then
              ga=sqrt(d(im,j)/hmax)
              uf(im,j,k)=ga*(.25d0*u(imm1,j-1,k)+.5d0*u(imm1,j,k)
     $                       +.25d0*u(imm1,j+1,k))
     $                    +(1.d0-ga)*(.25d0*ube(j-1,k)+.5d0*ube(j,k)
     $                      +.25d0*ube(j+1,k))
              vf(im,j,k)=vbe(j,k)
            end if
            ! west
            if(n_west.eq.-1) then
              ga=sqrt(d(1,j)/hmax)
              uf(2,j,k)=ga*(.25d0*u(3,j-1,k)+.5d0*u(3,j,k)
     $                      +.25d0*u(3,j+1,k))
     $                   +(1.d0-ga)*(.25d0*ubw(j-1,k)+.5d0*ubw(j,k)
     $                     +.25d0*ubw(j+1,k))
              uf(1,j,k)=uf(2,j,k)
              vf(1,j,k)=vbw(j,k)
            end if
          end do
        end do

        do k=1,kbm1
          do i=2,imm1
            ! south
            if(n_south.eq.-1) then
              ga=sqrt(d(i,1)/hmax)
              vf(i,2,k)=ga*(.25d0*v(i-1,3,k)+.5d0*v(i,3,k)
     $                      +.25d0*v(i+1,3,k))
     $                   +(1.d0-ga)*(.25d0*vbs(i-1,k)+.5d0*vbs(i,k)
     $                     +.25d0*vbs(i+1,k))
              vf(i,1,k)=vf(i,2,k)
              uf(i,1,k)=ubs(i,k)
            end if
            ! north
            if(n_north.eq.-1) then
              ga=sqrt(d(i,jm)/hmax)
              vf(i,jm,k)=ga*(.25d0*v(i-1,jmm1,k)+.5d0*v(i,jmm1,k)
     $                       +.25d0*v(i+1,jmm1,k))
     $                    +(1.d0-ga)*(.25d0*vbn(i-1,k)+.5d0*vbn(i,k)
     $                      +.25d0*vbn(i+1,k))
              uf(i,jm,k)=ubn(i,k)
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.4) then

! temperature and salinity boundary conditions (using uf and vf,
! respectively)
        do k=1,kbm1
          do j=1,jm
            ! east
            if(n_east.eq.-1) then
              u1=2.d0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
              if(u1.le.0.d0) then
                uf(im,j,k)=t(im,j,k)-u1*(tbe(j,k)-t(im,j,k))
                vf(im,j,k)=s(im,j,k)-u1*(sbe(j,k)-s(im,j,k))
              else
                uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k))
                vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5d0*(w(imm1,j,k)+w(imm1,j,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(imm1,j))
                  uf(im,j,k)=uf(im,j,k)-wm*(t(imm1,j,k-1)-t(imm1,j,k+1))
                  vf(im,j,k)=vf(im,j,k)-wm*(s(imm1,j,k-1)-s(imm1,j,k+1))
                endif
              end if
            end if

            ! west
            if(n_west.eq.-1) then
              u1=2.d0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
              if(u1.ge.0.d0) then
                uf(1,j,k)=t(1,j,k)-u1*(t(1,j,k)-tbw(j,k))
                vf(1,j,k)=s(1,j,k)-u1*(s(1,j,k)-sbw(j,k))
              else
                uf(1,j,k)=t(1,j,k)-u1*(t(2,j,k)-t(1,j,k))
                vf(1,j,k)=s(1,j,k)-u1*(s(2,j,k)-s(1,j,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5d0*(w(2,j,k)+w(2,j,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(2,j))
                  uf(1,j,k)=uf(1,j,k)-wm*(t(2,j,k-1)-t(2,j,k+1))
                  vf(1,j,k)=vf(1,j,k)-wm*(s(2,j,k-1)-s(2,j,k+1))
                end if
              end if
            end if
          end do

          do i=1,im
            ! south
            if(n_south.eq.-1) then
              u1=2.d0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
              if(u1.ge.0.d0) then
                uf(i,1,k)=t(i,1,k)-u1*(t(i,1,k)-tbs(i,k))
                vf(i,1,k)=s(i,1,k)-u1*(s(i,1,k)-sbs(i,k))
              else
                uf(i,1,k)=t(i,1,k)-u1*(t(i,2,k)-t(i,1,k))
                vf(i,1,k)=s(i,1,k)-u1*(s(i,2,k)-s(i,1,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5d0*(w(i,2,k)+w(i,2,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(i,2))
                  uf(i,1,k)=uf(i,1,k)-wm*(t(i,2,k-1)-t(i,2,k+1))
                  vf(i,1,k)=vf(i,1,k)-wm*(s(i,2,k-1)-s(i,2,k+1))
                end if
              end if
            end if

            ! north
            if(n_north.eq.-1) then
              u1=2.d0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
              if(u1.le.0.d0) then
                uf(i,jm,k)=t(i,jm,k)-u1*(tbn(i,k)-t(i,jm,k))
                vf(i,jm,k)=s(i,jm,k)-u1*(sbn(i,k)-s(i,jm,k))
              else
                uf(i,jm,k)=t(i,jm,k)-u1*(t(i,jm,k)-t(i,jmm1,k))
                vf(i,jm,k)=s(i,jm,k)-u1*(s(i,jm,k)-s(i,jmm1,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5d0*(w(i,jmm1,k)+w(i,jmm1,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(i,jmm1))
                  uf(i,jm,k)=uf(i,jm,k)-wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1))
                  vf(i,jm,k)=vf(i,jm,k)-wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1))
                end if
              end if
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.5) then

! vertical velocity boundary conditions
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.6) then

! q2 and q2l boundary conditions

        do k=1,kb
          do j=1,jm
            ! west
            if(n_west.eq.-1) then
              u1=2.d0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
              if(u1.ge.0.d0) then
                uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
              else
                uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
              end if
            end if

            ! east
            if(n_east.eq.-1) then
              u1=2.d0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
              if(u1.le.0.d0) then
                uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
                vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
              else
                uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
                vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
              end if
            end if
          end do

          do i=1,im
            ! south
            if(n_south.eq.-1) then
              u1=2.d0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
              if(u1.ge.0.d0) then
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
              else
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
              end if
            end if

            ! north
            if(n_north.eq.-1) then
              u1=2.d0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
              if(u1.le.0.d0) then
                uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
                vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
              else
                uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
                vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
              end if
            end if
          end do
        end do

        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.d-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.d-10
            end do
          end do
        end do

        return

      endif

      end

!_______________________________________________________________________
      subroutine bcondorl(idx)
! this is an optional subroutine replacing  bcond and using Orlanski's
! scheme (J. Comp. Phys. 21, 251-269, 1976), specialized for the
! seamount problem
      implicit none
      include 'pom.h'
      integer idx
      double precision cl,denom
!      double precision ar,eps
      integer i,j,k

      if(idx.eq.1) then

! external (2-D) elevation boundary conditions
        if(n_west.eq.-1) elf(1,:)=elf(2,:)
        if(n_east.eq.-1) elf(im,:)=elf(imm1,:)

        elf = elf*fsm

        return

      else if(idx.eq.2) then

! external (2-D) velocity  boundary conditions
        do j=2,jmm1
          ! east
          if(n_east.eq.-1) then
            denom=(uaf(im-1,j)+uab(im-1,j)-2.d0*ua(im-2,j))
            if(denom.eq.0.0d0)denom=0.01d0
            cl=(uab(im-1,j)-uaf(im-1,j))/denom
            if(cl.gt.1.d0) cl=1.d0
            if(cl.lt.0.d0) cl=0.d0
            uaf(im,j)=(uab(im,j)*(1.d0-cl)+2.d0*cl*ua(im-1,j))
     $                  /(1.d0+cl)
            vaf(im,j)=0.d0
          end if

          ! west
          if(n_west.eq.-1) then
            denom=(uaf(3,j)+uab(3,j)-2.d0*ua(4,j))
            if(denom.eq.0.0d0)denom=0.01d0
            cl=(uab(3,j)-uaf(3,j))/denom
            if(cl.gt.1.d0) cl=1.d0
            if(cl.lt.0.d0) cl=0.d0
            uaf(2,j)=(uab(2,j)*(1.d0-cl)+2.d0*cl*ua(3,j))
     $                 /(1.d0+cl)
            uaf(1,j)=uaf(2,j)
            vaf(1,j)=0.d0
          end if
        end do

        do i=2,imm1
          ! south
          if(n_south.eq.-1) then
            denom=(vaf(i,3)+vab(i,3)-2.d0*va(i,4))
            if(denom.eq.0.0d0)denom=0.01d0
            cl=(vab(i,3)-vaf(i,3))/denom
            if(cl.gt.1.d0) cl=1.d0
            if(cl.lt.0.d0) cl=0.d0
            vaf(i,2)=(vab(i,2)*(1.d0-cl)+2.d0*cl*va(i,3))
     $                 /(1.d0+cl)
            vaf(i,1)=vaf(i,2)
            uaf(i,1)=0.d0
          end if

          ! north
          if(n_north.eq.-1) then
            denom=(vaf(i,jm-1)+vab(i,jm-1)-2.d0*va(i,jm-2))
            if(denom.eq.0.0d0)denom=0.01d0
            cl=(vab(i,jm-1)-vaf(i,jm-1))/denom
            if(cl.gt.1.d0) cl=1.d0
            if(cl.lt.0.d0) cl=0.d0
            vaf(i,jm)=(vab(i,jm)*(1.d0-cl)+2.d0*cl*va(i,jm-1))
     $                  /(1.d0+cl)
            uaf(i,jm)=0.d0
          end if
        end do

        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do

        return

      else if(idx.eq.3) then

! internal (3-D) velocity boundary conditions

        do k=1,kbm1
          do j=2,jmm1
            ! east
            if(n_east.eq.-1) then
              denom=(uf(im-1,j,k)+ub(im-1,j,k)-2.d0*u(im-2,j,k))
              if(denom.eq.0.d0)denom=0.01d0
              cl=(ub(im-1,j,k)-uf(im-1,j,k))/denom
              if(cl.gt.1.d0) cl=1.d0
              if(cl.lt.0.d0) cl=0.d0
              uf(im,j,k)=(ub(im,j,k)*(1.d0-cl)+2.d0*cl*u(im-1,j,k))
     $                    /(1.d0+cl)
              vf(im,j,k)=0.d0
            end if

            ! west
            if(n_west.eq.-1) then
              denom=(uf(3,j,k)+ub(3,j,k)-2.d0*u(4,j,k))
              if(denom.eq.0.d0)denom=0.01d0
              cl=(ub(3,j,k)-uf(3,j,k))/denom
              if(cl.gt.1.d0) cl=1.d0
              if(cl.lt.0.d0) cl=0.d0
              uf(2,j,k)=(ub(2,j,k)*(1.d0-cl)+2.d0*cl*u(3,j,k))
     $                   /(1.d0+cl)
              uf(1,j,k)=uf(2,j,k)
              vf(1,j,k)=0.d0
            end if
          end do

          do i=2,imm1
            ! south
            if(n_south.eq.-1) then
              denom=(vf(i,3,k)+vb(i,3,k)-2.d0*v(i,4,k))
              if(abs(denom).eq.0.0d0)denom=0.01d0
              cl=(vb(i,3,k)-vf(i,3,k))/denom
              if(cl.gt.1.d0) cl=1.d0
              if(cl.lt.0.d0) cl=0.d0
              vf(i,2,k)=(vb(i,2,k)*(1.d0-cl)+2.d0*cl*v(i,3,k))
     $                   /(1.d0+cl)
              vf(i,1,k)=vf(i,2,k)
              uf(i,1,k)=0.d0
            end if

            ! north
            if(n_north.eq.-1) then
              denom=(vf(i,jm-1,k)+vb(i,jm-1,k)-2.d0*v(i,jm-2,k))
              if(abs(denom).eq.0.0d0)denom=0.01d0
              cl=(vb(i,jm-1,k)-vf(i,jm-1,k))/denom
              if(cl.gt.1.d0) cl=1.d0
              if(cl.lt.0.d0) cl=0.d0
              vf(i,jm,k)=(vb(i,jm,k)*(1.d0-cl)+2.d0*cl*v(i,jm-1,k))
     $                    /(1.d0+cl)
              uf(i,jm,k)=0.d0
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.4) then

! temperature and salinity boundary conditions (using uf and vf,
! respectively)
        do k=1,kbm1
          do j=1,jm
            ! east
            if(n_east.eq.-1) then
              ube(j,k)=ub(im,j,k)
              denom=(uf(im-1,j,k)+tb(im-1,j,k)-2.d0*t(im-2,j,k))
              if(denom.eq.0.d0) denom=0.01d0
              cl=(tb(im-1,j,k)-uf(im-1,j,k))/denom
              if(cl.gt.1.d0) cl=1.d0
              if(cl.lt.0.d0) cl=0.d0
              uf(im,j,k)=(tb(im,j,k)*(1.d0-cl)+2.d0*cl*t(im-1,j,k))
     $                    /(1.d0+cl)
              if(cl.eq.0.d0.and.ube(j,k).le.0.d0) uf(im,j,k)=tbe(j,k)

              denom=(vf(im-1,j,k)+sb(im-1,j,k)-2.d0*s(im-2,j,k))
              if(denom.eq.0.d0) denom=0.01d0
              cl=(sb(im-1,j,k)-vf(im-1,j,k))/denom
              if(cl.gt.1.d0) cl=1.d0
              if(cl.lt.0.d0) cl=0.d0
              vf(im,j,k)=(sb(im,j,k)*(1.d0-cl)+2.d0*cl*s(im-1,j,k))
     $                    /(1.d0+cl)
              if(cl.eq.0.d0.and.ube(j,k).le.0.d0) vf(im,j,k)=sbe(j,k)
            end if

            ! west
            if(n_west.eq.-1) then
              ubw(j,k)=ub(2,j,k)
              denom=(uf(2,j,k)+tb(2,j,k)-2.d0*t(3,j,k))
              if(denom.eq.0.d0) denom=0.01d0
              cl=(tb(2,j,k)-uf(2,j,k))/denom
              if(cl.gt.1.d0) cl=1.d0
              if(cl.lt.0.d0) cl=0.d0
              uf(1,j,k)=(tb(1,j,k)*(1.d0-cl)+2.d0*cl*t(2,j,k))/(1.d0+cl)
              if(cl.eq.0.d0.and.ubw(j,k).ge.0.d0) uf(1,j,k)=tbw(j,k)

              denom=(vf(2,j,k)+sb(2,j,k)-2.d0*s(3,j,k))
              if(denom.eq.0.d0) denom=0.01d0
              cl=(sb(2,j,k)-vf(2,j,k))/denom
              if(cl.gt.1.d0) cl=1.d0
              if(cl.lt.0.d0) cl=0.d0
              vf(1,j,k)=(sb(1,j,k)*(1.d0-cl)+2.d0*cl*s(2,j,k))/(1.d0+cl)
              if(cl.eq.0.d0.and.ubw(j,k).ge.0.d0) vf(1,j,k)=sbw(j,k)
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.5) then

! vertical velocity boundary conditions
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.6) then

! q2 and q2l boundary conditions
        do k=1,kb
          do j=1,jm
            if(n_east.eq.-1) then
              uf(im,j,k)=1.d-10
              vf(im,j,k)=1.d-10
            end if
            if(n_west.eq.-1) then
              uf(1,j,k)=1.d-10
              vf(1,j,k)=1.d-10
            end if
          end do

          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      endif

      end

! _____________________________________________________________________
      subroutine lateral_bc
! create variable lateral boundary conditions
      implicit none
      include 'pom.h'
      integer nz
      parameter(nz=40)
      integer i,j,k,ntime,ibc
      double precision tbc,fold,fnew
!      double precision z0(nz),hs(im,jm)
!      double precision t_w(jm,nz),s_w(jm,nz),t_e(jm,nz),s_e(jm,nz),q,
!     $       t_n(im,nz),s_n(im,nz),t_s(im,nz),s_s(im,nz)
!      double precision ts1(im,jm,nz),ss1(im,jm,nz)
!      double precision ts2(im,jm,kb),ss2(im,jm,kb)

      tbc=1./24. ! time between bc files (days)
      ibc=int(tbc*86400.d0/dti)
      ntime=int(time/tbc)
! read bc data
      ! read initial bc file
      if (iint.eq.1) then
        call read_boundary_conditions_pnetcdf((iint+cont_bry)/ibc+1,kb
     $                         ,tbwf,sbwf,ubwf,vbwf,tbef,sbef,ubef,vbef
     $                         ,tbnf,sbnf,vbnf,ubnf,tbsf,sbsf,vbsf,ubsf
     $                         ,elw,ele,eln,els)
! integrate by depth
        uabwf = 0.
        vabwf = 0.
        uabef = 0.
        vabef = 0.
        uabnf = 0.
        vabnf = 0.
        uabsf = 0.
        vabsf = 0.
        do k = 1,kb
          uabwf(:) = uabwf(:) + ubwf(:,k)*dz(k)
          vabwf(:) = vabwf(:) + vbwf(:,k)*dz(k)
          uabef(:) = uabef(:) + ubef(:,k)*dz(k)
          vabef(:) = vabef(:) + vbef(:,k)*dz(k)
          uabnf(:) = uabnf(:) + ubnf(:,k)*dz(k)
          vabnf(:) = vabnf(:) + vbnf(:,k)*dz(k)
          uabsf(:) = uabsf(:) + ubsf(:,k)*dz(k)
          vabsf(:) = vabsf(:) + vbsf(:,k)*dz(k)
        end do
!  south
!        do i=1,im
!          do j=1,jm
!            hs(i,j)=h(i,2)
!            do k=1,nz
!              ts1(i,j,k)=t_s(i,k)
!              ss1(i,j,k)=s_s(i,k)
!            end do
!          end do
!        end do
!        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        do i=1,im
!          do k=1,kb
!            !write(52,'(2(i4),e16.7)') i,k,ss2(i,jm/2,k)
!            tbsf(i,k)=ts2(i,jm/2,k)
!            sbsf(i,k)=ss2(i,jm/2,k)
!          end do
!        end do
! north
!        do i=1,im
!          do j=1,jm
!            hs(i,j)=h(i,jmm1)
!            do k=1,nz
!              ts1(i,j,k)=t_n(i,k)
!              ss1(i,j,k)=s_n(i,k)
!            end do
!          end do
!        end do
!        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        do i=1,im
!          do k=1,kb
!            tbnf(i,k)=ts2(i,jm/2,k)
!            sbnf(i,k)=ss2(i,jm/2,k)
!          end do
!        end do
! east
!        do i=1,im
!          do j=1,jm
!            hs(i,j)=h(imm1,j)
!            do k=1,nz
!              ts1(i,j,k)=t_e(j,k)
!              ss1(i,j,k)=s_e(j,k)
!            end do
!          end do
!        end do
!        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        do j=1,jm
!          do k=1,kb
!            tbef(j,k)=ts2(im/2,j,k)
!            sbef(j,k)=ss2(im/2,j,k)
!          end do
!        end do
! west
!        do i=1,im
!          do j=1,jm
!            hs(i,j)=h(2,j)
!            do k=1,nz
!              ts1(i,j,k)=t_w(j,k)
!              ss1(i,j,k)=s_w(j,k)
!            end do
!          end do
!        end do
!        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        do j=1,jm
!          do k=1,kb
!            tbwf(j,k)=ts2(im/2,j,k)
!            sbwf(j,k)=ss2(im/2,j,k)
!          end do
!        end do

      end if
      ! read bc file corresponding to next time
      if (iint.eq.1 .or. mod(iint+cont_bry,ibc).eq.0.) then
        tbwb = tbwf
        sbwb = sbwf
        tbeb = tbef
        sbeb = sbef
        tbnb = tbnf
        sbnb = sbnf
        tbsb = tbsf
        sbsb = sbsf
        uabwb = uabwf
        uabeb = uabef
        vabnb = vabnf
        vabsb = vabsf
        if (iint.ne.iend) then
          call read_boundary_conditions_pnetcdf(
     $                                     (iint+cont_bry+ibc)/ibc+1,kb
     $     ,tbwf,sbwf,ubwf,vbwf,tbef,sbef,ubef,vbef,tbnf,sbnf,vbnf,ubnf
     $                             ,tbsf,sbsf,vbsf,ubsf,elw,ele,eln,els)
! integrate by depth
          uabwf = 0.
          vabwf = 0.
          uabef = 0.
          vabef = 0.
          uabnf = 0.
          vabnf = 0.
          uabsf = 0.
          vabsf = 0.
          do k=1,kb
            uabwf(:) = uabwf(:) + ubwf(:,k)*dz(k)
            vabwf(:) = vabwf(:) + vbwf(:,k)*dz(k)
            uabef(:) = uabef(:) + ubef(:,k)*dz(k)
            vabef(:) = vabef(:) + vbef(:,k)*dz(k)
            uabnf(:) = uabnf(:) + ubnf(:,k)*dz(k)
            vabnf(:) = vabnf(:) + vbnf(:,k)*dz(k)
            uabsf(:) = uabsf(:) + ubsf(:,k)*dz(k)
            vabsf(:) = vabsf(:) + vbsf(:,k)*dz(k)
          end do
! south
!        do i=1,im
!          do j=1,jm
!            hs(i,j)=h(i,2)
!            do k=1,nz
!              ts1(i,j,k)=t_s(i,k)
!              ss1(i,j,k)=s_s(i,k)
!            end do
!          end do
!        end do
!        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        do i=1,im
!          do k=1,kb
!            tbsf(i,k)=ts2(i,jm/2,k)
!            sbsf(i,k)=ss2(i,jm/2,k)
!          end do
!        end do
! north
!        do i=1,im
!          do j=1,jm
!            hs(i,j)=h(i,jmm1)
!            do k=1,nz
!              ts1(i,j,k)=t_n(i,k)
!              ss1(i,j,k)=s_n(i,k)
!            end do
!          end do
!        end do
!        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        do i=1,im
!          do k=1,kb
!            tbnf(i,k)=ts2(i,jm/2,k)
!            sbnf(i,k)=ss2(i,jm/2,k)
!          end do
!        end do
! east
!        do i=1,im
!          do j=1,jm
!            hs(i,j)=h(imm1,j)
!            do k=1,nz
!              ts1(i,j,k)=t_e(j,k)
!              ss1(i,j,k)=s_e(j,k)
!            end do
!          end do
!        end do
!        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        do j=1,jm
!          do k=1,kb
!            tbef(j,k)=ts2(im/2,j,k)
!            sbef(j,k)=ss2(im/2,j,k)
!          end do
!        end do
! west
!        do i=1,im
!          do j=1,jm
!            hs(i,j)=h(2,j)
!            do k=1,nz
!              ts1(i,j,k)=t_w(j,k)
!              ss1(i,j,k)=s_w(j,k)
!            end do
!          end do
!        end do
!        call ztosig(z0,ts1,zz,hs,ts2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        call ztosig(z0,ss1,zz,hs,ss2,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        do j=1,jm
!          do k=1,kb
!            tbwf(j,k)=ts2(im/2,j,k)
!            sbwf(j,k)=ss2(im/2,j,k)
!          end do
!        end do
! end
        end if
      end if

! linear interpolation in time
      fnew=time/tbc-real(ntime, 8)
      fold=1.-fnew
      tbw(1:jm,1:kb) = fold*tbwb(1:jm,1:kb)+fnew*tbwf(1:jm,1:kb)
      sbw(1:jm,1:kb) = fold*sbwb(1:jm,1:kb)+fnew*sbwf(1:jm,1:kb)
      ubw(1:jm,1:kb) = fold*ubwb(1:jm,1:kb)+fnew*ubwf(1:jm,1:kb)
      tbe(1:jm,1:kb) = fold*tbeb(1:jm,1:kb)+fnew*tbef(1:jm,1:kb)
      sbe(1:jm,1:kb) = fold*sbeb(1:jm,1:kb)+fnew*sbef(1:jm,1:kb)
      ube(1:jm,1:kb) = fold*ubeb(1:jm,1:kb)+fnew*ubef(1:jm,1:kb)
      tbn(1:im,1:kb) = fold*tbnb(1:im,1:kb)+fnew*tbnf(1:im,1:kb)
      sbn(1:im,1:kb) = fold*sbnb(1:im,1:kb)+fnew*sbnf(1:im,1:kb)
      vbn(1:im,1:kb) = fold*vbnb(1:im,1:kb)+fnew*vbnf(1:im,1:kb)
      tbs(1:im,1:kb) = fold*tbsb(1:im,1:kb)+fnew*tbsf(1:im,1:kb)
      sbs(1:im,1:kb) = fold*sbsb(1:im,1:kb)+fnew*sbsf(1:im,1:kb)
      vbs(1:im,1:kb) = fold*vbsb(1:im,1:kb)+fnew*vbsf(1:im,1:kb)
      uabe = 0.
      uabw = 0.
      vabn = 0.
      vabs = 0.
      do k=1,kb
        uabe(1:jm) = uabe(1:jm) + ube(1:jm,k)*dz(k)
        uabw(1:jm) = uabw(1:jm) + ubw(1:jm,k)*dz(k)
        vabn(1:im) = vabn(1:im) + vbn(1:im,k)*dz(k)
        vabs(1:im) = vabs(1:im) + vbs(1:im,k)*dz(k)
      end do

      return
      end

!_______________________________________________________________________
      subroutine wind
! read and interpolate (in time) wind stress
      implicit none
      include 'pom.h'
      integer i,j,ntime,iwind
      double precision twind,fold,fnew
      double precision, dimension(im,jm) :: wu, wv

      twind= .125!30. ! time between wind files (days)
      iwind=int(twind*86400.d0/dti)

! read wind stress data
      ! read initial wind file
      if (iint.eq.1) then
        call read_wind_pnetcdf((iint+cont_bry)/iwind+1,wu,wv)
        wusurff(1:im,1:jm) = wu
        wvsurff(1:im,1:jm) = wv
      end if
      ! read wind file corresponding to next time
      if (iint.eq.1 .or. mod(iint+cont_bry,iwind).eq.0.) then
        do i=1,im
          do j=1,jm
            wusurfb(i,j)=wusurff(i,j)
            wvsurfb(i,j)=wvsurff(i,j)
          end do
        end do
        if (iint.ne.iend) then
          call read_wind_pnetcdf((iint+cont_bry+iwind)/iwind+1,wu,wv)
          wusurff(1:im,1:jm) = wu
          wvsurff(1:im,1:jm) = wv
        end if
      end if

! linear interpolation in time
      ntime=int(time/twind)
      fnew=time/twind-ntime
      fold=1.-fnew
      wusurf(1:im,1:jm)=fold*wusurfb(1:im,1:jm)+fnew*wusurff(1:im,1:jm)
      wvsurf(1:im,1:jm)=fold*wvsurfb(1:im,1:jm)+fnew*wvsurff(1:im,1:jm)

      return
      end

!_______________________________________________________________________
      subroutine heat
! read and interpolate heat flux in time
      implicit none
      include 'pom.h'
      integer i,j,ntime,iheat
      double precision theat,fold,fnew
      double precision, dimension(im,jm) :: shf, swr

      theat=.125 ! time between heat forcing (days)
      iheat=int(theat*86400.d0/dti)

! read heat stress data
      ! read initial heat file
      if (iint.eq.1) then
        call read_heat_pnetcdf((iint+cont_bry)/iheat+1,shf,swr)
        wtsurff(1:im,1:jm) = shf
        swradf(1:im,1:jm) = swr
      end if
      ! read heat forcing corresponding to next theat
      if (iint.eq.1 .or. mod(iint+cont_bry,iheat).eq.0.) then
        do i=1,im
          do j=1,jm
            wtsurfb(i,j)=wtsurff(i,j)
            swradb(i,j)=swradf(i,j)
          end do
        end do
        if (iint/=iend) then
          call read_heat_pnetcdf((iint+cont_bry+iheat)/iheat+1,shf,swr)
          wtsurff(1:im,1:jm) = shf
          swradf(1:im,1:jm) = swr
        end if
      end if

! linear interpolation in time
      ntime=int(time/theat)
      fnew=time/theat-ntime
      fold=1.-fnew
      do i=1,im
        do j=1,jm
          wtsurf(i,j)=fold*wtsurfb(i,j)+fnew*wtsurff(i,j)
          swrad(i,j)=fold*swradb(i,j)+fnew*swradf(i,j)
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine surface
! read and interpolate heat flux in time
      implicit none
      include 'pom.h'
      integer i,j,ntime,isrf
      double precision tsrf!,fold,fnew
      double precision, dimension(im,jm) :: sst, sss

      tsrf=.125 ! time between heat forcing (days)
      isrf=int(tsrf*86400.d0/dti)

! read heat stress data
      ! read surface fields (no interpolation)
      if (iint.eq.1 .or. mod(iint+cont_bry,isrf).eq.0.) then
        call read_surface_pnetcdf((iint+cont_bry)/isrf+1,sst,sss)
        tsurf(1:im,1:jm) = sst
!        ssurf(1:im,1:jm) = sss ! Do not overwrite SSS yet.
      end if

      return
      end

!_______________________________________________________________________
      subroutine water
! read and interpolate water flux in time
      implicit none
      include 'pom.h'
      integer i,j,ntime,iwater
      double precision twater,fold,fnew

      twater=1. ! time between wind forcing (days)
      iwater=int(twater*86400.d0/dti)

! read wind stress data
      ! read initial water file
      if (iint.eq.1) call read_water_pnetcdf(iint/iwater,wssurff)
      ! read heat forcing corresponding to next twind
      if (iint.eq.1 .or. mod(iint,iwater).eq.0.) then
        do i=1,im
          do j=1,jm
            wssurfb(i,j)=wssurff(i,j)
          end do
        end do
        call read_water_pnetcdf((iint+iwater)/iwater,wssurff)
      end if

! linear interpolation in time
      ntime=int(time/twater)
      fnew=time/twater-ntime
      fold=1.-fnew
      do i=1,im
        do j=1,jm
          wssurf(i,j)=fold*wssurfb(i,j)+fnew*wssurff(i,j)
        end do
      end do

      return
      end

! _____________________________________________________________________
      subroutine restore_interior
! read, interpolate (in time) and apply restore interior data
      implicit none
      include 'pom.h'
      integer nz
      parameter(nz=40)
      double precision z0(nz),tr(im,jm,kb),sr(im,jm,kb)
      integer i,j,k,ntime,irst
      double precision trst,fold,fnew

      trst=30. ! time between restore files (days)
      irst=int(trst*86400.d0/dti)
      ntime=int(time/trst)

! read restore data
      ! read initial restore file
      if (iint.eq.2) then
        call read_restore_ts_interior_pnetcdf((iint/irst)+1,kb,tr,sr)
        trstrf(1:im,1:jm,:) = tr
        srstrf(1:im,1:jm,:) = sr
        taurstrf = 1./trst
!        call ztosig(z0,f0,zz,h,trstrf,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        call read_restore_s_interior_pnetcdf(iint/irst,nz,z0,f0)
!        call ztosig(z0,f0,zz,h,srstrf,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!        call read_restore_tau_interior_pnetcdf(iint/irst,nz,z0,f0)
!        call ztosig(z0,f0,zz,h,taurstrf,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
      end if
      ! read restore file corresponding to next time
      if (iint.eq.2 .or. mod(iint,irst).eq.0.) then
        do k=1,kbm1
          do i=1,im
            do j=1,jm
              trstrb(i,j,k)=trstrf(i,j,k)
              srstrb(i,j,k)=srstrf(i,j,k)
              taurstrb(i,j,k)=taurstrf(i,j,k)
            end do
          end do
        end do
        if (iint.ne.iend) then
          call read_restore_ts_interior_pnetcdf((iint+irst)/irst+1,kb,
     $                                                            tr,sr)
          trstrf(1:im,1:jm,:) = tr
          srstrf(1:im,1:jm,:) = sr
          taurstrf = 1./trst
!          call ztosig(z0,f0,zz,h,trstrf,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!          call read_restore_s_interior_pnetcdf((iint+irst)/irst,nz,z0,
!     $                                                               f0)
!          call ztosig(z0,f0,zz,h,srstrf,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
!          call read_restore_tau_interior_pnetcdf((iint+irst)/irst,nz,
!     $                                                            z0,f0)
!          call ztosig(z0,f0,zz,h,taurstrf,im,jm,nz,kb,
!     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
        end if
      end if

! linear interpolation in time
      fnew=time/trst-ntime
      fold=1.-fnew
      do k=1,kbm1
        do i=1,im
          do j=1,jm
            trstr(i,j,k)=fold*trstrb(i,j,k)+fnew*trstrf(i,j,k)
            srstr(i,j,k)=fold*srstrb(i,j,k)+fnew*srstrf(i,j,k)
            taurstr(i,j,k)=fold*taurstrb(i,j,k)+fnew*taurstrf(i,j,k)
          end do
        end do
      end do

! restore
      do k=1,kbm1
        do i=1,im
          do j=1,jm
            t(i,j,k)=t(i,j,k)+2.*dti/86400.*taurstr(i,j,k)*
     $                                           (trstr(i,j,k)-t(i,j,k))
            tb(i,j,k)=tb(i,j,k)+2.*dti/86400.*taurstr(i,j,k)*
     $                                          (trstr(i,j,k)-tb(i,j,k))
            s(i,j,k)=s(i,j,k)+2.*dti/86400.*taurstr(i,j,k)*
     $                                           (srstr(i,j,k)-s(i,j,k))
            sb(i,j,k)=sb(i,j,k)+2.*dti/86400.*taurstr(i,j,k)*
     $                                          (srstr(i,j,k)-sb(i,j,k))
          end do
        end do
      end do

! mask
      do k=1,kbm1
        t(:,:,k)  = t(:,:,k)*fsm
        tb(:,:,k) = tb(:,:,k)*fsm
        s(:,:,k)  = s(:,:,k)*fsm
        sb(:,:,k) = sb(:,:,k)*fsm
      end do

      return
      end
