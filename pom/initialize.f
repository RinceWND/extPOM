! initialize.f

! initialize POM: define constant, read initial values, grid and initial
! conditions

!_______________________________________________________________________
      subroutine initialize
! initialize POM
      implicit none
      include 'pom.h'
!      integer i,j,k

! initialize the MPI execution environment and create communicator for
!internal POM communications
      call initialize_mpi

! distribute the model domain across processors
      call distribute_mpi

! read input values and define constants
      call read_input

! initialize arrays for safety (may be overwritten)
      call initialize_arrays

! read in grid data
      call read_grid

! read initial and lateral boundary conditions
      call initial_conditions

! update initial conditions and set the remaining conditions
      call update_initial

! calculate the bottom friction coefficient
      call bottom_friction

! read restart data from a previous run
      if(nread_rst.ne.0) call read_restart_pnetcdf
      
!      call write_aux_pnetcdf
!      call finalize_mpi
!      stop

!      call write_aux_pnetcdf
!      call finalize_mpi
!      stop

! check for errors
      call sum0d_mpi(error_status,master_task)
      call bcast0d_mpi(error_status,master_task)
      if(error_status.ne.0) then
        if(my_task.eq.master_task) write(*,'(/a)')
     $                                       'POM terminated with error'
        call finalize_mpi
        stop
      end if

! write grid and initial conditions
      if(netcdf_file.ne.'nonetcdf') call write_output_pnetcdf
      if(my_task.eq.master_task) write(*,'(/a)') 'End of initialization'

      return
      end

!_______________________________________________________________________
      subroutine read_input
! read input values and defines constants
      implicit none
      include 'pom.h'
      namelist/pom_nml/ title,wrk_pth,netcdf_file,mode,nadv,nitera,
     $                  sw,npg,dte,isplit,time_start,nread_rst,
     $                  read_rst_file,cont_bry,write_rst,write_rst_file,
     $                  days,prtd1,prtd2,swtch,ntp,nbct,nbcs

! Input of filenames and constants

! Logical for inertial ramp (.true. if inertial ramp to be applied
! to wind stress and baroclinic forcing, otherwise .false.)
      lramp=.false.

! Reference density (recommended values: 1025 for seawater,
! 1000 for freswater; S.I. units):
      rhoref=1025.d0

! Temperature bias (deg. C)
      tbias=0.d0

! Salinity bias
      sbias=0.d0

! gravity constant (S.I. units)
      grav=9.806d0

! von Karman's constant
      kappa=0.4d0

! Bottom roughness (metres)
      z0b=.01d0

! Minimum bottom friction coeff.
      cbcmin=.0025d0

! Maximum bottom friction coeff.
      cbcmax=1.d0

! Smagorinsky diffusivity coeff.
      horcon=0.1d0

! Inverse horizontal turbulent Prandtl number (ah/am; dimensionless):
! NOTE that tprni=0.d0 yields zero horizontal diffusivity!
      tprni=.1d0

! Background viscosity used in subroutines profq, proft, profu and
! profv (S.I. units):
      umol=1.d-6

! Maximum magnitude of vaf (used in check that essentially tests
! for CFL violation):
      vmaxl=100.d0

! Maximum allowable value of:
!   <difference of depths>/<sum of depths>
! for two adjacent cells (dimensionless). This is used in subroutine
! slpmax. If >= 1, then slpmax is not applied:
      slmax=2.d0

! Water type, used in subroutine proft.
!    ntp    Jerlov water type
!     1            i
!     2            ia
!     3            ib
!     4            ii
!     5            iii
      ntp=2

! Surface temperature boundary condition, used in subroutine proft:
!    nbct   prescribed    prescribed   short wave
!           temperature      flux      penetration
!     1        no           yes           no
!     2        no           yes           yes
!     3        yes          no            no
!     4        yes          no            yes
      nbct=1

! Surface salinity boundary condition, used in subroutine proft:
!    nbcs   prescribed    prescribed
!            salinity      flux
!     1        no           yes
!     3        yes          no
! NOTE that only 1 and 3 are allowed for salinity.
      nbcs=1

! Step interval during which external (2-D) mode advective terms are
! not updated (dimensionless):
      ispadv=1

! Constant in temporal filter used to prevent solution splitting
! (dimensionless):
      smoth=0.10d0

! Weight used for surface slope term in external (2-D) dynamic
! equation (a value of alpha = 0.d0 is perfectly acceptable, but the
! value, alpha=.225d0 permits a longer time step):
      alpha=0.d0

! Initial value of aam:
      aam_init=0.d0

! End of input of constants

! read input namelist (overwrites vars)
      open(73,file='pom.nml',status='old')
      read(73,nml=pom_nml)
      close(73)

! calculate some constants
      small=1.d-9           ! Small value
      pi=atan(1.d0)*4.d0    ! PI

      dti=dte*float(isplit)
      dte2=dte*2
      dti2=dti*2

      iend=max0(nint(days*24.d0*3600.d0/dti),2)
      iprint=nint(prtd1*24.d0*3600.d0/dti)
      iswtch=nint(swtch*24.d0*3600.d0/dti)
      irestart=nint(write_rst*24.d0*3600.d0/dti)

      ispi=1.d0/float(isplit)
      isp2i=1.d0/(2.d0*float(isplit))

! initialise time
      time0=0.d0
      time=0.d0

! check cont_bry
      if (nread_rst == 0) cont_bry = 0

! print initial summary
      if(my_task.eq.master_task) then
        write(6,'(/'' title      = '',a40)') title
        write(6,'(/'' mode       = '',i10)') mode
        write(6,'('' nadv       = '',i10)') nadv
        write(6,'('' nitera     = '',i10)') nitera
        write(6,'('' sw         = '',f10.4)') sw
        write(6,'('' nread_rst  = '',i10)') nread_rst
        write(6,'('' write_rst  = '',f10.4)') write_rst
        write(6,'('' irestart   = '',i10)') irestart
        write(6,'('' dte        = '',f10.2)') dte
        write(6,'('' dti        = '',f10.1)') dti
        write(6,'('' isplit     = '',i10)') isplit
        write(6,'('' time_start = '',a26)') time_start
        write(6,'('' days       = '',f10.4)') days
        write(6,'('' iend       = '',i10)') iend
        write(6,'('' prtd1      = '',f10.4)') prtd1
        write(6,'('' iprint     = '',i10)') iprint
        write(6,'('' prtd2      = '',f10.4)') prtd2
        write(6,'('' swtch      = '',f10.2)') swtch
        write(6,'('' iswtch     = '',i10)') iswtch
        write(6,'('' lramp      = '',l10)') lramp
        write(6,'('' rhoref     = '',f10.3)') rhoref
        write(6,'('' tbias      = '',f10.3)') tbias
        write(6,'('' sbias      = '',f10.3)') sbias
        write(6,'('' grav       = '',f10.4)') grav
        write(6,'('' kappa      = '',f10.4)') kappa
        write(6,'('' z0b        = '',f10.6)') z0b
        write(6,'('' cbcmin     = '',f10.6)') cbcmin
        write(6,'('' cbcmax     = '',f10.6)') cbcmax
        write(6,'('' horcon     = '',f10.3)') horcon
        write(6,'('' tprni      = '',f10.4)') tprni
        write(6,'('' umol       = '',f10.4)') umol
        write(6,'('' vmaxl      = '',f10.4)') vmaxl
        write(6,'('' slmax      = '',f10.4)') slmax
        write(6,'('' ntp        = '',i10)') ntp
        write(6,'('' nbct       = '',i10)') nbct
        write(6,'('' nbcs       = '',i10)') nbcs
        write(6,'('' ispadv     = '',i10)') ispadv
        write(6,'('' smoth      = '',f10.4)') smoth
        write(6,'('' alpha      = '',f10.4)') alpha
      end if

      return
      end

!_______________________________________________________________________
      subroutine initialize_arrays
! initialize arrays for safety
      implicit none
      include 'pom.h'
!      integer i,j,k

! boundary arrays
      vabn = 0.d0
      vabs = 0.d0
      uabe = 0.d0
      uabw = 0.d0
      
      eln  = 0.d0
      els  = 0.d0
      ele  = 0.d0
      elw  = 0.d0
      
      vbn  = 0.d0
      vbs  = 0.d0
      tbn  = 0.d0
      tbs  = 0.d0
      sbn  = 0.d0
      sbs  = 0.d0

      ube  = 0.d0
      ubw  = 0.d0
      tbe  = 0.d0
      tbw  = 0.d0
      sbe  = 0.d0
      sbw  = 0.d0

      fluxua = 0.d0
      fluxva = 0.d0
! 2-D arrays
      uab     = 0.d0
      vab     = 0.d0
      elb     = 0.d0
      etb     = 0.d0
      e_atmos = 0.d0
      vfluxb  = 0.d0
      vfluxf  = 0.d0
      wusurf  = 0.d0
      wusurff = 0.d0
      wusurfb = 0.d0
      wvsurf  = 0.d0
      wvsurff = 0.d0
      wvsurfb = 0.d0
      wtsurf  = 0.d0
      wtsurff = 0.d0
      wtsurfb = 0.d0
      wssurf  = 0.d0
      wssurff = 0.d0
      wssurfb = 0.d0
      swrad   = 0.d0
      drx2d   = 0.d0
      dry2d   = 0.d0
! 3-D arrays
      ub = 0.d0
      vb = 0.d0
      
      drhox = 0.d0
      drhoy = 0.d0

      return
      end

!_______________________________________________________________________
      subroutine read_grid
! set up vertical and horizontal grid, topography, areas and masks
      implicit none
      include 'pom.h'
      integer i,j,k
      double precision deg2rad

! degrees to radians
      deg2rad=pi/180.

! read grid
      call read_grid_pnetcdf

! derived vertical grid variables
      do k=1,kb-1
        dz(k)=z(k)-z(k+1)
        dzz(k)=zz(k)-zz(k+1)
      end do
      dz(kb)=0.
      dzz(kb)=0.

! print vertical grid information
      if(my_task.eq.master_task) then
        write(6,'(/2x,a,7x,a,9x,a,9x,a,9x,a)') 'k','z','zz','dz','dzz'
        do k=1,kb
          write(6,'(1x,i5,4f10.3)') k,z(k),zz(k),dz(k),dzz(k)
        end do
      end if

! set up Coriolis parameter
        do j=1,jm
          do i=1,im
            cor(i,j)=2.*7.29d-5*sin(north_e(i,j)*deg2rad)
          end do
        end do

! inertial period for temporal filter
      if (cor(im/2,jm/2) == 0.) then
        stop 'Coriolis problem. Division by zero coriolis.'
      else
        period=(2.d0*pi)/abs(cor(im/2,jm/2))/86400.d0
      end if

! calculate areas of "t" and "s" cells
      art = dx*dy

! calculate areas of "u" and "v" cells
      do j=2,jm
        do i=2,im
          aru(i,j)=.25d0*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j))
          arv(i,j)=.25d0*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1))
        end do
      end do
      call exchange2d_mpi(aru,im_local,jm_local)
      call exchange2d_mpi(arv,im_local,jm_local)

      if (n_west.eq.-1) then
        aru(1,:)=aru(2,:)
        arv(1,:)=arv(2,:)
      end if

      if (n_south.eq.-1) then
        aru(:,1)=aru(:,2)
        arv(:,1)=arv(:,2)
      end if

      do i=1,im
        do j=1,jm
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
        end do
      end do

      call check_cflmin_mpi

      return
      end

!_______________________________________________________________________
      subroutine initial_conditions
! set up initial and lateral boundary conditions
      implicit none
      include 'pom.h'
      integer nz
      parameter(nz=40)
      integer i,j,k
!      double precision sum1,sum2
!      double precision tb0(im,jm,nz),sb0(im,jm,nz)
!      double precision tb2(im,jm,nz),sb2(im,jm,nz)
!      double precision p1(im,jm,kb),p2(im,jm,kb)
!      double precision z2(nz)

! read initial temperature and salinity from ic file
      call read_initial_ts_pnetcdf(kb,tb,sb)
      call read_clim_ts_pnetcdf(kb,10,tclim,sclim)

! map onto sigma coordinate
!      call ztosig(z2,tb0,zz,h,tclim,im,jm,nz,kb,
!     $                                    n_west,n_east,n_south,n_north)
!      call ztosig(z2,sb0,zz,h,sclim,im,jm,nz,kb,
!     $                                    n_west,n_east,n_south,n_north)

! mean density
      call dens(sclim,tclim,rmean)

! map onto sigma coordinate
!      call ztosig(z2,tb2,zz,h,tb,im,jm,nz,kb,
!     $                                    n_west,n_east,n_south,n_north)
!      call ztosig(z2,sb2,zz,h,sb,im,jm,nz,kb,
!     $                                    n_west,n_east,n_south,n_north)

! density
      call dens(sb,tb,rho)

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
          end do
        end do
      end do

! iniital heat and water fluxes (see subroutine surface_forcing)
      do i=1,im
        do j=1,jm
          tsurf(i,j)=tb(i,j,1)
          ssurf(i,j)=sb(i,j,1)
        end do
      end do

! lateral boundary conditions
! boundary conditions are variable (see subroutine lateral_bc)
      rfe=1.d0
      rfw=1.d0
      rfn=1.d0
      rfs=1.d0

      do k=1,kbm1
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine update_initial
! update the initial conditions and set the remaining initial conditions
      implicit none
      include 'pom.h'
      integer i,j,k

      do i=1,im
        do j=1,jm
          ua(i,j)=uab(i,j)
          va(i,j)=vab(i,j)
          el(i,j)=elb(i,j)
          et(i,j)=etb(i,j)
          etf(i,j)=et(i,j)
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
          w(i,j,1)=vfluxf(i,j)
        end do
      end do

      do k=1,kb
        do j=1,jm
          do i=1,im
            l(i,j,k)=0.1*dt(i,j)
            q2b(i,j,k)=small
            q2lb(i,j,k)=l(i,j,k)*q2b(i,j,k)
            kh(i,j,k)=l(i,j,k)*sqrt(q2b(i,j,k))
            km(i,j,k)=kh(i,j,k)
            kq(i,j,k)=kh(i,j,k)
            aam(i,j,k)=aam_init
          end do
        end do
      end do

      do k=1,kbm1
        do i=1,im
          do j=1,jm
            q2(i,j,k)=q2b(i,j,k)
            q2l(i,j,k)=q2lb(i,j,k)
            t(i,j,k)=tb(i,j,k)
            s(i,j,k)=sb(i,j,k)
            u(i,j,k)=ub(i,j,k)
            v(i,j,k)=vb(i,j,k)
          end do
        end do
      end do

      if (npg.eq.1) then
        call baropg
      else if (npg.eq.2) then
        call baropg_mcc
      else
        error_status=1
        write(6,'(/''Error: invalid value for npg'')')
      end if

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
            dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine bottom_friction
! calculate the bottom friction coefficient
      implicit none
      include 'pom.h'
      integer i,j
!      double precision TCB(im,jm)
      
! calculate bottom friction
      !TCB(1:im,1:jm) = h(1:im,1:jm)
      
      do i=1,im
        do j=1,jm
          cbc(i,j)=(kappa/log((1.+zz(kbm1))*h(i,j)/z0b))**2
          cbc(i,j)=max(cbcmin,cbc(i,j))
! if the following is invoked, then it is probable that the wrong
! choice of z0b or vertical spacing has been made:
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do
      return
      end

!_______________________________________________________________________
      subroutine ztosig(zs,tb,zz,h,t,im,jm,ks,kb,
     $                  im_local,jm_local,n_west,n_east,n_south,n_north)
! interpolate vertically
      implicit none
      integer im,jm,ks,kb,im_local,jm_local
      double precision zs(ks),tb(im,jm,ks),zz(kb),h(im,jm),t(im,jm,kb),
     $                                          tin(ks),tout(kb),zzh(kb)
      integer n_west,n_east,n_south,n_north
      double precision tmax
      integer i,j,k
      
      t = 0.d0

      do i=2,im-1
      do j=2,jm-1
        if (h(i,j).gt.1.0) then
! special interp on z-lev for cases of no data because h smoothing
          do k=1,ks
            tin(k)=tb(i,j,k)
            if (zs(k).le.h(i,j) .and. tin(k).lt.0.01) then
              tmax=amax1(tb(i-1,j,k),tb(i+1,j,k),
     $                   tb(i,j-1,k),tb(i,j+1,k))
              tin(k)=tmax
            endif
            if (tin(k).lt.0.01 .and. k.ne.1) tin(k)=tin(k-1)
          end do

          do k=1,kb
            zzh(k)=-zz(k)*h(i,j)
          end do

! vertical spline interp
          call splinc(zs,tin,ks,2.d30,2.d30,zzh,tout,kb)

          do k=1,kb
              t(i,j,k)=tout(k)
          end do

        end if
      end do
      end do
      call exchange3d_mpi(t,im_local,jm_local,kb)

! boundaries
      do k=1,kb
        do j=1,jm
          if(n_west.eq.-1) t(1,j,k)=t(2,j,k)
          if(n_east.eq.-1) t(im,j,k)=t(im-1,j,k)
        end do
        do i=1,im
          if(n_south.eq.-1) t(i,1,k)=t(i,2,k)
          if(n_north.eq.-1) t(i,jm,k)=t(i,jm-1,k)
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine splinc(x,y,n,yp1,ypn,xnew,ynew,m)
! interpolate using splines
      double precision :: x, y, yp1, u, sig, p, y2, un, xnew, ynew, ypn
      parameter (nmax=300)
      dimension x(n),y(n),y2(nmax),u(nmax),xnew(m),ynew(m)

      if (yp1.gt..99d30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     $      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do

      if (ypn.gt..99d30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do

      do i=1,m
        call splint(x,y,y2,n,xnew(i),ynew(i))
      end do

      return
      end

!_______________________________________________________________________
      subroutine splint(xa,ya,y2a,n,x,y)
      
      double precision :: xa, ya, y2a, x, y, h, a, b
      dimension xa(n),ya(n),y2a(n)

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
        goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.)  then
        error_staus=1
        write(6,'(/a)') 'Error: bad xa input in splint'
      end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     $      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end
