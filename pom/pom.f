! pom.f

! main program

      program pom
      implicit none
      include 'pom.h'
      integer date_time(8)
      character*10 tmp(3)

      call date_and_time(tmp(1),tmp(2),tmp(3),date_time)

! initialize model
      call initialize

! main loop
      do iint=1,iend
        call advance ! advance model
      end do

! final step: print computing time
      if(my_task.eq.master_task) then
        write (*,'(/a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')
     $                                        'Initial date and time: ',
     $              date_time(1),'/',date_time(2),'/',date_time(3),'  ',
     $                    date_time(5),':',date_time(6),':',date_time(7)
        call date_and_time(tmp(1),tmp(2),tmp(3),date_time)
        write (*,'(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')
     $                                        'Final date and time:   ',
     $              date_time(1),'/',date_time(2),'/',date_time(3),'  ',
     $                    date_time(5),':',date_time(6),':',date_time(7)
        write (*,'(/a/)') 'POM terminated successfully '
      end if

! finalize mpi
      call finalize_mpi

      stop
      end
