module Date_Utility
  implicit none
  private
  public :: Date_since, Days_in_between, get_Month, Number_of_Days
  public :: get_Month_Int, get_Month_IntDays, get_Month_IntFactor
  public :: Days_to_Stamp, Days_to_TimeStamp, Days_plus, Date_to_Day_of_Year
  public :: day_of_year
  public :: T_TimeStamp, T_Zone, T_Date, T_Time
  integer*1, parameter :: N_MONTHS = 13
  integer*2, parameter :: MONTH_DAYS(N_MONTHS) = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 /)
!  integer*2, parameter :: YEAR_OFFSET = 1900
  type T_Date
  sequence
    integer*2 :: year  = 2000
    integer*1 :: month = 1
    integer*1 :: day   = 1
  end type T_Date
  
  type T_TimeShort
  sequence
    integer*1 :: hour   = 0
    integer*1 :: minute = 0
  end type T_TimeShort
  
  type T_Time
  sequence
    integer*1        :: hour    = 0
    integer*1        :: minute  = 0
    double precision :: seconds = 0.d0
  end type T_Time
  
  type T_Zone
  sequence
    logical   :: negative = .false.
    integer*1 :: hour     = 0
    integer*1 :: minute   = 0
  end type T_Zone
  
  type T_DateTime
  sequence
    type (T_Date) :: date
    type (T_Time) :: time
  end type T_DateTime
  
  type T_DateTimeShort
  sequence
    type (T_Date)      :: date
    type (T_TimeShort) :: time
  end type T_DateTimeShort
  
  type T_TimeStampShort
  sequence
    type (T_Date)      :: date
    type (T_TimeShort) :: time
    type (T_Zone)      :: zone
  end type T_TimeStampShort
  
  type T_TimeStamp
  sequence
    type (T_Date) :: date
    type (T_Time) :: time
    type (T_Zone) :: zone
  end type T_TimeStamp
  
  contains
!___________________________________________________________________
!
! Check for leap year
!
    logical function Is_Leap_Year(Year)
    
      integer*2, intent(in) :: Year
      Is_Leap_Year = .false.
      if ( ( mod( Year, 4 ) == 0 .and. mod( Year, 100 ) /= 0 ) .or. mod( Year, 400 ) == 0 ) Is_Leap_Year = .true.
      
    end function Is_Leap_Year
!___________________________________________________________________
!
! Get the number of a day in a year
!
!    function Date_to_Day_of_Year( Day_of_Month, Month, Year ) result(Day_Of_Year)
!      
!      integer*1, intent(in) :: Day_of_Month
!      integer*1, intent(in) :: Month
!      integer*2, intent(in) :: Year
!      integer*2 :: Day_of_Year
!      integer*2 :: Days_of_Month(N_MONTHS)
!
!      Day_Of_Year = -1
!! Compute days per month
!      Days_of_Month = MONTH_DAYS
!! Error checking
!      if ( Year < 1 ) return
!      if ( Month < 1 .or. Month > N_MONTHS ) return
!      if ( Day_of_Month > Days_of_Month(Month) ) return
!! Compute day of year
!      Day_of_Year = Days_of_Month(Month) + Day_of_Month
!      if ( Is_Leap_Year(Year) ) Day_of_Year = Day_of_Year + 1
!      
!    end function Date_to_Day_of_Year
!___________________________________________________________________
!
! Parse datetime string 'YYYY-MM-DD hh:ii:ss {+|-}HH:II'
!
!    double precision function get_Day_of_Year(Time)
!      
!      double precision, intent(in) :: Time
!      type(T_TimeStamp)            :: TimeStamp
!      type(T_Zone)                 :: Zone
!      double precision             :: TimeOfTheDay
!      
!      TimeStamp = Days_to_TimeStamp(Time, Zone)
!      TimeOfTheDay = Time - floor(Time)
!      
!      get_Day_Of_Year = -1.
!! Error checking
!      if ( TimeStamp%date%year < 1 ) return
!      if ( TimeStamp%date%month < 1 .or. TimeStamp%date%month > N_MONTHS-1 ) return
!! Compute day of year
!      get_Day_of_Year = MONTH_DAYS(TimeStamp%date%month) + TimeStamp%date%day
!      if ( Is_Leap_Year(TimeStamp%date%year) .and. TimeStamp%date%month > 2) get_Day_of_Year = get_Day_of_Year + 1
!      get_Day_of_Year = get_Day_of_Year + TimeOfTheDay
!    
!    end function get_Day_of_Year
!___________________________________________________________________
!
! Parse datetime string 'YYYY-MM-DD hh:ii:ss {+|-}HH:II'
!
    function TimeStamp_Parse(TimeStampString) result(TimeStamp)
        
      character(len=*), intent(in) :: TimeStampString
      type (T_TimeStamp) :: TimeStamp
      integer*1 :: offset
      integer   :: sec
      character :: sign
      
! Error checks
      if (len(TimeStampString)<10) then
        write(*,*) "[date_utility] /invalid timestamp format/"
        write(*,*) "    The TimeStamp string is too short."
        write(*,*) ""
        return
      end if
      
      read(TimeStampString, '(i4,1x,i2,1x,i2)') TimeStamp%date%year,  &
                                                TimeStamp%date%month, &
                                                TimeStamp%date%day
      offset = 18
      
      read(TimeStampString(offset-6:offset-2), '(i2,1x,i2)') TimeStamp%time%hour, TimeStamp%time%minute
! Detect a format possibility
      if (TimeStampString(offset-1:offset-1).eq.":") then
        read(TimeStampString(offset:offset+1), '(i2)') sec
        TimeStamp%time%seconds = real(sec, 8)
        offset = int(offset+3, 1)
      end if
      if (len(TimeStampString)>19) then
        read(TimeStampString(offset:offset+5), '(a1,i2,1x,i2)') sign, TimeStamp%zone%hour, TimeStamp%zone%minute
        if (sign=='-') TimeStamp%zone%negative = .true.
      end if
        
    end function TimeStamp_Parse
!___________________________________________________________________
!
! Build datetime string (with +00:00 timezone)
!
    character(len=26) function TimeStamp_Build(TimeStamp, mode)
        
      type (T_TimeStamp), intent(in) :: TimeStamp
      character*1 :: mode
      character*1 :: sign
        
      sign = "+"
      if ( TimeStamp%zone%negative ) sign = "-"
      if ( mode=="s" ) then
        write(TimeStamp_Build,                                      &
         '(i4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,a1,i0.2,a,i0.2,a3)')    &
                            TimeStamp%date%year,    "-",            &
                            TimeStamp%date%month,   "-",            &
                            TimeStamp%date%day,     " ",            &
                            TimeStamp%time%hour,    ":",            &
                            TimeStamp%time%minute,  " ",            &
                            sign, TimeStamp%zone%hour, ":",         &
                                  TimeStamp%zone%minute, "   "
      elseif ( mode=="d" ) then
        write(TimeStamp_Build, '(i4,i0.2,i0.2)')     &
                            TimeStamp%date%year,    &
                            TimeStamp%date%month,   &
                            TimeStamp%date%day
      elseif ( mode=="f" ) then
        write(TimeStamp_Build,                                      &
         '(i4,i0.2,i0.2,a,i0.2,i0.2,i0.2)')                         &
                            TimeStamp%date%year,                    &
                            TimeStamp%date%month,                   &
                            TimeStamp%date%day,     "_",            &
                            TimeStamp%time%hour,                    &
                            TimeStamp%time%minute,                  &
                            int(TimeStamp%time%seconds)
      else
        write(TimeStamp_Build,                                      &
         '(i4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,a1,i0.2,a,i0.2)')&
                            TimeStamp%date%year,    "-",            &
                            TimeStamp%date%month,   "-",            &
                            TimeStamp%date%day,     " ",            &
                            TimeStamp%time%hour,    ":",            &
                            TimeStamp%time%minute,  ":",            &
                            int(TimeStamp%time%seconds), " ",       &
                            sign, TimeStamp%zone%hour, ":",         &
                                  TimeStamp%zone%minute
      end if
        
    end function TimeStamp_Build
!___________________________________________________________________
!
! Build datetime string from days
!
    character(len=26) function TimeStamp(Days, TimeZone)
        
      double precision, intent(in) :: Days
      type (T_Zone) TimeZone
        
      TimeStamp = TimeStamp_Build(Days_to_TimeStamp(Days, TimeZone), "s")
        
    end function TimeStamp
!___________________________________________________________________
!
! Get number of days (derive from TimeStamp)
!
    double precision function days_number(TimeStamp)
    
      type (T_TimeStamp), intent(in) :: TimeStamp
      type (T_TimeStamp) :: TimeStampTmp
      integer :: era, yoe, doy, doe
      
      TimeStampTmp = TimeStamp
      if (TimeStampTmp%date%month <= 2) TimeStampTmp%date%year = int(TimeStampTmp%date%year-1, 2)
      if (TimeStampTmp%date%year < 0) TimeStampTmp%date%year = int(TimeStampTmp%date%year-399, 2)
      era = TimeStampTmp%date%year / 400
      yoe = TimeStampTmp%date%year - era*400
      
      if (TimeStampTmp%date%month > 2) then
        TimeStampTmp%date%month = int(TimeStampTmp%date%month - 3, 1)
      else
        TimeStampTmp%date%month = int(TimeStampTmp%date%month + 9, 1)
      end if
      doy = (153*(TimeStampTmp%date%month)+2)/5 + TimeStampTmp%date%day-1
      doe = yoe*365 + yoe/4 - yoe/100 + doy
      days_number = era*146097+doe
      
      days_number = days_number + time_to_days(TimeStampTmp%time)
      
      return
      
    end function days_number
!___________________________________________________________________
!
!
!
    double precision function Number_of_Days(TimeStampString)
    
      character(len=*), intent(in) :: TimeStampString
      type (T_TimeStamp) :: TimeStamp
      
      TimeStamp = TimeStamp_Parse(TimeStampString)

      Number_of_Days = days_number(TimeStamp)
      
    end function Number_of_Days
!
!    
    double precision function time_to_seconds(Time)
      type (T_Time), intent(in) :: Time
      
      time_to_seconds = Time%seconds + 60.*(Time%minute + Time%hour*60.)
      
    end function time_to_seconds
    
    double precision function time_to_days(Time)
      type (T_Time), intent(in) :: Time
      
      time_to_days = ((Time%seconds/60. + Time%minute)/60. + Time%hour)/24.
      
    end function time_to_days
    
    type (T_Time) function days_to_time(Days)
      double precision, intent(in) :: Days
      
      days_to_time%hour    = int(Days*24., 1)
      days_to_time%minute  = int((Days*24.-days_to_time%hour)*60., 1)
      days_to_time%seconds = ((Days*24.-days_to_time%hour)*60.-days_to_time%minute)*60.
      
    end function days_to_time
!
    function Days_to_TimeStamp(Days, TimeZone) result(TimeStamp)
      double precision, intent(in) :: Days
      type (T_Zone), intent(in) :: TimeZone
      type (T_TimeStamp) :: TimeStamp
!      double precision :: Time
!      integer*1 :: leap
      
      integer t, era, doe, yoe, doy
      
      t = int(Days)
      if (t < 0) t = t-146096
      era = t/146097
      doe = t - era*146097
      yoe = (doe - doe/1460 + doe/36524 - doe/146096)/365
      TimeStamp%date%year = int(yoe + era*400, 2)
      doy = doe - (365*yoe + yoe/4 - yoe/100)
      TimeStamp%date%month = int((5*doy + 2)/153, 1)
      TimeStamp%date%day = int(doy - (153*TimeStamp%date%month + 2)/5 + 1, 1)
      if (TimeStamp%date%month < 10) then
        TimeStamp%date%month = int(TimeStamp%date%month + 3, 1)
      else
        TimeStamp%date%month = int(TimeStamp%date%month - 9, 1)
      end if
      if (TimeStamp%date%month <= 2) TimeStamp%date%year = int(TimeStamp%date%year + 1, 2)
      
      TimeStamp%time%seconds = Days - floor(Days)
      TimeStamp%time%hour    = int(floor(TimeStamp%time%seconds*24.), 1)
      TimeStamp%time%minute  = int(floor((TimeStamp%time%seconds*24.-TimeStamp%time%hour)*60.), 1)
      TimeStamp%time%seconds = ((TimeStamp%time%seconds*24.-TimeStamp%time%hour)*60.-TimeStamp%time%minute)*60.
      
      TimeStamp%zone = TimeZone
      
    end function Days_to_TimeStamp
!______________________________________________________________________________
!
    character(len=26) function Days_to_Stamp(Days, form)
      
      double precision, intent(in) :: Days
      character(len=1), intent(in) :: form
      type (T_TimeStamp) :: TimeStamp
      type (T_Zone)      :: TimeZone
      
      TimeStamp = Days_to_TimeStamp(Days, TimeZone)

      Days_to_Stamp = TimeStamp_Build(TimeStamp, form)
      
    end function
!______________________________________________________________________________
!
    function Days_since_to_TimeStamp(TimeStartString, Time) result(TimeStamp)
      character(len=*), intent(in) :: TimeStartString
      double precision, intent(in) :: Time
      double precision :: NewTime
      type (T_TimeStamp) :: TimeStamp
      
      TimeStamp = TimeStamp_Parse(TimeStartString)
      NewTime = days_number(TimeStamp)+Time
      TimeStamp = Days_to_TimeStamp(NewTime, TimeStamp%zone)
      
    end function Days_since_to_TimeStamp
!
! [ Exported function ]
! 
!  Returns timestamp string of date plus specified amount of days.
!
    character(len=26) function Date_since(TimeStartString, Time, Mode) result(TimeEndString)
      character(len=*), intent(in) :: TimeStartString
      double precision, intent(in) :: Time
      character,        intent(in) :: Mode
      type (T_TimeStamp) :: TimeStamp
      
      TimeStamp = Days_since_to_TimeStamp(TimeStartString, Time)
      TimeEndString = TimeStamp_Build(TimeStamp, Mode)
    end function Date_since
!    
    double precision function Days_in_Between(time_start, time_end)
      character(len=*), intent(in) :: time_start, time_end
      type (T_TimeStamp) :: TimeStart, TimeEnd
      
      TimeStart = TimeStamp_Parse(time_start)
      TimeEnd   = TimeStamp_Parse(time_end)
      
      Days_in_Between = Number_of_Days(time_end) - Number_of_Days(time_start)
      
    end function Days_in_Between
!
!
!    character(len=26) function TimeStamp_plus(TimeStartString, PlusString) result(TimeEndString)
!    
!      character(len=*), intent(in) :: TimeStartString, PlusString
!      double precision :: days
!      integer :: DD, hh, ii, ss
!      type (T_TimeStamp) :: TimeStart, Plus
!      
!      TimeStart = TimeStamp_Parse(TimeStartString)
!      Plus      = TimeStamp_Parse(PlusString)
!      
!      ! Add year and month first, to avoid leap year inconsistencies
!      TimeStart%date%year  = TimeStart%date%year  + Plus%date%year
!      TimeStart%date%month = TimeStart%date%month + Plus%date%month
!      DD = Plus%date%day
!      hh = TimeStart%time%hour   + Plus%time%hour
!      ii = TimeStart%time%minute + Plus%time%minute
!      ss = TimeStart%time%seconds+ Plus%time%seconds
!      
!      TimeStart%date%year = TimeStart%date%year + TimeStart%date%month/12
!      TimeStart%date%month= modulo(TimeStart%date%month, 12)
!      
!      days = real(DD, 8) + ( ( real(ss, 8)/60. + real(ii, 8) )/60. + real(hh, 8) )/24.
!      
!      days = days_number(TimeStart)+days
!      TimeStart = Days_to_TimeStamp(days, TimeStart%zone)
!      
!      TimeEndString = TimeStamp_Build(TimeStart, 'l')
!      
!      return 
!    
!    end function TimeStamp_plus
!
!
    double precision function Days_plus(TimeStartString, PlusString)
    
      character(len=*), intent(in) :: TimeStartString, PlusString
      double precision :: days_start, days
      integer :: DD, hh, ii, ss
      type (T_TimeStamp) :: TimeStart, Plus
      
      TimeStart = TimeStamp_Parse(TimeStartString)
      Plus      = TimeStamp_Parse(PlusString)
      !write(*,*) "[ö] TS: ", TimeStart, Plus
      
      days_start = days_number(TimeStart)
      !write(*,*) "[ö] DS: ", days_start
      ! Add year and month first, to avoid leap year inconsistencies
      TimeStart%date%year  = TimeStart%date%year  + Plus%date%year
      TimeStart%date%month = TimeStart%date%month + Plus%date%month
      DD = Plus%date%day
      hh = Plus%time%hour
      ii = Plus%time%minute
      ss = int(Plus%time%seconds)
      
      TimeStart%date%year = TimeStart%date%year + int(TimeStart%date%month/12, 2)
      TimeStart%date%month= int(modulo(TimeStart%date%month, 12), 1)
      
      days = real(DD, 8) + ( ( real(ss, 8)/60. + real(ii, 8) )/60. + real(hh, 8) )/24.
      
      days = days_number(TimeStart) + days
      !write(*,*) "[ö] D: ", days
      
      Days_plus = days - days_start
      
      return 
    
    end function Days_plus
!___________________________________________________________________________________
!
    double precision function Date_to_Day_of_Year(TimeStampString)
    
      character(len=*), intent(in) :: TimeStampString
      type(T_TimeStamp) :: TimeStamp
      
      TimeStamp = TimeStamp_Parse(TimeStampString)
                
      Date_to_Day_of_Year = day_of_year(TimeStamp)
      
    end function Date_to_Day_of_Year
!___________________________________________________________________________________
!
    double precision function day_of_year(TimeStamp)
    
      type(T_TimeStamp), intent(in) :: TimeStamp
      type(T_TimeStamp)             :: TimeStampTmp
      double precision  :: tim
      integer :: doy
      
      TimeStampTmp = TimeStamp
      if (TimeStampTmp%date%month > 2) then
        TimeStampTmp%date%month = int(TimeStampTmp%date%month - 3, 1)
      else
        TimeStampTmp%date%month = int(TimeStampTmp%date%month + 9, 1)
      end if
      doy = (153*(TimeStampTmp%date%month)+2)/5 + TimeStampTmp%date%day-1
      if (doy>=306) then
        doy = doy-305
      else
        doy = doy+60
      end if
      if (Is_Leap_Year(TimeStampTmp%date%year).and.TimeStampTmp%date%month<10) doy=doy+1
      
      tim = ( ( TimeStampTmp%time%seconds/60. + &
                TimeStampTmp%time%minute )/60. + &
                TimeStampTmp%time%hour )/24.
                
      day_of_year = real(doy, 8) + tim
      
    end function day_of_year
!___________________________________________________________________________________
!
    integer function get_Month(Time)
      
      double precision, intent(in) :: Time
      type (T_TimeStamp)           :: TimeStamp
      type (T_Zone)                :: Zone
      
      TimeStamp = Days_to_TimeStamp(Time, Zone)
      
      get_Month = TimeStamp%date%month
      
    end function get_Month
!
!
!
    double precision function get_Month_IntFactor(TimeStampString)
      
      character(len=*), intent(in) :: TimeStampString
      type (T_TimeStamp)           :: TimeStamp
      double precision             :: Days, lb, rb
      integer*1 :: lbi, rbi
      integer*2 :: dom(N_MONTHS), year
      double precision :: mom(N_MONTHS+1)
      
      dom = MONTH_DAYS
      
      read(TimeStampString, '(i4)') year
      Days = Date_to_Day_of_Year(TimeStampString)
      
      if (Is_Leap_Year(year)) then
        dom(3:N_MONTHS) = dom(3:N_MONTHS)+1
      end if
      mom(2:N_MONTHS) = .5*(dom(1:N_MONTHS-1)+dom(2:N_MONTHS))
      mom(1) = mom(2)-31.
      mom(N_MONTHS+1) = mom(N_MONTHS)+31.
      
      rbi = modulo(minloc(mom, 1, mom>=Days)-1, 12)+1
      lbi = modulo(maxloc(mom, 1, mom<Days)-1, 12)+1
      
      rb  = mom(rbi)
      lb  = mom(lbi)
      
      get_Month_IntFactor = (Days-lb)/(rb-lb)
      
    end function
    
    integer*1 function get_Month_Int(TimeStampString)
      
      character(len=*), intent(in) :: TimeStampString
      double precision             :: Days
      integer*2 :: dom(N_MONTHS), year
      double precision :: mom(N_MONTHS+1)
      
      dom = MONTH_DAYS
      
      read(TimeStampString, '(i4)') year
      Days = Date_to_Day_of_Year(TimeStampString)
      
      if (Is_Leap_Year(year)) then
        dom(3:N_MONTHS) = dom(3:N_MONTHS)+1
      end if
      mom(2:N_MONTHS) = .5*(dom(1:N_MONTHS-1)+dom(2:N_MONTHS))
      mom(1) = mom(2)-31.
      mom(N_MONTHS+1) = mom(N_MONTHS)+31.
      
      get_Month_Int = modulo(maxloc(mom, 1, mom<=Days)-1, 12)+1
      
    end function
    
    double precision function get_Month_IntDays(TimeStampString)
      
      character(len=*), intent(in) :: TimeStampString
      character(len=26)            :: TimeString
      type (T_TimeStamp)           :: TimeStamp, TS1
      double precision             :: Days, Time
      integer*1 :: ind
      integer*2 :: dom(N_MONTHS), year
      double precision :: mom(N_MONTHS+1)
      
      dom = MONTH_DAYS
      
      TimeStamp = TimeStamp_Parse(TimeStampString)
!      write(*,*) "01: ", TimeStamp
      TS1 = TimeStamp
      Days = Date_to_Day_of_Year(TimeStampString)
!      write(*,*) "02: ", Days
      
      if (Is_Leap_Year(TimeStamp%date%year)) then
        dom(3:N_MONTHS) = dom(3:N_MONTHS)+1_2
      end if
!      write(*,*) "03: ", dom
      mom(2:N_MONTHS) = .5*(dom(1:N_MONTHS-1)+dom(2:N_MONTHS))
      mom(1) = mom(2)-31.
      mom(N_MONTHS+1) = mom(N_MONTHS)+31.

!      write(*,*) "04: ", mom
      
      ind = modulo(maxloc(mom, 1, mom<Days)-1, 12)+1
!      write(*,*) "05: ", ind
      
      get_Month_IntDays = mom(ind+1)-Days
!      write(*,*) "25: ", get_Month_IntDays
      
    end function
    
!    integer function Days_in_Between_I_S(time_start, time_end)
!      double precision, intent(in) :: time_start
!      character(len=*), intent(in) :: time_end
!      
!      Days_in_Between = abs(Number_of_Days(time_end) - time_start)
!      
!    end function Days_in_Between_I_S

end module Date_Utility