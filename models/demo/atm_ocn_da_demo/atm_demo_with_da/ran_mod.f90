Module ran_mod

  Implicit None
! ran return a uniform random number between 0-1  
! norma return a normal distribution  
contains 
  function ran()   !returns random number between 0 - 1  
    implicit none 
    integer , save :: flag = 0
    real*4 :: ran 
    if(flag==0) then 
      call random_seed()
      !flag = 1 
    endif 
    call random_number(ran)    
  end function ran
  
  function normal(mean,sigma) 
    implicit none 
    integer :: flag 
    real*4, parameter :: pi = 3.141592653589793239  
    real*4 :: u1, u2, y1, y2, normal, mean, sigma 
    save flag 
    data flag /0/ 
    u1 = ran(); u2 = ran() 
    if (flag.eq.0) then 
      y1 = sqrt(-2.0*log(u1))*cos(2.0*pi*u2) 
      normal = mean + sigma*y1 
      flag = 1 
    else 
      y2 = sqrt(-2.0*log(u1))*sin(2.0*pi*u2) 
      normal = mean + sigma*y2 
      flag = 0 
    endif  
  end function normal 

End Module ran_mod