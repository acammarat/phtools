subroutine show_logo
  use var, only: version

  write(*,'(a)')     "         _                      _                " 
  write(*,'(a)')     "   _ __ | |__   ___  _ __   ___| |__   __ _ _ __ "
  write(*,'(a)')     "  | '_ \| '_ \ / _ \| '_ \ / __| '_ \ / _` | '__|"
  write(*,'(a)')     "  | |_) | | | | (_) | | | | (__| | | | (_| | |   "
  write(*,'(a)')     "  | .__/|_| |_|\___/|_| |_|\___|_| |_|\__,_|_|   "
  write(*,'(a)')     "  |_|                                            "
  write(*,'(*(a))')  "                  ",version
  write(*,*)
         
  
  return
end subroutine show_logo
