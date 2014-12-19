module mydefs

!    integer, parameter :: rk = selected_real_kind(15,307)
    integer, parameter :: sp =kind(1.0)
    integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))
    integer, parameter :: rk =selected_real_kind(2*precision(1.0_dp))
    integer, parameter :: ik = selected_int_kind(5)
    integer, parameter :: iw1=15, iw2=16, iw3=17, iw4=18
    real(rk), parameter :: pi = 3.141592653589793238462643383279502884197_rk
!    integer, parameter :: maxatoms=3, maxmodes=3, step=50
    real(rk), parameter :: autoaa = 0.52917721092_rk!, ssize = 0.05000_rk
    

end module mydefs
