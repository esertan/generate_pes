program main
    
    use mydefs
    use pes_init
    use pes_calc

    implicit none
    
    real(rk), dimension(1:maxatoms,1:3) :: xyz0
    real(rk) :: r01, r02, a0
    real(rk), dimension(1:maxmodes,1:maxatoms,1:3) :: dxyz
    real(rk), dimension(1:maxmodes,1:3) :: drtheta
    real(rk), dimension(1:maxmodes,1:step) :: r1, r2, a
    real(rk), dimension(1:maxmodes,1:maxatoms,1:3,1:step) :: xyz 
   
    call init_geo(xyz0)

    call init_internal_geo(xyz0, r01, r02, a0)

    call read_cartd(dxyz)

    call read_interd(drtheta)

    call make_internal(r01, r02, a0, r1, r2, a, drtheta)

    call make_cartesian(xyz, xyz0, dxyz)

    call write_grid

    call write_cartesian(xyz)

    call write_internal(r1,r2,a)

end program main
