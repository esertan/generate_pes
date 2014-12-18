program main
    
    use mydefs
    use pes_init
    use pes_calc
    use pes_write

    implicit none
    
    integer(ik) :: maxatoms, maxmodes, step
    real(rk) :: xmin, ssize!, ar1, ar2, br1, amin, angstep, angsize, angle, rmin, rstep, rsize
    real(rk) :: r01, r02, a0
    real(rk), allocatable :: xyz0(:,:)
    real(rk), allocatable :: dxyz(:,:,:)
    real(rk), allocatable :: drtheta(:,:)
    real(rk), allocatable :: xyz(:,:,:,:)
    real(rk), allocatable :: r1(:,:)
    real(rk), allocatable :: r2(:,:)
    real(rk), allocatable :: a(:,:)
    character(len=3), allocatable :: aname(:)
    character(len=5) :: grid='GRIDI', mol='MOLIN', geo='GEOME', disp='DISPL', inte='INTER'!, ang='ANGLE', bond='BONDS' 
    integer :: hflag, iflag

    iflag = read_key(mol)

    if(iflag==1) then

        call get_molinfo(maxatoms, maxmodes)

    end if

    iflag = read_key(grid)
    
    if(iflag==1) then

        call get_gridinfo(xmin, step, ssize)

    end if

    hflag = read_key(geo)

    if(hflag==1) then

        call allocate_arrays(xyz0, dxyz, xyz, drtheta, r1, r2, a, maxmodes, maxatoms, step, aname)        

        call init_geo(xyz0, maxatoms, aname)
        
        call init_internal_geo(xyz0, r01, r02, a0, maxatoms)

        iflag = read_key(disp)    
        
        if(iflag==1) then
        
            call read_cartd(dxyz, maxatoms, maxmodes)
        
            call make_cartesian(xyz, xyz0, dxyz, maxatoms, maxmodes, step, xmin, ssize)
        
            call write_cartesian(xyz, maxatoms, maxmodes, step, aname) 
        
        end if

        iflag = read_key(inte)
    
        if(iflag==1) then
        
            call read_interd(drtheta, maxmodes)

            call make_internal(r01, r02, a0, r1, r2, a, drtheta, maxatoms, maxmodes, step, xmin, ssize)

            call write_internal(r1,r2,a, maxatoms, maxmodes, step,aname)
        
        end if

        call write_grid(step, ssize, xmin)

        call deallocate_arrays(xyz0,dxyz, xyz, drtheta, r1, r2, a, maxmodes, maxatoms, step, aname)

    end if

!    iflag = read_key(ang)

!    if(iflag==1) then
!        call angle_bending(ar1, ar2, amin, angstep, angsize)
!        call make_angle_bending(amin, angstep, angsize, acoord)
!        call write_angle_bending(acoord)
!    end if

!    iflag = read_key(bond)
!
!    if(iflag==1) then
!        call bond_stretching(br1, angle, rmin, rstep, rsize)
!        call make_bond_stretching(rmin, rstep, rsize, rcoord)
!        call write_bond_stretching(rcoord)
!    end if

end program main
