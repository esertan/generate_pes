module pes_init

    use mydefs

    implicit none

contains

    integer function read_key(checkkey)

    character(len=*), intent(in) :: checkkey
    character(len=5) :: key
!    integer, intent(inout) :: flag=0
    
    read(*,*) key
    write(*,*) key
    
    if(key==checkkey) then
        read_key=1
    else
        read_key=0
    end if

!    write(*,*) flag
        
    end function read_key

    subroutine get_molinfo(maxatoms, maxmodes)

    ! Reads number of atoms and number of vibration modes from input

    integer(ik), intent(inout) :: maxatoms, maxmodes

    read(*,*) maxatoms, maxmodes

    end subroutine get_molinfo

    subroutine get_gridinfo(xmin, step, ssize)

    ! Reads parameters xmin, #steps, step size for the distortion gird from input.

    integer(ik), intent(inout) :: step
    real(rk), intent(inout) :: ssize, xmin

    read(*,*) xmin, step, ssize

    end subroutine get_gridinfo
    subroutine allocate_arrays(xyz0, dxyz, xyz, drtheta, r1, r2, a, maxmodes, maxatoms, step)

    ! Subroutine allocates the arrays for atom positions and dispacements.

    integer :: istat
    integer(ik) :: maxmodes, maxatoms, step
    real(rk), allocatable :: xyz0(:,:)
    real(rk), allocatable :: dxyz(:,:,:)
    real(rk), allocatable :: drtheta(:,:)
    real(rk), allocatable :: xyz(:,:,:,:)
    real(rk), allocatable :: r1(:,:)
    real(rk), allocatable :: r2(:,:)
    real(rk), allocatable :: a(:,:)

    allocate(xyz0(1:maxatoms,1:3), STAT=istat)
    allocate(dxyz(1:maxmodes,1:maxatoms,1:3), STAT=istat)
    allocate(drtheta(1:maxmodes,1:3), STAT=istat)
    allocate(xyz(1:maxmodes,1:maxatoms,1:3,1:step), STAT=istat)
    allocate(r1(1:maxmodes,1:step), STAT=istat)
    allocate(r2(1:maxmodes,1:step), STAT=istat)
    allocate(a(1:maxmodes,1:step), STAT=istat)

    end subroutine allocate_arrays

    subroutine deallocate_arrays(xyz0, dxyz, xyz, drtheta, r1, r2, a, maxmodes, maxatoms, step)

    ! Subroutine deallocates arrays that were allocated in the allocate_arrays subroutine.

    integer :: istat
    integer(ik) :: maxmodes, maxatoms, step
    real(rk), allocatable :: xyz0(:,:)
    real(rk), allocatable :: dxyz(:,:,:)
    real(rk), allocatable :: drtheta(:,:)
    real(rk), allocatable :: xyz(:,:,:,:)
    real(rk), allocatable :: r1(:,:)
    real(rk), allocatable :: r2(:,:)
    real(rk), allocatable :: a(:,:)

    deallocate(xyz0, STAT=istat)
    deallocate(dxyz, STAT=istat)
    deallocate(drtheta, STAT=istat)
    deallocate(xyz, STAT=istat)
    deallocate(r1, STAT=istat)
    deallocate(r2, STAT=istat)
    deallocate(a, STAT=istat)    


    end subroutine deallocate_arrays    
  
    subroutine init_geo(xyz0, maxatoms)
    
    ! Subroutine reads in initial geometry from stdin in cartesian coordinates (Angstrom).
    integer(ik) :: maxatoms
    real(rk), dimension(1:maxatoms, 1:3), intent(inout) :: xyz0
    integer :: i,j,k

!    write(*,*) 'Reading input geometry (a.u)'
    do i=1, maxatoms
        read(*,*) (xyz0(i,k), k=1,3)
        xyz0(i,:) = xyz0(i,:)/autoaa
!        write(*,*) (xyz0(i,k), k=1,3)
    end do

    
    end subroutine init_geo

    subroutine init_internal_geo(xyz0, r01, r02, a0, maxatoms)

    ! Converts the initial geometry to internal coordinates    

    integer(ik) :: maxatoms
    real(rk), dimension(1:maxatoms, 1:3), intent(in) :: xyz0
    real(rk), intent(inout) :: r01, r02, a0

    r01 = dsqrt(sum((xyz0(1,:)-xyz0(2,:))**2))
    r02 = dsqrt(sum((xyz0(3,:)-xyz0(2,:))**2))
    a0 = dacos(sum((xyz0(3,:)-xyz0(2,:))*(xyz0(1,:)-xyz0(2,:)))/r01/r02)

 !   write(*,*) 'Input geo in internal coordinates'
 !   write(*,*) r01, r02, a0
    end subroutine init_internal_geo

    subroutine read_cartd(dxyz, maxatoms, maxmodes)

    ! Subroutine reads in displacement vectors from stdin in atomic units.
    integer(ik) :: maxatoms, maxmodes
    real(rk), dimension(1:maxmodes, 1:maxatoms, 1:3), intent(inout) :: dxyz
    integer :: i,j,k

!    write(*,*) 'Reading cartesian displacement'
    do i=1, maxatoms
        do k=1,3
            read(*,*) (dxyz(j,i,k), j=1,maxmodes)
            do j=1, maxmodes
                select case(j)
                    case(2)
                        dxyz(j,i,k) = -dxyz(j,i,k)
                end select
            end do
!            write(*,*) (dxyz(j,i,k), j=1,maxmodes)
        end do
    end do    
    end subroutine read_cartd

    subroutine read_interd(drtheta, maxmodes)

    ! Subroutine reads in displacement in internal coordinates (atomic units).

    integer(ik) :: maxmodes
    real(rk), dimension(1:maxmodes, 1:3), intent(inout) :: drtheta
    integer :: i, j, k

!    write(*,*) 'Reading displacement internal coordinates'
    do j=1,maxmodes
        read(*,*) (drtheta(j,k), k=1,3)
        select case(j)
            case(2)
                do k=1,3
                    drtheta(j,k) = -drtheta(j,k)
                end do
            end select
!        write(*,*) (drtheta(j,k), k=1,3)
    end do
    
    end subroutine read_interd

end module pes_init
