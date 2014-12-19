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

    ! Reads parameters xmin, #steps, step size for the distortion grid from input.

    integer(ik), intent(inout) :: step
    real(rk), intent(inout) :: ssize, xmin

    read(*,*) xmin, step, ssize

    end subroutine get_gridinfo
    subroutine allocate_arrays(xyz0, dxyz, xyz, drtheta, r1, r2, a, maxmodes, maxatoms, step, aname)

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
    character(len=3), allocatable :: aname(:)

    allocate(xyz0(1:maxatoms,1:3), STAT=istat)
    allocate(dxyz(1:maxmodes,1:maxatoms,1:3), STAT=istat)
    allocate(drtheta(1:maxmodes,1:3), STAT=istat)
    allocate(xyz(1:maxmodes,1:maxatoms,1:3,1:step), STAT=istat)
    allocate(r1(1:maxmodes,1:step), STAT=istat)
    allocate(r2(1:maxmodes,1:step), STAT=istat)
    allocate(a(1:maxmodes,1:step), STAT=istat)
    allocate(aname(1:maxatoms), STAT=istat)

    end subroutine allocate_arrays

!    subroutine allocate_2d_arrays(step2d,r2d1, r2d2, a2d)
    subroutine allocate_2d_arrays(maxatoms, step2d, xyz2d, r2d1, r2d2, a2d)
    ! Subroutine allocates the arrays for atom positions and dispacements.

    integer :: istat
    integer(ik) :: step2d,maxatoms
    real(rk), allocatable :: xyz2d(:,:,:,:)
    real(rk), allocatable :: r2d1(:,:)
    real(rk), allocatable :: r2d2(:,:)
    real(rk), allocatable :: a2d(:,:)

    allocate(xyz2d(1:maxatoms,1:3,1:step2d,1:step2d), STAT=istat)
    if(istat /= 0)write(*,*)"allocation error"
    allocate(r2d1(1:step2d,1:step2d), STAT=istat)
    if(istat /= 0)write(*,*)"allocation error"

    allocate(r2d2(1:step2d,1:step2d), STAT=istat)
    if(istat /= 0)write(*,*)"allocation error"

    allocate(a2d(1:step2d,1:step2d), STAT=istat)

    end subroutine allocate_2d_arrays

    subroutine deallocate_arrays(xyz0, dxyz, xyz, drtheta, r1, r2, a, aname)

    ! Subroutine deallocates arrays that were allocated in the allocate_arrays subroutine.

    integer :: istat
    real(rk), allocatable :: xyz0(:,:)
    real(rk), allocatable :: dxyz(:,:,:)
    real(rk), allocatable :: drtheta(:,:)
    real(rk), allocatable :: xyz(:,:,:,:)
    real(rk), allocatable :: r1(:,:)
    real(rk), allocatable :: r2(:,:)
    real(rk), allocatable :: a(:,:)
    character(len=3), allocatable :: aname(:)

    deallocate(xyz0, STAT=istat)
    deallocate(dxyz, STAT=istat)
    deallocate(drtheta, STAT=istat)
    deallocate(xyz, STAT=istat)
    deallocate(r1, STAT=istat)
    deallocate(r2, STAT=istat)
    deallocate(a, STAT=istat)      
    deallocate(aname, STAT=istat)
    ! TODO: Check istat !!

    end subroutine deallocate_arrays    

    subroutine deallocate_2d_arrays(xyz2d, r2d1, r2d2, a2d)
    !subroutine deallocate_2d_arrays(r2d1, r2d2, a2d)
    ! Subroutine deallocates arrays that were allocated in the allocate_arrays subroutine.

    integer :: istat
    real(rk), allocatable :: xyz2d(:,:,:,:)
    real(rk), allocatable :: r2d1(:,:)
    real(rk), allocatable :: r2d2(:,:)
    real(rk), allocatable :: a2d(:,:)

    deallocate(xyz2d, STAT=istat)
    if(istat /= 0)write(*,*)"deallocation error"
    deallocate(r2d1, STAT=istat)
    if(istat /= 0)write(*,*)"deallocation error"
    deallocate(r2d2, STAT=istat)
    if(istat /= 0)write(*,*)"deallocation error"
    deallocate(a2d, STAT=istat)   
    if(istat /= 0)write(*,*)"deallocation error"
    
    ! TODO: Check istat !!

    end subroutine deallocate_2d_arrays  
    subroutine init_geo(xyz0, maxatoms, aname)
    
    ! Subroutine reads in atom names and initial geometry from stdin in cartesian coordinates (Angstrom).
    integer(ik) :: maxatoms
    real(rk), allocatable, intent(inout) :: xyz0(:,:)
    integer :: i,k
    character(len=3), allocatable :: aname(:)

!    write(*,*) 'Reading input geometry (a.u)'

    read(*,*) (aname(i), i=1,maxatoms)
!    write(*,*) (aname(i), i=1,maxatoms)
    do i=1, maxatoms
        read(*,*) (xyz0(i,k), k=1,3)
        xyz0(i,:) = xyz0(i,:)/autoaa
!        write(*,*) (xyz0(i,k), k=1,3)
    end do

    
    end subroutine init_geo

    subroutine init_internal_geo(xyz0, r01, r02, a0)

    ! Converts the initial geometry to internal coordinates    

    real(rk), allocatable, intent(in) :: xyz0(:,:)
    real(rk), intent(inout) :: r01, r02, a0

    r01 = sqrt(sum((xyz0(1,:)-xyz0(2,:))**2))
    r02 = sqrt(sum((xyz0(3,:)-xyz0(2,:))**2))
    a0 = acos(sum((xyz0(3,:)-xyz0(2,:))*(xyz0(1,:)-xyz0(2,:)))/r01/r02)

 !   write(*,*) 'Input geo in internal coordinates'
 !   write(*,*) r01, r02, a0
    end subroutine init_internal_geo

    subroutine read_cartd(dxyz, maxatoms, maxmodes)

    ! Subroutine reads in displacement vectors from stdin in atomic units.
    integer(ik) :: maxatoms, maxmodes
    real(rk), allocatable, intent(inout) :: dxyz(:,:,:)
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
    real(rk), allocatable, intent(inout) :: drtheta(:,:)
    integer :: j, k

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

    subroutine read_2d(mode1, mode2, xmin2d, step2d, ssize2d)

    ! Reads which modes that will be combined to make (2d) displacement

    integer(ik), intent(inout) :: mode1, mode2, step2d
    real(rk), intent(inout) :: xmin2d, ssize2d

    read(*,*) mode1, mode2
    write(*,*) mode1, mode2
    read(*,*) xmin2d, step2d, ssize2d
    write(*,*) xmin2d, step2d, ssize2d

    end subroutine read_2d

!    subroutine angle_bending(ar1, ar2, amin, angstep, angsize)

    ! Subroutine reads in two distances and angle min and angle max for angle bending with fixed bond distances.

!    real(rk), intent(inout) :: ar1, ar2, amin, angstep, angsize

!    read(*,*) ar1, ar2, amin, angstep, angsize
!    write(*,*) ar1, ar2, amin, angstep, angsize

!    end subroutine angle_bending

!    subroutine bond_stretching(br1, angle, rmin, rstep, rsize)

    ! Subroutine reads in one distance and one angle and a min distance and max distance for bond stretching with fixed angle and fixed other bond.

!    real(rk), intent(inout) :: br1, angle, rmin, rstep, rsize

!    read(*,*) br1, angle, rmin, rstep, rsize
!    write(*,*) br1, angle, rmin, rstep, rsize

!    end subroutine bond_stretching
end module pes_init
