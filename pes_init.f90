module pes_init

    use mydefs

    implicit none

contains
  
    subroutine init_geo(xyz0)
    
    ! Subroutine reads in initial geometry from stdin in cartesian coordinates (Angstrom).

    real(rk), dimension(1:maxatoms, 1:3), intent(inout) :: xyz0
    integer :: i,j,k

!    write(*,*) 'Reading input geometry (a.u)'
    do i=1, maxatoms
        read(*,*) (xyz0(i,k), k=1,3)
        xyz0(i,:) = xyz0(i,:)/autoaa
!        write(*,*) (xyz0(i,k), k=1,3)
    end do

    
    end subroutine init_geo

    subroutine init_internal_geo(xyz0, r01, r02, a0)

    real(rk), dimension(1:maxatoms, 1:3), intent(in) :: xyz0
    real(rk), intent(inout) :: r01, r02, a0

    r01 = dsqrt(sum((xyz0(1,:)-xyz0(2,:))**2))
    r02 = dsqrt(sum((xyz0(3,:)-xyz0(2,:))**2))
    a0 = dacos(sum((xyz0(3,:)-xyz0(2,:))*(xyz0(1,:)-xyz0(2,:)))/r01/r02)

 !   write(*,*) 'Input geo in internal coordinates'
 !   write(*,*) r01, r02, a0
    end subroutine init_internal_geo

    subroutine read_cartd(dxyz)

    ! Subroutine reads in displacement vectors from stdin in atomic units.
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

    subroutine read_interd(drtheta)

    ! Subroutine reads in displacement in internal coordinates (atomic units).

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
