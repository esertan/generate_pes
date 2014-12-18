module pes_calc

    use mydefs

    implicit none

contains
  
    subroutine make_internal(r01, r02, a0, r1, r2, a, drtheta, maxatoms, maxmodes, step, xmin, ssize)

    ! Make the displacement in internal coordinates.
    
    integer(ik) :: maxatoms, maxmodes, step
    real(rk), intent(in) :: r01, r02, a0
    real(rk), allocatable, intent(in) :: drtheta(:,:)
    real(rk), allocatable, intent(inout) :: r1(:,:), r2(:,:), a(:,:)
    integer :: l, j
    real(rk) :: s, ssize, xmin

    write(*,*) 'Displaced internal coordinates'
    do l=1,step
        s=xmin+l*ssize
        do j=1,maxmodes
            r1(j,l) = r01+s*drtheta(j,1)
            r2(j,l) = r02+s*drtheta(j,2)
            a(j,l) = a0+s*drtheta(j,3)
!            write(*,*) drtheta(j,1), drtheta(j,2), drtheta(j,3)
!            write(*,*) s, r1(j,l)*autoaa, r2(j,l)*autoaa, a(j,l)*180.d0/pi
        end do
    end do

    end subroutine make_internal

    subroutine make_cartesian(xyz, xyz0, dxyz, maxatoms, maxmodes, step, xmin, ssize)

    ! Make displacement in cartesian coordinates

    integer(ik) :: maxatoms, maxmodes, step
    real(rk), allocatable, intent(in) :: xyz0(:,:)
    real(rk), allocatable, intent(in) :: dxyz(:,:,:)
    real(rk), allocatable, intent(inout) :: xyz(:,:,:,:)
    integer :: l, j, i, k
    real(rk) :: s, ssize, xmin
    
    do l = 1, step
        s=xmin+l*ssize
        do j=1,maxmodes
            do i=1, maxatoms
                do k=1,3
                    xyz(j,i,k,l)=xyz0(i,k)+s*dxyz(j,i,k)
                end do
            end do
        end do
    end do

    end subroutine make_cartesian

!    subroutine make_angle_bending(amin, angstep, angsize, acoord)

    ! Subroutine makes displacement along angle interval specified in input. Bond distances are held fixed.

!    real(rk), intent(in) :: amin, angstep, angsize
!    real(rk), allocatable, intent(inout) :: acoord
!    integer(ik) :: i, istat
!
!    allocate(acoord(1:angstep), STAT=istat)

!    do i=1,astep
!        acoord(i) = amin+i*angsize
!    end do

!    end subroutine make_angle_bending

!    subroutine make_bond_stretching(rmin, rstep, rsize, rcoord)

    ! Subroutine makes displacement along one bond with interval specified in input. Other bond and angle fixed.

!    real(rk), intent(in) :: rmin, rstep, rsize
!    real(rk), allocatable, intention(inout) :: rcoord(:)
!    integer(ik) :: i, istat

!    allocate(rcoord(1:rstep), STAT=istat)
!    do i=1,rstep
!        rcoord(i) = rmin+i*rsize
!    end do

!    end subroutine make_bond_stretching
    
!    subroutine write_angle_bending(acoord)

!    real(rk), allocatable, intent(in) :: acoord(:)
!    integer(ik) :: istat

!    deallocate(acoord, STAT=istat)

!    end subroutine write_angle_bending

!    subroutine write_bond_stretching(rcoord)

!    real(rk), allocatable, intent(in) :: rcoord(:)
!    integer(ik) :: istat

!    deallocate(rcoord, STAT=istat)

!    end subroutine write_bond_stretching

end module pes_calc
