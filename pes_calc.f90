module pes_calc

    use mydefs

    implicit none

contains
  
    subroutine make_internal(r01, r02, a0, r1, r2, a, drtheta, maxmodes, step, xmin, ssize)

    ! Make the displacement in internal coordinates.
    
    integer(ik) :: maxmodes, step
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

    subroutine make_2d_internal(r01, r02, a0, r2d1, r2d2, a2d, drtheta, step2d, xmin2d, ssize2d, mode1, mode2)

    ! Make the displacement in internal coordinates.
    
    integer(ik), intent(in) :: step2d, mode1, mode2
    real(rk), intent(in) :: r01, r02, a0, ssize2d, xmin2d
    real(rk), allocatable, intent(in) :: drtheta(:,:)
    real(rk), allocatable, intent(inout) :: r2d1(:,:), r2d2(:,:), a2d(:,:)
    integer :: l, m!, counter=0
    real(rk) :: s, t

    write(*,*) 'Displaced internal coordinates'
    do m=1,step2d
        t=xmin2d+m*ssize2d
        write(*,*) "t=", t
        do l=1,step2d
            s=xmin2d+l*ssize2d
            write(*,*) "s=", s
!            counter=counter+1
            r2d1(l,m) = r01+t*drtheta(mode1,1)+s*drtheta(mode2,1)
            r2d2(l,m) = r02+t*drtheta(mode1,2)+s*drtheta(mode2,2)
            a2d(l,m) = a0+t*drtheta(mode1,3)+s*drtheta(mode2,3)
            write(*,*) 0.d0, r2d1(l,m)*sin(0.5_rk*a2d(l,m))*autoaa, -r2d1(l,m)*cos(0.5_rk*a2d(l,m))*autoaa
        end do
    end do
!    write(*,*) counter
    end subroutine make_2d_internal

    subroutine make_2d_cartesian(xyz2d, xyz0, dxyz, maxatoms, step2d, xmin2d, ssize2d, mode1, mode2)

   ! Make displacement in cartesian coordinates using two the displacement vectors of normal modes (~2d)

    integer(ik), intent(in) :: maxatoms, step2d, mode1, mode2
    real(rk), allocatable, intent(in) :: xyz0(:,:)
    real(rk), allocatable, intent(in) :: dxyz(:,:,:)
    real(rk), allocatable, intent(inout) :: xyz2d(:,:,:,:)
    integer :: l, i, k, m
    real(rk) :: s,t, ssize2d, xmin2d
    
    do m=1,step2d
        t=xmin2d+m*ssize2d
        write(*,*) "t=", t
        do l = 1, step2d
            s=xmin2d+l*ssize2d
            write(*,*) "s=", s
            do i=1, maxatoms
                do k=1,3
                    xyz2d(i,k,l,m)=xyz0(i,k)+t*dxyz(mode1,i,k)+s*dxyz(mode2,i,k)
                end do
            end do
            write(*,*) ((xyz2d(1,k, l, m)-xyz2d(2,k, l, m))*autoaa, k=1,3)

        end do
    end do

    end subroutine make_2d_cartesian




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
