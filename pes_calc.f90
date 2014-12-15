module pes_calc

    use mydefs

    implicit none

contains
  
    subroutine make_internal(r01, r02, a0, r1, r2, a, drtheta, maxatoms, maxmodes, step, xmin, ssize)

    ! Make the displacement in internal coordinates.
    
    integer(ik) :: maxatoms, maxmodes, step
    real(rk), intent(in) :: r01, r02, a0
    real(rk), dimension(1:maxmodes, 1:3), intent(in) :: drtheta
    real(rk), dimension(1:maxmodes, 1:step), intent(inout) :: r1, r2, a
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
    real(rk), dimension(1:maxatoms, 1:3), intent(in) :: xyz0
    real(rk), dimension(1:maxmodes, 1:maxatoms, 1:3) :: dxyz
    real(rk), dimension(1:maxmodes, 1:maxatoms, 1:3, 1:step) :: xyz
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
    
    subroutine write_grid(step, ssize, xmin)

    ! Writes the distortion grid to file in atomic units
    
    integer(ik) :: step
    integer :: l
    real(rk) :: s, ssize, xmin
    
    open(unit=15,form='formatted',status='unknown',file='Q')

    do l=1,step
        s=xmin+l*ssize
        write(15,'(i3XXXF6.3)') l,s
    end do
    close(15)
    
    end subroutine write_grid
    subroutine write_cartesian(xyz, maxatoms, maxmodes, step)
    
    !Writes displaced cartesian coordinates to specified output file. Maximum number of output geometries are 999.

    integer(ik) :: step, maxatoms, maxmodes
    real(rk), dimension(1:maxmodes, 1:maxatoms, 1:3, 1:step), intent(in) :: xyz
    character(len=10) :: cfile, dfile, efile
    integer :: i, imol, l
    real(rk) :: s

    do imol=1,step
!        write(*,*) imol
            select case(imol)
                case(1:9)
                    write(cfile, '(a,i1)') 'H2OBEND00', imol
                    write(dfile, '(a,i1)') 'H2OSYMM00', imol
                    write(efile, '(a,i1)') 'H2OASYM00', imol
                case(10:99)
                    write(cfile, '(a,i2)') 'H2OBEND0', imol
                    write(dfile, '(a,i2)') 'H2OSYMM0', imol
                    write(efile, '(a,i2)') 'H2OASYM0', imol
                case(100:999)
                    write(cfile, '(a,i3)') 'H2OBEND', imol
                    write(dfile, '(a,i3)') 'H2OSYMM', imol
                    write(efile, '(a,i3)') 'H2OASYM', imol
            end select
            !write(*,*) 'Generating', cfile, dfile, efile

            open(unit=16, form='formatted',status='unknown',file=cfile)
            open(unit=17, form='formatted',status='unknown',file=dfile)
            open(unit=18, form='formatted',status='unknown',file=efile)

            do i=1, maxatoms
                select case(i)
                    case(1)             !write H2 positions in cartesian
                    write(16,'(A2XXF9.6XXF9.6XXF9.6)') "H2", xyz(1,i,1,imol)*autoaa, xyz(1,i,2,imol)*autoaa,&
                         xyz(1,i,3,imol)*autoaa
                    write(17,'(A2XXF9.6XXF9.6XXF9.6)') "H2", xyz(2,i,1,imol)*autoaa, xyz(2,i,2,imol)*autoaa,&
                         xyz(2,i,3,imol)*autoaa
                    write(18,'(A2XXF9.6XXF9.6XXF9.6)') "H2", xyz(3,i,1,imol)*autoaa, xyz(3,i,2,imol)*autoaa,&
                         xyz(3,i,3,imol)*autoaa
                case(2)                                                                                         !write O1 positions in cartesian
                    write(16,'(A2XXF9.6XXF9.6XXF9.6)') "O1", xyz(1,i,1,imol)*autoaa, xyz(1,i,2,imol)*autoaa,&
                         xyz(1,i,3,imol)*autoaa
                    write(17,'(A2XXF9.6XXF9.6XXF9.6)') "O1", xyz(2,i,1,imol)*autoaa, xyz(2,i,2,imol)*autoaa,&
                         xyz(2,i,3,imol)*autoaa
                    write(18,'(A2XXF9.6XXF9.6XXF9.6)') "O1", xyz(3,i,1,imol)*autoaa, xyz(3,i,2,imol)*autoaa,&
                         xyz(3,i,3,imol)*autoaa
                case(3)                                                                                         !write H3 positions in cartesian
                    write(16,'(A2XXF9.6XXF9.6XXF9.6)') "H3", xyz(1,i,1,imol)*autoaa, xyz(1,i,2,imol)*autoaa,&
                         xyz(1,i,3,imol)*autoaa
                    write(17,'(A2XXF9.6XXF9.6XXF9.6)') "H3", xyz(2,i,1,imol)*autoaa, xyz(2,i,2,imol)*autoaa,&
                         xyz(2,i,3,imol)*autoaa
                    write(18,'(A2XXF9.6XXF9.6XXF9.6)') "H3", xyz(3,i,1,imol)*autoaa, xyz(3,i,2,imol)*autoaa,&
                         xyz(3,i,3,imol)*autoaa
                end select
        end do                                                                                                  !write X positions in cartesian
            write(16,'(A1XXXF9.6XXF9.6XXF9.6)') "X", xyz(1,2,1,imol)*autoaa, xyz(1,2,2,imol)*autoaa,&
             (xyz(1,2,3,imol)+0.00001D0)*autoaa
            write(17,'(A1XXXF9.6XXF9.6XXF9.6)') "X", xyz(2,2,1,imol)*autoaa, xyz(2,2,2,imol)*autoaa,&
             (xyz(2,2,3,imol)+0.00001D0)*autoaa
            write(18,'(A1XXXF9.6XXF9.6XXF9.6)') "X", xyz(3,2,1,imol)*autoaa, xyz(3,2,2,imol)*autoaa,&
             (xyz(3,2,3,imol)+0.00001D0)*autoaa
    
        close(16)
        close(17)
        close(18)

        open(unit=19,form='formatted',status='unknown',file='TRAJBEND',position='append')
        open(unit=20,form='formatted',status='unknown',file='TRAJSYMM',position='append')
        open(unit=21,form='formatted',status='unknown',file='TRAJASYM',position='append')
        write(19,'(i1)') maxatoms
        write(19,'(i1)')
        write(19,'(AXXF9.6XXF9.6XXF9.6)') "H", xyz(1,1,1,imol)*autoaa, xyz(1,1,2,imol)*autoaa,&
                         xyz(1,1,3,imol)*autoaa
        write(19,'(AXXF9.6XXF9.6XXF9.6)') "O", xyz(1,2,1,imol)*autoaa, xyz(1,2,2,imol)*autoaa,&
                         xyz(1,2,3,imol)*autoaa
        write(19,'(AXXF9.6XXF9.6XXF9.6)') "H", xyz(1,3,1,imol)*autoaa, xyz(1,3,2,imol)*autoaa,&
                         xyz(1,3,3,imol)*autoaa
        write(20,'(i1)') maxatoms
        write(20,*)
        write(20,'(AXXF9.6XXF9.6XXF9.6)') "H", xyz(2,1,1,imol)*autoaa, xyz(2,1,2,imol)*autoaa,&
                         xyz(2,1,3,imol)*autoaa
        write(20,'(AXXF9.6XXF9.6XXF9.6)') "O", xyz(2,2,1,imol)*autoaa, xyz(2,2,2,imol)*autoaa,&
                         xyz(2,2,3,imol)*autoaa
        write(20,'(AXXF9.6XXF9.6XXF9.6)') "H", xyz(2,3,1,imol)*autoaa, xyz(2,3,2,imol)*autoaa,&
                         xyz(2,3,3,imol)*autoaa

        write(21,'(I1)') maxatoms
        write(21,*)
        write(21,'(AXXF9.6XXF9.6XXF9.6)') "H", xyz(3,1,1,imol)*autoaa, xyz(3,1,2,imol)*autoaa,&
                         xyz(3,1,3,imol)*autoaa
        write(21,'(AXXF9.6XXF9.6XXF9.6)') "O", xyz(3,2,1,imol)*autoaa, xyz(3,2,2,imol)*autoaa,&
                         xyz(3,2,3,imol)*autoaa
        write(21,'(AXXF9.6XXF9.6XXF9.6)') "H", xyz(3,3,1,imol)*autoaa, xyz(3,3,2,imol)*autoaa,&
                         xyz(3,3,3,imol)*autoaa

        close(19)
        close(20)
        close(21)
    end do

    end subroutine write_cartesian

    subroutine write_internal(r1, r2, a, maxatoms, maxmodes, step)

    ! Writes displaced internal coordinates to file. Output in cartesian format.
    
    integer(ik) :: maxatoms, maxmodes, step
    real(rk), dimension(1:maxmodes, 1:step), intent(in) :: r1, r2, a
    character(len=13) :: cfile, dfile, efile
    integer :: i, imol, l
    real(rk) :: s

    do imol=1,step
!        write(*,*) imol
        select case(imol)
            case(1:9)
                write(cfile, '(a,i1)') 'H2OBENDINT00', imol
                write(dfile, '(a,i1)') 'H2OSYMMINT00', imol
                write(efile, '(a,i1)') 'H2OASYMINT00', imol
            case(10:99)
                write(cfile, '(a,i2)') 'H2OBENDINT0', imol
                write(dfile, '(a,i2)') 'H2OSYMMINT0', imol
                write(efile, '(a,i2)') 'H2OASYMINT0', imol
            case(100:999)
                write(cfile, '(a,i3)') 'H2OBENDINT', imol
                write(dfile, '(a,i3)') 'H2OSYMMINT', imol
                write(efile, '(a,i3)') 'H2OASYMINT', imol
        end select
        !write(*,*) 'Generating files', cfile, dfile, efile

        open(unit=16,form='formatted',status='unknown',file=cfile)
        open(unit=17,form='formatted',status='unknown',file=dfile)
        open(unit=18,form='formatted',status='unknown',file=efile)

        ! Write coordinates to file converted to cartesian format.

        write(16,'(A2XXF9.6XXF9.6XXF9.6)') "H2", 0.d0, r1(1,imol)*dsin(0.5d0*a(1,imol))*autoaa, &
                -r1(1,imol)*dcos(0.5d0*a(1,imol))*autoaa
        write(17,'(A2XXF9.6XXF9.6XXF9.6)') "H2", 0.d0, r1(2,imol)*dsin(0.5d0*a(2,imol))*autoaa, &
                -r1(2,imol)*dcos(0.5d0*a(2,imol))*autoaa
        write(18,'(A2XXF9.6XXF9.6XXF9.6)') "H2", 0.d0, r1(3,imol)*dsin(0.5d0*a(3,imol))*autoaa, &
                -r1(3,imol)*dcos(0.5d0*a(3,imol))*autoaa

        write(16,'(A2XXF9.6XXF9.6XXF9.6)') "O1", 0d0, 0d0, 0d0
        write(17,'(A2XXF9.6XXF9.6XXF9.6)') "O1", 0d0, 0d0, 0d0
        write(18,'(A2XXF9.6XXF9.6XXF9.6)') "O1", 0d0, 0d0, 0d0

        write(16,'(A2XXF9.6XXF9.6XXF9.6)') "H3", 0.d0, -r2(1,imol)*dsin(0.5d0*a(1,imol))*autoaa,&
                    -r2(1,imol)*dcos(0.5d0*a(1,imol))*autoaa
        write(17,'(A2XXF9.6XXF9.6XXF9.6)') "H3", 0.d0, -r2(2,imol)*dsin(0.5d0*a(2,imol))*autoaa,&
                    -r2(2,imol)*dcos(0.5d0*a(2,imol))*autoaa
        write(18,'(A2XXF9.6XXF9.6XXF9.6)') "H3", 0.d0, -r2(3,imol)*dsin(0.5d0*a(3,imol))*autoaa,&
                    -r2(3,imol)*dcos(0.5d0*a(3,imol))*autoaa

        write(16,'(A1XXXF9.6XXF9.6XXF9.6)') "X", 0.d0, 0.d0, 0.00005d0
        write(17,'(A1XXXF9.6XXF9.6XXF9.6)') "X", 0.d0, 0.d0, 0.00005d0
        write(18,'(A1XXXF9.6XXF9.6XXF9.6)') "X", 0.d0, 0.d0, 0.00005d0

        close(16)
        close(17)
        close(18)

        ! Write coordinates to trajectory file (for append to work remove these files before running the code)

        open(unit=19,form='formatted',status='unknown',file='TRAJBENDINT',position='append')
        open(unit=20,form='formatted',status='unknown',file='TRAJSYMMINT',position='append')
        open(unit=21,form='formatted',status='unknown',file='TRAJASYMINT',position='append')

        write(19,'(i1)') maxatoms         
        write(19,'(i1)')
        write(19,'(AXXF9.6XXF9.6XXF9.6)') "H", 0.d0, r1(1,imol)*dsin(0.5d0*a(1,imol))*autoaa, &
                -r1(1,imol)*dcos(0.5d0*a(1,imol))*autoaa
        write(19,'(AXXF9.6XXF9.6XXF9.6)') "O", 0.d0, 0.d0, 0.d0
        write(19,'(AXXF9.6XXF9.6XXF9.6)') "H", 0.d0, -r2(1,imol)*dsin(0.5d0*a(1,imol))*autoaa,&
                    -r2(1,imol)*dcos(0.5d0*a(1,imol))*autoaa

        write(20,'(i1)') maxatoms         
        write(20,'(i1)')
        write(20,'(AXXF9.6XXF9.6XXF9.6)') "H", 0.d0, r1(2,imol)*dsin(0.5d0*a(2,imol))*autoaa, &
                -r1(2,imol)*dcos(0.5d0*a(2,imol))*autoaa
        write(20,'(AXXF9.6XXF9.6XXF9.6)') "O", 0.d0, 0.d0, 0.d0
        write(20,'(AXXF9.6XXF9.6XXF9.6)') "H", 0.d0, -r2(2,imol)*dsin(0.5d0*a(2,imol))*autoaa,&
                    -r2(2,imol)*dcos(0.5d0*a(2,imol))*autoaa

        write(21,'(i1)') maxatoms         
        write(21,'(i1)')
        write(21,'(AXXF9.6XXF9.6XXF9.6)') "H", 0.d0, r1(3,imol)*dsin(0.5d0*a(3,imol))*autoaa, &
                -r1(3,imol)*dcos(0.5d0*a(3,imol))*autoaa
        write(21,'(AXXF9.6XXF9.6XXF9.6)') "O", 0.d0, 0.d0, 0.d0
        write(21,'(AXXF9.6XXF9.6XXF9.6)') "H", 0.d0, -r2(3,imol)*dsin(0.5d0*a(3,imol))*autoaa,&
                    -r2(3,imol)*dcos(0.5d0*a(3,imol))*autoaa

        close(19)
        close(20)
        close(21)
    end do 

    end subroutine write_internal

end module pes_calc
