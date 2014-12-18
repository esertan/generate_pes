module pes_write

    use mydefs

    implicit none

contains

    subroutine write_grid(step, ssize, xmin)

    ! Writes the distortion grid to file in atomic units
    
    integer(ik) :: step
    integer :: l
    real(rk) :: s, ssize, xmin
    
    open(unit=iw1,form='formatted',status='unknown',file='Q')

    do l=1,step
        s=xmin+l*ssize
        write(iw1,'(i3XXXF6.3)') l,s
    end do
    close(iw1)
    
    end subroutine write_grid

    subroutine write_cartesian(xyz, maxatoms, maxmodes, step, aname)
    
    !Writes displaced cartesian coordinates to specified output file. Maximum number of output geometries are 999.

    integer(ik) :: step, maxatoms, maxmodes
    real(rk), allocatable, intent(in) :: xyz(:,:,:,:)
    character(len=3), allocatable, intent(in) :: aname(:)
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

            open(unit=iw2, form='formatted',status='unknown',file=cfile)
            open(unit=iw3, form='formatted',status='unknown',file=dfile)
            open(unit=iw4, form='formatted',status='unknown',file=efile)

            do i=1, maxatoms
                select case(i)
                    case(1)             !write H2 positions in cartesian
                    write(iw2,'(A2XXF9.6XXF9.6XXF9.6)') aname(i), xyz(1,i,1,imol)*autoaa, xyz(1,i,2,imol)*autoaa,&
                         xyz(1,i,3,imol)*autoaa
                    write(iw3,'(A2XXF9.6XXF9.6XXF9.6)') aname(i), xyz(2,i,1,imol)*autoaa, xyz(2,i,2,imol)*autoaa,&
                         xyz(2,i,3,imol)*autoaa
                    write(iw4,'(A2XXF9.6XXF9.6XXF9.6)') aname(i), xyz(3,i,1,imol)*autoaa, xyz(3,i,2,imol)*autoaa,&
                         xyz(3,i,3,imol)*autoaa
                case(2)                                                                                         !write O1 positions in cartesian
                    write(iw2,'(A2XXF9.6XXF9.6XXF9.6)') aname(i), xyz(1,i,1,imol)*autoaa, xyz(1,i,2,imol)*autoaa,&
                         xyz(1,i,3,imol)*autoaa
                    write(iw3,'(A2XXF9.6XXF9.6XXF9.6)') aname(i), xyz(2,i,1,imol)*autoaa, xyz(2,i,2,imol)*autoaa,&
                         xyz(2,i,3,imol)*autoaa
                    write(iw4,'(A2XXF9.6XXF9.6XXF9.6)') aname(i), xyz(3,i,1,imol)*autoaa, xyz(3,i,2,imol)*autoaa,&
                         xyz(3,i,3,imol)*autoaa
                case(3)                                                                                         !write H3 positions in cartesian
                    write(iw2,'(A2XXF9.6XXF9.6XXF9.6)') aname(i), xyz(1,i,1,imol)*autoaa, xyz(1,i,2,imol)*autoaa,&
                         xyz(1,i,3,imol)*autoaa
                    write(iw3,'(A2XXF9.6XXF9.6XXF9.6)') aname(i), xyz(2,i,1,imol)*autoaa, xyz(2,i,2,imol)*autoaa,&
                         xyz(2,i,3,imol)*autoaa
                    write(iw4,'(A2XXF9.6XXF9.6XXF9.6)') aname(i), xyz(3,i,1,imol)*autoaa, xyz(3,i,2,imol)*autoaa,&
                         xyz(3,i,3,imol)*autoaa
                end select
        end do                                                                                                  !write X positions in cartesian
            write(iw2,'(A1XXXF9.6XXF9.6XXF9.6)') "X", xyz(1,2,1,imol)*autoaa, xyz(1,2,2,imol)*autoaa,&
             (xyz(1,2,3,imol)+0.00001D0)*autoaa
            write(iw3,'(A1XXXF9.6XXF9.6XXF9.6)') "X", xyz(2,2,1,imol)*autoaa, xyz(2,2,2,imol)*autoaa,&
             (xyz(2,2,3,imol)+0.00001D0)*autoaa
            write(iw4,'(A1XXXF9.6XXF9.6XXF9.6)') "X", xyz(3,2,1,imol)*autoaa, xyz(3,2,2,imol)*autoaa,&
             (xyz(3,2,3,imol)+0.00001D0)*autoaa
    
        close(iw2)
        close(iw3)
        close(iw4)

        open(unit=iw2,form='formatted',status='unknown',file='TRAJBEND',position='append')
        open(unit=iw3,form='formatted',status='unknown',file='TRAJSYMM',position='append')
        open(unit=iw4,form='formatted',status='unknown',file='TRAJASYM',position='append')
        write(iw2,'(i1)') maxatoms
        write(iw2,'(i1)')
        write(iw2,'(A1XXF9.6XXF9.6XXF9.6)') aname(1), xyz(1,1,1,imol)*autoaa, xyz(1,1,2,imol)*autoaa,&
                         xyz(1,1,3,imol)*autoaa
        write(iw2,'(A1XXF9.6XXF9.6XXF9.6)') aname(2), xyz(1,2,1,imol)*autoaa, xyz(1,2,2,imol)*autoaa,&
                         xyz(1,2,3,imol)*autoaa
        write(iw2,'(A1XXF9.6XXF9.6XXF9.6)') aname(3), xyz(1,3,1,imol)*autoaa, xyz(1,3,2,imol)*autoaa,&
                         xyz(1,3,3,imol)*autoaa
        write(iw3,'(i1)') maxatoms
        write(iw3,*)
        write(iw3,'(A1XXF9.6XXF9.6XXF9.6)') aname(1), xyz(2,1,1,imol)*autoaa, xyz(2,1,2,imol)*autoaa,&
                         xyz(2,1,3,imol)*autoaa
        write(iw3,'(A1XXF9.6XXF9.6XXF9.6)') aname(2), xyz(2,2,1,imol)*autoaa, xyz(2,2,2,imol)*autoaa,&
                         xyz(2,2,3,imol)*autoaa
        write(iw3,'(A1XXF9.6XXF9.6XXF9.6)') aname(3), xyz(2,3,1,imol)*autoaa, xyz(2,3,2,imol)*autoaa,&
                         xyz(2,3,3,imol)*autoaa

        write(iw4,'(I1)') maxatoms
        write(iw4,*)
        write(iw4,'(A1XXF9.6XXF9.6XXF9.6)') aname(1), xyz(3,1,1,imol)*autoaa, xyz(3,1,2,imol)*autoaa,&
                         xyz(3,1,3,imol)*autoaa
        write(iw4,'(A1XXF9.6XXF9.6XXF9.6)') aname(2), xyz(3,2,1,imol)*autoaa, xyz(3,2,2,imol)*autoaa,&
                         xyz(3,2,3,imol)*autoaa
        write(iw4,'(A1XXF9.6XXF9.6XXF9.6)') aname(3), xyz(3,3,1,imol)*autoaa, xyz(3,3,2,imol)*autoaa,&
                         xyz(3,3,3,imol)*autoaa

        close(iw2)
        close(iw3)
        close(iw4)
    end do

    end subroutine write_cartesian

    subroutine write_internal(r1, r2, a, maxatoms, maxmodes, step, aname)

    ! Writes displaced internal coordinates to file. Output in cartesian format.
    
    integer(ik) :: maxatoms, maxmodes, step
    real(rk), allocatable, intent(in) :: r1(:,:), r2(:,:), a(:,:)
    character(len=3), allocatable, intent(in) :: aname(:)
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

        open(unit=iw2,form='formatted',status='unknown',file=cfile)
        open(unit=iw3,form='formatted',status='unknown',file=dfile)
        open(unit=iw4,form='formatted',status='unknown',file=efile)

        ! Write coordinates to file converted to cartesian format.

        write(iw2,'(A2XXF9.6XXF9.6XXF9.6)') aname(1), 0.d0, r1(1,imol)*dsin(0.5d0*a(1,imol))*autoaa, &
                -r1(1,imol)*dcos(0.5d0*a(1,imol))*autoaa
        write(iw3,'(A2XXF9.6XXF9.6XXF9.6)') aname(1), 0.d0, r1(2,imol)*dsin(0.5d0*a(2,imol))*autoaa, &
                -r1(2,imol)*dcos(0.5d0*a(2,imol))*autoaa
        write(iw4,'(A2XXF9.6XXF9.6XXF9.6)') aname(1), 0.d0, r1(3,imol)*dsin(0.5d0*a(3,imol))*autoaa, &
                -r1(3,imol)*dcos(0.5d0*a(3,imol))*autoaa

        write(iw2,'(A2XXF9.6XXF9.6XXF9.6)') aname(2), 0d0, 0d0, 0d0
        write(iw3,'(A2XXF9.6XXF9.6XXF9.6)') aname(2), 0d0, 0d0, 0d0
        write(iw4,'(A2XXF9.6XXF9.6XXF9.6)') aname(2), 0d0, 0d0, 0d0

        write(iw2,'(A2XXF9.6XXF9.6XXF9.6)') aname(3), 0.d0, -r2(1,imol)*dsin(0.5d0*a(1,imol))*autoaa,&
                    -r2(1,imol)*dcos(0.5d0*a(1,imol))*autoaa
        write(iw3,'(A2XXF9.6XXF9.6XXF9.6)') aname(3), 0.d0, -r2(2,imol)*dsin(0.5d0*a(2,imol))*autoaa,&
                    -r2(2,imol)*dcos(0.5d0*a(2,imol))*autoaa
        write(iw4,'(A2XXF9.6XXF9.6XXF9.6)') aname(3), 0.d0, -r2(3,imol)*dsin(0.5d0*a(3,imol))*autoaa,&
                    -r2(3,imol)*dcos(0.5d0*a(3,imol))*autoaa

        write(iw2,'(A1XXXF9.6XXF9.6XXF9.6)') "X", 0.d0, 0.d0, 0.00005d0
        write(iw3,'(A1XXXF9.6XXF9.6XXF9.6)') "X", 0.d0, 0.d0, 0.00005d0
        write(iw4,'(A1XXXF9.6XXF9.6XXF9.6)') "X", 0.d0, 0.d0, 0.00005d0

        close(iw2)
        close(iw3)
        close(iw4)

        ! Write coordinates to trajectory file (for append to work remove these files before running the code)

        open(unit=iw2,form='formatted',status='unknown',file='TRAJBENDINT',position='append')
        open(unit=iw3,form='formatted',status='unknown',file='TRAJSYMMINT',position='append')
        open(unit=iw4,form='formatted',status='unknown',file='TRAJASYMINT',position='append')

        write(iw2,'(i1)') maxatoms         
        write(iw2,'(i1)')
        write(iw2,'(A1XXF9.6XXF9.6XXF9.6)') aname(1), 0.d0, r1(1,imol)*dsin(0.5d0*a(1,imol))*autoaa, &
                -r1(1,imol)*dcos(0.5d0*a(1,imol))*autoaa
        write(iw2,'(A1XXF9.6XXF9.6XXF9.6)') aname(2), 0.d0, 0.d0, 0.d0
        write(iw2,'(A1XXF9.6XXF9.6XXF9.6)') aname(3), 0.d0, -r2(1,imol)*dsin(0.5d0*a(1,imol))*autoaa,&
                    -r2(1,imol)*dcos(0.5d0*a(1,imol))*autoaa

        write(iw3,'(i1)') maxatoms         
        write(iw3,'(i1)')
        write(iw3,'(A1XXF9.6XXF9.6XXF9.6)') aname(1), 0.d0, r1(2,imol)*dsin(0.5d0*a(2,imol))*autoaa, &
                -r1(2,imol)*dcos(0.5d0*a(2,imol))*autoaa
        write(iw3,'(A1XXF9.6XXF9.6XXF9.6)') aname(2), 0.d0, 0.d0, 0.d0
        write(iw3,'(A1XXF9.6XXF9.6XXF9.6)') aname(3), 0.d0, -r2(2,imol)*dsin(0.5d0*a(2,imol))*autoaa,&
                    -r2(2,imol)*dcos(0.5d0*a(2,imol))*autoaa

        write(iw4,'(i1)') maxatoms         
        write(iw4,'(i1)')
        write(iw4,'(A1XXF9.6XXF9.6XXF9.6)') aname(1), 0.d0, r1(3,imol)*dsin(0.5d0*a(3,imol))*autoaa, &
                -r1(3,imol)*dcos(0.5d0*a(3,imol))*autoaa
        write(iw4,'(A1XXF9.6XXF9.6XXF9.6)') aname(2), 0.d0, 0.d0, 0.d0
        write(iw4,'(A1XXF9.6XXF9.6XXF9.6)') aname(3), 0.d0, -r2(3,imol)*dsin(0.5d0*a(3,imol))*autoaa,&
                    -r2(3,imol)*dcos(0.5d0*a(3,imol))*autoaa

        close(iw2)
        close(iw3)
        close(iw4)
    end do 

    end subroutine write_internal

end module pes_write
