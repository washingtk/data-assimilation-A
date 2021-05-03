program main
    use global_variable
    use model_functions
    use other_subroutines
    implicit none
    integer :: t_idx, n_idx
    real(8) :: XX(1:NN)
    real(8) :: XX_NEX(1:NN)
    character(500) :: dirname, filename, hfile
    character(10) :: header(1:NN+1)

    !
    !check
    !
    write(*,"(a20,i9)") 'nubmber of variable: ', NN
    write(*,"(a20,i9)") 'time max steps: ', t_max
    write(*,"(a20,f9.4)") 'dt:', dt
    write(*,"(a20,f9.4)") 'F:', FF

    !
    !make directory
    !
    write(dirname, '("output_model_data/")')
    call makedirs(dirname)

    !
    !make infodata file
    !
    write(filename, '("result_N", i2.2, "t_max", i8.8, "FF", i2.2, f0.3  ".txt")') &
        & NN, t_max, int(FF), FF-int(FF)
    hfile    = trim(dirname)//"header_"//trim(filename)
    filename = trim(dirname)//trim(filename)
    open(200, file=hfile, status="replace")
    write(200, '(a20,a)')     'datafile ', trim(filename)
    write(200, '(a20,i8)')    'NN ', NN
    write(200, '(a20,f8.5)')  'F ', FF
    write(200, '(a20,f8.5)')  'dt ', dt
    write(200, '(a20,i8)')    'max_time_step ', t_max
    close(200)

    !
    !make data file
    !
    write(*,*) 'save data in >> ', trim(filename)
    write(*,*) 'info data in >> ', trim(hfile)
    open(100, file=filename, status="replace")


    !
    !header in data file
    !
    write(header(1), '("time")')
    do n_idx = 1, NN
        write(header(n_idx+1), '("X", i2.2)') n_idx
    end do
    write(100,*) header

    !
    !initialize values 1
    !
    ! do n_idx = 1, NN
    !     XX(n_idx) = sin(1d0*n_idx/NN*2d0*pi)
    ! end do

    !
    !initialize values 2
    !
    XX(:) = FF
    XX(20) = FF + 0.008d0

    !
    !time step in model
    !
    write(100,*) 0d0, XX
    do t_idx = 1, t_max
        call runge_kutta_method(NN, XX, XX_NEX)
        XX = XX_NEX
        write(100,*) dt*t_idx, XX
    end do

    close(100)

end program main
