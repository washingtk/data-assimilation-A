!===========================================
!リアプノフ指数の推定プログラム
!===========================================
program main
    use global_variable, only: NN, t_1year, t_1day, dt
    use model_functions
    use other_subroutines
    implicit none
    integer, parameter :: sampling_num = 10000
    integer, parameter :: t_dataset = 3*t_1year
    integer, parameter :: t_lyapunov = 14*t_1day
    real(8), parameter :: sigma = 5.0d0
    real(8) :: XX(1:NN)
    real(8) :: XX_NEX(1:NN)
    real(8) :: XX_ATRC(1:t_dataset, 1:NN)
    real(8) :: ERROR_GROW(1:sampling_num,1:t_lyapunov)
    real(8) :: ERROR(1:NN)
    real(8) :: random, error_norm, error_norm_next
    integer :: t_idx, n_idx
    integer :: sampling_idx, data_idx, seedsize
    integer, allocatable :: seed(:)
    character(500) :: dirname, filename, hfile
    character(10) :: header(1:NN+1)

    !
    !check
    !
    write(*,*) "--------------------------------------"
    write(*,"(a30,i9)") 'nubmber of variable: ', NN
    write(*,"(a30,i9)") 'time step for 1 year: ', t_1year
    write(*,"(a30,i9)") 'time step for 1 day: ', t_1day
    write(*,"(a30,f9.4)") 'dt: ', dt
    write(*,"(a30,f9.4)") 'F: ', FF
    write(*,"(a30,f9.4)") 'sigma, sampling: ', sigma
    write(*,"(a30,i9)") 'sampling number: ', sampling_num
    write(*,*) "--------------------------------------"

    !
    !make directory
    !
    write(dirname, '("output_lyapnov_data/")')
    call makedirs(dirname)

    !
    !make infodata file
    !
    write(filename, '("result_N", i2.2, "FF", i2.2, f0.3, "std", i2.2, f0.2, ".txt")') &
        & NN, int(FF), FF-int(FF), int(sigma), sigma-int(sigma)
    hfile    = trim(dirname)//"header_"//trim(filename)
    filename = trim(dirname)//trim(filename)
    open(200, file=hfile, status="replace")
    write(200, '(a30,a)')     'datafile ', trim(filename)
    write(200, '(a30,i8)')    'NN ', NN
    write(200, '(a30,f8.5)')  'F ', FF
    write(200, '(a30,f8.5)')  'dt ', dt
    write(200, '(a30,i9)') 'time_step_for_1_year: ', t_1year
    write(200, '(a30,i9)') 'time_step_for_1_day: ', t_1day
    write(200, '(a30,f9.4)') 'sigma_sampling: ', sigma
    write(200, '(a30,i9)') 'sampling_number: ', sampling_num
    close(200)

    !
    !initialize values
    !
    XX(:) = FF
    XX(20) = FF + 0.008d0

    !
    !initial run for 1 year .when. F=8
    !
    do t_idx = 1, t_1year
        call runge_kutta_method(NN, XX(1:NN), XX_NEX(1:NN))
        XX(:) = XX_NEX(:)
    end do

    !
    !attractor for 3 years
    !
    XX_ATRC(1,1:NN) = XX(1:NN)
    do t_idx = 2, t_dataset
        call runge_kutta_method(NN, XX_ATRC(t_idx-1,1:NN), XX_ATRC(t_idx,1:NN))
    end do

    !
    !attractor for 3 years database
    !
    open(300, file="lyapunov.dat", status="replace")
    do t_idx = 1, t_dataset
        write(300, *) t_idx/t_1day, XX_ATRC(t_idx,:)
    end do
    close(300)

    !
    !一様乱数
    !
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    call random_seed(get=seed(:))
    seed = 1
    call random_seed(put=seed(:))

    !
    !lyapunov sampling
    !
    do sampling_idx = 1, sampling_num
        call random_number(random)
        data_idx = floor((t_dataset-t_lyapunov)*random)
        !
        !誤差成長率の調査
        !
        call box_muller(NN, ERROR(1:NN), 0d0, sigma, sampling_idx)
        XX(1:NN) = XX_ATRC(data_idx,1:NN) + ERROR(1:NN)
        error_norm = norm2(ERROR(1:NN))
        !
        !t_lyapunov 分の時間ステップにおける誤差発達率
        !
        do t_idx = 1, t_lyapunov
            call runge_kutta_method(NN, XX(1:NN), XX_NEX(1:NN))
            XX(:) = XX_NEX(:)
            ERROR(:) = XX(:) - XX_ATRC(data_idx+t_idx,:)
            error_norm_next = norm2(ERROR(:))
            ERROR_GROW(sampling_idx,t_idx) = error_norm_next/error_norm
        end do
    end do

    !
    !make data file
    !
    write(*,*) 'save data in >> ', trim(filename)
    write(*,*) 'info data in >> ', trim(hfile)
    open(100, file=filename, status="replace")
    do t_idx = 1, 7*t_1day
        write(100,*) 1d0*t_idx/t_1day, sum(ERROR_GROW(:,t_idx))/sampling_num
    end do
    close(100)

end program
