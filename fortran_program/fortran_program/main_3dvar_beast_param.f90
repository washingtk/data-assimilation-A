!===========================================
!Kalman Filter の実装
!===========================================
program main
    use global_variable, only: NN, NM, t_1year
    use model_functions
    use other_subroutines
    use data_assimilation
    implicit none
    integer :: w_idx, t_idx, n_idx, sigma_idx, i, j
    integer, parameter :: window = 5
    real(8) :: XX(1:NN)
    real(8) :: XX_NEX(1:NN)
    real(8) :: XF(1:NN)
    real(8) :: XA(1:NN)
    real(8) :: YY(1:NM)
    real(8) :: ER(1:NN)
    !---------------
    !KF parameter
    !---------------
    !MM:    接線型行列
    !BB:    状態変数の誤差分散行列
    !HH:    観測行列
    !RR:    観測誤差行列
    !KK:    ゲイン
    !DD:    D値
    !sigma: 観測誤差
    !---------------
    real(8) :: MM(1:NN,1:NN)
    real(8) :: BB(1:NN,1:NN)
    real(8) :: HH(1:NM,1:NN)
    real(8) :: RR(1:NM,1:NM)
    real(8) :: KK(1:NN,1:NM)
    real(8) :: DD(1:NN)
    real(8) :: rmse, trace_pf, trace_pa
    real(8) :: sigma_b
    real(8), parameter :: sigma = 1d0
    character(500) :: dirname, filename

    !
    !観測モデル行列H
    !
    if (mod(NN,NM) .ne. 0) then
        print*, "Error! NN/NM is not 0"
        return
    else
        j = int(NN/NM)
        HH(:,:) = 0d0
        write(*,*) "--------------------------------------"
        write(*,*) "matrix H"
        do i = 1, NM
            HH(i,1+j*(i-1)) = 1d0
            write(*,"(40i3)"), int(HH(i,:))
        end do
    end if


    do sigma_idx = 1, 40
        !
        !set alpha
        !
        sigma_b = 0.1 + 0.1*(sigma_idx-1)

        !
        !行列の初期化
        !
        BB(:,:) = 0d0
        do i = 1, NN
            BB(i,i) = sigma_b * sigma_b
        end do
        RR(:,:) = 0d0
        do i = 1, NM
            RR(i,i) = sigma
        end do

        !
        !check
        !
        write(*,*) "--------------------------------------"
        write(*,"(a30,i9)") 'nubmber of variable: ', NN
        write(*,"(a30,i9)") 'nubmber of obervation: ', NM
        write(*,"(a30,i9)") 'NN/NM: ', int(NN/NM)
        write(*,"(a30,i9)") 'window: ', window
        write(*,"(a30,f9.4)") 'dt: ', dt
        write(*,"(a30,f9.4)") 'F: ', FF
        write(*,"(a30,f9.4)") 'sigma, sampling: ', sigma
        write(*,"(a30,f9.4)") 'sigma, matrix B:', sigma_b
        write(*,*) "--------------------------------------"

        !
        !make directory
        !
        write(dirname, '("output_3dvar_data/")')
        call makedirs(dirname)
        write(filename, '("3dvar_out", "window", i3.3, "NM", i2.2, "sigma_b", i2.2, f0.3, ".txt")') &
            & window, NM, int(sigma_b), sigma_b-int(sigma_b)
        filename = trim(dirname)//trim(filename)
        write(*,*) 'save data in >> ', trim(filename)
        open(200, file=filename, status="replace")

        !
        !initialize values 1
        !
        XX(:) = FF
        XX(20) = FF + 0.008d0
        XA(:) = FF
        XA(20) = FF + 0.004d0
        !
        !time step in model
        !
        do t_idx = 1, t_1year
            call runge_kutta_method(NN, XX, XX_NEX)
            XX = XX_NEX
            call runge_kutta_method(NN, XA, XX_NEX)
            XA = XX_NEX
        end do

        !
        !header
        !
        ! write(*,'(8a10)') "sycle", "obs", "forcast", "analysis", "true", "rmse", "trace_pf", "trace_pa"

        !
        !time step in model
        !
        do w_idx = 1, 1000
            XF = XA
            do t_idx = 1, window
                !真値
                call runge_kutta_method(NN, XX, XX_NEX)
                XX = XX_NEX
                !観測値
                call runge_kutta_method(NN, XF, XX_NEX)
                XF = XX_NEX
            end do
            !
            !観測値
            !
            call box_muller(NN, ER(1:NN), 0d0, sigma, w_idx)
            YY(1:NM) = matmul( HH(1:NM,1:NN), XX(1:NN)+ER(1:NN) )
            !
            !ゲイン
            !
            call KF_matrix_get_K(NN, NM, KK(1:NN,1:NM), BB(1:NN,1:NN), &
                & HH(1:NM,1:NN), RR(1:NM,1:NM))
            !
            !解析値の計算
            !
            call KF_matrix_get_XA(NN, NM, XA(1:NN), XF(1:NN), YY(1:NM), &
                & KK(1:NN,1:NM), HH(1:NM,1:NN))
            !
            !RMSE
            !
            call get_rmse(NN, rmse, XA(1:NN), XX(1:NN))
            write(200,*) w_idx, YY(1), XF(1), XA(1), XX(1), rmse
            ! write(*,'(i10, 7f10.5)') w_idx, YY(1), XF(1), XA(1), XX(1), rmse, sqrt(trace_pf/dble(NN)), sqrt(trace_pa/dble(NN))
        end do

        close(200)

    end do
end program
