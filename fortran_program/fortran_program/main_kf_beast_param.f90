!===========================================
!Kalman Filter の実装
!===========================================
program main
    use global_variable, only: NN, NM, t_1year
    use model_functions
    use other_subroutines
    use data_assimilation
    implicit none
    integer :: w_idx, t_idx, n_idx, alpha_idx, i, j
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
    !PF:    予報誤差行列
    !PA:    解析誤差行列
    !HH:    観測行列
    !RR:    観測誤差行列
    !KK:    カルマンゲイン
    !sigma: 観測誤差
    !---------------
    real(8) :: MM(1:NN,1:NN)
    real(8) :: PF(1:NN,1:NN)
    real(8) :: PA(1:NN,1:NN)
    real(8) :: HH(1:NM,1:NN)
    real(8) :: RR(1:NM,1:NM)
    real(8) :: KK(1:NN,1:NM)
    real(8) :: rmse, trace_pf, trace_pa
    real(8) :: alpha
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


    do alpha_idx = 1, 20
        !
        !set alpha
        !
        alpha = 1d0 + 0.01*(alpha_idx-1)

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
        write(*,"(a30,f9.4)") 'inflation alpha:', alpha
        write(*,*) "--------------------------------------"

        !
        !make directory
        !
        write(dirname, '("output_kf_data/")')
        call makedirs(dirname)
        write(filename, '("kf_out", "window", i3.3, "NM", i2.2, "alpha", i2.2, f0.3, ".txt")') &
            & window, NM, int(alpha), alpha-int(alpha)
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
        !行列の初期化
        !
        RR(:,:) = 0d0
        PF(:,:) = 0d0
        PA(:,:) = 0d0
        do i = 1, NN
            PF(i,i) = 100d0
            PA(i,i) = 100d0
        end do
        do i = 1, NM
            RR(i,i) = sigma
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
            !接線型行列, 1d-n(n:integer) = 10^-n in double precision
            !
            call KF_matrix_get_M(NN, MM(1:NN,1:NN), XA(1:NN), window, 1d-8)
            !
            !予報誤差共分散行列
            !
            call KF_matrix_get_PF(NN, PF(1:NN,1:NN), PA(1:NN,1:NN), MM(1:NN,1:NN), alpha)
            !
            !カルマンゲイン
            !
            call KF_matrix_get_K(NN, NM, KK(1:NN,1:NM), PF(1:NN,1:NN), &
                & HH(1:NM,1:NN), RR(1:NM,1:NM))
            !
            !解析誤差共分散行列の計算
            !
            call KF_matrix_get_PA(NN, NM, PA(1:NN,1:NN), KK(1:NN,1:NM), &
                & HH(1:NM,1:NN), PF(1:NN,1:NN))
            !
            !解析値の計算
            !
            call KF_matrix_get_XA(NN, NM, XA(1:NN), XF(1:NN), YY(1:NM), &
                & KK(1:NN,1:NM), HH(1:NM,1:NN))
            !
            !RMSE,Trace(PF),Trace(PA)の計算
            !
            call get_rmse(NN, rmse, XA(1:NN), XX(1:NN))
            call get_trace(NN, PF(1:NN,1:NN), trace_pf)
            call get_trace(NN, PA(1:NN,1:NN), trace_pa)
            write(200,*) w_idx, YY(1), XF(1), XA(1), XX(1), rmse, sqrt(trace_pf/dble(NN)), sqrt(trace_pa/dble(NN))
            ! write(*,'(i10, 7f10.5)') w_idx, YY(1), XF(1), XA(1), XX(1), rmse, sqrt(trace_pf/dble(NN)), sqrt(trace_pa/dble(NN))
        end do

        close(200)

    end do
end program
