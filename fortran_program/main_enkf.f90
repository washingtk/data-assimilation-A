!===========================================
!Ensemble Kalman Filter の実装
!===========================================
program main
    use global_variable, only: NN, NM, t_1year
    use model_functions
    use other_subroutines
    use data_assimilation
    implicit none
    integer :: w_idx, t_idx, n_idx, ens_idx, i, j
    integer, parameter :: window = 5
    integer, parameter :: ENS_SIZE = 1000
    real(8) :: XX(1:NN)
    real(8) :: XX_NEX(1:NN)
    real(8) :: XF(1:ENS_SIZE,1:NN)
    real(8) :: XA(1:ENS_SIZE,1:NN)
    real(8) :: YY(1:NM)
    real(8) :: ER(1:NN)
    !---------------
    !KF parameter
    !---------------
    !PF:    予報誤差行列
    !PA:    解析誤差行列
    !HH:    観測行列
    !RR:    観測誤差行列
    !KK:    カルマンゲイン
    !WW:    共分散局所化行列
    !sigma: 観測誤差
    !dist:  局所化距離
    !---------------
    real(8) :: PF(1:NN,1:NN)
    real(8) :: PA(1:NN,1:NN)
    real(8) :: HH(1:NM,1:NN)
    real(8) :: RR(1:NM,1:NM)
    real(8) :: KK(1:NN,1:NM)
    real(8) :: WW(1:NN,1:NN)
    real(8) :: rmse, rmse_xf, trace_pf, trace_pa
    real(8), parameter :: sigma = 1d0
    real(8), parameter :: dist = 20d0
    character(500) :: dirname, filename1, filename2

    !
    !make directory
    !
    write(dirname, '("output_enkf_data/")')
    call makedirs(dirname)
    write(filename1, '("enkf_out", "window", i3.3, "NM", i2.2, "ens_num", i5.5, "dist", i2.2, ".txt")') &
        & window, NM, ENS_SIZE, int(dist)
    write(filename2, '("enkf_forcast", "window", i3.3, "NM", i2.2, "ens_num", i5.5, ".txt")') &
        & window, NM, ENS_SIZE
    filename1 = trim(dirname)//trim(filename1)
    filename2 = trim(dirname)//trim(filename2)
    write(*,*) 'save data in >> ', trim(filename1)
    write(*,*) 'save data in >> ', trim(filename2)
    open(200, file=filename1, status="replace")
    open(300, file=filename2, status="replace")

    !
    !check
    !
    write(*,*) "--------------------------------------"
    write(*,"(a30,i9)") 'nubmber of variable: ', NN
    write(*,"(a30,i9)") 'nubmber of obervation: ', NM
    write(*,"(a30,i9)") 'ensemble size: ', ENS_SIZE
    write(*,"(a30,i9)") 'NN/NM: ', int(NN/NM)
    write(*,"(a30,i9)") 'window: ', window
    write(*,"(a30,f9.4)") 'dt: ', dt
    write(*,"(a30,f9.4)") 'F: ', FF
    write(*,*) "--------------------------------------"

    !
    !initialize values 1
    !
    XX(:) = FF
    XX(20) = FF + 0.008d0
    !
    !time step in model
    !
    do t_idx = 1, t_1year
        call runge_kutta_method(NN, XX, XX_NEX)
        XX = XX_NEX
    end do
    !
    !行列の初期化
    !
    RR(:,:) = 0d0
    PF(:,:) = 0d0
    PA(:,:) = 0d0
    do i = 1, NN
        PF(i,i) = 1d0
        PA(i,i) = 1d0
    end do
    do i = 1, NM
        RR(i,i) = sigma
    end do
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

    !
    !initial ensemble
    !真値から平均2、分散１に従う正規乱数を誤差に加えてアンサンブルを作成
    do ens_idx = 1, ENS_SIZE
        call box_muller(NM, ER(1:NN), 0d0, 10d0, ens_idx+100000)
        XA(ens_idx,1:NN) = XX(1:NN) + ER(1:NN)
    end do

    !
    !予報共分散の局所化行列
    !
    call EnKF_matrix_get_WW(NN, dist, WW(1:NN,1:NN))

    !
    !time step in model
    !
    do w_idx = 1, 1000
        do t_idx = 1, window
            call runge_kutta_method(NN, XX, XX_NEX)
            XX(:) = XX_NEX(:)
            do ens_idx = 1, ENS_SIZE
                call runge_kutta_method(NN, XA(ens_idx,1:NN), XF(ens_idx,1:NN))
            end do
            XA(:,:) = XF(:,:)
        end do
        !
        !観測値
        !
        call box_muller(NM, ER(1:NN), 0d0, sigma, w_idx)
        YY(1:NM) = matmul( HH(1:NM,1:NN), XX(1:NN)+ER(1:NN) )
        !
        !予報誤差共分散行列
        !
        call EnKF_matrix_get_PF(NN, ENS_SIZE, PF(1:NN,1:NN), XF(1:ENS_SIZE,1:NN), 1.0d0)
        !
        !予報共分散の局所化
        !
        PF(1:NN,1:NN) = PF(:,:) * WW(:,:)
        !
        !カルマンゲイン
        !
        call KF_matrix_get_K(NN, NM, KK(1:NN,1:NM), PF(1:NN,1:NN), &
            & HH(1:NM,1:NN), RR(1:NM,1:NM))
        !
        !解析値の計算
        !
        call EnKF_matrix_get_XA(NN, NM, ENS_SIZE, sigma, XA(1:ENS_SIZE,1:NN), XF(1:ENS_SIZE,1:NN), &
            & YY(1:NM), KK(1:NN,1:NM), HH(1:NM,1:NN))
        !
        !解析誤差共分散行列の計算
        !
        call EnKF_matrix_get_PA(NN, ENS_SIZE, PA(1:NN,1:NN), XA(1:ENS_SIZE,1:NN))
        !
        !RMSEの計算
        !
        call get_rmse(NN, rmse, sum(XA(1:ENS_SIZE,1:NN), dim=1)/dble(ENS_SIZE), XX(1:NN))
        call get_rmse(NN, rmse_xf, sum(XF(1:ENS_SIZE,1:NN), dim=1)/dble(ENS_SIZE), XX(1:NN))
        call get_trace(NN, PF(1:NN,1:NN), trace_pf)
        call get_trace(NN, PA(1:NN,1:NN), trace_pa)
        write(200,*) w_idx, YY(1), XF(1,1), sum(XA(1:ENS_SIZE,1))/dble(ENS_SIZE), XX(1), &
            & rmse, rmse_xf, sqrt(trace_pf/NN), sqrt(trace_pa/NN)
        write(300,*) XF(1:ENS_SIZE,1)
        ! write(*,*) w_idx, YY(1), XF(1,1), XA(1,1), XX(1), rmse, trace
    end do

    close(200)

end program
