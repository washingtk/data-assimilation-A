!
! データ同化を実行するプログラムモジュール
! KF, EnKF, 4DVARを作成予定
!
module data_assimilation
    implicit none
contains

    !===========================================
    !extended Karman Filter matrix M
    !接線型行列Mを求めるサブルーチン
    !===========================================
    !n      :変数Xの配列サイズ
    !MM     :接線型行列
    !X_INIT :同化１ステップ前の解析値
    !t_step :同化間隔のタイムステップ
    !alpha  :摂動の大きさ
    subroutine KF_matrix_get_M(n, MM, X_INIT, t_step, alpha)
        use model_functions
        implicit none
        integer, intent(in) :: n, t_step
        real(8), intent(in) :: X_INIT(1:n)
        real(8), intent(out) :: MM(1:n,1:n)
        real(8), intent(in) :: alpha
        real(8) :: X_IN(1:n), X_OUT(1:n), X1(1:n), X2(1:n)
        integer :: idx_tstep, i, j
        !
        !摂動なしの場合
        !
        X_IN(:) = X_INIT(:)
        do idx_tstep = 1, t_step
            call runge_kutta_method(n, X_IN, X_OUT)
            X_IN(:) = X_OUT(:)
        end do
        X1(:) = X_OUT(:)

        !
        !接線型行列 M, (i:行のインデックス, j:列のインデックス)
        !
        do j = 1, n
            !
            !j成分に摂動を与えた場合
            !
            X_IN(:) = X_INIT(:)
            X_IN(j) = X_INIT(j) + 1d0*alpha
            do idx_tstep = 1, t_step
                call runge_kutta_method(n, X_IN, X_OUT)
                X_IN(:) = X_OUT(:)
            end do
            X2(:) = X_OUT(:)
            do i = 1, n
                MM(i,j) = ( X2(i) - X1(i) ) / alpha
            end do
        end do
    end subroutine


    !===========================================
    !extended Kalman Filter get PF
    !予報誤差共分散行列を求めるサブルーチン
    !===========================================
    !n      :変数Xの配列サイズ
    !PF     :予報誤差共分散行列
    !PA     :解析誤差共分散行列
    !MM     :接線型行列
    subroutine KF_matrix_get_PF(n, PF, PA, MM)
        implicit none
        integer, intent(in) :: n
        real(8), intent(out) :: MM(1:n,1:n), PF(1:n,1:n), PA(1:n,1:n)
        PF(1:n,1:n) = matmul( MM(1:n,1:n), matmul( PA(1:n,1:n), transpose(MM(1:n,1:n)) ) )
    end subroutine


    !===========================================
    !extended Kalman Filter get Kalman Gain, K
    !カルマンゲインKを求めるサブルーチン
    !===========================================
    !n      :変数Xの配列サイズ
    subroutine KF_matrix_get_K(n)
        implicit none
        integer, intent(in) :: n
    end subroutine


end module
