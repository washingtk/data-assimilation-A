module lorenz_model
    implicit none
    real(8), parameter :: FF = 8d0
    real(8), parameter :: dt = 0.01d0

contains

    !===========================================
    !力学モデル右辺の計算
    !===========================================
    !n:     配列サイズ
    !P_IN:  入力
    !P_OUT: 出力
    subroutine model_right_term(n, P_IN, P_OUT)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: P_IN(1:n)
        real(8), intent(out) :: P_OUT(1:n)
        real(8) :: P_TMP(-1:n+1)
        integer :: j
        !
        !initialize & boundary condition
        !
        P_TMP(1:n) = P_IN(1:n)
        P_TMP(-1)  = P_IN(n-1)
        P_TMP(0)   = P_IN(n)
        P_TMP(n+1) = P_IN(1)
        !
        !output, right side hand of equation
        !
        do j = 1, n
            P_OUT(j) = ( P_TMP(j+1) - P_TMP(j-2) ) * P_TMP(j-1) - P_TMP(j) + FF
        end do
    end subroutine model_right_term


    !===========================================
    !4次のRunge-Kutta 法による時間発展
    !===========================================
    !n:     配列サイズ
    !X_IN:  入力
    !X_OUT: 出力
    subroutine runge_kutta_method(n, X_IN, X_OUT)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: X_IN(1:n)
        real(8), intent(out) :: X_OUT(1:n)
        real(8) :: P(1:n), K1(1:n), K2(1:n), K3(1:n), K4(1:n)
        !
        !runge-kutta step1
        !
        P = X_IN
        call model_right_term(n, P, K1)
        !
        !runge-kutta step2
        !
        P = X_IN + 0.5d0 * dt * K1
        call model_right_term(n, P, K2)
        !
        !runge-kutta step3
        !
        P = X_IN + 0.5d0 * dt * K2
        call model_right_term(n, P, K3)
        !
        !runge-kutta step4
        !
        P = X_IN + dt * K3
        call model_right_term(n, P, K4)
        !
        !runge-kutta step5
        !
        X_OUT = X_IN + dt/6d0 * ( K1 + 2d0*K2 + 2d0*K3 + K4  )
    end subroutine runge_kutta_method


    !===========================================
    !extended Karman Filter matrix M
    !接線型行列Mを求めるサブルーチン
    !===========================================
    !n          :変数Xの配列サイズ
    !MM         :接線型行列
    !X_INIT     :同化１ステップ前の解析値
    !assim_step :同化間隔のタイムステップ(1同化間隔における時間積分回数)
    !alpha      :摂動の大きさ
    subroutine KF_matrix_get_M(n, MM, X_INIT, assim_step, alpha)
        implicit none
        integer, intent(in) :: n, assim_step
        real(8), intent(in) :: X_INIT(1:n)
        real(8), intent(out) :: MM(1:n,1:n)
        real(8), intent(in) :: alpha
        real(8) :: X_IN(1:n), X_OUT(1:n), X1(1:n), X2(1:n)
        integer :: idx_tstep, i, j
        !
        !摂動なしの場合
        !
        X_IN(:) = X_INIT(:)
        do idx_tstep = 1, assim_step
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
            do idx_tstep = 1, assim_step
                call runge_kutta_method(n, X_IN, X_OUT)
                X_IN(:) = X_OUT(:)
            end do
            X2(:) = X_OUT(:)
            do i = 1, n
                MM(i,j) = ( X2(i) - X1(i) ) / alpha
            end do
        end do
    end subroutine

end module
