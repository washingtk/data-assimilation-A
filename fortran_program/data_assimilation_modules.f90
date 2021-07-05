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
    !n          :変数Xの配列サイズ
    !MM         :接線型行列
    !X_INIT     :同化１ステップ前の解析値
    !assim_step :同化間隔のタイムステップ(1同化間隔における時間積分回数)
    !alpha      :摂動の大きさ
    subroutine KF_matrix_get_M(n, MM, X_INIT, assim_step, alpha)
        use model_functions
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


    !===========================================
    !extended Kalman Filter get PF
    !予報誤差共分散行列を求めるサブルーチン
    !===========================================
    !n          :変数Xの配列サイズ
    !PF         :予報誤差共分散行列
    !PA         :解析誤差共分散行列
    !MM         :接線型行列
    !inflation  :インフレーションファクター
    subroutine KF_matrix_get_PF(n, PF, PA, MM, inflation)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: MM(1:n,1:n), PA(1:n,1:n), inflation
        real(8), intent(out) :: PF(1:n,1:n)
        integer :: i
        PF(1:n,1:n) = matmul( MM(1:n,1:n), matmul( PA(1:n,1:n), transpose(MM(1:n,1:n)) ) ) * inflation
    end subroutine


    !===========================================
    !extended Kalman Filter get Kalman Gain, K
    !カルマンゲインKを求めるサブルーチン
    !lapack を利用して逆行列を計算する
    !===========================================
    !n:     変数Xの配列サイズ
    !m:     観測変数Yの配列サイズ
    !KK:    カルマンゲイン
    !PF:    予報誤差行列
    !HH:    観測行列
    !RR:    観測誤差行列
    subroutine KF_matrix_get_K(n, m, KK, PF, HH, RR)
        implicit none
        integer, intent(in) :: n, m
        real(8), intent(in) :: PF(1:n,1:n)
        real(8), intent(in) :: HH(1:m,1:n)
        real(8), intent(in) :: RR(1:m,1:m)
        real(8), intent(out) :: KK(1:n,1:m)
        integer :: INFO, i
        real(8) :: TEMP(1:m,1:m)
        real(8) :: TEMP2(1:m,1:m)
        TEMP = RR + matmul( HH, matmul( PF, transpose(HH) ) )
        TEMP2 = TEMP
        ! LAPACKのルーチンを呼び出して行列をコレスキー分解します。
        ! K には上三角行列にコレスキー分解した結果が格納されます。
        ! K を書き換えられたくない場合は注意が必要です。
        call DPOTRF( 'U', m, TEMP(1:m,1:m), m, INFO )
        ! LAPACKのルーチンを直接呼び出して
        ! DPOTRFでコレスキー分解した結果を用いてもとの行列の逆行列を求めます。
        ! もとの行列 K の逆行列が上三角要素に格納されています。
        call DPOTRI( 'U', m, TEMP(1:m,1:m), m, INFO )
        ! 下半分の要素を埋めます。
        do i = 1, m
            TEMP(i+1:m,i) = TEMP(i,i+1:m)
        end do
        !
        !カルマンゲインの計算
        !
        KK = matmul( PF, matmul( transpose(HH), TEMP ) )
        !
        !逆行列のチェック
        !
        ! TEMP = abs(matmul(TEMP, TEMP2))
        ! print*, 'TEMP'
        ! write(*,'(40F5.1)') ( TEMP(i,:), i=1,m )
        ! print*, 'KK'
        ! write(*,'(40F5.1)') ( KK(i,:), i=1,n )
    end subroutine


    !===========================================
    !extended Kalman Filter get XA
    !解析値XAを求めるサブルーチン
    !===========================================
    !n:     変数Xの配列サイズ
    !m:     観測変数Yの配列サイズ
    !XA:    解析値
    !XF:    予報値
    !YY:    観測値
    !KK:    カルマンゲイン
    !HH:    観測行列
    subroutine KF_matrix_get_XA(n, m, XA, XF, YY, KK, HH)
        implicit none
        integer, intent(in)  :: n, m
        real(8), intent(out) :: XA(1:n)
        real(8), intent(in)  :: XF(1:n)
        real(8), intent(in)  :: YY(1:m)
        real(8), intent(in)  :: KK(1:n,1:m)
        real(8), intent(in)  :: HH(1:m,1:n)
        real(8) :: INOVATION(1:m)
        INOVATION(:) = YY(:) - matmul( HH(:,:), XF(:) )
        XA(:) = XF(:) + matmul( KK(:,:), INOVATION(:) )
    end subroutine


    !===========================================
    !extended Kalman Filter get PA
    !予報誤差共分散行列を求めるサブルーチン
    !===========================================
    !n:     変数Xの配列サイズ
    !m:     観測変数Yの配列サイズ
    !PA:    解析誤差共分散行列
    !KK:    カルマンゲイン
    !HH:    観測行列
    !PF:    予報誤差共分散行列
    subroutine KF_matrix_get_PA(n, m, PA, KK, HH, PF)
        implicit none
        integer, intent(in)  :: n, m
        real(8), intent(out) :: PA(1:n,1:n)
        real(8), intent(in)  :: KK(1:n,1:m)
        real(8), intent(in)  :: HH(1:m,1:n)
        real(8), intent(in)  :: PF(1:n,1:n)
        PA(:,:) = PF(:,:) - matmul( KK(:,:), matmul( HH(:,:), PF(:,:) ) )
    end subroutine


    !===========================================
    !Ensenble Kalman Filter get PF
    !重みつき行列の計算
    !===========================================
    !n          :変数Xの配列サイズ
    !dist       :局所化距離
    !WW         :重み行列
    subroutine EnKF_matrix_get_WW(n, dist, WW)
        implicit none
        integer, intent(in)  :: n
        real(8), intent(in)  :: dist
        real(8), intent(out) :: WW(1:n,1:n)
        integer :: i, j
        real(8) :: temp
        do i = 1, n
            do j = 1, n
                temp = abs(i-j)
                if (temp>n/2) then
                    temp = 20 - temp
                end if
                WW(i,j) = exp( -temp**2 / (2d0*dist) )
            end do
        end do
    end subroutine



    !===========================================
    !Ensenble Kalman Filter get PF
    !予報誤差共分散行列を求めるサブルーチン
    !===========================================
    !n          :変数Xの配列サイズ
    !XF         :予報アンサンブル
    !PF         :予報誤差共分散行列
    !inflation  :共分散膨張係数
    subroutine EnKF_matrix_get_PF(n, ens_size, PF, XF, inflation)
        implicit none
        integer, intent(in) :: n, ens_size
        real(8), intent(in) :: inflation
        real(8), intent(inout) :: XF(1:ens_size,1:n)
        real(8), intent(out) :: PF(1:n,1:n)
        integer :: i_idx, j_idx, ens_idx
        real(8) :: i_temp, j_temp
        real(8) :: xf_avr(1:n)
        !
        !予報平均ベクトル
        !
        xf_avr(1:n) = sum( XF(1:ens_size,1:n), dim=1 ) / dble(ens_size)
        !
        !共分散膨張
        !
        do ens_idx = 1, ens_size
            XF(ens_idx,1:n) = xf_avr(1:n) + inflation * ( XF(ens_idx,1:n) - xf_avr(1:n) )
        end do
        !
        !共分散誤差行列
        !
        PF(:,:) = 0d0
        do ens_idx = 1, ENS_SIZE
            do i_idx = 1, n
                i_temp = XF(ens_idx,i_idx) - xf_avr(i_idx)
                do j_idx = 1, n
                    j_temp = XF(ens_idx,j_idx) - xf_avr(j_idx)
                    PF(i_idx,j_idx) = PF(i_idx,j_idx) + i_temp * j_temp
                end do
            end do
        end do
        PF(:,:) = PF(:,:) / dble(ENS_SIZE-1)
    end subroutine


    !===========================================
    !Ensenble Kalman Filter get PA
    !解析誤差共分散行列を求めるサブルーチン
    !===========================================
    !n      :変数Xの配列サイズ
    !XA     :解析アンサンブル
    !PA     :予報誤差共分散行列
    subroutine EnKF_matrix_get_PA(n, ens_size, PA, XA)
        implicit none
        integer, intent(in) :: n, ens_size
        real(8), intent(in) :: XA(1:ens_size,1:n)
        real(8), intent(out) :: PA(1:n,1:n)
        integer :: i_idx, j_idx, ens_idx
        real(8) :: i_temp, j_temp
        real(8) :: xa_avr(1:n)
        xa_avr(1:n) = sum( XA(1:ens_size,1:n), dim=1 ) / (1d0*ens_size)
        PA(:,:) = 0d0
        do ens_idx = 1, ENS_SIZE
            do i_idx = 1, n
                i_temp = XA(ens_idx,i_idx) - xa_avr(i_idx)
                do j_idx = 1, n
                    i_temp = XA(ens_idx,j_idx) - xa_avr(j_idx)
                    PA(i_idx,j_idx) = PA(i_idx,j_idx) + i_temp * j_temp
                end do
            end do
        end do
        PA(:,:) = PA(:,:) / dble(ENS_SIZE-1)
    end subroutine


    !===========================================
    !EnKF Kalman Filter get XA
    !解析値XAを求めるサブルーチン
    !摂動観測法を利用
    !===========================================
    !n:         変数Xの配列サイズ
    !m:         観測変数Yの配列サイズ
    !ens_size:  アンサンブルサイズ
    !sigma:     観測誤差
    !XA:        解析値
    !XF:        予報値
    !YY:        観測値
    !KK:        カルマンゲイン
    !HH:        観測行列
    subroutine EnKF_matrix_get_XA(n, m, ens_size, sigma, XA, XF, YY, KK, HH)
        use other_subroutines
        implicit none
        integer, intent(in)  :: n, m, ens_size
        real(8), intent(out) :: XA(1:ens_size,1:n)
        real(8), intent(in)  :: sigma
        real(8), intent(in)  :: XF(1:ens_size,1:n)
        real(8), intent(in)  :: YY(1:m)
        real(8), intent(in)  :: KK(1:n,1:m)
        real(8), intent(in)  :: HH(1:m,1:n)
        integer :: ens_idx
        real(8) :: INOVATION(1:m), YY_TMP(1:m)
        do ens_idx = 1, ens_size
            call box_muller_without_seedindex(m, YY_TMP(1:m), 0d0, sigma)
            YY_TMP(1:m) = YY(1:m) + YY_TMP(1:m)
            INOVATION(:) = YY_TMP(:) - matmul( HH(:,:), XF(ens_idx,:) )
            XA(ens_idx,:) = XF(ens_idx,:) + matmul( KK(:,:), INOVATION(:) )
        end do
    end subroutine


    !===========================================
    !LETKF
    !DXF
    !===========================================
    !n:         変数Xの配列サイズ
    !ens_size:  アンサンブルサイズ
    !XF:        予報値アサンブル
    !D_XF:      予報アンサンブル摂動
    subroutine LETKF_get_DXF(n, ens_size, XF, D_XF)
        implicit none
        integer, intent(in)  :: n, ens_size
        real(8), intent(in)  :: XF(1:n, 1:ens_size)
        real(8), intent(out) :: D_XF(1:n, 1:ens_size)
        real(8) :: XF_AVR(1:n)
        integer :: idx
        ! XF_AVR(1:n) = sum(XF(:,:), dim=2)/ens_size
        ! do idx = 1, ens_size
        !     D_XF(1:n,idx) = XF(1:n) - XF_AVR(1:n)
        ! end do
    end subroutine



    !===========================================
    !LETKF
    !DYF
    !===========================================
    !n:         変数Xの配列サイズ
    !m:         観測変数Yの配列サイズ
    !ens_size:  アンサンブルサイズ
    !HH:        観測行列
    !D_XF:      予報アンサンブル摂動
    !D_YF:      予報アンサンブル摂動
    subroutine LETKF_get_DYF(n, m, ens_size, HH, D_XF, D_YF)
        implicit none
        integer, intent(in)  :: n, m, ens_size
        real(8), intent(in) :: HH(1:m,1:n)
        real(8), intent(in)  :: D_XF(1:n, 1:ens_size)
        real(8), intent(out) :: D_YF(1:m, 1:ens_size)
        integer :: idx
        ! do ens_idx = 1, ens_size
        !     D_YF(1:m,ens_idx) = matmul( HH(1:m,1:n), D_XF(1:n,ens_idx) )
        ! end do
    end subroutine


    !===========================================
    !LETKF
    !固有値分解
    !RR は対角行列
    !===========================================
    !n:         変数Xの配列サイズ
    !m:         観測変数Yの配列サイズ
    !ens_size:  アンサンブルサイズ
    !D_XF:      予報アンサンブル摂動
    !D_YF:      予報アンサンブル摂動
    !HH:        観測行列
    !RR:        観測誤差共分散行列
    !UU:        固有ベクトルの正規直行行列
    !DD:        固有値の対角行列
    subroutine LETKF_eigen_decomp(n, m, ens_size, D_XF, D_YF, HH, RR, UU, DD)
        implicit none
        !
        !引数
        !
        integer, intent(in)  :: n, m, ens_size
        real(8), intent(in)  :: D_XF(1:n,1:ens_size)
        real(8), intent(in)  :: D_YF(1:m,1:ens_size)
        real(8), intent(in)  :: HH(1:m,1:n)
        real(8), intent(in)  :: RR(1:m,1:m)
        real(8), intent(out) :: UU(1:ens_size,1:ens_size)
        real(8), intent(out) :: DD(1:ens_size,1:ens_size)
        !
        !サブルーチン中
        !
        real(8) :: RR_inv(1:m,1:m)
        real(8) :: II(1:ens_size,1:ens_size)
        real(8) :: TEMP(1:ens_size,1:ens_size)
        integer :: i, j, ens_idx
        !
        !固有値分解
        !
        real(8) :: WORK(1:ens_size), W(1:ens_size)
        integer :: LWORK, INFO
        !
        !単位行列
        !
        II(:,:) = 0d0
        do i = 1, ens_size
            II(i,i) = 1d0
        end do
        !
        !観測誤差共分散行列
        !
        RR_inv(:,:) = 0d0
        do i = 1, m
            RR_inv(i,i) = 1d0 / RR(i,i)
        end do
        !
        !固有値分解する行列
        !
        TEMP(1:ens_size,1:ens_size) = (ens_size-1) * II(:,:) &
            & + matmul( transpose(D_YF(1:m,1:ens_size)), matmul(RR_inv(1:m,1:m), D_YF(1:m,1:ens_size)) )
        !
        !実対称行列の固有値問題
        !
        LWORK = 3*ens_size - 1
        UU(:,:) = TEMP(:,:)
        call DSYEV('V', 'U', ens_size, UU(1:ens_size,1:ens_size), ens_size, &
            & W(1:ens_size), WORK(1:ens_size), LWORK, INFO)
        if (INFO==0) then
            DD(:,:) = 0d0
            do ens_idx = 1, ens_size
                DD(ens_idx,ens_idx) = W(ens_idx)
            end do
        else
            print*, 'eigenvalue decomposition failed'
            call exit(1)
        end if
    end subroutine



    !===========================================
    !LETKF
    !予報値の修正量
    !===========================================
    !n:         変数Xの配列サイズ
    !m:         観測変数Yの配列サイズ
    !ens_size:  アンサンブルサイズ
    !D_XF:      予報アンサンブル摂動
    !D_YF:      予報アンサンブル摂動
    !HH:        観測行列
    !RR:        観測誤差共分散行列
    !UU:        固有ベクトルの正規直行行列
    !DD:        固有値の対角行列
    subroutine LETKF_get_INOV(n, m, ens_size, XF, D_XF, D_YF, HH, RR, UU, DD, INOV)
        implicit none
        !
        !引数
        !
        integer, intent(in)  :: n, m, ens_size
        real(8), intent(in)  :: XF(1:n, 1:ens_size)
        real(8), intent(in)  :: D_XF(1:n,1:ens_size)
        real(8), intent(in)  :: D_YF(1:m,1:ens_size)
        real(8), intent(in)  :: HH(1:m,1:n)
        real(8), intent(in)  :: RR(1:m,1:m)
        real(8), intent(in)  :: UU(1:ens_size,1:ens_size)
        real(8), intent(in)  :: DD(1:ens_size,1:ens_size)
        real(8), intent(in)  :: INOV(1:ens_size,1:ens_size)
    end subroutine

end module
