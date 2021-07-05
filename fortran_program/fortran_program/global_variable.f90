!
!グローバル変数の定義
!
module global_variable
    implicit none
    !NN:        モデルの変数の数
    !MM:        観測点数
    !t_max:     時間積分回数
    !t_1year:   1年相当の時間積分回数(F=8)
    !t_1day:    1日相当の時間積分回数(F=8)
    !dt:        モデルのタイムステップ
    !FF:        モデルの強制項
    !pi:        円周率
    !XX:        モデル変数
    integer, parameter :: NN = 40
    integer, parameter :: NM = 40
    integer, parameter :: t_max = 5000
    integer, parameter :: t_1year = 7300
    integer, parameter :: t_1day = 20
    real(8), parameter :: dt = 0.01d0
    real(8), parameter :: FF = 8d0
    real(8), parameter :: pi = 3.14159265358979323846264338327950288d0

end module global_variable
