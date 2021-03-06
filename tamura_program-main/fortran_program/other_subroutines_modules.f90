module other_subroutines
    implicit none

contains

    !
    !ガウス関数
    !
    function gaussian(x, mean, var) result(output)
        implicit none
        real(8), intent(in) :: x, mean, var
        real(8) :: output
        output = exp( -(x-mean)**2/(2.0*var**2) )
        return
    end function gaussian


    !
    !boxmuller method
    !正規乱数の生成
    !
    subroutine box_muller(nsample, rnd, mean, std, idx_seed)
        use global_variable, only: pi
        implicit none
        integer, intent(in) :: nsample, idx_seed
        real(8), intent(in) :: mean, std
        real(8), intent(out) :: rnd(nsample)
        real(8), allocatable :: x(:), y(:), r(:), t(:)
        integer :: i, seedsize
        integer, allocatable :: seed(:)
        !
        !allocation
        !
        allocate(x(nsample), y(nsample), r(nsample), t(nsample))
        !
        !random setting
        !
        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        call random_seed(get=seed(:))
        seed = idx_seed
        call random_seed(put=seed(:))
        !
        !generate randam variable
        !
        call random_number(x)
        call random_number(y)
        !
        !normal randam variable N(0,1)
        !
        r = sqrt(-2d0 * log(x))
        t = 2d0 * pi * y
        rnd = r * cos(t)
        !
        !requred randam variable ~ N(mean, std)
        !
        do i = 1, nsample
            rnd(i) = rnd(i)*std + mean
        end do
        !
        !deallocation
        !
        deallocate(seed)
    end subroutine box_muller


    !
    !directory を作成するサブルーチン
    !
    subroutine makedirs(outdir)
        character(len=*), intent(in) :: outdir
        character(len=900) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        ! write(*, *) trim(command)
        call system(command)
    end subroutine makedirs

end module
