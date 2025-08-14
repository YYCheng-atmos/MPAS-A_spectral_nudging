!==================================================================================================
module diffusion_filter_coefficients

use mpas_kind_types
use mpas_derived_types
use mpas_pool_routines

use mpas_atmphys_constants
use mpas_atmphys_utilities
use mpas_atmphys_vars


implicit none

!==================================================================================================
!> \brief MPAS grid diffusion-based filter.
!> 
!> \author Yiyuan Cheng (send comments to yiyuan_cheng@foxmail.com).
!> \date 2025-05-06.
!>
!> \details
!> Subroutines in mpas_atmphys_fdda:
!> --------------------------------
!> 1. compute_filter_coefficients:
!>    - Main subroutine for computing the `p` vector (called from subroutine apply_fdda).
!>    - Double precision results match the GCM-filter. Single precision results are biased by ~1e-5.
!>
!> Development Log:
!> ----------------
!> 2025.05:
!> - Implemented basic GCM-filter coefficient calculations.
!> - Achieved accurate results under Gaussian grid.
!> - Noted potential need for second-order precision due to the absence of LAPACK or similar libraries.
!> - Identified significant issues with taper coefficient calculations.
!>
!> 2025.06:
!> - Modified interpolation methods to align with PchipInterpolator in GCM-filter's Python implementation.
!> - Defined a new series of functions for coefficient calculations.
!> - Noted ~1e-4 differences in TAPER coefficients compared to Python implementation.
!==================================================================================================


! PCHIP 插值类型
type :: pchip_interpolator
    real(kind=RKIND), allocatable :: x(:)      ! 插值节点
    real(kind=RKIND), allocatable :: y(:)      ! 节点值
    real(kind=RKIND), allocatable :: d(:)      ! 节点导数
    integer, allocatable :: n                      ! 节点数量
end type pchip_interpolator

contains


    ! 端点导数计算
    real(kind=RKIND) function pchip_endpoint(h1, h2, del1, del2)
        real(kind=RKIND), intent(in) :: h1, h2, del1, del2
        real(kind=RKIND) :: d, d_prop
        
        d = ((2.0d0*h1 + h2)*del1 - h1*del2) / (h1 + h2)
        
        ! 单调性约束
        if (d * del1 < 0.0d0) then
            d = 0.0d0
        else if (abs(d) > 3.0d0*abs(del1)) then
            d = 3.0d0 * del1
        end if
        
        pchip_endpoint = d
    end function pchip_endpoint

    ! PCHIP 插值计算
    real(kind=RKIND) function pchip_interpolate(interp, x_val)
        type(pchip_interpolator), intent(inout) :: interp
        real(kind=RKIND), intent(in) :: x_val
        integer :: i, k
        real(kind=RKIND) :: h, t, c0, c1, c2, c3

        ! 查找区间
        k = 1
        do i = 1, interp%n-1
            if (x_val >= interp%x(i) .and. x_val <= interp%x(i+1)) then
                k = i
                exit
            end if
        end do
        
        ! 三次 Hermite 插值

        h = interp%x(k+1) - interp%x(k)
        t = (x_val - interp%x(k)) / h
        
        c0 = interp%y(k)
        c1 = interp%d(k)

        c2 = (3.0d0*(interp%y(k+1)-interp%y(k))/h - 2.0d0*interp%d(k) - interp%d(k+1))
        c3 = (2.0d0*(interp%y(k)-interp%y(k+1))/h + interp%d(k) + interp%d(k+1))
        
        pchip_interpolate = c0 + c1*h*t + c2*h*t**2 + c3*h*t**3
    end function pchip_interpolate

! 主计算函数
subroutine compute_filter_coefficients(n_steps, filter_scale, dx_min, filter_shape, p)

    integer, intent(in) :: n_steps
    real(kind=RKIND), intent(in) :: filter_scale, dx_min
    character(*), intent(in) :: filter_shape
    real(kind=RKIND), intent(out) :: p(n_steps + 1)

    real(kind=RKIND), parameter :: pi = 3.141592653589793238462643383279502884197_RKIND
    integer :: n, i, k, j, n_points, order
    real(kind=RKIND) :: s_max, F_1, x_val, t, xk, term, F_val, integral
    real(kind=RKIND), allocatable :: M(:, :), b(:), c_hat(:), points(:), weights(:), tmp(:), phi(:)
    real(kind=RKIND), allocatable :: T_vals(:)  ! 用于存储 Chebyshev 多项式
    character :: filter_type
    type(pchip_interpolator) :: taper_interp
    real(kind=RKIND) :: interp_x(4), interp_y(4), transition_width

    n = n_steps - 1  ! c_hat 向量的长度
    n_points = n_steps + 1  ! 积分点数
    transition_width = pi
    s_max = 2.0d0 * (2.0d0 / dx_min)**2  ! 对于 2D 情况

    filter_type = 'g'
    if (index(filter_shape, 'taper') > 0) filter_type = 't'

    ! 为 Taper 滤波器创建 PCHIP 插值器
    if (filter_type == 't') then
        interp_x = [0.0d0, &
                   2.0d0 * pi / (transition_width * filter_scale), &
                   2.0d0 * pi / filter_scale, &
                   8.0d0 * sqrt(s_max)]
        interp_y = [1.0d0, 1.0d0, 0.0d0, 0.0d0]
        call create_pchip_interpolator(interp_x, interp_y, taper_interp)
    end if

    allocate(M(n, n), b(n), c_hat(n), points(n_points), weights(n_points), &
             tmp(n_points), phi(n_points), T_vals(0:n_points-1))

    ! 构造质量矩阵 M
    M = 0.0d0
    do i = 1, n
        do j = 1, n
            if (i == j) then
                M(i, j) = pi
                if (i == 1) M(i, j) = 3.0d0 * pi / 2.0d0
            else if (abs(i - j) == 2) then
                M(i, j) = -pi / 2.0d0
            end if
        end do
    end do

    ! 计算 F(1.0)
    if (filter_type == 'g') then
        F_1 = exp(-s_max * (1.0d0 + 1.0d0) / 2.0d0 * filter_scale**2 / 24.0d0)
    else
        x_val = sqrt(s_max)  ! 在 t=1 时计算 x_val
        F_1 = pchip_interpolate(taper_interp, x_val)
    end if

    ! 生成 Chebyshev-Gauss 积分点和权重
    call chebgauss(n_points, points, weights)

    ! 计算右端项 b
    b = 0.0d0
    do i = 1, n
        ! 基函数: phi = T_{i-1} - T_{i+1} (阶数从0开始)
        tmp = 0.0d0
        order = i - 1  ! 对应 Python 中的 i (阶数偏移)
        if (order < n_points) tmp(order + 1) = 1.0d0
        if (order + 2 < n_points) tmp(order + 3) = -1.0d0

        do k = 1, n_points
            xk = points(k)
            call chebyshev_polynomials(n_points-1, xk, T_vals)
            phi(k) = dot_product(tmp(1:n_points), T_vals(0:n_points-1))
        end do

        integral = 0.0d0
        do k = 1, n_points
            t = points(k)
            if (filter_type == 'g') then
                F_val = exp(-s_max * (t + 1.0d0) / 2.0d0 * filter_scale**2 / 24.0d0)
            else
                x_val = sqrt((t + 1.0d0) * s_max / 2.0d0)
                F_val = pchip_interpolate(taper_interp, x_val)
            end if
            term = weights(k) * phi(k) * (F_val - ((1.0d0 - t)/2.0d0 + F_1*(t + 1.0d0)/2.0d0))
            integral = integral + term
        end do
        b(i) = integral
    end do

    ! 调用自定义线性方程组求解器
    call solve_linear_system(M, b, c_hat, n)

    ! 转换为 Chebyshev 基系数 p
    p = 0.0d0
    p(1) = c_hat(1) + (1.0d0 + F_1)/2.0d0
    
    if (n >= 1) p(2) = c_hat(2) - (1.0d0 - F_1)/2.0d0

    do i = 3, n_steps - 1
        if (i - 2 >= 1 .and. i <= n) then
            p(i) = c_hat(i) - c_hat(i-2)
        end if
    end do

    ! 修正最后两项
    if (n >= 2) p(n_steps) = -c_hat(n-1)
    p(n_steps + 1) = -c_hat(n)

    deallocate(M, b, c_hat, points, weights, tmp, phi, T_vals)
    if (filter_type == 't') then
        deallocate(taper_interp%x, taper_interp%y, taper_interp%d, taper_interp%n)
    end if

contains

    subroutine chebgauss(n, x, w)
        integer, intent(in) :: n
        real(kind=RKIND), intent(out) :: x(n), w(n)
        integer :: k
        real(kind=RKIND) :: theta
        do k = 1, n
            theta = (2*k - 1) * pi / (2.0d0 * n)
            x(k) = cos(theta)
            w(k) = pi / n
        end do
    end subroutine chebgauss

    subroutine chebyshev_polynomials(n, x, T)
        integer, intent(in) :: n
        real(kind=RKIND), intent(in) :: x
        real(kind=RKIND), intent(out) :: T(0:n)
        integer :: k
        T(0) = 1.0d0
        if (n >= 1) T(1) = x
        do k = 2, n
            T(k) = 2.0d0 * x * T(k-1) - T(k-2)
        end do
    end subroutine chebyshev_polynomials

    ! 带部分主元选择的高斯消元法
    subroutine solve_linear_system(A, b, x, n)
        integer, intent(in) :: n
        real(kind=RKIND), intent(in) :: A(n, n), b(n)
        real(kind=RKIND), intent(out) :: x(n)
        real(kind=RKIND) :: Aug(n, n+1), factor, max_val, temp
        integer :: i, j, k, max_row

        ! 设置增广矩阵
        Aug(1:n, 1:n) = A
        Aug(1:n, n+1) = b

        ! 高斯消元过程
        do k = 1, n-1
            ! 部分主元选择
            max_val = abs(Aug(k, k))
            max_row = k
            do i = k+1, n
                if (abs(Aug(i, k)) > max_val) then
                    max_val = abs(Aug(i, k))
                    max_row = i
                end if
            end do
            
            ! 交换行
            if (max_row /= k) then
                do j = k, n+1
                    temp = Aug(k, j)
                    Aug(k, j) = Aug(max_row, j)
                    Aug(max_row, j) = temp
                end do
            end if
            
            ! 消元
            do i = k+1, n
                factor = Aug(i, k) / Aug(k, k)
                do j = k+1, n+1
                    Aug(i, j) = Aug(i, j) - factor * Aug(k, j)
                end do
            end do
        end do
        
        x(n) = Aug(n, n+1) / Aug(n, n)
        do i = n-1, 1, -1
            x(i) = Aug(i, n+1)
            do j = i+1, n
                x(i) = x(i) - Aug(i, j) * x(j)
            end do
            x(i) = x(i) / Aug(i, i)
        end do
        
    end subroutine solve_linear_system
    
    ! 创建 PCHIP 插值器
    subroutine create_pchip_interpolator(x, y, interp)
        real(kind=RKIND), intent(in) :: x(:), y(:)
        type(pchip_interpolator), intent(inout) :: interp
        integer :: n, i
        real(kind=RKIND), allocatable :: h(:), delta(:)
        real(kind=RKIND) :: w1, w2, d_max, d_min

        n = size(x)
        allocate(interp%x(n), interp%y(n), interp%d(n))
        interp%n = n
        interp%x = x
        interp%y = y
        
        ! 计算区间长度和斜率
        allocate(h(n-1), delta(n-1))
        do i = 1, n-1
            h(i) = x(i+1) - x(i)
            delta(i) = (y(i+1) - y(i)) / h(i)
        end do
        
        ! 端点处理
        interp%d(1) = pchip_endpoint(h(1), h(2), delta(1), delta(2))
        interp%d(n) = pchip_endpoint(h(n-1), h(n-2), delta(n-1), delta(n-2))
        
        ! 内部点处理
        do i = 2, n-1
            if (delta(i-1) * delta(i) <= 0.0d0) then
                interp%d(i) = 0.0d0
            else
                w1 = 2.0d0*h(i) + h(i-1)
                w2 = h(i) + 2.0d0*h(i-1)
                interp%d(i) = (w1 + w2) / (w1/delta(i-1) + w2/delta(i))
            end if
        end do
        
        deallocate(h, delta)
    end subroutine create_pchip_interpolator


end subroutine compute_filter_coefficients

end module diffusion_filter_coefficients