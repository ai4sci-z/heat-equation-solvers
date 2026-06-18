! 1D Heat equation: u_t = alpha*u_xx
! Domain: [0,1], BC: u=0, IC: sin(pi*x)
! Exact: exp(-alpha*pi^2*t)*sin(pi*x)
! Method: explicit FTCS
! Output: heat1d_result.csv

program heat1d
    implicit none

    real(8), parameter :: alpha = 1.0d0
    integer, parameter :: N     = 50
    real(8), parameter :: T_END = 0.1d0
    real(8), parameter :: R     = 0.4d0   ! CFL: alpha*dt/dx^2

    real(8) :: dx, dt, t_final, l2, e
    integer :: M, steps, n, i
    real(8), allocatable :: x(:), u(:), u_new(:), u_ex(:)

    real(8), parameter :: PI = acos(-1.0d0)

    M  = N + 2
    dx = 1.0d0 / (N + 1)
    dt = R * dx * dx / alpha
    steps = int(T_END / dt) + 1

    allocate(x(M), u(M), u_new(M), u_ex(M))

    do i = 1, M
        x(i) = (i-1) * dx
        u(i) = sin(PI * x(i))
    end do
    u(1) = 0.0d0; u(M) = 0.0d0

    ! time march
    do n = 1, steps
        do i = 2, M-1
            u_new(i) = u(i) + R*(u(i+1) - 2.0d0*u(i) + u(i-1))
        end do
        u_new(1) = 0.0d0; u_new(M) = 0.0d0
        u = u_new
    end do

    ! exact and L2 error
    t_final = (steps - 1) * dt
    l2 = 0.0d0
    do i = 1, M
        u_ex(i) = exp(-alpha * PI**2 * t_final) * sin(PI * x(i))
        e = u(i) - u_ex(i)
        l2 = l2 + e*e * dx
    end do
    l2 = sqrt(l2)
    write(*,'(A,I4,A,F6.3,A,ES10.3)') &
        '1D Heat  N=', N, '  T=', T_END, '  L2 error = ', l2

    ! save CSV
    open(unit=10, file='heat1d_result.csv', status='replace')
    write(10,'(A)') 'x,u_numerical,u_exact,error'
    do i = 1, M
        write(10,'(4(ES14.6,","))') x(i), u(i), u_ex(i), u(i)-u_ex(i)
    end do
    close(10)
    write(*,*) 'Saved heat1d_result.csv'

    deallocate(x, u, u_new, u_ex)
end program heat1d
