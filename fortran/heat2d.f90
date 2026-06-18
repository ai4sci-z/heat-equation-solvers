! 2D Heat equation: u_t = alpha*(u_xx+u_yy)
! Domain: [0,1]^2, Dirichlet BC=0, IC: sin(pi*x)*sin(pi*y)
! Exact: exp(-2*alpha*pi^2*t)*sin(pi*x)*sin(pi*y)
! Method: explicit FTCS
! Output: heat2d_result.csv

program heat2d
    implicit none

    real(8), parameter :: alpha = 1.0d0
    integer, parameter :: N     = 50
    real(8), parameter :: T_END = 0.1d0
    real(8), parameter :: R     = 0.2d0

    real(8), parameter :: PI = acos(-1.0d0)

    real(8) :: dx, dt, t_final, l2, e, decay, ex
    integer :: M, steps, n, i, j
    real(8), allocatable :: x(:), u(:,:), u_new(:,:)

    M  = N + 2
    dx = 1.0d0 / (N + 1)
    dt = R * dx * dx / alpha
    steps = int(T_END / dt) + 1

    allocate(x(M), u(M,M), u_new(M,M))

    do i = 1, M
        x(i) = (i-1) * dx
    end do

    ! initial condition
    do j = 1, M
        do i = 1, M
            u(i,j) = sin(PI*x(i)) * sin(PI*x(j))
        end do
    end do
    u(1,:) = 0.0d0; u(M,:) = 0.0d0
    u(:,1) = 0.0d0; u(:,M) = 0.0d0

    ! time march
    do n = 1, steps
        do j = 2, M-1
            do i = 2, M-1
                u_new(i,j) = u(i,j) + R*(u(i+1,j)+u(i-1,j) &
                                        +u(i,j+1)+u(i,j-1)-4.0d0*u(i,j))
            end do
        end do
        u_new(1,:) = 0.0d0; u_new(M,:) = 0.0d0
        u_new(:,1) = 0.0d0; u_new(:,M) = 0.0d0
        u = u_new
    end do

    ! exact & L2
    t_final = (steps - 1) * dt
    decay   = exp(-2.0d0 * alpha * PI**2 * t_final)
    l2 = 0.0d0
    do j = 2, M-1
        do i = 2, M-1
            ex = decay * sin(PI*x(i)) * sin(PI*x(j))
            e  = u(i,j) - ex
            l2 = l2 + e*e * dx*dx
        end do
    end do
    l2 = sqrt(l2)
    write(*,'(A,I4,A,F6.3,A,ES10.3)') &
        '2D Heat  N=', N, '  T=', T_END, '  L2 error = ', l2

    ! save CSV
    open(unit=10, file='heat2d_result.csv', status='replace')
    write(10,'(A)') 'x,y,u_numerical,u_exact,error'
    do j = 1, M
        do i = 1, M
            ex = decay * sin(PI*x(i)) * sin(PI*x(j))
            write(10,'(5(ES14.6,","))') x(i), x(j), u(i,j), ex, u(i,j)-ex
        end do
    end do
    close(10)
    write(*,*) 'Saved heat2d_result.csv'

    deallocate(x, u, u_new)
end program heat2d
