! 3D Heat equation: u_t = alpha*(u_xx+u_yy+u_zz)
! Domain: [0,1]^3, Dirichlet BC=0, IC: sin(pi*x)*sin(pi*y)*sin(pi*z)
! Exact: exp(-3*alpha*pi^2*t)*sin(pi*x)*sin(pi*y)*sin(pi*z)
! Method: explicit FTCS
! Output: heat3d_midplane.csv (z=0.5 slice)

program heat3d
    implicit none

    real(8), parameter :: alpha = 1.0d0
    integer, parameter :: N     = 20
    real(8), parameter :: T_END = 0.05d0
    real(8), parameter :: R     = 0.13d0

    real(8), parameter :: PI = acos(-1.0d0)

    real(8) :: dx, dt, t_final, l2, e, decay, ex
    integer :: M, steps, n, i, j, k, kmid
    real(8), allocatable :: x(:), u(:,:,:), u_new(:,:,:)

    M  = N + 2
    dx = 1.0d0 / (N + 1)
    dt = R * dx * dx / alpha
    steps = int(T_END / dt) + 1

    allocate(x(M), u(M,M,M), u_new(M,M,M))

    do i = 1, M
        x(i) = (i-1) * dx
    end do

    ! initial condition
    do k = 1, M
        do j = 1, M
            do i = 1, M
                u(i,j,k) = sin(PI*x(i))*sin(PI*x(j))*sin(PI*x(k))
            end do
        end do
    end do
    u(1,:,:)=0.0d0; u(M,:,:)=0.0d0
    u(:,1,:)=0.0d0; u(:,M,:)=0.0d0
    u(:,:,1)=0.0d0; u(:,:,M)=0.0d0

    ! time march
    do n = 1, steps
        do k = 2, M-1
            do j = 2, M-1
                do i = 2, M-1
                    u_new(i,j,k) = u(i,j,k) + R*( &
                        u(i+1,j,k)+u(i-1,j,k) &
                       +u(i,j+1,k)+u(i,j-1,k) &
                       +u(i,j,k+1)+u(i,j,k-1) &
                       -6.0d0*u(i,j,k))
                end do
            end do
        end do
        u_new(1,:,:)=0.0d0; u_new(M,:,:)=0.0d0
        u_new(:,1,:)=0.0d0; u_new(:,M,:)=0.0d0
        u_new(:,:,1)=0.0d0; u_new(:,:,M)=0.0d0
        u = u_new
    end do

    ! exact & L2
    t_final = (steps - 1) * dt
    decay   = exp(-3.0d0 * alpha * PI**2 * t_final)
    l2 = 0.0d0
    do k = 2, M-1
        do j = 2, M-1
            do i = 2, M-1
                ex = decay*sin(PI*x(i))*sin(PI*x(j))*sin(PI*x(k))
                e  = u(i,j,k) - ex
                l2 = l2 + e*e * dx**3
            end do
        end do
    end do
    l2 = sqrt(l2)
    write(*,'(A,I4,A,F7.4,A,ES10.3)') &
        '3D Heat  N=', N, '  T=', T_END, '  L2 error = ', l2

    ! save midplane z=0.5
    kmid = N/2 + 1
    open(unit=10, file='heat3d_midplane.csv', status='replace')
    write(10,'(A)') 'x,y,u_numerical,u_exact,error'
    do j = 1, M
        do i = 1, M
            ex = decay*sin(PI*x(i))*sin(PI*x(j))*sin(PI*x(kmid))
            write(10,'(5(ES14.6,","))') x(i), x(j), u(i,j,kmid), ex, u(i,j,kmid)-ex
        end do
    end do
    close(10)
    write(*,*) 'Saved heat3d_midplane.csv'

    deallocate(x, u, u_new)
end program heat3d
