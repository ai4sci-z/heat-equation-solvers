/*
 * 3D Heat equation: u_t = alpha*(u_xx+u_yy+u_zz)
 * Domain: [0,1]^3, Dirichlet BC=0, IC: sin(pi*x)*sin(pi*y)*sin(pi*z)
 * Exact: exp(-3*alpha*pi^2*t)*sin(pi*x)*sin(pi*y)*sin(pi*z)
 * Method: explicit FTCS
 * Output: heat3d_midplane.csv (z=0.5 slice)
 */

#include <cmath>
#include <cstdio>
#include <vector>
#include <fstream>

static const double PI    = std::acos(-1.0);
static const double ALPHA = 1.0;
static const int    N     = 20;
static const double T_END = 0.05;
static const double R     = 0.13;   // alpha*dt/dx^2 < 1/6

inline int idx(int i, int j, int k) {
    return i*(N+2)*(N+2) + j*(N+2) + k;
}

int main()
{
    int    M  = N + 2;
    double dx = 1.0 / (N + 1);
    double dt = R * dx * dx / ALPHA;
    int    steps = static_cast<int>(T_END / dt) + 1;

    std::vector<double> x(M);
    for (int i = 0; i < M; ++i) x[i] = i * dx;

    long long total = (long long)M * M * M;
    std::vector<double> u(total, 0.0), u_new(total, 0.0);

    // initial condition
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            for (int k = 0; k < M; ++k)
                u[idx(i,j,k)] = std::sin(PI*x[i])*std::sin(PI*x[j])*std::sin(PI*x[k]);

    // time march
    for (int n = 0; n < steps; ++n) {
        for (int i = 1; i <= N; ++i)
            for (int j = 1; j <= N; ++j)
                for (int k = 1; k <= N; ++k) {
                    double lap = (u[idx(i+1,j,k)]+u[idx(i-1,j,k)]
                                 +u[idx(i,j+1,k)]+u[idx(i,j-1,k)]
                                 +u[idx(i,j,k+1)]+u[idx(i,j,k-1)]
                                 -6.0*u[idx(i,j,k)]) / (dx*dx);
                    u_new[idx(i,j,k)] = u[idx(i,j,k)] + ALPHA*dt*lap;
                }
        for (int i = 1; i <= N; ++i)
            for (int j = 1; j <= N; ++j)
                for (int k = 1; k <= N; ++k)
                    u[idx(i,j,k)] = u_new[idx(i,j,k)];
    }

    // exact & error
    double t_final = (steps-1)*dt;
    double decay   = std::exp(-3.0*ALPHA*PI*PI*t_final);
    double l2 = 0.0;
    for (int i = 1; i <= N; ++i)
        for (int j = 1; j <= N; ++j)
            for (int k = 1; k <= N; ++k) {
                double ex = decay*std::sin(PI*x[i])*std::sin(PI*x[j])*std::sin(PI*x[k]);
                double e  = u[idx(i,j,k)] - ex;
                l2 += e*e*dx*dx*dx;
            }
    l2 = std::sqrt(l2);
    std::printf("3D Heat  N=%d  T=%.3f  L2 error = %.3e\n", N, T_END, l2);

    // save midplane z=0.5
    int kk = N/2 + 1;
    std::ofstream f("heat3d_midplane.csv");
    f << "x,y,u_numerical,u_exact,error\n";
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) {
            double ex = decay*std::sin(PI*x[i])*std::sin(PI*x[j])*std::sin(PI*x[kk]);
            f << x[i] << "," << x[j] << ","
              << u[idx(i,j,kk)] << "," << ex << ","
              << (u[idx(i,j,kk)]-ex) << "\n";
        }
    std::printf("Saved heat3d_midplane.csv\n");
    return 0;
}
