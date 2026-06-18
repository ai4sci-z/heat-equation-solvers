/*
 * 2D Heat equation: u_t = alpha*(u_xx + u_yy)
 * Domain: [0,1]^2, Dirichlet BC=0, IC: sin(pi*x)*sin(pi*y)
 * Exact: exp(-2*alpha*pi^2*t)*sin(pi*x)*sin(pi*y)
 * Method: explicit FTCS
 * Output: heat2d_result.csv (x, y, u_numerical, u_exact, error)
 */

#include <cmath>
#include <cstdio>
#include <vector>
#include <fstream>
#include <string>

static const double PI    = std::acos(-1.0);
static const double ALPHA = 1.0;
static const int    N     = 50;
static const double T_END = 0.1;
static const double R     = 0.2;    // alpha*dt/dx^2

int idx(int i, int j) { return i * (N+2) + j; }

int main()
{
    int    M  = N + 2;
    double dx = 1.0 / (N + 1);
    double dt = R * dx * dx / ALPHA;
    int    steps = static_cast<int>(T_END / dt) + 1;

    std::vector<double> x(M);
    for (int i = 0; i < M; ++i) x[i] = i * dx;

    int total = M * M;
    std::vector<double> u(total, 0.0), u_new(total, 0.0);

    // initial condition
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            u[idx(i,j)] = std::sin(PI*x[i]) * std::sin(PI*x[j]);
    // zero boundary (already 0 at i=0,M-1,j=0,M-1)

    // time march
    for (int n = 0; n < steps; ++n) {
        for (int i = 1; i <= N; ++i)
            for (int j = 1; j <= N; ++j) {
                double lap = (u[idx(i+1,j)] + u[idx(i-1,j)]
                            + u[idx(i,j+1)] + u[idx(i,j-1)]
                            - 4.0*u[idx(i,j)]) / (dx*dx);
                u_new[idx(i,j)] = u[idx(i,j)] + ALPHA*dt*lap;
            }
        // copy interior, keep boundary=0
        for (int i = 1; i <= N; ++i)
            for (int j = 1; j <= N; ++j)
                u[idx(i,j)] = u_new[idx(i,j)];
    }

    // exact & error
    double t_final = (steps-1)*dt;
    double decay   = std::exp(-2.0*ALPHA*PI*PI*t_final);
    double l2 = 0.0;
    for (int i = 1; i <= N; ++i)
        for (int j = 1; j <= N; ++j) {
            double ex = decay * std::sin(PI*x[i]) * std::sin(PI*x[j]);
            double e  = u[idx(i,j)] - ex;
            l2 += e*e * dx*dx;
        }
    l2 = std::sqrt(l2);
    std::printf("2D Heat  N=%d  T=%.2f  L2 error = %.3e\n", N, T_END, l2);

    std::ofstream f("heat2d_result.csv");
    f << "x,y,u_numerical,u_exact,error\n";
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) {
            double ex = decay * std::sin(PI*x[i]) * std::sin(PI*x[j]);
            f << x[i] << "," << x[j] << ","
              << u[idx(i,j)] << "," << ex << ","
              << (u[idx(i,j)]-ex) << "\n";
        }
    std::printf("Saved heat2d_result.csv\n");
    return 0;
}
