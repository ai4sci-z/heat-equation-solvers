/*
 * 1D Heat equation: u_t = alpha * u_xx
 * Domain: [0,1], BC: u=0, IC: sin(pi*x)
 * Exact: exp(-alpha*pi^2*t)*sin(pi*x)
 * Method: explicit FTCS
 * Output: heat1d_result.csv (columns: x, u_numerical, u_exact, error)
 */

#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <fstream>

static const double PI    = std::acos(-1.0);
static const double ALPHA = 1.0;
static const int    N     = 50;
static const double T_END = 0.1;
static const double R     = 0.4;     // CFL: alpha*dt/dx^2

void save_csv(const std::string& fname,
              const std::vector<double>& x,
              const std::vector<double>& u_num,
              const std::vector<double>& u_ex)
{
    std::ofstream f(fname);
    f << "x,u_numerical,u_exact,error\n";
    for (int i = 0; i <= N+1; ++i)
        f << x[i] << "," << u_num[i] << "," << u_ex[i]
          << "," << (u_num[i]-u_ex[i]) << "\n";
}

int main()
{
    int M = N + 2;
    std::vector<double> x(M), u(M), u_new(M);

    double dx = 1.0 / (N + 1);
    for (int i = 0; i < M; ++i) x[i] = i * dx;

    double dt    = R * dx * dx / ALPHA;
    int    steps = static_cast<int>(T_END / dt) + 1;

    // initial condition
    for (int i = 0; i < M; ++i)
        u[i] = std::sin(PI * x[i]);
    u[0] = u[M-1] = 0.0;

    // time march
    for (int n = 0; n < steps; ++n) {
        for (int i = 1; i <= N; ++i)
            u_new[i] = u[i] + R * (u[i+1] - 2.0*u[i] + u[i-1]);
        u_new[0] = u_new[M-1] = 0.0;
        std::swap(u, u_new);
    }

    // exact solution and L2 error
    double t_final = (steps - 1) * dt;
    std::vector<double> u_ex(M);
    double l2 = 0.0;
    for (int i = 0; i < M; ++i) {
        u_ex[i] = std::exp(-ALPHA * PI * PI * t_final) * std::sin(PI * x[i]);
        double e = u[i] - u_ex[i];
        l2 += e * e * dx;
    }
    l2 = std::sqrt(l2);
    std::printf("1D Heat  N=%d  T=%.2f  L2 error = %.3e\n", N, T_END, l2);

    save_csv("heat1d_result.csv", x, u, u_ex);
    std::printf("Saved heat1d_result.csv\n");
    return 0;
}
