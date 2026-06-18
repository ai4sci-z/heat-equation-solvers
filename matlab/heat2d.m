% 2D Heat equation: u_t = alpha*(u_xx + u_yy)
% Domain: [0,1]^2, Dirichlet BC, IC: sin(pi*x)*sin(pi*y)
% Exact: exp(-2*alpha*pi^2*t)*sin(pi*x)*sin(pi*y)
% Method: explicit FTCS

clear; close all;

alpha = 1.0;
N     = 50;
T     = 0.1;
x     = linspace(0, 1, N+2);
dx    = x(2) - x(1);
r     = 0.2;
dt    = r * dx^2 / alpha;
steps = floor(T / dt) + 1;

[X, Y] = meshgrid(x, x);   % X(i,j)=x_j, Y(i,j)=y_i (MATLAB row=y)

u = sin(pi*X) .* sin(pi*Y);
u(1,:)=0; u(end,:)=0; u(:,1)=0; u(:,end)=0;

exact2d = @(t) exp(-2*alpha*pi^2*t) .* sin(pi*X) .* sin(pi*Y);

snap_steps = [1, floor(steps/2), steps];
snaps = containers.Map('KeyType','double','ValueType','any');

for n = 1:steps
    if ismember(n, snap_steps)
        snaps(n*dt) = u;
    end
    lap = (u(3:end,2:end-1) - 2*u(2:end-1,2:end-1) + u(1:end-2,2:end-1)) ...
        + (u(2:end-1,3:end) - 2*u(2:end-1,2:end-1) + u(2:end-1,1:end-2));
    lap = lap / dx^2;
    u_new = u;
    u_new(2:end-1,2:end-1) = u(2:end-1,2:end-1) + alpha*dt*lap;
    u_new(1,:)=0; u_new(end,:)=0; u_new(:,1)=0; u_new(:,end)=0;
    u = u_new;
end

t_final = (steps-1)*dt;
l2_err  = sqrt(dx^2 * sum(sum((u - exact2d(t_final)).^2)));
fprintf('2D Heat  N=%d  T=%.2f  L2 error = %.3e\n', N, T, l2_err);

% ── plot ──────────────────────────────────────────────────────────────────────
if ~exist('../results', 'dir'), mkdir('../results'); end
t_keys = sort(cell2mat(keys(snaps)));

figure('Position', [100 100 1200 700]);
for k = 1:length(t_keys)
    t_val = t_keys(k);
    u_s   = snaps(t_val);
    err   = u_s - exact2d(t_val);

    subplot(2, 3, k);
    contourf(X, Y, u_s, 20, 'LineColor', 'none');
    colorbar; title(sprintf('Numerical t=%.3f', t_val));
    xlabel('x'); ylabel('y');

    subplot(2, 3, k+3);
    contourf(X, Y, err, 20, 'LineColor', 'none');
    colorbar; colormap(subplot(2,3,k+3), 'RdBu');
    title(sprintf('Error t=%.3f', t_val));
    xlabel('x'); ylabel('y');
end
sgtitle(sprintf('MATLAB — 2D Heat Equation (FTCS)  N=%d  L2=%.2e', N, l2_err));
saveas(gcf, '../results/matlab_heat2d.png');
fprintf('Saved ../results/matlab_heat2d.png\n');
