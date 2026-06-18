% 1D Heat equation: u_t = alpha*u_xx
% Domain: [0,1], BC: u(0)=u(1)=0, IC: sin(pi*x)
% Exact: exp(-alpha*pi^2*t)*sin(pi*x)
% Method: explicit FTCS

clear; close all;

alpha = 1.0;
N     = 50;
T     = 0.1;
x     = linspace(0, 1, N+2)';
dx    = x(2) - x(1);
r     = 0.4;
dt    = r * dx^2 / alpha;
steps = floor(T / dt) + 1;

u = sin(pi * x);
u(1) = 0; u(end) = 0;

snap_steps = [1, floor(steps/2), steps];
snaps = containers.Map('KeyType','double','ValueType','any');

for n = 1:steps
    if ismember(n, snap_steps)
        snaps(n * dt) = u;
    end
    u_new = u;
    u_new(2:end-1) = u(2:end-1) + r*(u(3:end) - 2*u(2:end-1) + u(1:end-2));
    u_new(1) = 0; u_new(end) = 0;
    u = u_new;
end

exact = @(t) exp(-alpha*pi^2*t) .* sin(pi*x);
t_final = (steps-1)*dt;
l2_err = sqrt(dx * sum((u - exact(t_final)).^2));
fprintf('1D Heat  N=%d  T=%.2f  L2 error = %.3e\n', N, T, l2_err);

% ── plot ──────────────────────────────────────────────────────────────────────
if ~exist('../results', 'dir'), mkdir('../results'); end

figure('Position', [100 100 1000 400]);
subplot(1,2,1); hold on;
colors = lines(length(snap_steps));
t_keys = sort(cell2mat(keys(snaps)));
for k = 1:length(t_keys)
    t_val = t_keys(k);
    plot(x, snaps(t_val), '-',  'Color', colors(k,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('num t=%.3f', t_val));
    plot(x, exact(t_val), '--', 'Color', colors(k,:), 'LineWidth', 1);
end
xlabel('x'); ylabel('u'); title('1D Heat (solid=num, dashed=exact)');
legend('show', 'Location', 'best', 'FontSize', 7); grid on;

subplot(1,2,2);
plot(x, u - exact(t_final), 'r-', 'LineWidth', 1.5);
xlabel('x'); ylabel('error');
title(sprintf('Error at t=%.3f  L2=%.2e', t_final, l2_err)); grid on;

sgtitle(sprintf('MATLAB — 1D Heat Equation (FTCS)  N=%d', N));
saveas(gcf, '../results/matlab_heat1d.png');
fprintf('Saved ../results/matlab_heat1d.png\n');
