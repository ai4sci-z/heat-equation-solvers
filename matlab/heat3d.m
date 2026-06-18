% 3D Heat equation: u_t = alpha*(u_xx+u_yy+u_zz)
% Domain: [0,1]^3, Dirichlet BC, IC: sin(pi*x)*sin(pi*y)*sin(pi*z)
% Exact: exp(-3*alpha*pi^2*t)*sin(pi*x)*sin(pi*y)*sin(pi*z)
% Method: explicit FTCS

clear; close all;

alpha = 1.0;
N     = 20;
T     = 0.05;
x     = linspace(0, 1, N+2);
dx    = x(2) - x(1);
r     = 0.13;
dt    = r * dx^2 / alpha;
steps = floor(T / dt) + 1;

[X, Y, Z] = meshgrid(x, x, x);

u = sin(pi*X) .* sin(pi*Y) .* sin(pi*Z);
u(1,:,:)=0; u(end,:,:)=0;
u(:,1,:)=0; u(:,end,:)=0;
u(:,:,1)=0; u(:,:,end)=0;

exact3d = @(t) exp(-3*alpha*pi^2*t) .* sin(pi*X) .* sin(pi*Y) .* sin(pi*Z);

snap_steps = [1, floor(steps/2), steps];
snaps = containers.Map('KeyType','double','ValueType','any');

for n = 1:steps
    if ismember(n, snap_steps)
        snaps(n*dt) = u;
    end
    i = 2:N+1; im = 1:N; ip = 3:N+2;
    lap = ( u(ip,i,i) - 2*u(i,i,i) + u(im,i,i) ...
          + u(i,ip,i) - 2*u(i,i,i) + u(i,im,i) ...
          + u(i,i,ip) - 2*u(i,i,i) + u(i,i,im) ) / dx^2;
    u_new = u;
    u_new(i,i,i) = u(i,i,i) + alpha*dt*lap;
    u_new(1,:,:)=0; u_new(end,:,:)=0;
    u_new(:,1,:)=0; u_new(:,end,:)=0;
    u_new(:,:,1)=0; u_new(:,:,end)=0;
    u = u_new;
end

t_final = (steps-1)*dt;
l2_err  = sqrt(dx^3 * sum(sum(sum((u - exact3d(t_final)).^2))));
fprintf('3D Heat  N=%d  T=%.3f  L2 error = %.3e\n', N, T, l2_err);

% ── plot: midplane z=0.5 ──────────────────────────────────────────────────────
if ~exist('../results', 'dir'), mkdir('../results'); end
mid = N/2 + 1;
t_keys = sort(cell2mat(keys(snaps)));

figure('Position', [100 100 1200 700]);
for k = 1:length(t_keys)
    t_val  = t_keys(k);
    u_s    = snaps(t_val);
    u_ex_s = exact3d(t_val);
    slice_num = squeeze(u_s(:,:,mid));
    slice_err = squeeze(u_s(:,:,mid) - u_ex_s(:,:,mid));
    Xs = squeeze(X(:,:,mid)); Ys = squeeze(Y(:,:,mid));

    subplot(2, 3, k);
    contourf(Xs, Ys, slice_num, 20, 'LineColor', 'none');
    colorbar; title(sprintf('Numerical z=0.5  t=%.3f', t_val));
    xlabel('x'); ylabel('y');

    subplot(2, 3, k+3);
    contourf(Xs, Ys, slice_err, 20, 'LineColor', 'none');
    colorbar; colormap(subplot(2,3,k+3), 'RdBu');
    title(sprintf('Error z=0.5  t=%.3f', t_val));
    xlabel('x'); ylabel('y');
end
sgtitle(sprintf('MATLAB — 3D Heat Equation (FTCS)  N=%d  L2=%.2e', N, l2_err));
saveas(gcf, '../results/matlab_heat3d.png');
fprintf('Saved ../results/matlab_heat3d.png\n');
