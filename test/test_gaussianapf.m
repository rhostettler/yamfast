% Gaussian APF Test
%
% test_gaussianapf.m -- 2017-02-20
% Roland Hostettler <roland.hostettler@aalto.fi>

% Housekeeping
clear variables;

% Add library
addpath ../src

%% 
K = 100;   % No. of MC simulations
M = 250;  % No. of particles
N = 100;  % No. of time samples
Ts = 1; % Sampling time

t = (1:N)*Ts;

%% 
xs = zeros(1, N, K);
xhat_bpf = xs;
xhat_epf = xs;
xhat_upf = xs;
xhat_plupf = xs;

y = zeros(1, N, K);

u = zeros(1, N);
h = pbar(K);
parfor k = 1:K
% for k = 1:K
    %% Initialize Filter & Model
    model = UNGModel([], 0.1);
    bpf = BootstrapFilter(model, M);
    epf = GaussianAPF(model, M); % EKF approximation
    upf = GaussianAPF(model, M, UnscentedTransform());
    plupf = GaussianAPF(model, M, UnscentedTransform());
    plupf.J = 1;

    %%
    x = model.px0_rand(1);
    for n = 1:N
        %% Simulation
        Q = model.Q(x, t(n), u(n));
        q = chol(Q).'*randn(size(Q, 1), 1);
        x = model.f(x, q, t(n), u(n));
        R = model.R(x, t(n), u(n));
        r = chol(R).'*randn(size(R, 1), 1);
        y(:, n, k) = model.g(x, r, t(n), u(n));
        xs(:, n, k) = x;

        %% Filter
        xhat_bpf(:, n, k) = bpf.update(y(:, n, k), t(n), u(n));
        xhat_epf(:, n, k) = epf.update(y(:, n, k), t(n), u(n));
        xhat_upf(:, n, k) = upf.update(y(:, n, k), t(n), u(n));
        xhat_plupf(:, n, k) = plupf.update(y(:, n, k), t(n), u(n));

        %% Update
        model.t = t(n);
    end
    
    %% Progress Update
    pbar(k, h);
end
pbar(0, h);

%% 
e_bpf = xs-xhat_bpf;
e_epf = xs-xhat_epf;
e_upf = xs-xhat_upf;
e_plupf = xs-xhat_plupf;

trms = @(e) mean(sqrt(sum(e.^2, 1)), 2);
mean(trms(e_bpf))
mean(trms(e_epf))
mean(trms(e_upf))
mean(trms(e_plupf))

%%
figure(1); clf();
plot(t, xs(:, :, 1)); hold on;
plot(t, xhat_bpf(:, :, 1));
plot(t, xhat_epf(:, :, 1));
plot(t, xhat_upf(:, :, 1));
plot(t, xhat_plupf(:, :, 1));
legend('State', 'Bootstrap', 'E-APF', 'U-APF', 'PL-APF');

figure(2); clf();
plot(t, y(:, :, 1));


