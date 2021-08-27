% RBCKF Model Test
%
% <long description>
%
% Stochasitc oscillator with random walk frequency
%
% 
% 
% xn = 
%
% test_clgssmodel.m -- 2016-08-25
% Roland Hostettler <roland.hostettler@aalto.fi>

% Housekeeping
clear variables;

% Add library
addpath ../src

%% Parameters
% Sampling time
Ts = 0.05;

% No. of samples
T = 200;

% Process noise
var_q = 0.2^2;
B = [0, 0, 1]';
Q = B*var_q*B';

% Measurement noise
R = 0.5^2;

% Initial state
m0 = [0, 1, 2*pi*1].';
P0 = eye(3).^2;

%% System Description
il = (1:2).';
in = 3;
fn = @(xn, u, t) xn;
An = @(xn, u, t) zeros(1, 2);
fl = @(xn, u, t) zeros(size(xn));
Al = @(xn, u, t) [
    cos(xn*Ts), -sin(xn*Ts);
    sin(xn*Ts),  cos(xn*Ts);
];
h = @(xn, u, t) 0;
C = @(xn, u, t) [1, 0];
Qf = @(xn, u, t) Q;
Rf = @(xn, u, t) R;
model = GenericCLGSSModel(fn, An, fl, Al, Qf, h, C, Rf, m0, P0, in, il);

%% Filter
filter = RBGF(model);

%% Simulate
xhat = zeros(3, T);
x = zeros(3, T);
y = zeros(1, T);
xt = m0 + chol(P0)'*randn(3, 1);

for t = 1:T
    %% System Simulation
    % Propagate the state
    qt = B*sqrt(var_q)*randn(1);
    xt = model.f(xt, [], qt, t);
    x(:, t) = xt;
    
    % Measurement
    rt = sqrt(R)*randn(1);
    y(:, t) = model.g(xt, [], rt, t);
    
    %% Filter Update
    filter.update([], y(:, t), t);
    xhat(:, t) = filter.m;
end


%% Visualize
figure(1); clf();
for i = 1:3
    subplot(3, 1, i);
    plot(x(i, :)); hold on;
    plot(xhat(i, :));
end

figure(2); clf();
plot(y);
