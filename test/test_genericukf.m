% Simple test case for the GenericUKF
%
% test.m -- 2016-03-02
% Roland Hostettler <roland.hostettler@aalto.fi>

% Housekeeping
clear variables;

addpath ..

%% Parameters
% No. of samples
K = 100;

% Model parameters
F = 1;
G = 0.5;
Q = 0.1;
R = 0.5;
m0 = 0;
P0 = 2;

%% Initialize
x = m0;
P = P0;
model = LGSSModel(F, G, Q, R, m0, P0);
filter = GenericUKF(model);

%% Simulate
x_save = zeros(1, K);
y_save = zeros(1, K);
m_save = zeros(1, K);

for k = 1:K
    %% Generate data
    x = F*x + sqrt(Q)*randn(1);
    y = G*x + sqrt(R)*randn(1);        
    
    %% Update the filter
    filter.update(0, y, k);
    
    %% Store the values
    y_save(:, k) = y;
    x_save(:, k) = x;
    m_save(:, k) = filter.m;
end

%% Visualize the result
figure(1); clf();
plot(x_save); hold on;
plot(y_save, 'x');
plot(m_save);
legend('x', 'y', 'm');
fprintf('RMSE w/ inversion: %.2f\n', rms(x_save-2*y_save));
fprintf('RMSE w/ UKF: %.2f\n', rms(x_save-m_save));
