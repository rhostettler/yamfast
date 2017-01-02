% Simple test case for the IMMUKF
%
% A one-dimensional tracking problem that alters between a constant
% position and constant velocity motion.
%
% test.m -- 2016-04-14
% Roland Hostettler <roland.hostettler@aalto.fi>

% Housekeeping
clear variables;

addpath ..

%% Parameters
% No. of samples
K = 100;

% No. of switching points
Ns = 2;

% Sampling frequency
Ts = 0.5;

% Stationary model
F0 = [
    1, 0; 
    0, 0;
];
Q0 = [
    0.1, 0;
      0, 0;
];

% Const. velocity model
F1 = [
    1, Ts;
    0,  1;
];
Q1 = [
    Ts, 0;
     0, 1;
];

% Observation model
G = [0.5, 0];
R = 0.5;

% Initial states
m00 = 0;
P00 = 2;
m01 = zeros(2, 1);
P01 = diag([2, 1]);

%% Model & Filter Definition
imm = NonlinearIMM();
model0 = LGSSModel(1, G(1), Q0(1, 1), R, m00, P00);
imm.addModel(model0, 1);
model1 = LGSSModel(F1, G, Q1, R, m01, P01);
imm.addModel(model1, (1:2).');

% Transition probabilities
imm.p_ij = [
    0.9, 0.1;
    0.1, 0.9;
];
    
% Initial model probabilities
imm.mu_i0 = 0.5*ones(2, 1); 

% Initialize the filter
filter = IMMUKF(imm);

%% Initialize the System
% Determine the switching points
s = randi(K, [Ns, 1]);

% Initial mode
mode = round(rand(1));

% Initial state
if mode == 0
    x = m00 + sqrt(P00)*randn(1, 1);
    x = [x; 0];
else
    x = m01 + sqrt(P01)*randn(2, 1);
end

%% Simulate
x_save = zeros(2, K);
mode_save = zeros(1, K);
y_save = zeros(1, K);
m_save = zeros(2, K);
mu_save = zeros(2, K);

for k = 1:K
    %% System Simulation
    % Determine if we have to switch and do so
    if ~isempty(find(k == s, 1))
        mode = 1-mode;
    end
    
    % Propagate the state
    if mode == 0
        % Stationary
        x = F0*x + sqrt(Q0)*randn(2, 1);
    else
        % Moving w/ constant velocity
        x = F1*x + sqrt(Q1)*randn(2, 1);
    end
    
    % Measurement
    y = G*x + sqrt(R)*randn(1);        
    
    %% Update the filter
    filter.update(0, y, k);
    
    %% Store the values
    y_save(:, k) = y;
    x_save(:, k) = x;
    mode_save(:, k) = mode;
    m_save(:, k) = filter.m;
    mu_save(:, k) = filter.mu_i;
end

%% Visualize the result
figure(1); clf();
subplot(211);
plot([x_save(1, :).', m_save(1, :).']);
legend('True', 'Estimated');
subplot(212);
plot([x_save(2, :).', m_save(2, :).']);
legend('True', 'Estimated');

figure(2); clf();
plot([mu_save(1, :).', mu_save(2, :).', mode_save.']);
legend('Pr(M0)', 'Pr(M1)', 'True');
