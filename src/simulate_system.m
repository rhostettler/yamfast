function [x, y] = simulate_system(model, t, u)

N = length(t);
xt = model.px0_rand(1);
x = zeros(size(xt, 1), N);

for n = 1:N
    xt = model.px_rand(xt, t(n), u(:, n));
    
    % TODO: Replace with model.py_rand();
    R = model.R(xt, t(n), u(:, n));
    r = chol(R).'*randn(size(R, 1), 1);
    y(:, n) = model.g(xt, r, t(n), u(:, n));
    x(:, n) = xt;
end
