function [x_opt, f_opt, k, x_history] = bfgs_optimizer_new(fun, x0, tol, max_iter)
    if nargin < 4
        max_iter = 1000;
    end
    if nargin < 3
        tol = 1e-6;
    end

    x = x0(:);
    n = length(x);
    A = eye(n);
    x_history = x;
    k = 0;

    f_val = fun(x);
    g = numerical_gradient(fun, x);
    norm_g = norm(g);

    while norm_g > tol && k < max_iter
        d = -A * g;

        % 检查搜索方向是否为下降方向
        if g' * d >= 0
            warning('非下降方向，终止迭代');
            break;
        end

        % Armijo 回溯线搜索
        phi = @(alpha) fun(x + alpha * d);
        alpha = armijo_search(phi, f_val, g, d, 1, 1e-4, 0.5, 20);

        x_new = x + alpha * d;
        g_new = numerical_gradient(fun, x_new);
        delta_x = x_new - x;
        delta_g = g_new - g;

        % BFGS 更新条件
        if delta_g' * delta_x > 1e-10
            rho = 1 / (delta_g' * delta_x);
            V = A * delta_g;
            A = A + (delta_x * delta_x') * rho - (V * V') / (delta_g' * V);
        else
            warning('BFGS 更新被跳过');
        end

        x = x_new;
        g = g_new;
        norm_g = norm(g);
        k = k + 1;
        x_history = [x_history, x];
        f_val = fun(x);
        fprintf('Iter %3d: f = %.6f, |grad| = %.6f\n', k, f_val, norm_g);
    end

    x_opt = x;
    f_opt = f_val;
end

function g = numerical_gradient(f, x)
    h = 1e-6;
    n = length(x);
    g = zeros(n, 1);
    f0 = f(x);  % 初始函数值缓存
    for i = 1:n
        e = zeros(n, 1);
        e(i) = 1;
        h_scaled = h * (1 + abs(x(i)));  % 动态步长
        g(i) = (f(x + h_scaled * e) - f0) / h_scaled;  % 前向差分
    end
end

function alpha = armijo_search(phi, phi0, grad, d, alpha_init, c1, rho, max_steps)
    alpha = alpha_init;
    dphi0 = grad' * d;
    for i = 1:max_steps
        phi_alpha = phi(alpha);
        if phi_alpha <= phi0 + c1 * alpha * dphi0
            break;
        else
            alpha = rho * alpha;
        end
    end
end