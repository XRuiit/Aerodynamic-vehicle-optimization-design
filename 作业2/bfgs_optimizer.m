function [x_opt, f_opt, k, x_history] = bfgs_optimizer(fun, x0, tol, max_iter)
% BFGS 准牛顿法优化算法
% 输入：
%   fun      - 目标函数句柄
%   x0       - 初始点（列向量）
%   tol      - 收敛容忍度（默认1e-6）
%   max_iter - 最大迭代次数（默认1000）
% 输出：
%   x_opt     - 最优解
%   f_opt     - 最优目标函数值
%   k         - 实际迭代次数
%   x_history - 每一步的变量记录（列向量）

    if nargin < 4
        max_iter = 1000;
    end
    if nargin < 3
        tol = 1e-6;
    end

    x = x0(:);              % 转列向量
    n = length(x);
    A = eye(n);             % 初始Hessian逆矩阵
    x_history = x;          % 记录历史
    k = 0;

    f_val = fun(x);         % 初始函数值
    g = numerical_gradient(fun, x);
    norm_g = norm(g);

    % 主迭代循环
    while norm_g > tol && k < max_iter
        d = -A * g;

        % 黄金分割法线搜索（动态边界）
        alpha = golden_section(@(alpha) fun(x + alpha * d), 0, 1, 1e-6);

        % 更新变量
        x_new = x + alpha * d;
        g_new = numerical_gradient(fun, x_new);

        delta_x = x_new - x;
        delta_g = g_new - g;

        % BFGS 正定性判断
        if delta_g' * delta_x <= 1e-10
            warning('BFGS 更新被跳过：delta_g^T * delta_x <= 0');
            break;
        end

        % BFGS 两项更新公式
        rho = 1 / (delta_g' * delta_x);
        V = A * delta_g;
        A = A + (delta_x * delta_x') * rho - (V * V') / (delta_g' * V);

        % 更新状态
        x = x_new;
        g = g_new;
        norm_g = norm(g);
        k = k + 1;
        x_history = [x_history, x];  % 添加列

        f_val = fun(x);  % 缓存函数值
        fprintf('Iter %3d: f = %.6f, |grad| = %.6f\n', k, f_val, norm_g);
    end

    x_opt = x;
    f_opt = f_val;
end

% ================== 辅助函数 ==================== %

% 数值梯度（中心差分）
function g = numerical_gradient(f, x)
    h = 1e-6;
    n = length(x);
    g = zeros(n, 1);
    for i = 1:n
        e = zeros(n, 1);
        e(i) = 1;
        g(i) = (f(x + h * e) - f(x - h * e)) / (2 * h);
    end
end

% 黄金分割法（用于一维线搜索）
function alpha = golden_section(phi, a, b, tol)
    rho = (sqrt(5) - 1) / 2;
    c = b - rho * (b - a);
    d = a + rho * (b - a);

    for iter = 1:100
        if phi(c) < phi(d)
            b = d;
        else
            a = c;
        end
        if abs(b - a) < tol
            break;
        end
        c = b - rho * (b - a);
        d = a + rho * (b - a);
    end
    alpha = (a + b) / 2;
end
