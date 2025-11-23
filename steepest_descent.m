%最速下降法
function [x_opt, f_opt, k, x_history] = steepest_descent(fun, x0, tol, max_iter)
% 最速下降法优化算法（数值梯度 + 黄金分割法）
% 输入：
%   fun      - 目标函数（函数句柄）
%   x0       - 初始点（列向量）
%   tol      - 容忍度（默认1e-6）
%   max_iter - 最大迭代次数（默认1000）
% 输出：
%   x_opt    - 最优解
%   f_opt    - 最优目标函数值
%   k        - 实际迭代次数
%   x_history- 每次迭代的点（列向量集合）

% 处理输入参数
if nargin < 4
    max_iter = 1000;
end
if nargin < 3
    tol = 1e-6;
end

% 初始化
x = x0(:); % 确保为列向量
k = 0;
g = numerical_gradient(fun, x);
norm_g = norm(g);
x_history = x; % 存储迭代历史

% 主循环
while (norm_g > tol) && (k < max_iter)
    d = -g; % 最速下降方向
    
    % 定义方向上的目标函数
    phi = @(alpha) fun(x + alpha*d);
    
    % 寻找包含极小值的区间
    [a, b] = bracket_minimum(phi);
    
    % 黄金分割法线搜索
    alpha = golden_section_search(phi, a, b, 1e-6);
    
    % 更新点
    x = x + alpha*d;
    g = numerical_gradient(fun, x);
    norm_g = norm(g);
    k = k + 1;
    x_history = [x_history, x]; % 记录历史点
end

x_opt = x;
f_opt = fun(x);

% ======== 嵌套辅助函数 ========
function g = numerical_gradient(f, x)
    % 中心差分计算数值梯度
    h = 1e-6;
    n = length(x);
    g = zeros(n,1);
    for i = 1:n
        e = zeros(n,1);
        e(i) = 1;
        g(i) = (f(x + h*e) - f(x - h*e))/(2*h);
    end
end

function [a, b] = bracket_minimum(phi)
    % 扩展区间找到包含极小点的[a, b]
    a = 0;
    step = 0.1;
    b = a + step;
    fa = phi(a);
    fb = phi(b);
    
    % 寻找函数值上升的边界
    while fb < fa
        a = b;
        fa = fb;
        step = 2*step;
        b = a + step;
        fb = phi(b);
    end
    
    % 回退一步确保极小点在区间内
    if a > 0
        a = a - step/2;
        b = b - step/2;
    end
end

function alpha = golden_section_search(phi, a, b, tol)
    % 黄金分割法一维搜索
    rho = (sqrt(5)-1)/2; % 0.618
    c = b - rho*(b - a);
    d = a + rho*(b - a);
    
    for iter = 1:100
        if phi(c) < phi(d)
            b = d;
        else
            a = c;
        end
        if abs(b - a) < tol
            break;
        end
        % 更新中间点
        c = b - rho*(b - a);
        d = a + rho*(b - a);
    end
    alpha = (a + b)/2;
end

end