
x = linspace(0, 1, 200);
a_upper = [0.1, 0.4, 0.9, 0.5]; % 4个系数（3阶Bernstein多项式）
a_lower = [-0.2, -0.4, -0.4, -0.1];

[z_upper, z_lower] = cst_airfoil(x, a_upper, a_lower);

% 绘图
figure;
plot(x, z_upper, 'r', 'LineWidth', 1.5, 'DisplayName', 'Upper Surface');
hold on;
plot(x, z_lower, 'b', 'LineWidth', 1.5, 'DisplayName', 'Lower Surface');
xlabel('Chord (x)');
ylabel('Thickness (z)');
title('CST Parameterized Airfoil');
axis equal;
legend;
grid on;


function [z_upper, z_lower] = cst_airfoil(x, a_upper, a_lower, N1, N2, delta_te)
    % CST 方法生成翼型
    % 输入：
    %     x         : 弦长坐标 (0~1 的数组)
    %     a_upper   : 上表面 Bernstein 多项式系数 [a0, a1, ..., an]
    %     a_lower   : 下表面 Bernstein 多项式系数 [b0, b1, ..., bn]
    %     N1, N2    : 类别函数指数（默认圆头尖尾）
    %     delta_te  : 后缘厚度修正（默认0，自动闭合）
    % 返回：
    %     z_upper : 上表面z坐标
    %     z_lower : 下表面z坐标

    % 默认参数
    if nargin < 4 || isempty(N1)
        N1 = 0.5;
    end
    if nargin < 5 || isempty(N2)
        N2 = 1.0;
    end
    if nargin < 6 || isempty(delta_te)
        delta_te = 0.0;
    end

    % 类别函数
    C = x.^N1 .* (1 - x).^N2;

    % 计算形状函数 S(x)（Bernstein 多项式）
    n_upper = length(a_upper) - 1;
    S_upper = zeros(size(x));
    for i = 0:n_upper
        S_upper = S_upper + a_upper(i+1) * bernstein_poly(n_upper, i, x);
    end

    n_lower = length(a_lower) - 1;
    S_lower = zeros(size(x));
    for i = 0:n_lower
        S_lower = S_lower + a_lower(i+1) * bernstein_poly(n_lower, i, x);
    end

    % 计算上下表面坐标
    z_upper = C .* S_upper + x * delta_te;
    z_lower = C .* S_lower - x * delta_te; % 下表面符号调整

    % 后缘闭合修正（可选）
    z_upper(end) = (z_upper(end) + z_lower(end)) / 2;
    z_lower(end) = z_upper(end);

end

function B = bernstein_poly(n, i, x)
    % 计算 Bernstein 多项式
    coeff = nchoosek(n, i);
    B = coeff * x.^i .* (1 - x).^(n - i);
end

% 示例：生成对称翼型（类似 NACA0012）
