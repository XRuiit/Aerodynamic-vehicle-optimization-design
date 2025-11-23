%反优化
clc;clear;
close all;
clear functions ; % 重置 persistent 变量 iter


%目标翼型
% 读取翼型数据
filename = 'target_airfoil.dat'; % 替换为你的文件名
% filename = '0012.dat';
data = dlmread(filename); % 读取数据
x_t = data(:, 1); % 提取 x 坐标
y_t = data(:, 2); % 提取 y 坐标
%目标翼型x数据
x_low = flip(x_t(1:59));
x_up = x_t(59:end);

%%CST参数化 初始翼型
N1 = 0.5;
N2 = 1;
a_up = [0.1,0.2,0.3,0.3,0.3,0.2]; % 系数
a_low = -[0.1,0.2,0.3,0.3,0.3,0.2];
% a_up = [0.12, 0.14, 0.19, 0.16,0.2,0.2]; % 系数
% a_low = [-0.12, -0.16, -0.19, -0.1,-0.15,0.1];
a0 = [a_up,a_low]';
%上下翼面x不同，故分开计算
[y_up,~] = cst_airfoil(x_up, a_up, a_low, N1, N2);
[~,y_low] = cst_airfoil(x_low, a_up, a_low, N1, N2);
yt_low = flip(y_t(1:59));%下翼面反序，从x=0开始
yt_up = y_t(59:end); 

target.y = y_t;
target.xu = x_up;
target.xl = x_low;
cst.yu = y_up;
cst.yl = y_low;
cst.a = a0;

f = @(a) Target(a,target);

%% 最速下降
% [a_opt, f_opt, k, x_history] = steepest_descent(f, a0, 1e-6, 50);
% % 显示结果
% fprintf('最优解：(%f, %f)\n', a_opt(1), a_opt(2));
% fprintf('目标函数值：%f\n', f_opt);
% fprintf('迭代次数：%d\n', k);


%% 拟牛顿法
% [a_opt, f_opt, k, x_history]=bfgs_optimizer(f, a0, 1e-4, 500);
% % 显示结果
% fprintf('最优解：(%f, %f)\n', a_opt(1), a_opt(2));
% fprintf('目标函数值：%f\n', f_opt);
% fprintf('迭代次数：%d\n', k);


%% 遗传算法
ub = [1;1;1;1;1;1;1;1;1;1;1;1];
lb = -[1;1;1;1;1;1;1;1;1;1;1;1];
params = struct(...
    'PopulationSize', 500, ...    % 种群大小
    'MaxGenerations', 1000, ...    % 最大迭代次数
    'CrossoverProb', 0.9, ...     % 交叉概率
    'MutationProb', 0.01, ...     % 变异概率
    'BitsPerVar', 10, ...         % 每变量二进制位数
    'Tol', 1e-6, ...              % 收敛阈值
    'Display', 'iter' ...         % 显示模式
);
[a_opt, best_fitness, convergence] = genetic_algorithm(f, 12, lb, ub,params);
figure;
plot(convergence);
    title('遗传算法收敛曲线');
    xlabel('迭代次数');
    ylabel('最优适应度值');
    grid on;
%% 绘制曲线

    [y_up,~] = cst_airfoil(x_up,a_opt(1:6), a_opt(7:12), N1, N2);
    [~,y_low] = cst_airfoil(x_low, a_opt(1:6),a_opt(7:12), N1, N2);
    figure;
    plot(x_up, y_up, 'b', 'LineWidth', 1.5, 'DisplayName', 'Upper Surface');
    hold on;
    plot(x_low, y_low, 'b', 'LineWidth', 1.5, 'DisplayName', 'Lower Surface');
    hold on;
    plot(x_t, y_t,'r', 'LineWidth', 1.5); % 绘制翼型轮廓
    axis equal;

%% 目标函数
function f = Target(a, target)
    persistent iter
    if isempty(iter)
        iter = 1;
    end

    a_up = a(1:6);
    a_low = a(7:12);
    yt_low = flip(target.y(1:59));  % 目标下翼面（反序后第一点为头部）
    yt_up = target.y(59:end); 

    N1 = 0.5;
    N2 = 1;
    
    % 生成原始下翼面坐标
    [y_up, y_low_raw] = cst_airfoil(target.xu, a_up, a_low, N1, N2);
    
    % 强制约束：下翼面头部必须严格小于0
    if y_low_raw(1) >= 0
        % 取反下翼面参数，确保头部向下
        a_low_corrected = -a_low;  % 参数取反
        [~, y_low] = cst_airfoil(target.xl, a_up, a_low_corrected, N1, N2);
        % 如果取反后头部仍为0，则强制添加微小负值
        if abs(y_low(1)) < 1e-10  % 检测浮点零值
            y_low(1) = -1e-6;     % 强制设为微小负值
        end
        % 更新参数以反映修正（可选，根据需求）
        a(7:12) = a_low_corrected;
    else
        y_low = y_low_raw;
    end

    % 尾部连接点强制闭合
    y_up(end) = (y_up(end) + y_low(end)) / 2;
    y_low(end) = y_up(end);
    
    % 计算误差
    f1 = sum((y_up - yt_up).^2);
    f2 = sum((y_low - yt_low).^2);
    f = f1 + f2;

    % ---------- 可视化 ----------
    if iter == 500
        figure(100); clf;
        hold on; grid on;
        plot(target.xu, yt_up, 'r--', 'LineWidth', 1.5);   % 目标上翼面
        plot(target.xl, yt_low, 'r--', 'LineWidth', 1.5);  % 目标下翼面
        plot(target.xu, y_up, 'b-', 'LineWidth', 1.5);     % 当前上翼面
        plot(target.xl, y_low, 'b-', 'LineWidth', 1.5);    % 修正后的下翼面
        title(['Iteration ', num2str(iter), ', Error = ', num2str(f, '%.6f')]);
        legend('Target Upper', 'Target Lower', 'Current Upper', 'Current Lower');
        axis equal;
        xlabel('x'); ylabel('y');
        iter = 1;
        drawnow;
    end

    iter = iter + 1;
end