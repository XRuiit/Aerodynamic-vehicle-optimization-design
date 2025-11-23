clc;clear;
%% RAE2822标准翼型
% 上翼面数据
x_upper = [0.000000, 0.000602, 0.002408, 0.005412, 0.009607, 0.014984, 0.021530, 0.029228, 0.038060, 0.048005, 0.059039, 0.071136, 0.084265, 0.098396, 0.113495, 0.129524, 0.146447, 0.164221, 0.182803, 0.202150, 0.222215, 0.242949, 0.264302, 0.286222, 0.308658, 0.331555, 0.354858, 0.378510, 0.402455, 0.426635, 0.450991, 0.475466, 0.500000, 0.524534, 0.549009, 0.573365, 0.597545, 0.621490, 0.645142, 0.668445, 0.691342, 0.713778, 0.735698, 0.757051, 0.777785, 0.797850, 0.817197, 0.835779, 0.853553, 0.870476, 0.886505, 0.901604, 0.915735, 0.928864, 0.940961, 0.951995, 0.961940, 0.970772, 0.978470, 0.985016, 0.990393, 0.994588, 0.997592, 0.999398, 1.000000];
y_upper = [0.000000, 0.003165, 0.006306, 0.009416, 0.012480, 0.015489, 0.018441, 0.021348, 0.024219, 0.027062, 0.029874, 0.032644, 0.035360, 0.038011, 0.040585, 0.043071, 0.045457, 0.047729, 0.049874, 0.051885, 0.053753, 0.055470, 0.057026, 0.058414, 0.059629, 0.060660, 0.061497, 0.062133, 0.062562, 0.062779, 0.062774, 0.062530, 0.062029, 0.061254, 0.060194, 0.058845, 0.057218, 0.055344, 0.053258, 0.050993, 0.048575, 0.046029, 0.043377, 0.040641, 0.037847, 0.035017, 0.032176, 0.029347, 0.026554, 0.023817, 0.021153, 0.018580, 0.016113, 0.013769, 0.011562, 0.009508, 0.007622, 0.005915, 0.004401, 0.003092, 0.002001, 0.001137, 0.000510, 0.000128, 0.000000];

% 下翼面数据
x_lower = [0.000000, 0.000602, 0.002408, 0.005412, 0.009607, 0.014984, 0.021530, 0.029228, 0.038060, 0.048005, 0.059039, 0.071136, 0.084265, 0.098396, 0.113495, 0.129524, 0.146447, 0.164221, 0.182803, 0.202150, 0.222215, 0.242949, 0.264302, 0.286222, 0.308658, 0.331555, 0.354858, 0.378510, 0.402455, 0.426635, 0.450991, 0.475466, 0.500000, 0.524534, 0.549009, 0.573365, 0.597545, 0.621490, 0.645142, 0.668445, 0.691342, 0.713778, 0.735698, 0.757051, 0.777785, 0.797850, 0.817197, 0.835779, 0.853553, 0.870476, 0.886505, 0.901604, 0.915735, 0.928864, 0.940961, 0.951995, 0.961940, 0.970772, 0.978470, 0.985016, 0.990393, 0.994588, 0.997592, 0.999398, 1.000000];
y_lower = [0.000000, -0.003160, -0.006308, -0.009443, -0.012559, -0.015649, -0.018707, -0.021722, -0.024685, -0.027586, -0.030416, -0.033170, -0.035843, -0.038431, -0.040929, -0.043326, -0.045610, -0.047773, -0.049805, -0.051694, -0.053427, -0.054994, -0.056376, -0.057547, -0.058459, -0.059046, -0.059236, -0.058974, -0.058224, -0.056979, -0.055257, -0.053099, -0.050563, -0.047719, -0.044642, -0.041397, -0.038043, -0.034631, -0.031207, -0.027814, -0.024495, -0.021289, -0.018232, -0.015357, -0.012690, -0.010244, -0.008027, -0.006048, -0.004314, -0.002829, -0.001592, -0.000600, 0.000157, 0.000694, 0.001033, 0.001197, 0.001212, 0.001112, 0.000935, 0.000719, 0.000497, 0.000296, 0.000137, 0.000035, 0.000000];

%% 输入参数
leading_edge = [0, 0]; % 前缘顶点坐标
trailing_edge = [1, 0]; % 后缘顶点坐标
upper_leading_slope = tan(deg2rad(32)); % 上表面前缘斜率
upper_trailing_slope = tan(deg2rad(-8));% 上表面后缘斜率
lower_leading_slope = tan(deg2rad(-36)); % 下表面前缘斜率（新增）
lower_trailing_slope = tan(deg2rad(-4)); % 下表面后缘斜率（新增）
upper_peak = [0.38, 0.063]; % 上表面极值点坐标
lower_peak = [0.293, -0.058]; % 下表面极值点坐标

%% 生成上下表面贝塞尔曲线控制点
[upper_ctrl, lower_ctrl] = calculate_control_points(...
    leading_edge, trailing_edge,...
    upper_leading_slope, upper_trailing_slope,...
    lower_leading_slope, lower_trailing_slope,...
    upper_peak, lower_peak);

%% 绘制结果
t = linspace(0, 1, 100);
upper_curve = bezier_curve(upper_ctrl, t);
lower_curve = bezier_curve(lower_ctrl, t);
figure
hold on
plot(upper_curve(:,1), upper_curve(:,2), 'b-', 'LineWidth', 2)
plot(lower_curve(:,1), lower_curve(:,2), 'b-', 'LineWidth', 2)
plot(upper_ctrl(:,1), upper_ctrl(:,2), 'bo--')
plot(lower_ctrl(:,1), lower_ctrl(:,2), 'bo--')

plot(x_upper, y_upper, 'r', 'LineWidth', 2); % 上翼面
plot(x_lower, y_lower, 'r', 'LineWidth', 2); % 下翼面
axis equal
grid on
title('贝塞尔曲线翼型生成')
legend('生成上表面','生成下表面','控制点','RAE2822上表面','RAE2822下表面')

%% 计算误差
% 在相同x坐标位置比较生成翼型和RAE2822的y值
% 上表面误差计算
upper_err = zeros(length(x_upper),1);
for i = 1:length(x_upper)
    % 找到生成曲线上最接近当前x坐标的点
    [~, idx] = min(abs(upper_curve(:,1) - x_upper(i)));
    upper_err(i) = abs(upper_curve(idx,2) - y_upper(i));
end

% 下表面误差计算
lower_err = zeros(length(x_lower),1);
for i = 1:length(x_lower)
    [~, idx] = min(abs(lower_curve(:,1) - x_lower(i)));
    lower_err(i) = abs(lower_curve(idx,2) - y_lower(i));
end

% 计算统计误差
max_upper_err = max(upper_err);
mean_upper_err = mean(upper_err);
rmse_upper = sqrt(mean(upper_err.^2));

max_lower_err = max(lower_err);
mean_lower_err = mean(lower_err);
rmse_lower = sqrt(mean(lower_err.^2));

total_err = [upper_err; lower_err];
max_total_err = max(total_err);
mean_total_err = mean(total_err);
rmse_total = sqrt(mean(total_err.^2));

% 显示误差结果
disp('========== 误差分析结果 ==========');
fprintf('上表面最大误差: %.6f (%.2f%% 相对误差)\n', max_upper_err, 100*max_upper_err/max(y_upper));
fprintf('上表面平均误差: %.6f\n', mean_upper_err);
fprintf('上表面均方根误差(RMSE): %.6f\n', rmse_upper);
fprintf('--------------------------------\n');
fprintf('下表面最大误差: %.6f (%.2f%% 相对误差)\n', max_lower_err, 100*max_lower_err/abs(min(y_lower)));
fprintf('下表面平均误差: %.6f\n', mean_lower_err);
fprintf('下表面均方根误差(RMSE): %.6f\n', rmse_lower);
fprintf('--------------------------------\n');
fprintf('整体最大误差: %.6f\n', max_total_err);
fprintf('整体平均误差: %.6f\n', mean_total_err);
fprintf('整体均方根误差(RMSE): %.6f\n', rmse_total);


%% 辅助函数
% 控制点函数
function [upper_ctrl, lower_ctrl] = calculate_control_points(...
    le, te, ...
    upper_lead_slope, upper_trail_slope, ...
    lower_lead_slope, lower_trail_slope, ...
    upper_peak, lower_peak)
    
    n_ctrl = 4; 
    upper_ctrl = zeros(n_ctrl, 2);
    lower_ctrl = zeros(n_ctrl, 2);
    
    % 设置端点
    upper_ctrl(1,:) = le;
    upper_ctrl(end,:) = te;
    lower_ctrl(1,:) = le;
    lower_ctrl(end,:) = te;
    
    % 前缘控制点
    delta = 0.2*norm(te-le);
    upper_ctrl(2,:) = le + delta*[1, upper_lead_slope];
    lower_ctrl(2,:) = le + delta*[1, lower_lead_slope];
    
    % 后缘控制点
    upper_ctrl(end-1,:) = te - delta*[1, upper_trail_slope];
    lower_ctrl(end-1,:) = te - delta*[1, lower_trail_slope];
    
    % 极值点优化
    options = optimoptions('fmincon', 'Display', 'off');
    upper_ctrl(3,:) = fmincon(@(p) peak_error(p, upper_ctrl, upper_peak),...
        upper_ctrl(3,:), [],[],[],[],[],[], @(p) peak_constraint(p, upper_ctrl, upper_peak), options);
    
    lower_ctrl(3,:) = fmincon(@(p) peak_error(p, lower_ctrl, lower_peak),...
        lower_ctrl(3,:), [],[],[],[],[],[], @(p) peak_constraint(p, lower_ctrl, lower_peak), options);
end

function curve = bezier_curve(ctrl_pts, t)
    % 计算贝塞尔曲线
    n = length(t);
    curve = zeros(n,2);
    order = size(ctrl_pts,1)-1;
    for i = 0:order
        weight = nchoosek(order,i) .* (1-t).^(order-i) .* t.^i;
        curve = curve + weight' .* ctrl_pts(i+1,:);
    end
end
% 极值点误差函数
function err = peak_error(p, ctrl_pts, target)
    
    ctrl_pts(3,:) = p;
    t_peak = fzero(@(t) derivative_zero(t, ctrl_pts), 0.5);
    current_peak = bezier_curve(ctrl_pts, t_peak);
    err = norm(current_peak - target);
end
% 极值点约束条件
function [c, ceq] = peak_constraint(p, ctrl_pts, target)
    
    ctrl_pts(3,:) = p;
    t_peak = fzero(@(t) derivative_zero(t, ctrl_pts), 0.5);
    current_peak = bezier_curve(ctrl_pts, t_peak);
    ceq = norm(current_peak - target);
    c = [];
end
% 计算导数为零的条件
function dy = derivative_zero(t, ctrl_pts)
    
    order = size(ctrl_pts,1)-1;
    deriv = zeros(1,2);
    for i = 0:order-1
        weight = nchoosek(order-1,i) * (1-t).^(order-1-i) .* t.^i;
        deriv = deriv + order*(ctrl_pts(i+2,:)-ctrl_pts(i+1,:)) * weight;
    end
    dy = deriv(2); % y方向导数为零
end