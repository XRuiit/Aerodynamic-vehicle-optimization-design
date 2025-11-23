clc;clear;

%目标翼型
% 读取翼型数据
filename = 'rae2822.dat'; % 替换为你的文件名
% filename = '0012.dat';
data = dlmread(filename); % 读取数据
% data = readtable(filename); % 读取数据
x = data(:, 1); % 提取 x 坐标
y = data(:, 2); % 提取 y 坐标

% 绘制翼型
figure; % 创建新图形窗口
% plot(x, y, 'r','LineWidth', 1.5); % 绘制翼型轮廓
grid on; % 添加网格
axis equal; % 保持比例
xlabel('x'); % x 轴标签
ylabel('y'); % y 轴标签
title('翼型轮廓'); % 添加标题
hold on;

% a_up = [0.1, 0.1, 0.13, 0.15,0.1,0.18]; % low系数
% a_low = [-0.1, -0.1, -0.13, -0.1,-0.1,0.1];
% a_up = [0.16, 0.19, 0.2, 0.24,0.25,0.23]; % up系数
% a_low = [-0.17, -0.18, -0.23, -0.17,-0.15,0.08];
a_up = [0.133489736070381	0.146627565982405	0.187956989247312	0.223372434017595	0.152639296187683	0.211671554252199]; 
a_low = [-0.138044965786901	-0.148328445747801	-0.228142717497556	-0.139755620723363	-0.138954056695992	0.0999413489736070];

x_cst = linspace(0, 1, 200);
[y_up_cst,y_low_cst] = cst_airfoil(x_cst, a_up, a_low, 0.5, 1);
% figure; % 创建新图形窗口
plot(x_cst, y_up_cst, 'b','LineWidth', 1.5); % 绘制翼型轮廓
plot(x_cst, y_low_cst, 'b','LineWidth', 1.5); % 绘制翼型轮廓
grid on; % 添加网格
axis equal; % 保持比例
xlabel('x'); % x 轴标签
ylabel('y'); % y 轴标签



