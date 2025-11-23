%% 调用exe
clc;clear;

% addpath('C:\Users\许景瑞\Desktop\main');
% system('main.exe');

dos('E:\00\main\main.exe');

%% 核心函数实现
function fitness = constrained_fitness(weights, n)
global target_data constraints current_individual_idx
% 1. 生成翼型并保存到geo_example.txt
x = linspace(0,1,100)';
[y_upper, y_lower] = cst_airfoil(x, weights, @(x) x.^0.5.*(1-x), n);
save_airfoil_to_geo_example(x, y_upper, y_lower, current_individual_idx);
% 2. 执行外部程序main.exe
if current_individual_idx == 1 % 只需运行一次
if ~exist('main.exe', 'file')
error('main.exe 未找到');
end
[status,~] = system('main.exe');
if status ~= 0
error('main.exe 执行失败');
end
end
% 3. 从result/1/读取气动数据
[Cl, Cd] = read_aero_data(current_individual_idx);
% 4. 计算几何参数
[thickness, max_camber] = calculate_airfoil_geometry(...
[x; x], [y_upper; y_lower]);
% 5. 约束检查与罚函数
penalty = 0;
if Cd <= 0 || isnan(Cd)
penalty = penalty + constraints.penalty_weights.Cd_negative;
else
if Cd > constraints.max_Cd
penalty = penalty + constraints.penalty_weights.Cd_violation * ...
(Cd - constraints.max_Cd);
end
if Cl < constraints.min_Cl
penalty = penalty + constraints.penalty_weights.Cl_violation * ...
(constraints.min_Cl - Cl);
end
end
if thickness < constraints.min_thickness
penalty = penalty + constraints.penalty_weights.thickness * ...
(constraints.min_thickness - thickness);
end
if max_camber > constraints.max_camber
penalty = penalty + constraints.penalty_weights.camber * ...
(max_camber - constraints.max_camber);
end
% 6. 计算目标函数
shape_error = calculate_shape_error(x, y_upper, y_lower, target_data);
fitness = Cd + 0.1*shape_error + penalty;
end