function [best_solution, best_fitness, convergence] = genetic_algorithm(fitness_func, m, lb, ub, params)
% 遗传算法实现（二进制编码）
% 输入：
%   fitness_func : 目标函数（极小化问题）
%   m           : 设计变量个数
%   lb          : 变量下界（标量或向量）
%   ub          : 变量上界（标量或向量）
%   params      : 算法参数结构体
% 输出：
%   best_solution : 最优解
%   best_fitness  : 最优适应度值
%   convergence   : 收敛曲线

% 默认参数设置
default_params = struct(...
    'PopulationSize', 100, ...    % 种群大小
    'MaxGenerations', 200, ...    % 最大迭代次数
    'CrossoverProb', 0.8, ...     % 交叉概率
    'MutationProb', 0.01, ...     % 变异概率
    'BitsPerVar', 20, ...         % 每变量二进制位数
    'Tol', 1e-6, ...              % 收敛阈值
    'Display', 'iter' ...         % 显示模式
);

% 合并参数
if nargin < 5
    params = default_params;
else
    params = merge_params(default_params, params);
end

% 参数初始化
n_p = params.PopulationSize;
max_gen = params.MaxGenerations;
p_c = params.CrossoverProb;
p_m = params.MutationProb;
q = params.BitsPerVar;
total_bits = m * q;

% 变量边界处理
if isscalar(lb), lb = repmat(lb, 1, m); end
if isscalar(ub), ub = repmat(ub, 1, m); end

% 初始化种群
population = randi([0 1], n_p, total_bits);

% 主循环
convergence = zeros(max_gen, 1);
for gen = 1:max_gen
    % 解码并计算适应度
    [fitness, decoded_pop] = evaluate_population(population, fitness_func, m, q, lb, ub);
    
    % 记录最优解
    [current_best, idx] = min(fitness);
    convergence(gen) = current_best;
    best_solution = decoded_pop(idx, :);
    best_fitness = current_best;
    
    % 显示迭代信息
    if strcmp(params.Display, 'iter') && mod(gen, 10) == 0
        fprintf('Generation %d: Best Fitness = %.4f\n', gen, current_best);
    end
    
    % 收敛检查
    if std(fitness) < params.Tol
        break;
    end
    
    % 选择操作（轮盘赌）
    parents = roulette_selection(population, fitness);
    
    % 交叉操作（单点交叉）
    offspring = crossover(parents, p_c);
    
    % 变异操作（基本位变异）
    population = mutate(offspring, p_m);
end

% 最终解码
[best_fitness, decoded_pop] = evaluate_population(population, fitness_func, m, q, lb, ub);
[best_fitness, idx] = min(fitness);
best_solution = decoded_pop(idx, :);

% ================== 辅助函数 ==================
function params = merge_params(default, user)
    % 合并参数结构体
    params = default;
    if ~isempty(user)
        fields = fieldnames(user);
        for i = 1:length(fields)
            params.(fields{i}) = user.(fields{i});
        end
    end
end

function [fitness, decoded] = evaluate_population(pop, fun, m, q, lb, ub)
    % 解码二进制种群并计算适应度
    n = size(pop, 1);
    decoded = zeros(n, m);
    for i = 1:n
        decoded(i,:) = binary2real(pop(i,:), m, q, lb, ub);
    end
    fitness = arrayfun(@(k) fun(decoded(k,:)), 1:n)';
end

function real_vec = binary2real(bin_str, m, q, lb, ub)
    % 二进制转实数
    real_vec = zeros(1, m);
    for i = 1:m
        start = (i-1)*q + 1;
        end_idx = i*q;
        bin_part = bin_str(start:end_idx);
        int_val = bin2dec(num2str(bin_part));
        real_vec(i) = lb(i) + (ub(i)-lb(i)) * int_val/(2^q -1);
    end
end

function parents = roulette_selection(pop, fitness)
    % 轮盘赌选择
    n = size(pop, 1);
    prob = 1./(fitness + eps); % 适应度倒数（极小化问题）
    prob = prob / sum(prob);   % 归一化概率
    
    cum_prob = cumsum(prob);
    parents = zeros(size(pop));
    
    for i = 1:n
        r = rand();
        idx = find(r <= cum_prob, 1);
        parents(i,:) = pop(idx,:);
    end
end

function offspring = crossover(parents, p_c)
    % 单点交叉
    [n, len] = size(parents);
    offspring = parents;
    num_pairs = floor(n * p_c / 2);
    
    for i = 1:num_pairs
        % 选择父代
        idx = randperm(n, 2);
        p1 = parents(idx(1), :);
        p2 = parents(idx(2), :);
        
        % 随机交叉点
        cross_point = randi([1 len-1]);
        
        % 生成子代
        c1 = [p1(1:cross_point), p2(cross_point+1:end)];
        c2 = [p2(1:cross_point), p1(cross_point+1:end)];
        
        offspring(idx(1), :) = c1;
        offspring(idx(2), :) = c2;
    end
end

function mutated = mutate(pop, p_m)
    % 基本位变异
    [n, len] = size(pop);
    mutated = pop;
    for i = 1:n
        for j = 1:len
            if rand() < p_m
                mutated(i,j) = ~pop(i,j);
            end
        end
    end
end

end