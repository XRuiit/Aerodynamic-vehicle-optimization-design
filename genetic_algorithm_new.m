function [best_solution, best_fitness, convergence] = genetic_algorithm_new(fitness_func, m, lb, ub, params)
% 遗传算法实现（二进制编码）含并行计算和精英保留
% 输入参数说明同上，新增UseParallel和EliteCount参数

% 默认参数设置
default_params = struct(...
    'PopulationSize', 100, ...    
    'MaxGenerations', 200, ...    
    'CrossoverProb', 0.8, ...     
    'MutationProb', 0.01, ...     
    'BitsPerVar', 20, ...         
    'Tol', 1e-6, ...              
    'Display', 'iter', ...
    'UseParallel', false, ...     % 新增：是否启用并行计算
    'EliteCount', 2 ...           % 新增：精英保留数量
);

% 合并参数
% params = merge_params(default_params, params);
if nargin < 5
    params = default_params;
else
    params = merge_params(default_params, params);
end
% 初始化变量
n_p = params.PopulationSize;
max_gen = params.MaxGenerations;
p_c = params.CrossoverProb;
p_m = params.MutationProb;
q = params.BitsPerVar;
total_bits = m * q;
e = params.EliteCount;

% 变量边界处理
if isscalar(lb), lb = repmat(lb, 1, m); end
if isscalar(ub), ub = repmat(ub, 1, m); end

% 初始化种群
population = randi([0 1], n_p, total_bits);

% 主循环
convergence = zeros(max_gen, 1);
for gen = 1:max_gen
    % 解码并计算适应度（含并行支持）
    [fitness, decoded_pop] = evaluate_population(population, fitness_func, m, q, lb, ub, params.UseParallel);
    
    % 记录最优解
    [current_best, idx] = min(fitness);
    convergence(gen) = current_best;
    best_solution = decoded_pop(idx, :);
    best_fitness = current_best;
    
    % 显示信息
    if strcmp(params.Display, 'iter') && mod(gen,1) == 0
        fprintf('Generation %d: Best Fitness = %.4f\n', gen, current_best);
    end
    % 早停检查（新增）
    if gen > 20 && std(convergence(gen-20:gen)) < params.Tol
        fprintf('Converged at generation %d\n', gen);
        break;
    end
    
    % 精英保留
    [sorted_fit, sort_idx] = sort(fitness);
    elite_pop = population(sort_idx(1:e), :);
    elite_fitness = sorted_fit(1:e);
    
    % 选择父代（剩余数量）
    parents = roulette_selection(population, fitness, n_p - e);
    
    % 交叉（生成指定数量的子代）
    offspring = crossover(parents, p_c, n_p - e);
    
    % 向量化变异（替换循环）
    offspring = mutate_vectorized(offspring, p_m);
    
    % 合并新种群
    population = [offspring; elite_pop];
end

% 最终评估
[best_fitness, decoded_pop] = evaluate_population(population, fitness_func, m, q, lb, ub);
[best_fitness, idx] = min(fitness);
best_solution = decoded_pop(idx, :);

% ================== 优化后的辅助函数 ==================
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

function [fitness, decoded] = evaluate_population(pop, fun, m, q, lb, ub, use_parallel)
% 并行评估适应度
n = size(pop, 1);
decoded = binary2real_vectorized(pop, m, q, lb, ub);

if use_parallel
    fitness = zeros(n, 1);
    parfor k = 1:n
        fitness(k) = fun(decoded(k,:));
    end
else
    fitness = arrayfun(@(k) fun(decoded(k,:)), 1:n)';
end
end

function real_vec = binary2real_vectorized(bin_mat, m, q, lb, ub)
% 向量化解码
n = size(bin_mat, 1);
real_vec = zeros(n, m);

% 预计算权重矩阵
weights = 2.^(q-1:-1:0)';

for i = 1:m
    % 提取每变量的二进制位
    bits = bin_mat(:, (i-1)*q+1:i*q);
    % 矩阵乘法计算整数值
    int_vals = double(bits) * weights;
    % 线性映射到实数范围
    real_vec(:,i) = lb(i) + (ub(i)-lb(i)) * int_vals/(2^q -1);
end
end

function parents = roulette_selection(pop, fitness, n_parents)
% 高效轮盘赌选择
prob = 1./(fitness + eps);
prob = prob / sum(prob);
idx = randsample(size(pop,1), n_parents, true, prob);
parents = pop(idx, :);
end

function offspring = crossover(parents, p_c, n_offspring)
% 可控数量的交叉
[~, len] = size(parents);
offspring = zeros(n_offspring, len);
pair_count = ceil(n_offspring/2);

for i = 1:pair_count
    % 随机选择父代
    idx = randperm(size(parents,1), 2);
    p1 = parents(idx(1),:);
    p2 = parents(idx(2),:);
    
    if rand() < p_c
        cp = randi(len-1);
        c1 = [p1(1:cp) p2(cp+1:end)];
        c2 = [p2(1:cp) p1(cp+1:end)];
    else
        c1 = p1; c2 = p2;
    end
    
    pos = 2*i-1;
    if pos+1 <= n_offspring
        offspring(pos,:) = c1;
        offspring(pos+1,:) = c2;
    else
        offspring(pos,:) = c1;
    end
end
end

function mutated = mutate_vectorized(pop, p_m)
% 向量化变异
mutation_mask = rand(size(pop)) < p_m;
mutated = pop;
mutated(mutation_mask) = 1 - pop(mutation_mask);
end

% 其他辅助函数同原实现
end