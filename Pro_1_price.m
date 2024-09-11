clc;clear;
% 读取 Excel 文件
data = readmatrix('整理后数据.xlsx','Range','A2:L61');

% 初始化 P_min 和 P_max 矩阵
num_plots = 9;   % 地块数量
num_seasons = 3;  % 季度数量
num_crops = 41;   % 作物种类数量

P_min = zeros(num_plots, num_seasons, num_crops);  % 最低价格矩阵
P_max = zeros(num_plots, num_seasons, num_crops);  % 最高价格矩阵
C = zeros(num_plots, num_seasons, num_crops);  
D = zeros(num_plots, num_seasons, num_crops);  

% 遍历数据并填充 P_min 和 P_max
for k = 1:height(data)

    %读取下标
    n = data(k, 1);  % 地块编号
    i = data(k, 6);  % 季度编号
    j = data(k, 2);  % 作物编号

    price_min = data(k,11);%价格下限
    price_max = data(k,12);%价格上限

    %亩产矩阵
    muchan = data(k, 7);

    %成本矩阵
    chengben = data(k, 9);
    
    % 分割价格区间为下界和上界
    %prices = split(price_range, '-');
    P_min(n, i, j) = price_min;  % 最低价格
    P_max(n, i, j) = price_max;  % 最高价格
    C(n,i,j) = chengben;
    D(n,i,j) = muchan;
end

% 检查结果
disp(P_min(:, :, 1));  % 显示第1种作物的最低价格矩阵
disp(P_max(:, :, 1));  % 显示第1种作物的最高价格矩阵

%蒙特卡洛算法求最佳定价
num_iterations = 1000;  % 蒙特卡洛模拟次数
best_Z = -inf;  % 用于保存最大利润
best_P = zeros(num_plots, num_seasons, num_crops);  % 用于保存最佳定价矩阵

% 蒙特卡洛模拟
for iter = 1:num_iterations
    % 随机生成定价矩阵 P
    P = P_min + (P_max - P_min) .* rand(num_plots, num_seasons, num_crops);
    
    % 计算总利润 Z
    total_profit = 0;
    for n = 1:num_plots
        for i = 1:num_seasons
            for j = 1:num_crops
                % 计算每个作物的收益和成本
                revenue = P(n, i, j) * D(n, i, j);  % 收入，亩产等于销量
                cost = C(n, i, j);  % 成本
                profit = revenue - cost;  % 利润
                total_profit = total_profit + profit;  % 累加总利润
            end
        end
    end
    
    % 如果当前模拟的总利润更大，更新最优解
    if total_profit > best_Z
        best_Z = total_profit;
        best_P = P;
    end
end


% 准备数据导出到 Excel
[num_plots, num_seasons, num_crops] = size(best_P);
export_data = [];

for n = 1:num_plots
    for i = 1:num_seasons
        for j = 1:num_crops
            row = [n, i, j, best_P(n, i, j)];
            export_data = [export_data; row];  % 拼接数据
        end
    end
end

% 导出数据到 Excel 文件
filename = '最佳定价.xlsx';
writematrix(export_data, filename, 'Sheet', 1, 'Range', 'A1');
