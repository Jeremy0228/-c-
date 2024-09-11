clc;
clear;

% 读取地块面积
data2 = readmatrix("最最新地块面积编号.xlsx",'Range','A3:C84');

% 初始化地块面积数组
num_plots = 82;
A = zeros(num_plots, 1);

for u = 1:82
    bianhao = data2(u, 2);
    mianji = data2(u, 3);
    A(bianhao) = mianji;
end

num_crops = 41;
num_years = 7;

% 读取数据文件
file_sales = '预销售量.xlsx';
file_yield = '亩产.xlsx';
file_price = '销售单价.xlsx';
file_cost = '成本.xlsx';

% 读取销售量数据
sales_data = readmatrix(file_sales, 'Range', 'B2:AQ84'); % 根据实际数据范围调整范围
sales_plot = sales_data(2:end, 1); % 地块编号
sales_values = sales_data(2:end, 2:end); % 销量值

% 读取单位亩产数据
yield_data = readmatrix(file_yield, 'Range', 'B2:AQ84'); % 根据实际数据范围调整范围
yield_plot = yield_data(2:end, 1); % 地块编号
yield_values = yield_data(2:end, 2:end); % 单位亩产值

% 读取售价数据
price_data = readmatrix(file_price, 'Range', 'B2:AQ84'); % 根据实际数据范围调整范围
price_plot = price_data(2:end, 1); % 地块编号
price_values = price_data(2:end, 2:end); % 售价值

% 读取成本数据
cost_data = readmatrix(file_cost, 'Range', 'B2:AQ84'); % 根据实际数据范围调整范围
cost_plot = cost_data(2:end, 1); % 地块编号
cost_values = cost_data(2:end, 2:end); % 成本值

% 确保四个矩阵的地块编号一致
assert(isequal(sales_plot, yield_plot, price_plot, cost_plot), '地块编号不一致，请检查数据。');

% 初始化矩阵
D = zeros(num_plots, num_crops, num_years); % 销量矩阵
P = zeros(num_plots, num_crops, num_years); % 单位亩产矩阵
S = zeros(num_plots, num_crops, num_years); % 售价矩阵
C = zeros(num_plots, num_crops, num_years); % 成本矩阵
X = zeros(num_plots, num_crops, num_years); % 决策变量
% 将数据填入矩阵
D(:, :, 1) = sales_values; % 销量
P(:, :, 1) = yield_values; % 单位亩产
S(:, :, 1) = price_values; % 售价
C(:, :, 1) = cost_values; % 成本

% 复制到所有年份
for u = 1:num_years
    D(:,:,u) = D(:, :, 1); 
    P(:,:,u) = P(:, :, 1);
    S(:,:,u) = S(:, :, 1); 
    C(:,:,u) = C(:, :, 1);
end

%生成随机数波动量矩阵
a = zeros(num_crops, num_years);
b = zeros(num_crops, num_years);
d = zeros(num_crops, num_years);
e = zeros(num_crops, num_years);

%随机化
% 随机生成a、b、d、e的值
a = 0.05 + (0.1 - 0.05) * rand(num_crops, num_years);  % a在0.05到0.1之间
b = -0.05 + (0.05 + 0.05) * rand(num_crops, num_years); % b在-0.05到0.05之间
d = -0.1 + (0.1 + 0.1) * rand(num_crops, num_years);    % d在-0.1到0.1之间
e = 0.01 + (0.05 - 0.01) * rand(num_crops, num_years);  % e在0.01到0.05之间

%更新矩阵的值
for i = 1:num_plots
    for j = 1:num_crops
        for k = 1:6
            if j>=1&&j<=16
                S(i,j,k) = S(i,j,k);

            elseif j>=17&&j<=37
                S(i,j,k) = S(i,j,k).*((1+0.05).^k);

            elseif j>=38&&j<=40
                S(i,j,k) = S(i,j,k).*prod(1-e(j,k),2);

            elseif j == 41
                S(i,j,k) = S(i,j,k).*((1-0.05).^k);

            end

            P(i,j,k) = P(i,j,k).*prod(1+d(j,k),2);

            C(i,j,k) = C(i,j,k).*((1+0.05).^k);

            if ismember(j,[6,7])
                D(i,j,k) = D(i,j,k).*prod(1+a(j,k),2);
            else
                D(i,j,k) = D(i,j,k).*(1+b(j,k));
            end
        end
    end
end

%*******************************************************************************************

%X变成一维向量，
X_flat = X(:);
objectiveFunction = @(X_flat) obj(X_flat,P,S,C,A,D,num_plots, num_crops, num_years);


%设置函数
% 设置模拟退火算法的选项  
options = optimoptions('simulannealbnd', ...
            'MaxIterations', 1000, ...         % 最大迭代次数  
            'StallIterLim', 300, ...      % 停滞迭代次数限制（如果连续这么多次迭代没有改进，则停止）  
            'TolFun', 1e-4, ...          % 函数值容忍度（如果函数值的改进小于这个值，则视为收敛）  
            'InitialTemperature', 100, ... % 初始温度  
            'TemperatureFcn', @temperatureexp, ...%123
            'PlotFcn', @saplotbestf);

% 定义变量的上下界
lb = zeros(num_plots, num_crops, num_years); % 下界
ub = repmat(80, num_plots, num_crops, num_years); % 上界

% 将下界和上界展平为向量
lb = lb(:);
ub = ub(:);

%生成初始解x0即2023年种植方案
X0 = zeros(num_plots, num_crops, num_years);
data3 = readmatrix("种植面积.xlsx",'Range','B2:AQ84');
first_plot = data3(2:end, 1); % 地块编号
first_values = data3(2:end, 2:end); %种植方案

X0(:, :, 1) = first_values; % 
X0(:, :, 3) = first_values;
X0(:, :, 5) = first_values;
X0(:, :, 7) = first_values;

X0_flat = X0(:);


%退火函数

% 执行模拟退火算法  
[X_best_flat, fval] = simulannealbnd(objectiveFunction, X0_flat, lb, ub, options);

X_best = reshape(X_best_flat, num_plots, num_crops, num_years);

% 创建 Excel 文件名
filename = '最终结果3.xlsx';

% 将 X_best 导出为 Excel 文件
for k = 1:num_years
    % 提取第 k 年的数据
    data = X_best(:, :, k);
    
    % 创建一个工作表名
    sheetName = sprintf('Year%d', k);
    
    % 写入 Excel 文件
    writematrix(data, filename, 'Sheet', sheetName);
end

disp('优化结果已导出为 Excel 文件。');
disp("最大值：");
disp(-fval);





%*******************************************************************************************

function Z = obj(X_flat,P,S,C,A,D, num_plots, num_crops, num_years)
    Z = 0;
    punish = 0;%惩罚度
    punish_fact = 5500;%因子5500

    %将X矩阵重构
    X = reshape(X_flat, num_plots, num_crops, num_years);
    %
    
    for i = 1:num_plots
        for j = 1:num_crops
            for k = 1:num_years
                Z = Z + (P(i,j,k).*X(i,j,k).*S(i,j,k)+(P(i,j,k).*X(i,j,k)-D(i,j,k)).*0.5.*S(i,j,k)-C(i,j,k).*X(i,j,k));
            end    
        end
    end

    %新增约束条件，相关性互补性等
    for i = 1:num_plots
        for k = 1:num_years
            punish = punish_fact*(X(i,9,k).*P(i,9,k)+X(i,15,k).*P(i,15,k)-40000);
        end
    end

    for i = 1:num_plots
        for j = 1:num_crops
            for k = 1:num_years
                punish = punish_fact*(D(i,j,k)+0.204.*S(i,j,k)+0.2431.*C(i,j,k));
            end
        end
    end
    

    %1.三年种豆>0
    for i = 1:num_plots
        for j = 1:num_crops

            if ismember(j,[1,2,3,4,5,17,18,19])
                for k0 = 1:5
                    punish = punish + punish_fact*isequal(sum(X(i,j,k0:k0+2),3),0);
                end
            end
           
            
        end
    end

    %2.???????????????=0
    %mem10 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26];
    for k = 1:num_years
        for i = 1:26
            punish = punish + punish_fact*(sum(X(i,16:41,k),2));
        end
        %punish = punish + punish_fact*(sum(X(mem10,1:15 , k),2)-sum(X(mem10,1:41, k),2)).^2;
    end


   %3.水稻---456约束
   for k = 1:num_years
       for j = 1:num_crops
           for i = 1:num_plots

                if j==16&&(i>=27&&i<=34)
                    punish = punish + punish_fact*(X(i,j,k)==0||X(i+28,j,k)>0);
                %
                elseif (j>=17&&j<=34)&&(i>=27&&i<=34)
                    for u = 35:37
                        punish = punish + punish_fact*((X(i,j,k)==0)||X(i+28,u,k)==0);
                    end

                %
                elseif (j>=17&&j<=34)&&(i>=35&&i<=50)
                    for u = 38:41
                        punish = punish + punish_fact*((X(i,j,k)==0)||X(i+28,u,k)==0);
                    end


                elseif (j>=17&&j<=34)&&(i>=51&&i<=54)
                    punish = punish + punish_fact*((X(i,j,k)==0)||(X(i+28,j,k)==0));
                

                end

           end
       end
   end







    %%%%%%%%%7.=0
    for i = 1:num_plots
        for j = 1:num_crops
            for k = 1:6
                punish = punish + punish_fact*(X(i,j,k)*X(i,j,k+1));%%%%%%%%
            end
        end
    end
    %8.<=
    for i = 1:num_plots
        for k = 1:num_years
            punish = punish + punish_fact*max((sum(X(i,1:41,k),2)-A(i)),0);
        end
    end



   

    Z = -1*Z + punish;%反过来
end
