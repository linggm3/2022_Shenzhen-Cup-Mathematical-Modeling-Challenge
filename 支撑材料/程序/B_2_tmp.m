clear; clc;

global center2; global center2_load; global edge;
global i;  global dott; global center_load;

%% 

% 读取/生成数据部分
data = readmatrix("12q.xlsx");
data = data ./ 1000;
supply = [data(1, 1:2); data(2, 1:2)];
x = data(3:size(data, 1), 1)'; % 点的x坐标
y = data(3:size(data, 1), 2)'; % 点的y坐标
power = data(2:size(data, 1), 3); % 各个负荷的功率
len = length(x);

% data = readmatrix("电网数据.xlsx");
% supply = data(1:2, 1:2);
% x = data(4:size(data, 1), 1)'; % 点的x坐标
% y = data(4:size(data, 1), 2)'; % 点的y坐标
% power = data(4:size(data, 1), 3); % 各个负荷的功率
% len = length(x);



%% 

% 第一层聚类部分, idx为各负荷所属簇的下标，center为各簇中心坐标
cluster_num = 4; % 聚类簇数k
[idx, center] = kmeans(cat(1, x, y)', cluster_num, 'Replicates', 2000); % 聚类2000次，选SSD最小的结果
center_load = count_num(idx, cluster_num); % 每簇对应的负荷数
% gscatter(x, y, idx); % 根据聚类结果，用不同颜色画出各点


%%

% 将所有簇分成两组，分别对应0，1两个电源
% 用二进制串表示每个簇的归属，0表示此簇归属到电源0
idxx = 1:cluster_num;
belong = zeros(1, cluster_num); % 长度为cluster_num的二进制串
L = zeros(1, 2^cluster_num-2);

% 算出最有可能的组合情况
for n = 1:2^cluster_num-2
    num = n;
    for m = cluster_num:-1:1  % 利用反复除2的方法得到二进制串
        belong(1, m) = mod(num, 2);
        num = floor(num / 2);
    end
                      
    L0 = minspan(cat(1, center(belong==0, :), supply(1, :) ) );
    L1 = minspan(cat(1, center(belong==1, :), supply(2, :)) );

    L(n) = max(L1, L0);  % 等于 (abs(L0 -L1) + L0 + L1) / 2
%     L(n) = L0 + L1;  

%     if sum(belong == 1) < 2 || sum(belong == 0) < 2  % 簇数约束
%         L(n) = 99999;
%     end
end

[~, sorted_index] = sort(L, 'ascend');


%% 


% 对可能的组合情况运用第一问的算法
ttt = 2; % 双层规划次数
for n = 1:1
    num = sorted_index(n);
    for m = cluster_num:-1:1  % 利用反复除2的方法得到二进制串     
        belong(1, m) = mod(num, 2);
        num = floor(num / 2);
    end
    figure;
    idxidx = idxx(belong == 0);
    [cost0, info0] = B1_fun_tmp(x, y, supply(1, :), cluster_num-sum(belong), idx, center(belong==0, :), idxidx, power, ttt);
    idxidx = idxx(belong == 1);
    [cost1, info1] = B1_fun_tmp(x, y, supply(2, :), sum(belong), idx, center(belong==1, :), idxidx, power, ttt);
    total_cost(n) = cost0 + cost1;
    fprintf("第二问最终花费%f千元\n", total_cost(n) );

    capacity0 = 0; capacity1 = 0;
    % 化整，方便背包问题求解
    for k = 1:size(info0{7}, 2)
        capacity0 = capacity0 + sum(info0{7}{k});
        info0{7}{k} = round(info0{7}{k} * 100);
    end
    for k = 1:size(info1{7}, 2)
        capacity1 = capacity1 + sum(info1{7}{k});
        info1{7}{k} = round(info1{7}{k} * 100);
    end
    capacity0 = capacity0 * 1.1 * 0.5; 
    capacity1 = capacity1 * 1.1 * 0.5; 
end

figure;
plot(total_cost)
