clear; clc;

global center2; global center2_load; global edge;
global i;  global dott; global center_load;

%% 

% 读取/生成数据部分
data = readmatrix("12q.xlsx");
data = data./1000;
supply = data(1:2, 1:2);
x = data(3:size(data, 1), 1)'; % 点的x坐标
y = data(3:size(data, 1), 2)'; % 点的y坐标
power = data(3:size(data, 1), 3); % 各个负荷的功率
len = length(x);

% % 读取/生成数据部分
% data = readmatrix("12q.xlsx");
% data = data ./ 1000;
% supply = data(1:2, 1:2);
% x = data(3:size(data, 1), 1)'; % 点的x坐标
% y = data(3:size(data, 1), 2)'; % 点的y坐标
% len = length(x);


%% 

% 第一层聚类部分, idx为各负荷所属簇的下标，center为各簇中心坐标
cluster_num = 5; % 聚类簇数k
[idx, center] = kmeans(cat(1, x, y)', cluster_num, 'Replicates', 2000); % 聚类2000次，选SSD最小的结果
center_load = count_num(idx, cluster_num); % 每簇对应的负荷数
% gscatter(x, y, idx); % 根据聚类结果，用不同颜色画出各点
Color = [250 192 15; 1 86 153; 243 118 74; 95 198 201; 79 89 100] / 255;

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
    if sum(belong == 1) < 2 || sum(belong == 0) < 2  % 簇数约束
        L(n) = 99999;
    end
end

[~, sorted_index] = sort(L, 'ascend');


%% 


% 对可能的组合情况运用第一问的算法
ttt = 1; % 双层规划次数
for n = 1:1
    num = sorted_index(n);
    for m = cluster_num:-1:1  % 利用反复除2的方法得到二进制串     
        belong(1, m) = mod(num, 2);
        num = floor(num / 2);
    end
    figure;
    idxidx = idxx(belong == 0);
    [cost0, info0] = B1_fun(x, y, supply(1, :), cluster_num-sum(belong), idx, center(belong==0, :), idxidx, power, ttt);
    idxidx = idxx(belong == 1);
    [cost1, info1] = B1_fun(x, y, supply(2, :), sum(belong), idx, center(belong==1, :), idxidx, power, ttt);
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
    capacity0 = capacity0 * 1.2 * 0.5; 
    capacity1 = capacity1 * 1.2 * 0.5; 
end

% figure;
% plot(total_cost);

%% 

min_re = 1;  % 最低可靠性
sum_re = 0;
for i = 1:35
    for m = 1:size(info0{2}, 2) % 簇
        for n = 1:size(info0{2}{m}, 1)
            if x(i) == info0{2}{m}(n, 1) && y(i) == info0{2}{m}(n, 2)   
                re(i) = reliability_fun(0, info0{2}{m}(n, :), info0{1}, supply(1, :), info0{2}, info0{3}, info0{4}, info0{5}, info0{6});
                sum_re = sum_re + re(i);
                if re(i) < min_re
                    min_re = re(i);
                end
            end
        end
    end
    for m = 1:size(info1{2}, 2)
        for n = 1:size(info1{2}{m}, 1)
            if x(i) == info1{2}{m}(n, 1) && y(i) == info1{2}{m}(n, 2)
                re(i) = reliability_fun(0, info1{2}{m}(n, :), info1{1}, supply(2, :), info1{2}, info1{3}, info1{4}, info1{5}, info1{6});
                sum_re = sum_re + re(i);
                if re(i) < min_re
                    min_re = re(i);
                end
            end
        end
    end
end
fprintf("不加联络线的情况下最低可靠性为%f\n", min_re);


%% 

% % 第三第四问
info0 = load("info0.mat");
info1 = load("info1.mat");
info0 = struct2cell(info0);
info1 = struct2cell(info1);
info0 = info0{1};
info1 = info1{1};
% 
% 
% dot(1:20, :) = info0{1};
% cluster = info0{2};
% center = info0{3};
% cluster2 = info0{4};
% center2 = info0{5};
% center2_load = info0{6} ;
% dotx = [supply(1, 1); dot(:, 1)];
% doty = [supply(1, 2); dot(:, 2)];
% 
% hold on;
% scatter(x, y, 20, Color(3,:), "filled");
% hold on;
% % scatter(center(:, 1), center(:, 2), 200, "diamond");
% hold on;
% scatter(supply(1), supply(2), 700, Color(1, :), "filled");
% hold on;
% line([dotx(edge(:, 1)); dotx(edge(:, 2))], [doty(edge(:, 1)); doty(edge(:, 2))], 'Color', color(1, :), 'LineWidth', 4)
% hold on;
% for i = 1:cluster_num  % 对第一层的每个聚类中心
%     % 将第一层聚类的每个簇提取出来
%     cluster{i} = [x(idx == idxidx(i))', y(idx == idxidx(i))'];
% 
%     % 根据每簇的点个数 进行 第二层聚类
%     clear SilhouetteCoefficient;
%     for n = 1:ceil(center_load(idxidx(i) ) / 2 + 1) 
%         [idx2{i}, center2{i}] = kmeans(cluster{i}, n, 'Replicates', 200);
%         SilhouetteCoefficient(n) = mean(silhouette(cluster{i}, idx2{i}) );
%         if center_load(idxidx(i) ) == 1
%             break;
%         end
%     end
%     [~, cluster2_num] = max(SilhouetteCoefficient);
%     [idx2{i}, center2{i}] = kmeans(cluster{i}, cluster2_num, 'Replicates', 500);
%     center2_load{i} = count_num(idx2{i}, cluster2_num);
%     cluster2{i} = cell(1, cluster2_num);
%     for t = 1:cluster2_num
%         cluster2{i}{t} = cluster{i}(idx2{i} == t, :);
%     end
% 
% 
%     % 计算费用
%     % 计算一级支线费用
%     d = zeros(1, size(edge, 1));
%     for j = 1:size(edge, 1)
%         % 二级分叉点到总线每段的距离
%         [d(1, j), nearest(j, :)] = distance_fun(point(i, :), dot(edge(j, 1), :), dot(edge(j, 2), :) );
%     end
%     % 二级分叉点到总线的最短距离
%     [distance, min_index] = min(d);
%     hold on;
%     line([point(i, 1), nearest(min_index, 1)], [point(i, 2), nearest(min_index, 2)], 'Color', Color(4, :), 'LineWidth', 3);
% 
%     % 计算二级支线费用
%     for j = 1:cluster2_num % 对第二层聚类的每簇
%         % 第i簇里面的第j簇的负荷数
%         if center2_load{i}(j) >= 3
%             weight = 239.4; % 支线B
%         else
%             weight = 188.6; % 支线A
%         end
%         cost_2 = cost_2 + weight * dist(point(i, :), center2{i}(j, :)');
%         hold on;
%         line([point(i, 1), center2{i}(j, 1)], [point(i, 2), center2{i}(j, 2)], 'Color', 'b', 'LineWidth', 2);
%         for k = 1:center2_load{i}(j) % 对第二层每簇的每个点
%             box_tmp = cluster{i}(idx2{i} == j, :);
%             hold on; % 画图，连线
%             line([box_tmp(k, 1), center2{i}(j, 1)], [box_tmp(k, 2), center2{i}(j, 2)], 'Color', color(4, :), 'LineWidth', 1);
%         end
% 
%     end
% 
% end
% 
% 
% 
% dot(1:20, :) = info1{1};
% cluster = info1{2};
% center = info1{3};
% cluster2 = info1{4};
% center2 = info1{5};
% center2_load = info1{6} ;
% dotx = [supply(2, 1); dot(:, 1)];
% doty = [supply(2, 2); dot(:, 2)];
% 
% hold on;
% scatter(x, y, 20, Color(3,:), "filled");
% hold on;
% % scatter(center(:, 1), center(:, 2), 200, "diamond");
% hold on;
% scatter(supply(1), supply(2), 700, Color(1, :), "filled");
% hold on;
% line([dotx(edge(:, 1)); dotx(edge(:, 2))], [doty(edge(:, 1)); doty(edge(:, 2))], 'Color', color(1, :), 'LineWidth', 4)
% hold on;
% for i = 1:cluster_num  % 对第一层的每个聚类中心
%     % 将第一层聚类的每个簇提取出来
%     cluster{i} = [x(idx == idxidx(i))', y(idx == idxidx(i))'];
% 
%     % 根据每簇的点个数 进行 第二层聚类
%     clear SilhouetteCoefficient;
%     for n = 1:ceil(center_load(idxidx(i) ) / 2 + 1) 
%         [idx2{i}, center2{i}] = kmeans(cluster{i}, n, 'Replicates', 200);
%         SilhouetteCoefficient(n) = mean(silhouette(cluster{i}, idx2{i}) );
%         if center_load(idxidx(i) ) == 1
%             break;
%         end
%     end
%     [~, cluster2_num] = max(SilhouetteCoefficient);
%     [idx2{i}, center2{i}] = kmeans(cluster{i}, cluster2_num, 'Replicates', 500);
%     center2_load{i} = count_num(idx2{i}, cluster2_num);
%     cluster2{i} = cell(1, cluster2_num);
%     for t = 1:cluster2_num
%         cluster2{i}{t} = cluster{i}(idx2{i} == t, :);
%     end
% 
% 
%     % 计算费用
%     % 计算一级支线费用
%     d = zeros(1, size(edge, 1));
%     for j = 1:size(edge, 1)
%         % 二级分叉点到总线每段的距离
%         [d(1, j), nearest(j, :)] = distance_fun(point(i, :), dot(edge(j, 1), :), dot(edge(j, 2), :) );
%     end
%     % 二级分叉点到总线的最短距离
%     [distance, min_index] = min(d);
%     hold on;
%     line([point(i, 1), nearest(min_index, 1)], [point(i, 2), nearest(min_index, 2)], 'Color', Color(4, :), 'LineWidth', 3);
% 
%     % 计算二级支线费用
%     for j = 1:cluster2_num % 对第二层聚类的每簇
%         % 第i簇里面的第j簇的负荷数
%         if center2_load{i}(j) >= 3
%             weight = 239.4; % 支线B
%         else
%             weight = 188.6; % 支线A
%         end
%         cost_2 = cost_2 + weight * dist(point(i, :), center2{i}(j, :)');
%         hold on;
%         line([point(i, 1), center2{i}(j, 1)], [point(i, 2), center2{i}(j, 2)], 'Color', 'b', 'LineWidth', 2);
%         for k = 1:center2_load{i}(j) % 对第二层每簇的每个点
%             box_tmp = cluster{i}(idx2{i} == j, :);
%             hold on; % 画图，连线
%             line([box_tmp(k, 1), center2{i}(j, 1)], [box_tmp(k, 2), center2{i}(j, 2)], 'Color', color(4, :), 'LineWidth', 1);
%         end
% 
%     end
% 
% end

% % % 第三问 一条线的情况
% [val, min_reliability, len1, len2, len3, cost, contact, h] = reliability_fun_3(supply, info0, info1, 324, -1, 3, [1 0 0], 1);
% [val, min_reliability, len1, len2, len3, cost, contact, h] = reliability_fun_3(supply, info0, info1, 400, -1, 3, [1 0 0], 1);
% [val, min_reliability, len1, len2, len3, cost, contact, h] = reliability_fun_3(supply, info0, info1, 500, -1, 3, [1 0 0], 1);
% title("单联络线双供电网拓扑图");
% tmp = 0.5;
% for i = 1:68
%     tmp = cat(2, tmp, 0.5);
% end
% w = 1.5;
% tmp(1) = w; tmp(2) = w; tmp(3) = w; tmp(5) = w; tmp(6) = w; tmp(8) = w; tmp(9) = w;
% tmp(11) = w; tmp(14) = w; tmp(15) = w; tmp(16) = w; tmp(18) = w; tmp(19) = w; tmp(21) = w;
% tmp(22) = w; tmp(24) = w; tmp(25) = w; tmp(27) = w; tmp(12) = 4;
% h.LineWidth = tmp;
% tmp = "-";
% for i = 1:68
%     tmp = cat(2, tmp, "-");
% end
% tmp(12) = ":";
% h.LineStyle = tmp;
% tmp = [0 0 0];
% for i = 1:68
%     tmp = cat(1, tmp, [0 0 0]);
% end
% tmp(12, :) = [1 0 0];
% h.EdgeColor = tmp;
% 
% % % 第三问 两条线的情况
% % [val, min_reliability, len1, len2, len3, cost, contact, h] = reliability_fun_3(supply, info0, info1, 723, -1, 3, [1 0 0], 1);
% % [val, min_reliability, len1, len2, len3, cost, contact, h] = reliability_fun_3(supply, info0, info1, 765, -1, 3, [1 0 0], 1);
% % % 第三问 三条线的情况
% % [val, min_reliability, len1, len2, len3, cost, contact, h] = reliability_fun_3(supply, info0, info1, 1164, -1, 3, [1 0 0], 1);
% % % 第四问 一条线的情况
% % [val, min_reliability, len, cost, contact, h] = reliability_fun_1(supply, info0, info1, 10000, 0.979, 4, [1 0 0], 1);
% [val, min_reliability, len, cost, contact, h] = reliability_fun_1(supply, info0, info1, 10000, 0.9812528, 4, [1 0 0], 1);
% % % 第四问 两条线的情况
% [val, min_reliability, len, cost, contact, h] = reliability_fun_1(supply, info0, info1, 10000, 0.9812818, 4, [1 0 0], 1);
% % % 第四问 三条线的情况
[val, min_reliability, len, cost, contact, h] = reliability_fun_1(supply, info0, info1, 10000, 0.9812821, 4, [1 0 0], 1);
% 
% fprintf("最终总建设费用为%f千元\n", cost + total_cost(4) + 16.5 * (capacity0 + capacity1) );

% line([-25.62, -25.05], [-12.73, -13.59], 'Color', [1 0 0], 'LineWidth', 4);
% disp("不加联络线的情况下最低可靠性为0.979349")
% disp("一条联络线无法达到可靠性要求")
% disp("正在尝试连接两条联络线")
% disp("两条联络线的情况下满足了可靠性要求")
% disp("最低费用为764.176280")
