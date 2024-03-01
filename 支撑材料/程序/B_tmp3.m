clear; clc;
global center2; global center2_load; global edge;
global i;  global dott; global center_load; global idx;  global x;  global y; 

%% 

% 读取/生成数据部分
data = readmatrix("电网数据.xlsx");
supply = data(1, 1:2);
x = data(4:size(data, 1), 1)'; % 点的x坐标
y = data(4:size(data, 1), 2)'; % 点的y坐标
len = length(x);

% 这是随机生成负荷坐标的代码
% supply = [floor(1000 * rand(1)), floor(750 * rand(1))];
% x = floor(1000 * rand(1, len)); % 点的x坐标
% y = floor(750 * rand(1, len)); % 点的y坐标


%% 

% 第一层聚类部分, idx为各负荷所属簇的下标，center为各簇中心坐标
cluster_num = 7; % 聚类簇数k
[idx, center] = kmeans(cat(1, x, y)', cluster_num, 'Replicates', 1000); % 聚类1000次，选SSD最小的结果
center_load = count_num(idx, cluster_num); % 每簇对应的负荷数
% gscatter(x, y, idx); % 根据聚类结果，用不同颜色画出各点


%% 

% 画出 二级分叉点和电源 的外接矩形，生成均匀点集
box_tmp = cat(1, center, supply);
spacing = 0.01 * min( max(box_tmp) - min(box_tmp) ); % 生成点的间距为外接矩形的宽度的1%
[boxx, boxy] = bounding_box(box_tmp(:, 1), box_tmp(:, 2), spacing); % 画出外接矩形，生成均匀点集
clear box_tmp;
% hold on;
% scatter(boxx, boxy, 2);


%% 

times = 1;
color = [0,0,0; 0,0,255; 0,255,0; 255,0,0] / 255;
for time = 1:times % 双层规划次数
    % 遗传算法
    epoch = 400;
    pop_size = 30; % 群体大小
    individual_size = 20; % 个体长度
    pm = 0.1; % 变异概率
    
    if time == 1
        [best_individual, fit, run_time1, cooo] = Genetic_algorithm_tmp(epoch, pop_size, individual_size, pm, supply, center, center_load);
    end

%     % 这是不用模拟退火，单用遗传的代码
%     [dot, fit, run_time1] = Genetic_algorithm(1000, pop_size, individual_size, pm, supply, center, center_load);


    % 进行模拟退火
    temterature = 5;  % 初始温度
    innerloop = 100;    % 马尔科夫链长 （内层循环次数）
    dc = 0.99; % 退温系数Dewarming_coefficient
    spacing = 0.01 * min( max(center) - min(center) ); % 扰动步长为外接矩形的宽度的1%

    % 将遗传算法的最优解作为初始点的位置
    if time == 1
        dot_index = dsearchn([boxx', boxy'], best_individual);
        dot = [boxx(dot_index)', boxy(dot_index)'];
        previous_cost = cost_fun(cat(1, supply, dot), center, center_load);
    end

    [dot, cost, run_time2] = Simulated_annealing_tmp(temterature, innerloop, dc, dot, supply, center, center_load, spacing, boxx, boxy, cooo);

    close all;
    if time == 1
        figure;
        plot(run_time1, fit);
        hold on;
        plot(max(run_time1)+run_time2, cost);
        hold on;
        line([max(run_time1), max(run_time1)+run_time2(1)], [fit(epoch), cost(1)], 'Color', [1, 0, 0]);
        title("遗传模拟退火算法下费用的变动情况");
        xlabel('运行时间 / s');
        ylabel('费用 / 千元');
    end

    dott = dot;
    % 第二层规划
    cluster = cell(1, cluster_num); % 第一次聚类中的簇
    idx2 = cell(1, cluster_num); % 第二次聚类中各点的归属
    cluster2 = cell(1, cluster_num); % 第二次聚类中各簇
    center2 = cell(1, cluster_num); % 第二次聚类中各中心
    center2_load = cell(1, cluster_num); % 第二次聚类中各簇点数
    point = zeros(cluster_num, 2); % 优化后的二级分叉点坐标
    nearest = zeros(size(edge, 1), 2);

    cost_1 = 0; cost_2 = 0; cost_3 = 0; cost_main = 0;
    cost_4 = 2.6 * len + 56.8 * cluster_num;

    for i = 1:cluster_num  % 对第一层的每个聚类中心
        % 将第一层聚类的每个簇提取出来
        cluster{i} = [x(idx == i)', y(idx == i)'];

        % 根据每簇的点个数 进行 第二层聚类
        clear SilhouetteCoefficient;
        for n = 1:ceil(center_load(i) / 2 + 1)
            [idx2{i}, center2{i}] = kmeans(cluster{i}, n, 'Replicates', 200);
            SilhouetteCoefficient(n) = mean(silhouette(cluster{i}, idx2{i}) );
            if center_load(i) == 1
                break;
            end
        end
        [~, cluster2_num] = max(SilhouetteCoefficient);
        [idx2{i}, center2{i}] = kmeans(cluster{i}, cluster2_num, 'Replicates', 500);
        center2_load{i} = count_num(idx2{i}, cluster2_num);
        cluster2{i} = cell(1, cluster2_num);
        for t = 1:cluster2_num
            cluster2{i}{t} = cluster{i}(idx2{i} == t, :);
        end


        hold on;
        %     scatter(center2{i}(:, 1), center2{i}(:, 2), 100);


        % 以 一级支线的费用和二级支线的费用最小为目标 优化二级分叉点
        point(i, :) = fminsearch(@optimize_fun, rand(1, 2));
        for n = 1:99
            point(i, :) = point(i, :) + fminsearch(@optimize_fun, rand(1, 2));
        end
        point(i, :) = point(i, :) / 100;

        % 计算费用
        % 计算一级支线费用
        d = zeros(1, size(edge, 1));
        for j = 1:size(edge, 1)
            % 二级分叉点到总线每段的距离
            [d(1, j), nearest(j, :)] = distance_fun(point(i, :), dot(edge(j, 1), :), dot(edge(j, 2), :) );
        end
        % 二级分叉点到总线的最短距离
        [distance, min_index] = min(d);        
        if time == times
            hold on;
%             line([point(i, 1), nearest(min_index, 1)], [point(i, 2), nearest(min_index, 2)], 'Color', color(2, :), 'LineWidth', 3);
        end

        if center_load(i) >= 3
            weight = 239.4; % 支线B
        else
            weight = 188.6; % 支线A
        end
        % 加上一级支线的价格
        cost_1 = cost_1 + weight * distance;

        % 计算二级支线费用
        for j = 1:cluster2_num % 对第二层聚类的每簇
            % 第i簇里面的第j簇的负荷数
            if center2_load{i}(j) >= 3
                weight = 239.4; % 支线B
            else
                weight = 188.6; % 支线A
            end
            cost_2 = cost_2 + weight * dist(point(i, :), center2{i}(j, :)');
            if time == times
                hold on;
%                 line([point(i, 1), center2{i}(j, 1)], [point(i, 2), center2{i}(j, 2)], 'Color', color(3, :), 'LineWidth', 2);
            end

            % 计算三级支线费用
            weight = 188.6; % 支线A
            distance = dist(cluster{i}(idx2{i} == j, :), center2{i}(j, :)');
            cost_3 = cost_3 + weight * sum(distance);

            for k = 1:center2_load{i}(j) % 对第二层每簇的每个点
                box_tmp = cluster{i}(idx2{i} == j, :);
                if time == times
                    hold on; % 画图，连线
%                     line([box_tmp(k, 1), center2{i}(j, 1)], [box_tmp(k, 2), center2{i}(j, 2)], 'Color', color(4, :), 'LineWidth', 1);
                end
            end

        end

    end

    center = point;
%     [~, ~, main_line_length] = cost_fun(dot, center, center_load);
    cost_main = 325.7 * main_line_length;
    fprintf("最终花费%f千元\n", cost_1 + cost_2 + cost_3 + cost_4 + cost_main);
    total_cost(time) = cost_1 + cost_2 + cost_3 + cost_4 + cost_main;

    if time == times
        hold on;
        scatter(point(:, 1), point(:, 2), 500, '*');
        axis([-1.75, 2.25, -1.5, 2.5]);
    end

end
