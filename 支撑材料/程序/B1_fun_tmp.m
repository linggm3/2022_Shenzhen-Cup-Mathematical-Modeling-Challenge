function [cost, info] = B1_fun_tmp(x, y, supply, cluster_num, idx, center, idxidx, power, ttt)

global center2; global center2_load; global edge;
global i;  global dott; global center_load;
len = size(x, 1);
Color = [250 192 15; 1 86 153; 243 118 74; 95 198 201; 79 89 100] / 255;


% 画出外接矩形，生成均匀点集
box_tmp = cat(1, center, supply);
spacing = 0.01 * min( max(box_tmp) - min(box_tmp) ); % 生成点的间距为外接矩形的宽度的1%
[boxx, boxy] = bounding_box(box_tmp(:, 1), box_tmp(:, 2), spacing); % 画出外接矩形，生成均匀点集


color = [0 0 0; 0 0 255; 0 255 0; 95 198 201] / 255;
for kkk = 1:ttt
    % 遗传算法
    epoch = 200;
    pop_size = 30; % 群体大小
    individual_size = 20; % 个体长度
    pm = 0.1; % 变异概率
    if kkk == 1
        [best_individual, ~] = Genetic_algorithm(epoch, pop_size, individual_size, pm, supply, center, center_load);
    end


    % 进行模拟退火
    temterature = 5;  % 初始温度
    innerloop = 100;    % 马尔科夫链长 （内层循环次数）
    dc = 0.99; % 退温系数Dewarming_coefficient
    spacing = 0.01 * min( max(center) - min(center) ); % 扰动步长为外接矩形的宽度的1%

    % 将遗传算法的最优解作为初始点的位置
    if kkk == 1
        dot_index = dsearchn([boxx', boxy'], best_individual);
        dot = [boxx(dot_index)', boxy(dot_index)'];
        previous_cost = cost_fun(cat(1, supply, dot), center, center_load);
    end

    [dot, ~] = Simulated_annealing(temterature, innerloop, dc, dot, supply, center, center_load, spacing, boxx, boxy);



    % 画图检查结果
    hold on;
    % scatter(center(:, 1), center(:, 2), 125); % 画出各簇中心
    % hold on;
    % scatter(supply(1), supply(2), 200);

    dot = cat(1, dot, supply);
    dotx = dot(:, 1)';
    doty = dot(:, 2)';
    [~, edge] = cost_fun(dot, center, center_load);

    % hold on;
    % scatter(dotx, doty, 50);
    % hold on;
    % line([dotx(edge(:, 1)); dotx(edge(:, 2))], [doty(edge(:, 1)); doty(edge(:, 2))] );


    % 画图
    if kkk == ttt
        hold on;
        scatter(x, y, 20, Color(3,:), "filled");
        hold on;
        % scatter(center(:, 1), center(:, 2), 200, "diamond");
        hold on;
        scatter(supply(1), supply(2), 700, Color(1, :), "filled");
        hold on;
        line([dotx(edge(:, 1)); dotx(edge(:, 2))], [doty(edge(:, 1)); doty(edge(:, 2))], 'Color', color(1, :), 'LineWidth', 4)
        hold on;
    end


    dott = dot;
    % 第二层规划
    cluster = cell(1, cluster_num); % 第一次聚类中的簇
    idx2 = cell(1, cluster_num); % 第二次聚类中各点的归属
    center2 = cell(1, cluster_num); % 第二次聚类中各中心
    center2_load = cell(1, cluster_num); % 第二次聚类中各簇点数
    cluster2 = cell(1, cluster_num); % 第二次聚类中各簇
    point = zeros(cluster_num, 2); % 优化后的二级分叉点坐标
    nearest = zeros(size(edge, 1), 2);

    cost_1 = 0; cost_2 = 0; cost_3 = 0; cost_main = 0;
    cost_4 = 2.6 * len + 56.8 * cluster_num;

    for i = 1:cluster_num  % 对第一层的每个聚类中心
        % 将第一层聚类的每个簇提取出来
        cluster{i} = [x(idx == idxidx(i))', y(idx == idxidx(i))'];

        % 根据每簇的点个数 进行 第二层聚类
        clear SilhouetteCoefficient;
        for n = 1:ceil(center_load(idxidx(i) ) / 2 + 1) 
            [idx2{i}, center2{i}] = kmeans(cluster{i}, n, 'Replicates', 200);
            SilhouetteCoefficient(n) = mean(silhouette(cluster{i}, idx2{i}) );
            if center_load(idxidx(i) ) == 1
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
        hold on;
        if kkk == ttt
            line([point(i, 1), nearest(min_index, 1)], [point(i, 2), nearest(min_index, 2)], 'Color', Color(4, :), 'LineWidth', 3);
        end

        if center_load(idxidx(i) ) >= 3
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
            hold on;
            if kkk == ttt
%                 line([point(i, 1), center2{i}(j, 1)], [point(i, 2), center2{i}(j, 2)], 'Color', color(3, :), 'LineWidth', 2);
            end

            % 计算三级支线费用
            weight = 188.6; % 支线A
            distance = dist(cluster{i}(idx2{i} == j, :), center2{i}(j, :)');
            cost_3 = cost_3 + weight * sum(distance);

            for k = 1:center2_load{i}(j) % 对第二层每簇的每个点
                box_tmp = cluster{i}(idx2{i} == j, :);
                hold on; % 画图，连线
                if kkk == ttt
%                     line([box_tmp(k, 1), center2{i}(j, 1)], [box_tmp(k, 2), center2{i}(j, 2)], 'Color', color(4, :), 'LineWidth', 1);
                end
            end

        end

    end

    center = point;
    [~, ~, main_line_length] = cost_fun(dot, center, center_load);
    cost_main = 325.7 * main_line_length;
%     fprintf("最终花费%f千元\n", cost_1 + cost_2 + cost_3 + cost_4 + cost_main);
    cost = cost_1 + cost_2 + cost_3 + cost_4 + cost_main;
    
    hold on;
    scatter(center(:, 1), center(:, 2), 150, Color(2, :), "filled");
    axis([-1.75, 2.25, -1.5, 2.5]);

end


info = cell(1, 7);
info{1} = dot(1:individual_size, :);
info{2} = cluster;
info{3} = center;
info{4} = cluster2;
info{5} = center2;
info{6} = center2_load;
for i = 1:size(info{2}, 2)
    for j = 1:size(info{2}{i}, 1)
        for k = 1:size(x, 2)
            if x(k) == info{2}{i}(j, 1) && y(k) == info{2}{i}(j, 2)
                break;
            end
        end
        info{7}{i}(j) = power(k);
    end
end