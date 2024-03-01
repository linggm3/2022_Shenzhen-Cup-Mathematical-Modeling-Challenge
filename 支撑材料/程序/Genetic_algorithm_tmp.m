function [best_individual, fit, run_time, counter] = Genetic_algorithm_tmp(epochs, pop_size, individual_size, pm, supply, center, center_load, ttt)
global idx; global x;  global y;  global cluster_num;
counter = 1
color = [0,0,0; 0,0,255; 0,255,0; 255,0,0] / 255;
Color = [250 192 15; 1 86 153; 243 118 74; 95 198 201; 79 89 100] / 255;
box_tmp = cat(1, center, supply);
xrange = max(box_tmp) - min(box_tmp);
yrange = xrange(2);
xrange = xrange(1);
minimun = min(box_tmp);
cost_4 = 2.6 * 35 + 56.8 * cluster_num;
mkdir("D:\picture\pic" + num2str(ttt)); % 创建文件夹

pop(:, 1, :) = xrange * rand([individual_size, 1, pop_size]) + minimun(1); % 随机产生初始群体
pop(:, 2, :) = yrange * rand([individual_size, 1, pop_size]) + minimun(2); % 随机产生初始群体

tmp(:, :, 1) = supply;
for epoch = 1:2*pop_size-1
    tmp = cat(3, tmp, supply); % 计算适应度时要带上电源
end
fitvalue = fitvalue_fun(cat(1, tmp(:, :, 1:pop_size), pop), center, center_load);

% 
% scatter(x, y, 20, Color(3,:), "filled");
% hold on;
% axis([-1.50, 2.50, -1.5, 2.5]);
% scatter(supply(1), supply(2), 700, Color(1, :), "filled");
% axis([-1.50, 2.50, -1.5, 2.5]);
% title("双层规划效果图");
% [~, lgd] = legend("负荷", "电源");
% disp(lgd)
% lgd(3).Children.MarkerSize = 6;
% lgd(4).Children.MarkerSize = 12.5;
% h1 = gcf;
% saveas(h1, "D:\picture\pic" + num2str(ttt) + "\picture" + num2str(counter) + ".jpg");
% counter = counter + 1;
% 
% figure;
% scatter(x, y, 20, Color(3,:), "filled");
% hold on;
% scatter(center(:, 1), center(:, 2), 150, Color(2, :), "filled");
% axis([-1.50, 2.50, -1.5, 2.5]);
% scatter(supply(1), supply(2), 700, Color(1, :), "filled");
% title("双层规划效果图");
% [~, lgd] = legend("负荷", "二级分叉点", "电源");
% disp(lgd)
% lgd(4).Children.MarkerSize = 6;
% lgd(5).Children.MarkerSize = 8;
% lgd(6).Children.MarkerSize = 12.5;
% h1 = gcf;
% saveas(h1, "D:\picture\pic" + num2str(ttt) + "\picture" + num2str(counter) + ".jpg");
% counter = counter + 1;
% 
% 
% figure;
% [dot] = best(pop, fitvalue); % 求出群体中最优的个体及其适应值
% title("双层规划效果图");
% hold on;
% scatter(dot(:, 1), dot(:, 2), '*');
% axis([-1.50, 2.50, -1.5, 2.5]);
% scatter(x, y, 20, Color(3,:), "filled");
% hold on;
% scatter(center(:, 1), center(:, 2), 150, Color(2, :), "filled");
% axis([-1.50, 2.50, -1.5, 2.5]);
% scatter(supply(1), supply(2), 700, Color(1, :), "filled");
% title("双层规划效果图");
% [~, lgd] = legend("主线端点", "负荷", "聚类中心点", "电源");
% disp(lgd)
% lgd(5).Children.MarkerSize = 6;
% lgd(6).Children.MarkerSize = 6;
% lgd(7).Children.MarkerSize = 10;
% lgd(8).Children.MarkerSize = 12.5;
% h1 = gcf;
% saveas(h1, "D:\picture\pic" + num2str(ttt) + "\picture" + num2str(counter) + ".jpg");
% counter = counter + 1;
% 
% 
% dot = cat(1, dot, supply);
% dotx = dot(:, 1)';
% doty = dot(:, 2)';
% [~, edge] = minspan(dot);
% hold on;
% title("双层规划效果图");
% hold on;
% line([dotx(edge(:, 1)); dotx(edge(:, 2))], [doty(edge(:, 1)); doty(edge(:, 2))], 'Color', color(1, :), 'LineWidth', 4)
% hold on;
% h1 = gcf;
% saveas(h1, "D:\picture\pic" + num2str(ttt) + "\picture" + num2str(counter) + ".jpg");
% counter = counter + 1;
% 
% 
% dot = cat(1, dot, supply);
% dotx = dot(:, 1)';
% doty = dot(:, 2)';
% [~, edge] = minspan(dot);
% hold on;
% title("双层规划效果图");
% hold on;
% scatter(center(:, 1), center(:, 2), 150, Color(2, :), "filled");
% hold on;
% scatter(supply(1), supply(2), 700, Color(1, :), "filled");
% hold on;
% line([dotx(edge(:, 1)); dotx(edge(:, 2))], [doty(edge(:, 1)); doty(edge(:, 2))], 'Color', color(1, :), 'LineWidth', 4)
% hold on;

%  % 计算费用
% % 计算一级支线费用
% d = zeros(1, size(edge, 1));
% for ii = 1:7
%     for j = 1:size(edge, 1)
%         % 二级分叉点到总线每段的距离
%         [d(1, j), nearest(j, :)] = distance_fun(center(ii, :), dot(edge(j, 1), :), dot(edge(j, 2), :) );
%     end
%     [~, pos] = min(d);
%     res(ii, :) = nearest(pos, :);
%     hold on;
%     line([res(ii, 1); center(ii, 1)], [res(ii, 2); center(ii, 2)], 'Color', color(2, :), 'LineWidth', 3)
% end
% 
% cost_1 = 0; cost_2 = 0; cost_3 = 0; 
% for i = 1:7  % 对第一层的每个聚类中心
%     % 将第一层聚类的每个簇提取出来
%     cluster{i} = [x(idx == i)', y(idx == i)'];
% 
%     % 根据每簇的点个数 进行 第二层聚类
%     clear SilhouetteCoefficient;
%     for n = 1:ceil(center_load(i) / 2 + 1)
%         [idx2{i}, center2{i}] = kmeans(cluster{i}, n, 'Replicates', 200);
%         SilhouetteCoefficient(n) = mean(silhouette(cluster{i}, idx2{i}) );
%         if center_load(i) == 1
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
%     hold on;
%     %     scatter(center2{i}(:, 1), center2{i}(:, 2), 100);
% 
% 
%     % 以 一级支线的费用和二级支线的费用最小为目标 优化二级分叉点
%     point = center;
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
%     line([point(i, 1), nearest(min_index, 1)], [point(i, 2), nearest(min_index, 2)], 'Color', Color(4, :), 'LineWidth', 3);
% 
%     if center_load(i) >= 3
%         weight = 239.4; % 支线B
%     else
%         weight = 188.6; % 支线A
%     end
%     % 加上一级支线的价格
%     cost_1 = cost_1 + weight * distance;
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
% %         line([point(i, 1), center2{i}(j, 1)], [point(i, 2), center2{i}(j, 2)], 'Color', color(3, :), 'LineWidth', 2);
% 
%         % 计算三级支线费用
%         weight = 188.6; % 支线A
%         distance = dist(cluster{i}(idx2{i} == j, :), center2{i}(j, :)');
%         cost_3 = cost_3 + weight * sum(distance);
% 
%         for k = 1:center2_load{i}(j) % 对第二层每簇的每个点
%             box_tmp = cluster{i}(idx2{i} == j, :);
%             hold on; % 画图，连线
% %             line([box_tmp(k, 1), center2{i}(j, 1)], [box_tmp(k, 2), center2{i}(j, 2)], 'Color', color(4, :), 'LineWidth', 1);
%         end
%     end
% end
%  [~, ~, main_line_length] = cost_fun(dot, center, center_load);
% cost_main = 325.7 * main_line_length;
% fprintf("最终花费%f千元\n", cost_1 + cost_2 + cost_3 + cost_4 + cost_main);
% total_cost = cost_1 + cost_2 + cost_3 + cost_4 + cost_main;
% hold on;
% 
% hold on;
% scatter(x, y, 20, Color(3,:), "filled");
% axis([-1.50, 2.50, -1.5, 2.5]);
% 
% h1=gcf ;
% saveas(h1, "D:\picture\pic" + num2str(ttt) + "\picture" + num2str(counter) + ".jpg");
% counter = counter + 1;



tic;
for epoch = 1:epochs % 遗传代数

    newpop = crossover(pop, fitvalue); % 交叉
    newpop = mutation(newpop, pm, 0.05 * min(xrange, yrange) ); % 变异
    newpop = cat(3, pop, newpop);
    
    fitvalue = fitvalue_fun(cat(1, tmp, newpop), center, center_load); % 计算群体中每个个体的适应度
    [pop, fitvalue] = selection(newpop, fitvalue); % 自然选择

    [dot, bestfit] = best(pop, fitvalue); % 求出群体中最优的个体及其适应值
    fit(epoch) = bestfit;
%     if epoch > 100 && bestfit > 1800 
%         best_individual = dot;
%         return;
%     end
    run_time(epoch) = toc;
    if mod(epoch, 50) == 0
        fprintf("---遗传第%d轮迭代---\n", epoch);
        fprintf(" 目前花费为%.2f\n", fit(epoch) );
    end

    if epoch < 20 && mod(epoch, 10) == 1 || epoch >= 20 && mod(epoch, 80) == 1 
        dot = [cat(1, dot, supply)];
        dotx = dot(:, 1)';
        doty = dot(:, 2)';
        [~, edge] = minspan(dot);
        figure;
        title("双层规划效果图");
        scatter(x, y, 20, Color(3,:), "filled");
        hold on;
        scatter(center(:, 1), center(:, 2), 150, Color(2, :), "filled");
        hold on;
        scatter(supply(1), supply(2), 700, Color(1, :), "filled");
        hold on;
        line([dotx(edge(:, 1)); dotx(edge(:, 2))], [doty(edge(:, 1)); doty(edge(:, 2))], 'Color', color(1, :), 'LineWidth', 4)
        hold on;
        
         % 计算费用
        % 计算一级支线费用
        d = zeros(1, size(edge, 1));
        for ii = 1:cluster_num
            for j = 1:size(edge, 1)
                % 二级分叉点到总线每段的距离
                [d(1, j), nearest(j, :)] = distance_fun(center(ii, :), dot(edge(j, 1), :), dot(edge(j, 2), :) );
            end
            [~, pos] = min(d);
            res(ii, :) = nearest(pos, :);
            hold on;
            line([res(ii, 1); center(ii, 1)], [res(ii, 2); center(ii, 2)], 'Color', color(2, :), 'LineWidth', 3)
        end
        
        cost_1 = 0; cost_2 = 0; cost_3 = 0; 
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
            point = center;
    
            % 计算费用
            % 计算一级支线费用
            d = zeros(1, size(edge, 1));
            for j = 1:size(edge, 1)
                % 二级分叉点到总线每段的距离
                [d(1, j), nearest(j, :)] = distance_fun(point(i, :), dot(edge(j, 1), :), dot(edge(j, 2), :) );
            end
            % 二级分叉点到总线的最短距离
            [distance, min_index] = min(d);        
            line([point(i, 1), nearest(min_index, 1)], [point(i, 2), nearest(min_index, 2)], 'Color', Color(4, :), 'LineWidth', 3);
    
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
                hold on;
%                 line([point(i, 1), center2{i}(j, 1)], [point(i, 2), center2{i}(j, 2)], 'Color', color(3, :), 'LineWidth', 2);
    
                % 计算三级支线费用
                weight = 188.6; % 支线A
                distance = dist(cluster{i}(idx2{i} == j, :), center2{i}(j, :)');
                cost_3 = cost_3 + weight * sum(distance);
    
                for k = 1:center2_load{i}(j) % 对第二层每簇的每个点
                    box_tmp = cluster{i}(idx2{i} == j, :);
                    hold on; % 画图，连线
%                     line([box_tmp(k, 1), center2{i}(j, 1)], [box_tmp(k, 2), center2{i}(j, 2)], 'Color', color(4, :), 'LineWidth', 1);
                end
    
            end
    
        end
    
        [~, ~, main_line_length] = cost_fun(dot, center, center_load);
        cost_main = 325.7 * main_line_length;
        fprintf("最终花费%f千元\n", cost_1 + cost_2 + cost_3 + cost_4 + cost_main);
        total_cost = cost_1 + cost_2 + cost_3 + cost_4 + cost_main;
        hold on;
        text(-26, -13.8, '遗传算法', 'FontSize', 16)
        text(-26, -13.9, '上层费用', 'FontSize', 16)
        text(-26, -14.0, num2str(round(bestfit, 2)), 'FontSize', 16)
        text(-26, -14.1, '总体费用', 'FontSize', 16)
        text(-26, -14.2, num2str(round(total_cost, 2)), 'FontSize', 16)

        hold on;
        
%         axis([-1.50, 2.50, -1.5, 2.5]);


        h1=gcf ;
        saveas(h1, "D:\picture\pic" + num2str(ttt) + "\picture" + num2str(counter) + ".jpg");
        counter = counter + 1;
    end
end
toc;
% 挑出遗传算法结果中最好的解作为模拟退火的初始解
[best_individual, ~] = best(pop, fitvalue);

% figure;
% plot(fit);
% title('遗传算法下解的变动情况');
% xlabel('迭代次数');
% ylabel('费用');
