function [reliability, sub_line_reliability, main_line_reliability, p, length, sorted_nearest_dot, sorted_index] = reliability_fun(option, target, dot, supply, cluster, center, cluster2, center2, center2_load)

main_line = cat(1, dot, supply);  % 总线上的点
[~, endnode] = minspan(main_line);  % 总线
tmp_distance = zeros(1, size(endnode, 1) );  % 某个中心到总线每段的距离
tmp_nearest_dot = zeros(size(endnode, 1), 2);  % 某个中心到总线每段最近的点
distance = zeros(1, size(center, 1) );   % 某个中心到总线的最短距离
nearest_dot = zeros(size(center, 1), 2);  % 总线到某个中心的最近的点

% 计算每个中心到总线的最短距离和最近点（一级分叉点）
for i = 1:size(center, 1)
    for j = 1:size(endnode, 1)
        % 第i个中心到总线各段的最短距离
        [tmp_distance(j), tmp_nearest_dot(j, :)] = distance_fun(center(i, :), main_line(endnode(j, 1), :), main_line(endnode(j, 2), :) );
    end
    [distance(i), min_index] = min(tmp_distance); % 第i个中心到总线的最短距离
    nearest_dot(i, :) = tmp_nearest_dot(min_index, :); % 第i个中心对应的一级分叉点
end

dist_from_supply = dist(nearest_dot, supply');  % 一级分叉点到电源的欧式距离
[~, sorted_index] = sort(dist_from_supply, 'ascend');  % 距离从小到大排序
sorted_nearest_dot = nearest_dot(sorted_index, :);
length = zeros(1, size(sorted_index, 1) ); % 总线各段长度
length(1, 1) = dist(nearest_dot(sorted_index(1), :), supply'); % 第一段是最近的一级分叉点到电源的距离
for i = 2:size(sorted_index, 1)
    % 总线各段长度：各两个一级分叉点间的距离
    length(1, i) = dist(nearest_dot(sorted_index(i-1), :), nearest_dot(sorted_index(i), :)' );
end

possibility = cumprod(1 - 0.002 * length);  % 第i个一级分叉点处故障概率
for i = 1:size(possibility, 2)
    possibility(i) = possibility(i) * (1 - 0.002) ^ i; % 开关故障率0.002
    possibility(i) = possibility(i) * (1 - 0.005); % 电源故障率0.005
end


if option == 0
    for i = 1:size(cluster, 2)
        for j = 1:size(cluster2{i}, 2)
            for k = 1:center2_load{i}(j)
                if target(1) == cluster2{i}{j}(k, 1) && target(2) == cluster2{i}{j}(k, 2) 
                    length2 = dist(target, center2{i}(j, :)' ); % 三级支线
                    length2 = length2 + dist(center2{i}(j, :), center(i, :)' ); % 二级支线
                    length2 = length2 + dist(center(i, :), nearest_dot(i, :)' ); % 一级支线
                    % 负荷可靠性 = 总线部分可靠性 * 支线部分可靠性（开关故障率0.002 用户故障率0.005）
                    main_line_reliability = possibility(sorted_index == i);
                    sub_line_reliability = (1 - 0.002 * length2) * 0.998 * 0.995;
                    reliability = main_line_reliability * sub_line_reliability;
                    p = nearest_dot(i, :); % p是对应的一级分叉点
                    return;
                end
            end
        end
    end
end

if option == 1
    for i = 1:size(center, 1)
        if nearest_dot(sorted_index(i), :) == target
            reliability = possibility(i);
            main_line_reliability = possibility(i);
            sub_line_reliability = 1;
            p = target;
            return;
        end 
    end
end

reliability = -1; % 报错
return;