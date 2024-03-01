function [reliability, sub_line_reliability, main_line_reliability] = reliability_fun00(option, target, cluster, center, cluster2, center2, center2_load, nearest_dot, sorted_index, possibility)

if option == 0
    for i = 1:size(cluster, 2)
        for j = 1:size(cluster2{i}, 2)
            for k = 1:center2_load{i}(j)
                if target(1) == cluster2{i}{j}(k, 1) && target(2) == cluster2{i}{j}(k, 2) 
                    length2 = dist(target, center2{i}(j, :)' ); % 三级支线
                    length2 = length2 + dist(center2{i}(j, :), center(i, :)' ); % 二级支线
                    length2 = length2 + dist(center(i, :), nearest_dot(i, :)' ); % 一级支线
                    % 负荷可靠性 = 总线部分可靠性 * 支线部分可靠性（开关故障率0.002 用户故障率0.005）
                    main_line_reliability = possibility(i);
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
        if nearest_dot(i, :) == target
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