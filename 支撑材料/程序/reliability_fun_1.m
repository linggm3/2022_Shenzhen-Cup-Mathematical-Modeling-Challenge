function [val, min_reliability, len1, cost, contact, h, avg] = reliability_fun_1(supply, info0, info1, money, request_min_reliability, target, color, invoke_time)
% 计算双供网络，一条联络线情况下的可靠性


% 求电源0的各个一级分叉点
[~, ~, ~, ~, length0, nearest_dot0, sorted_index0] = reliability_fun(0, info0{2}{1}(1, :), info0{1}, supply(1, :), info0{2}, info0{3}, info0{4}, info0{5}, info0{6});

% 求电源1的各个一级分叉点
[~, ~, ~, ~, length1, nearest_dot1, sorted_index1] = reliability_fun(0, info1{2}{1}(1, :), info1{1}, supply(2, :), info1{2}, info1{3}, info1{4}, info1{5}, info1{6});


% 计算可供功率上限
capacity0 = 0; capacity1 = 0;
for k = 1:size(info0{7}, 2)
    capacity0 = capacity0 + sum(info0{7}{k});
end
for k = 1:size(info1{7}, 2)
    capacity1 = capacity1 + sum(info1{7}{k});
end
capacity0 = capacity0 * 1.2 * 0.5;
capacity1 = capacity1 * 1.2 * 0.5;
capacity = round([capacity0, capacity1]);


if invoke_time == 1
    % 换顺序
    info0{2} = info0{2}(sorted_index0);
    info1{2} = info1{2}(sorted_index1);
    info0{3} = info0{3}(sorted_index0, :);
    info1{3} = info1{3}(sorted_index1, :);
    for i = 4:7
        info0{i} = info0{i}(sorted_index0);
        info1{i} = info1{i}(sorted_index1);
    end
end



% 在每两个相邻的一级分叉点间插入点
nearest_dot0 = cat(1, supply(1, :), nearest_dot0);
nearest_dot1 = cat(1, supply(2, :), nearest_dot1);
tmp_var = zeros(2*size(nearest_dot0, 1)-1, 2);
tmp_var(1:2:2*size(nearest_dot0, 1)-1, :) = nearest_dot0;
nearest_dot0 = tmp_var;
for i = 2:2:size(nearest_dot0, 1)-1
    nearest_dot0(i, :) = (nearest_dot0(i-1, :) + nearest_dot0(i+1, :) ) / 2;
end
tmp_var = zeros(2*size(nearest_dot1, 1)-1, 2);
tmp_var(1:2:2*size(nearest_dot1, 1)-1, :) = nearest_dot1;
nearest_dot1 = tmp_var;
for i = 2:2:size(nearest_dot1, 1)-1
    nearest_dot1(i, :) = (nearest_dot1(i-1, :) + nearest_dot1(i+1, :) ) / 2;
end


% 故障发生时决策对哪些负荷供电0
possibility0 = rot90(tril(ones(size(info0{2}, 2) ) ) );  decesion_mat0_size2 = 0;
for m = 1:size(info0{2}, 2) % 簇
    if size(info0{2}{m}, 1) > decesion_mat0_size2
        decesion_mat0_size2 = size(info0{2}{m}, 1);
    end
end   % 故障发生时决策对哪些负荷供电
decision_mat0 = zeros(size(info0{2}, 2), decesion_mat0_size2, size(info0{2}, 2));
candidate0 = cell(1, size(info0{2}, 2) ); % 故障发生时待供电的负荷的功率
for e = 1:size(info0{2}, 2)  % 对每种故障情况

    candidate0{e} = info0{7}{size(info0{7}, 2)}; % 待供电的
    for time = 2:e
        % 待供电的负荷
        candidate0{e} = cat(2, info0{7}{size(info0{2}, 2) - time + 1}, candidate0{e});
    end

    for cluster = 1:size(info0{2}, 2)
        for time = 1:cluster
            % 算出对应的故障概率
            if time == 1 && time == size(info0{2}, 2) - e  + 1
                possibility0(cluster, e) = possibility0(cluster, e) *  (1 - 0.998 * (1 - 0.002 * length0(time) ) * 0.995 );
            elseif time == 1
                possibility0(cluster, e) = possibility0(cluster, e) * 0.998 * (1 - 0.002 * length0(time) ) * 0.995;
            elseif time == size(info0{2}, 2) - e  + 1
                possibility0(cluster, e) = possibility0(cluster, e) *  (1 - 0.998 * (1 - 0.002 * length0(time) ) );
            else
                possibility0(cluster, e) = possibility0(cluster, e) * 0.998 * (1 - 0.002 * length0(time) );
            end
        end
    end

    [max_volumn0(e), choices0{e}] = package_problem2(capacity(1), candidate0{e}, candidate0{e});
    tmp_counter = 1;
    for m = size(info0{2}, 2)-e+1:size(info0{2}, 2)
        for n = 1:size(info0{2}{m}, 1)
            if sum(tmp_counter == choices0{e}) == 1
                decision_mat0(m, n, e) = 1;
            end
            tmp_counter = tmp_counter + 1;
        end
    end
end



% 故障发生时决策对哪些负荷供电1
possibility1 = rot90(tril(ones(size(info1{2}, 2) ) ) );  decesion_mat1_size2 = 0;
for m = 1:size(info1{2}, 2) % 簇
    if size(info1{2}{m}, 1) > decesion_mat1_size2
        decesion_mat1_size2 = size(info1{2}{m}, 1);
    end
end   % 故障发生时决策对哪些负荷供电
decision_mat1 = zeros(size(info1{2}, 2), decesion_mat1_size2, size(info1{2}, 2));
candidate1 = cell(1, size(info1{2}, 2) ); % 故障发生时待供电的负荷的功率
for e = 1:size(info1{2}, 2)  % 对每种故障情况（最近的故障有e种情况）

    candidate1{e} = info1{7}{size(info1{7}, 2)}; % 待供电的
    for time = 2:e
        % 待供电的负荷
        candidate1{e} = cat(2, info1{7}{size(info1{2}, 2) - time + 1}, candidate1{e});
    end

    for cluster = 1:size(info1{2}, 2)
        for time = 1:cluster
            % 算出对应的故障概率
            if time == 1 && time == size(info1{2}, 2) - e  + 1
                possibility1(cluster, e) = possibility1(cluster, e) *  (1 - 0.998 * (1 - 0.002 * length1(time) ) * 0.995 );
            elseif time == 1
                possibility1(cluster, e) = possibility1(cluster, e) * 0.998 * (1 - 0.002 * length1(time) ) * 0.995;
            elseif time == size(info1{2}, 2) - e  + 1
                possibility1(cluster, e) = possibility1(cluster, e) *  (1 - 0.998 * (1 - 0.002 * length1(time) ) );
            else
                possibility1(cluster, e) = possibility1(cluster, e) * 0.998 * (1 - 0.002 * length1(time) );
            end
        end
    end

    [max_volumn1(e), choices1{e}] = package_problem2(capacity(2), candidate1{e}, candidate1{e});
    tmp_counter = 1;
    for m = size(info1{2}, 2)-e+1:size(info1{2}, 2)
        for n = 1:size(info1{2}{m}, 1)
            if sum(tmp_counter == choices1{e}) == 1
                decision_mat1(m, n, e) = 1;
            end
            tmp_counter = tmp_counter + 1;
        end
    end
end



sum_reliability = zeros(1, size(nearest_dot0, 1) * size(nearest_dot1, 1) );
min_reliability = zeros(1, size(nearest_dot0, 1) * size(nearest_dot1, 1) );
[nd0, si0, po0] = reliability_fun0(info0{1}, supply(1, :), info0{3} );
[nd1, si1, po1] = reliability_fun0(info1{1}, supply(2, :), info1{3} );
len1 = 99999 * ones(1, size(nearest_dot0, 1) * size(nearest_dot1, 1) );
% 电源0的第i个一级分叉点 连接到 电源1的第j个一级分叉点
for i = 2:size(nearest_dot0, 1)
    for j = 2:size(nearest_dot1, 1)
        contact_line_length = dist(nearest_dot0(i, :), nearest_dot1(j, :)' );
        main_counter = (i-1) * size(nearest_dot1, 1) + j;
        if 325.7 * contact_line_length  + 56.8 > money
            continue;
        end


        clear reliability0; counter = 1;
        for m = 1:size(info0{2}, 2) % 簇
            for n = 1:size(info0{2}{m}, 1)
                % 本树的主线可靠性和支线可靠性
                [~, sub_line_re0, main_line_re0] = reliability_fun00(0, info0{2}{m}(n, :), info0{2}, info0{3}, info0{4}, info0{5}, info0{6}, nd0, si0, po0);

                % 联络线在另一树的端点的可靠性
                if j == 1
                    contact_line_re0 = 0.995;
                elseif mod(j, 2) == 1
                    contact_line_re0 = reliability_fun00(1, nearest_dot1(j, :), info1{2}, info1{3}, info1{4}, info1{5}, info1{6}, nd1, si1, po1);
                else
                    contact_line_re0 = reliability_fun00(1, nearest_dot1(j+1, :), info1{2}, info1{3}, info1{4}, info1{5}, info1{6}, nd1, si1, po1);
                    contact_line_re0 = contact_line_re0 / (1 - 0.002 * 0.5 * length1(j/2) );
                end
                if contact_line_re0 == -1
                    disp("ERROR"); % 报错
                end


                % 联络线在本树的端点 到 每个一级分叉点 的距离
                tmp1 = 2*m+1;
                tmp2 = i;
                tmp_reliability = 1;
                % line1
                tmp_len = 0;
                for o = min(tmp1, tmp2):max(tmp1, tmp2)-1
                    tmp_len = tmp_len + 0.5 * length0(ceil(o / 2) );
                    if mod(o, 2) == 0
                        tmp_reliability = tmp_reliability * 0.998; % 开关
                        tmp_reliability = tmp_reliability * (1 - 0.002 * tmp_len);
                        tmp_len = 0;
                    elseif o == max(tmp1, tmp2)-1 && mod(o, 2) == 1
                        tmp_reliability = tmp_reliability * 0.998; % 开关
                        tmp_reliability = tmp_reliability * (1 - 0.002 * tmp_len);
                        tmp_len = 0;
                    end
                end

                % 算上联络线本身的可靠性
                contact_line_re0 = contact_line_re0 * 0.998 * (1 - 0.002 * contact_line_length);
                % 算上联络线到一级分叉点的可靠性
                contact_line_re0 = contact_line_re0 * tmp_reliability;


                % 负荷可靠性计算
                tmp_re = 0;
                for e = 1:size(info0{2}, 2)
                    tmp_contact_re = 1;
                    if i >= 2 * (size(info0{2}, 2) - e) + 1 % 故障发生在联络线1之前，可靠性要计算联络线1
                        tmp_contact_re = tmp_contact_re * (1 - contact_line_re0);
                    end
                    tmp_re = tmp_re + decision_mat0(m, n, e) * possibility0(m, e) * (1 - tmp_contact_re );
                end
                tmp_re = tmp_re + main_line_re0;
                reliability0(counter) = sub_line_re0 * tmp_re;
                counter = counter + 1;
            end
        end


        clear reliability1; counter = 1;
        for m = 1:size(info1{2}, 2)
            for n = 1:size(info1{2}{m}, 1)
                % 本树的主线可靠性和支线可靠性
                [~, sub_line_re1, main_line_re1] = reliability_fun00(0, info1{2}{m}(n, :), info1{2}, info1{3}, info1{4}, info1{5}, info1{6}, nd1, si1, po1);

                % 联络线在另一树的端点的可靠性
                if i == 1
                    contact_line_re1 = 0.995;
                elseif mod(i, 2) == 1
                    contact_line_re1 = reliability_fun00(1, nearest_dot0(i, :), info0{2}, info0{3}, info0{4}, info0{5}, info0{6}, nd0, si0, po0);
                else
                    contact_line_re1 = reliability_fun00(1, nearest_dot0(i+1, :), info0{2}, info0{3}, info0{4}, info0{5}, info0{6}, nd0, si0, po0);
                    contact_line_re1 = contact_line_re1 / (1 - 0.002 * 0.5 * length0(i/2) );
                end
                if contact_line_re1 == -1
                    disp("ERROR"); % 报错
                end


                % 联络线在本树的端点 到 每个一级分叉点 的距离
                tmp1 = 2*m+1;
                tmp2 = j;
                tmp_reliability = 1;
                % line1
                tmp_len = 0;
                for o = min(tmp1, tmp2):max(tmp1, tmp2)-1
                    tmp_len = tmp_len + 0.5 * length1(ceil(o / 2) );
                    if mod(o, 2) == 0
                        tmp_reliability = tmp_reliability * 0.998; % 开关
                        tmp_reliability = tmp_reliability * (1 - 0.002 * tmp_len);
                        tmp_len = 0;
                    elseif o == max(tmp1, tmp2)-1 && mod(o, 2) == 1
                        tmp_reliability = tmp_reliability * 0.998; % 开关
                        tmp_reliability = tmp_reliability * (1 - 0.002 * tmp_len);
                        tmp_len = 0;
                    end
                end
                % 算上联络线本身的可靠性
                contact_line_re1 = contact_line_re1 * 0.998 * (1 - 0.002 * contact_line_length);
                % 算上联络线到一级分叉点的可靠性
                contact_line_re1 = contact_line_re1 * tmp_reliability;


                % 负荷可靠性计算
                tmp_re = 0;
                for e = 1:size(info1{2}, 2)
                    tmp_contact_re = 1;
                    if j >= 2 * (size(info1{2}, 2) - e) + 1 % 故障发生在联络线1之前，可靠性要计算联络线1
                        tmp_contact_re = tmp_contact_re * (1 - contact_line_re1);
                    end
                    tmp_re = tmp_re + decision_mat1(m, n, e) * possibility1(m, e) * (1 - tmp_contact_re );
                end
                tmp_re = tmp_re + main_line_re1;
                reliability1(counter) = sub_line_re1 * tmp_re;
                counter = counter + 1;
            end
        end
        sum_reliability(main_counter) = sum(reliability0) + sum(reliability1);
        min_reliability(main_counter) = min(min(reliability0), min(reliability1) );
        len1(main_counter) = dist(nearest_dot0(i, :), nearest_dot1(j, :)' );

    end
end

if sum(min_reliability) == 0
    val = -1; % 费用不够建设任何一条联络线
    cost = -1; h = -1;
    disp("费用不够建设任何一条联络线");  
    return;
end

if target == 3 % 第3问最低可靠性最大的方案
    [val, max_index1] = max(min_reliability); % 最低可靠性最大
    avg = sum_reliability(max_index1) / 35;
    fprintf("一条联络线的情况下最大的最低可靠性是%f\n", val);
    fprintf("此方案用户平均可靠性是%f\n", avg);

    % 得到方案对应的index
    col1 = mod(max_index1, size(nearest_dot1, 1) );
    if col1 == 0
        col1 = size(nearest_dot1, 1);
    end
    row1 = ceil(max_index1 / size(nearest_dot1, 1) );


    len1 = len1(max_index1);
    cost = 56.8 + 325.7 * len1; % 联络线花费


    hold on; % 画图
    line([nearest_dot0(row1, 1), nearest_dot1(col1, 1)], [nearest_dot0(row1, 2), nearest_dot1(col1, 2)], 'Color', color, 'LineWidth', 4);
    contact = [row1, col1];
    figure;
    h = topology_plot2(info0{3}, info0{4}, info0{5}, info1{3}, info1{4}, info1{5}, contact);
    

elseif target == 4 % 第4问费用最低的方案
    while true
        [min_len, min_index] = min(len1); % 联络线最短的方案
        if min_len == 99999 
            % 如果连一条线怎么样都达不到可靠性的要求
            disp("一条联络线无法达到可靠性要求");
            disp("正在尝试连接两条联络线");
            [val, min_reliability, len1, ~, cost, contact, h, avg] = reliability_fun_2(supply, info0, info1, money, request_min_reliability, target, color, invoke_time+1);
            return;
        elseif min_reliability(min_index) < request_min_reliability % 如果不满足最低可靠性的要求
            len1(min_index) = 99999; % 则舍去这个方案
        else
            cost = 56.8 + 325.7 * min_len; % 总花费
            val = min_reliability(min_index); % 本方案最低可靠性
            avg = sum_reliability(min_index) / 35;
            disp("一条联络线的情况下满足了可靠性要求");
            fprintf("最低费用为%f\n", cost);

            % 得到方案对应的index
            col1 = mod(min_index, size(nearest_dot1, 1) );
            if col1 == 0
                col1 = size(nearest_dot1, 1);
            end
            row1 = ceil(min_index / size(nearest_dot1, 1) );
            contact = [row1, col1];


            hold on; % 画图
            line([nearest_dot0(row1, 1), nearest_dot1(col1, 1)], [nearest_dot0(row1, 2), nearest_dot1(col1, 2)], 'Color', color, 'LineWidth', 4);
            figure;
            h = topology_plot2(info0{3}, info0{4}, info0{5}, info1{3}, info1{4}, info1{5}, contact);
            return;
        end
    end
end

return;