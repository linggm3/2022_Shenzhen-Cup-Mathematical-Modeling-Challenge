function [dot, cost, run_time] = Simulated_annealing(T, innerloop, dc, dot, supply, center, center_load, spacing, boxx, boxy)
% 模拟退火算法

previous_cost = cost_fun(cat(1, supply, dot), center, center_load);
epoch = 0; % 外层循环计数器

% 开始退火
tic;
while T > 0.01
    for i = 1:innerloop
        % 扰动：
        new_dot = dot_disturb(dot, boxx, boxy, spacing);
        % 计算扰动后的cost
        new_cost = cost_fun(cat(1, supply, new_dot), center, center_load) ;
        delta_cost = new_cost - previous_cost;
        if delta_cost < 0
            dot = new_dot;
            previous_cost = new_cost;
        else
            if exp(-delta_cost/T)>rand()
               dot = new_dot;
               previous_cost = new_cost;
            end
        end
    end
    epoch = epoch + 1;
    
    cost(epoch) = previous_cost;
    run_time(epoch) = toc;
    T = T * dc; % 退温系数
    if mod(epoch, 50) == 0
        fprintf("---退火第%d轮迭代---\n", epoch);
        fprintf("  目前温度为%.2f\n", T);
        fprintf("  目前花费为%.2f\n", previous_cost);
    end
end
toc;

% figure;
% plot(cost);
% title('模拟退火算法下解的变动情况');
% xlabel('迭代次数');
% ylabel('费用');

return;