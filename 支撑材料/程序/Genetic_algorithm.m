function [best_individual, fit, run_time] = Genetic_algorithm(epoch, pop_size, individual_size, pm, supply, center, center_load)

box_tmp = cat(1, center, supply);
xrange = max(box_tmp) - min(box_tmp);
yrange = xrange(2);
xrange = xrange(1);
minimun = min(box_tmp);

pop(:, 1, :) = xrange * rand([individual_size, 1, pop_size]) + minimun(1); % 随机产生初始群体
pop(:, 2, :) = yrange * rand([individual_size, 1, pop_size]) + minimun(2); % 随机产生初始群体

tmp(:, :, 1) = supply;
for i = 1:2*pop_size-1
    tmp = cat(3, tmp, supply); % 计算适应度时要带上电源
end
fitvalue = fitvalue_fun(cat(1, tmp(:, :, 1:pop_size), pop), center, center_load);

tic;
for i = 1:epoch % 遗传代数

    newpop = crossover(pop, fitvalue); % 交叉
    newpop = mutation(newpop, pm, 0.05 * min(xrange, yrange) ); % 变异
    newpop = cat(3, pop, newpop);
    
    fitvalue = fitvalue_fun(cat(1, tmp, newpop), center, center_load); % 计算群体中每个个体的适应度
    [pop, fitvalue] = selection(newpop, fitvalue); % 自然选择

    [~, bestfit] = best(pop, fitvalue); % 求出群体中最优的个体及其适应值
    fit(i) = bestfit;
    run_time(i) = toc;
    if mod(i, 50) == 0
        fprintf("---遗传第%d轮迭代---\n", i);
        fprintf(" 目前花费为%.2f\n", fit(i) );
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