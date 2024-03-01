function [newpop] = crossover(pop, fitvalue) 
% 交叉，更好的个体应该优先享有繁衍的权利
% 交叉规则：1：更好的个体参与交叉的概率更大，
% 2：对两个参与交叉的个体，随机选择对应的一段染色体片段，对这段上的每一位对应的数字做出如下操作：
% 个体A的某个染色体数字a 对应着 个体B的某个染色体数字b
% 更新a和b 都等于 mean(a + b) + 1.5 * abs(a - b) * (rand - 0.5)
% 其中rand为[0, 1]中均匀分布的随机数

[sx, sy, sz] = size(pop);

newpop = ones(sx, sy, sz);
for i = 1:2:sz
    a = ceil(sz * rand(1) );
    b = ceil(sz * rand(1) );

    for j = 1:sy
        % 交叉的染色体长度
        cross_length = ceil(sx / 2 * rand(1) );
        % 交叉的染色体片段的起始索引
        cross_index = ceil(rand(1) * (sx - cross_length + 1) );
        % 两个个体在 交叉的染色体片段上的 每个染色体的均值
        mean_value = (pop(cross_index:cross_index+cross_length-1, j, a) + pop(cross_index:cross_index+cross_length-1, j, b) ) / 2;
        % 两个个体 交叉的染色体片段上的 每个染色体的变化大小
        value_length = abs(pop(cross_index:cross_index+cross_length-1, j, a) - pop(cross_index:cross_index+cross_length-1, j, b) );
        % 不交叉的位置直接复制
        newpop(:, :, i) = pop(:, :, a);
        newpop(:, :, i+1) = pop(:, :, b);
        % 对交叉方式进行试验
%         newpop(cross_index:cross_index+cross_length-1, j, i) = pop(cross_index:cross_index+cross_length-1, j, b);
%         newpop(cross_index:cross_index+cross_length-1, j, i+1) = pop(cross_index:cross_index+cross_length-1, j, a);
        % 交叉得出新染色体片段
        for k = cross_index:cross_index+cross_length-1
            newpop(k, j, i) = 1.5 * value_length(k - cross_index + 1) * (rand(1) - 0.5) + mean_value(k - cross_index + 1);
            newpop(k, j, i+1) = 1.5 * value_length(k - cross_index + 1) * (rand(1) - 0.5) + mean_value(k - cross_index + 1);
        end
    end
end

return;