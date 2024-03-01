function [max_val, choice, f] = package_problem2(capacity, v, w)

% n件物品（重量为w，价值为v），背包容量为m 的 01背包问题
num_of_load = length(v);
f = zeros(1, capacity); % f(i, j) j体积下前i个物品的最大价值
counter = 1;


for i = 1:num_of_load
    for j = 1:capacity
        if i == 1 && j >= v(i)
            f(i, j) = w(i);
            if j == capacity
                choice(counter) = i;
                counter = counter + 1;
            end
        elseif i == 1 && j < v(i)
            f(i, j) = 0;
        elseif j < v(i)
            f(i, j) = f(i-1, j);
        elseif j == v(i)
            f(i, j) = max(f(i-1, j), w(i) );
            if w(i) > f(i-1, j) && j == capacity
                choice(counter) = i;
                counter = counter + 1;
            end
        else
            f(i, j) = max(f(i-1, j), f(i-1, j-v(i) ) + w(i) );
            if f(i-1, j-v(i) ) + w(i) > f(i-1, j) && j == capacity
                choice(counter) = i;
                counter = counter + 1;
            end
        end
    end
end

max_val = f(num_of_load, capacity);
