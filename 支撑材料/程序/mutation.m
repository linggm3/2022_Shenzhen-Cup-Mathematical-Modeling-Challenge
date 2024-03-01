function [pop] = mutation(pop, pm, pace)
% 值传递
[sx, sy, sz] = size(pop);
for i = 1:sz
    for j = 1:sx
        for k = 1:sy
            if rand(1) < pm
                pop(j, k, i) = pop(j, k, i) + 2*pace * rand(1) - pace;
            end
        end
    end
end

return;