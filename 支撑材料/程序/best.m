function [best_individual, best_fitval] = best(pop, fitvalue)

[~, ~, sz] = size(pop);
best_individual = pop(:, :, 1);
best_fitval = fitvalue(1);
for i = 2:sz
    if fitvalue(i) < best_fitval
        best_individual = pop(:, :, i);
        best_fitval = fitvalue(i);
    end
end

return;