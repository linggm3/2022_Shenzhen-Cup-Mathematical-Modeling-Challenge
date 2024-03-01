function [newpop, fitvalue]=selection(pop, fitvalue)
% 自然选择优胜劣汰，选择最好的一半个体组成新群体

[~, ~, sz] = size(pop);
[~, idx] = sort(fitvalue, 'ascend'); % 从小到大排列

newpop = pop(:, :, idx(1:sz/2) );
fitvalue = fitvalue(idx(1:sz/2) );
    
return;