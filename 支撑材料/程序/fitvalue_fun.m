function fitvalue = fitvalue_fun(pop, center, center_load)
% 计算适应值,此适应值对应线路的花费，应越小越好

fitvalue = zeros(1, size(pop, 3) );
for i = 1:size(pop, 3)
    fitvalue(i) = cost_fun(pop(:, :, i), center, center_load);
end

return;