function [money, node, L] = cost_fun(dot, center, center_load)
% 此函数算出总线和一级支线的费用之和

% dot = [dot(1, :); [-25287.5,-13343.7]./1000; dot(2:size(dot,1), :)];
[L, node, ~, sp_num] = minspan(dot); % 用一条最短的折线把dot连接起来
% disp(sp_num)

D = zeros(1, size(center, 1) ); % 一级支线长度
distance = zeros(1, size(node, 1) ); % 某中心到总线每条线段的距离
for i = 1:size(center, 1)
    for j = 1:size(node, 1)
        distance(j) = distance_fun(center(i, :), dot(node(j, 1), :), dot(node(j, 2), :) ); % 第i个中心到所有线段的距离
    end
    D(i) = min(distance); % 第i个中心到总线的最短距离
end

money = 325.7 * L * (sp_num + 1);

for i = 1:size(center, 1)
    if center_load(i) >= 3
        money = money + 239.4 * D(i); 
%         money = money + 600 * D(i); 
    else
        money = money + 188.6 * D(i);
%         money = money + 450 * D(i); 
    end
end

return;
