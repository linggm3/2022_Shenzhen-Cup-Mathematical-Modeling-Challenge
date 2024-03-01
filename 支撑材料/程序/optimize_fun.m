function cost = optimize_fun(x)

global center2; global center2_load; global edge;
global i;  global dott; global center_load;

cost = 0;
d = zeros(1, size(edge, 1) );
for k = 1:size(edge, 1)
    % 点x到总线各段的距离
    d(i) = distance_fun(x, dott(edge(k, 1), :), dott(edge(k, 2), :) );
end
distance = min(d); % 点x到总线的最短距离

if center_load(i) >= 3
    weight = 239.4; % 支线B
else
    weight = 188.6; % 支线A
end

cost = weight * distance;
for j = 1:size(center2{i}, 1)
    if center2_load{i}(j) >= 3
        weight = 239.4; % 支线B
    else
        weight = 188.6; % 支线A
    end
    cost = cost + weight * sqrt( ( center2{i}(j, 1) - x(1) )^2 + ( center2{i}(j, 2) - x(2) )^2 );
end


return;