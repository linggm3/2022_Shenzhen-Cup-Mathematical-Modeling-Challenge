function dot = dot_disturb(dot, boxx, boxy, spacing)
% 此函数对点进行扰动

size = length(dot);
disturb_num = ceil(size * rand(1) );
disturb_index = ceil(size * rand(1, disturb_num) ); % 扰动的点的索引
% disturbed = ceil(length(boxx) * rand(1, disturb_num) ); % 扰动后点的box索引 % 一开始用的方法
for i = 1:disturb_num
    % 对每个要扰动的点，在x方向上走[-2, 2]个spacing
    dot(disturb_index(i), 1) = dot(disturb_index(i), 1) + (floor(5 * rand(1) ) - 2) * spacing;
    % 对每个要扰动的点，在y方向上走[-2, 2]个spacing
    dot(disturb_index(i), 2) = dot(disturb_index(i), 2) + (floor(5 * rand(1) ) - 2) * spacing;
end
% dot(disturb_index, 1:2) = [boxx(disturbed)', boxy(disturbed)']; % 一开始用的方法

% matlab函数默认值传递，改变传入形参对main函数的实参没有影响
return;