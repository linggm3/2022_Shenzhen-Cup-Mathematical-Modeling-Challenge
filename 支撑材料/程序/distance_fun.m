function [distance, nearest_dot] = distance_fun(p, a, b)
% 计算点p到一条线段的距离，这条线段以a和b为端点

if a == b
    distance = norm(p - a);
    nearest_dot = a;
    return;
end
ap = p - a; % 向量ap
ab = b - a; % 向量ab
projection = ap * ab' / norm(ab)^2 * ab; % ap向ab作投影
factor = projection / ab;
if factor <= 0
    distance = norm(p - a);
    nearest_dot = a;
elseif factor >= 1
    distance = norm(p - b);
    nearest_dot = b;
else
    distance = norm(ap - projection);
    nearest_dot = p - ap + projection;
end

return;