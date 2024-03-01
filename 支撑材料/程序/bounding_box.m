function [boxx, boxy] = bounding_box(x, y, spacing);
% 此函数画出点集X（坐标为x和y）的外接矩形，在矩形中均匀生成一系列点

xx = min(x) : spacing : max(x); % 为生成点做准备
yy = min(y) : spacing : max(y); % 为生成点做准备

xlen = length(xx); % box在x方向的长度（元素数量）
ylen = length(yy); % box在y方向的长度（元素数量）

boxx = zeros(1, xlen * ylen); % 预分配内存
boxy = zeros(1, xlen * ylen); % 预分配内存

for i = 1 : ylen % 对每一行
    for j = 1 : xlen % 对每一列
        boxx(xlen * (i-1) + j) = xx(j);
    end
end

for i = 1 : ylen % 对每一行
    for j = 1 : xlen % 对每一列
        boxy(xlen * (i-1) + j) = yy(i);
    end
end

return;

