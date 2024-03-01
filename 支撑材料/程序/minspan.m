function [length, endnode, weight, supply_num] = minspan(dot)
% 函数目的：用一条最短的折线把dot连接起来

len = size(dot, 1);

weight = zeros(1, len^2); % 预分配内存
for i = 1:len
    for j = 1:len
        % 计算权重，值为各点间的欧氏距离
        weight(len*(i-1) + j) = sqrt( (dot(i, 1)-dot(j, 1))^2 + (dot(i, 2)-dot(j, 2))^2 );
    end
end

a = zeros(1, len^2); % 预分配内存
for i = 1:len
    for j = 1:len
        % 为了方便图论计算，a为各边起点
        a(len*(i-1) + j) = i;
    end
end

b = zeros(1, len^2); % 预分配内存
for i = 1:len:len^2 -len+1
    % 为了方便图论计算，b为各边终点
    b(i:i+len-1) = 1:len;
end

G = graph(a, b, weight); % 根据 各边起点a 和 各边起点b 创建图论图
[T, pred] = minspantree(G); % 计算 最小生成树，目的是用一条最短的折线把dot连接起来

tmp1 = 0; tmp2 = 0; e = 0; first = true;
for i = 1:size(dot, 1)
    if pred(i) == 1
        if first
            tmp1 = i;
            first = false;
        else
            tmp2 = i;
        end
    end
end

if tmp2 == 0
    supply_num = 0;
else
    supply_num = 1;
end

while true
    e = 0;
    for i = 1:size(dot, 1) 
        if pred(i) == tmp1
            tmp1 = i;
            e = e + 1;
        end
        if pred(i) == tmp2
            tmp2 = i;
            e = e + 1;
        end
        
    end
%     disp(tmp1);
%     disp(tmp2);
    if e <= 1
        break;
    else
        supply_num = supply_num + 1;
    end
end
% disp(supply_num)
endnode = T.Edges.EndNodes;
weight = T.Edges.Weight;
length = sum(weight);

return;