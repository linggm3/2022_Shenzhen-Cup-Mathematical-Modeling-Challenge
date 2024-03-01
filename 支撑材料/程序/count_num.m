function [flag] = count_num(index, k)

flag = zeros(k, 1);
for i = 1:length(index)
    for j = 1:k
        if index(i) == j
            flag(j) = flag(j) + 1;
            break;
        end
    end
end

end

