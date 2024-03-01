function h = topology_plot2(center, cluster2, center2, center_2, cluster2_2, center2_2, contact)

Color = [250 192 15; 1 86 153; 243 118 74; 95 198 201; 79 89 100; 1 86 240; ] / 255;
Size = [10, 8, 7, 6, 4, 0.01];
main_line = ["电源1"];
name = ["电源1", "P", "C", "c", "n", "m"];
main_line_2 = ["电源2"];
name_2 = ["电源2", "P.", "C.", "c.", "n.", "m."];
NodeSize = Size(1);
Node_Color = Color(1, :);

for i = 1 : size(center, 1)
    % 一级分叉点
    main_line = cat(2, main_line, name(6)+num2str(i) );
    main_line = cat(2, main_line, name(2)+num2str(i) );
    NodeSize = cat(2, NodeSize, Size(6) );
    NodeSize = cat(2, NodeSize, Size(2) );
    Node_Color = cat(1, Node_Color, Color(6, :) );
    Node_Color = cat(1, Node_Color, Color(2, :) );
end
NodeSize = cat(2, NodeSize, Size(1) );
Node_Color = cat(1, Node_Color, Color(1, :) );
for i = 1 : size(center_2, 1)
    % 一级分叉点
    main_line_2 = cat(2, main_line_2, name_2(6)+num2str(i) );
    main_line_2 = cat(2, main_line_2, name_2(2)+num2str(i) );
    NodeSize = cat(2, NodeSize, Size(6) );
    NodeSize = cat(2, NodeSize, Size(2) );
    Node_Color = cat(1, Node_Color, Color(6, :) );
    Node_Color = cat(1, Node_Color, Color(2, :) );
end

tmp = diag(ones(1, 2*size(cat(1, center, center_2), 1)+1 ), 1);
tmp(1 + 2*size(center, 1), 2 + 2*size(center, 1) ) = 0;
G = graph(tmp, cat(2, main_line, main_line_2), "upper");

for i = 1:size(contact, 1)
    if mod(contact(i, 1), 2) == 0
        left = name(6) + num2str(contact(i, 1) / 2);
    else
        left = name(2) + num2str(floor(contact(i, 1) / 2) );
    end
    if mod(contact(i, 2), 2) == 0
        right = name_2(6) + num2str(contact(i, 2) / 2);
    else
        right = name_2(2) + num2str(floor(contact(i, 2) / 2) );
    end
    G = G.addedge(left, right, 1);
end


counter = 1;
for i = 1 : size(center, 1)
    G = G.addnode(name(3)+num2str(i) );
    G = G.addedge(name(2)+num2str(i), name(3)+num2str(i), 1);
    NodeSize = cat(2, NodeSize, Size(3) );
    Node_Color = cat(1, Node_Color, Color(3, :) );
    for j = 1 : size(center2{i}, 1)
        if size(cell2mat(cluster2{i}(j) ), 1) == 1
            G = G.addnode(name(5)+num2str(counter) );
            G = G.addedge(name(3)+num2str(i), name(5)+num2str(counter), 1);
            counter = counter + 1;
            NodeSize = cat(2, NodeSize, Size(5) );
            Node_Color = cat(1, Node_Color, Color(5, :) );
        else
            G = G.addnode(name(4)+num2str(i)+num2str(j) );
            G = G.addedge(name(3)+num2str(i), name(4)+num2str(i)+num2str(j), 1);
            NodeSize = cat(2, NodeSize, Size(4) );
            Node_Color = cat(1, Node_Color, Color(4, :) );
            for k = 1 : size(cell2mat(cluster2{i}(j) ), 1)
                G = G.addnode(name(5)+num2str(counter) );
                G = G.addedge(name(4)+num2str(i)+num2str(j), name(5)+num2str(counter), 1);
                counter = counter + 1;
                NodeSize = cat(2, NodeSize, Size(5) );
                Node_Color = cat(1, Node_Color, Color(5, :) );
            end
        end
    end
end


counter = 1;
for i = 1 : size(center_2, 1)
    G = G.addnode(name_2(3)+num2str(i) );
    G = G.addedge(name_2(2)+num2str(i), name_2(3)+num2str(i), 1);
    NodeSize = cat(2, NodeSize, Size(3) );
    Node_Color = cat(1, Node_Color, Color(3, :) );
    for j = 1 : size(center2_2{i}, 1)   
        if size(cell2mat(cluster2_2{i}(j) ), 1) == 1
            G = G.addnode(name_2(5)+num2str(counter) );
            G = G.addedge(name_2(3)+num2str(i), name_2(5)+num2str(counter), 1);
            counter = counter + 1;
            NodeSize = cat(2, NodeSize, Size(5) );
            Node_Color = cat(1, Node_Color, Color(5, :) );
        else
            G = G.addnode(name_2(4)+num2str(i)+num2str(j) );
            G = G.addedge(name_2(3)+num2str(i), name_2(4)+num2str(i)+num2str(j), 1);
            NodeSize = cat(2, NodeSize, Size(4) );
            Node_Color = cat(1, Node_Color, Color(4, :) );
            for k = 1 : size(cell2mat(cluster2_2{i}(j) ), 1)
                G = G.addnode(name_2(5)+num2str(counter) );
                G = G.addedge(name_2(4)+num2str(i)+num2str(j), name_2(5)+num2str(counter), 1);
                counter = counter + 1;
                NodeSize = cat(2, NodeSize, Size(5) );
                Node_Color = cat(1, Node_Color, Color(5, :) );
            end
        end
    end
end

hold on;
h = plot(G, "Layout", "force");
% h = plot(G, "Layout", "layered");

h.MarkerSize = NodeSize;
h.NodeColor = Node_Color;
layout(h,'force','UseGravity',true)


tmp = cell(1, size(h.MarkerSize, 2) );
% for i = 1:size(h.MarkerSize, 2)
%     if i <= 2*size(center, 1)+1 && mod(i, 2) == 1
%         tmp{i} = h.NodeLabel{i};
%     elseif i > 2*size(center, 1)+1 && i <= 2*size(cat(1, center, center_2), 1)+2 && mod(i, 2) == 0
%         tmp{i} = h.NodeLabel{i};
%     else
%         tmp{i} = ' ';
%     end
% end
for i = 1:size(h.MarkerSize, 2)
     tmp{i} = ' ';
end
h.NodeLabel = tmp;
h.NodeFontSize = 6;
h.EdgeColor = [0 0 0];
h.LineWidth = 0.5;