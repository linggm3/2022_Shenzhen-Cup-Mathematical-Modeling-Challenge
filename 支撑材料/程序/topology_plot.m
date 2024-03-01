function topology_plot(center, cluster2, center2)

Color = [250 192 15; 1 86 153; 243 118 74; 95 198 201; 79 89 100] / 255;
Size = [10, 8, 7, 6, 4];
main_line = ["supply"];
name = ["supply", "P", "C", "c", "n"];
NodeSize = Size(1);
Node_Color = Color(1, :);

for i = 1 : size(center, 1)
    % 一级分叉点
    main_line = cat(2, main_line, name(2)+num2str(i) );
    NodeSize = cat(2, NodeSize, Size(2) );
    Node_Color = cat(1, Node_Color, Color(2, :) );
end
G = graph(diag(ones(1, size(center, 1) ), 1), main_line, "upper");

counter = 1;
for i = 1 : size(center, 1)
    G = G.addnode(name(3)+num2str(i) );
    G = G.addedge(name(2)+num2str(i), name(3)+num2str(i), 1);
    Node_Color = cat(1, Node_Color, Color(3, :) );
    NodeSize = cat(2, NodeSize, Size(3) );
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
            Node_Color = cat(1, Node_Color, Color(4, :) );
            NodeSize = cat(2, NodeSize, Size(4) );
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

hold on;
h = plot(G, "Layout", "force");
h.MarkerSize = NodeSize;
layout(h,'force','UseGravity',true);
tmp = cell(1, size(h.MarkerSize, 2) );
for i = 1:size(h.MarkerSize, 2)
     tmp{i} = ' ';
end
h.NodeLabel = tmp;
h.NodeColor = Node_Color;
h.EdgeColor = [0 0 0];



