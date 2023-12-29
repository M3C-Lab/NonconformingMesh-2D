function PlotMesh(DataArrays, meshsize, nb, fork)
% To plot the mesh

max_x = max(DataArrays.Nodes(1, :)); min_x = min(DataArrays.Nodes(1, :));
max_y = max(DataArrays.Nodes(2, :)); min_y = min(DataArrays.Nodes(2, :));

if fork == 0
    IEN_v = DataArrays.IEN_v;
elseif fork == 1
    IEN_v = DataArrays.IEN_v1;
elseif fork == 2
    IEN_v = DataArrays.IEN_v2;
end

nodes = [IEN_v(1, :), IEN_v(2, :), IEN_v(3, :)];
nodes = unique(nodes);
X = zeros(length(nodes), 1); Y = zeros(length(nodes), 1);
temp = 1;
for ii = 1 : DataArrays.nNode
    for jj = 1 : length(nodes)
        if ii == nodes(jj)
            X(temp) = DataArrays.Nodes(1, ii);
            Y(temp) = DataArrays.Nodes(2, ii);
            temp = temp + 1;
        end
    end
end

MESH = alphaShape(X, Y, 0.75 * meshsize, 'HoleThreshold', 1e-6);

figure(nb)
plot(MESH);
title('MESH', 'fontsize', 16);
xlabel('X - axis', 'fontsize', 13); ylabel('Y - axis', 'fontsize', 13);
axis([min_x max_x min_y max_y]);
set(gcf, 'unit', 'centimeters', 'position', [24 20 20 17.5]);

end

% EOF
