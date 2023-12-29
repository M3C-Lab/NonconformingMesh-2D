function ID = modify_ID(input_ID, DiriNode, direction)
% Modify ID array with selected Dirichlet nodes and direction

ID = input_ID;
for nn = 1 : length(DiriNode)
    node = DiriNode(nn);
    ID(direction, node) = -ID(direction, node);
end

end

