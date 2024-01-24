function [K, F] = EssBC(K, F, total_dof, Nodes, BCs, u_exact)
% Essential BC

% Number of Dirichlet nodes
num_DN = length(BCs.Dirichlet.DN);
for nn = 1 : num_DN
    DiriNode = BCs.Dirichlet.DN(nn);
    node_xyz = Nodes(:, DiriNode);

    F(DiriNode) = u_exact(node_xyz(1), node_xyz(2));
    K(:, DiriNode) = zeros(total_dof, 1);
    K(DiriNode, :) = zeros(1, total_dof);
    K(DiriNode, DiriNode) = 1.0;
end

end

