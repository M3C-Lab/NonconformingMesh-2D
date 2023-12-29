function [K, F] = EssBC(K, F, total_dof, BCs)
% Essential BC

% Number of Dirichlet nodes
num_DN = length(BCs.Dirichlet.DN);
for nn = 1 : num_DN
    DiriNode = BCs.Dirichlet.DN(nn);

    F(DiriNode) = 0.0;
    K(:, DiriNode) = zeros(total_dof, 1);
    K(DiriNode, :) = zeros(1, total_dof);
    K(DiriNode, DiriNode) = 1.0;
end

end

