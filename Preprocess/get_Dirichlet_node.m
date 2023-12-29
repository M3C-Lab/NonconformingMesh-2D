function DiriNode = get_Dirichlet_node(msh, surface_name, ele_order)
% Get Dirichlet nodes' id

% Number of physical groups
num_PG = size(msh.PhyGrp, 1);

% Find surface id
for ii = 1 : num_PG
    if strcmp(surface_name, msh.PhyGrp{ii, 3})
        PG_id = msh.PhyGrp{ii, 2};
        break;
    end
end

if ele_order == 1
    num_total_sur_ele = msh.nbLines;
    sur_nLocBas = 2;
    sur_ele = msh.LINES;
elseif ele_order == 2
    num_total_sur_ele = msh.nbLines3;
    sur_nLocBas = 3;
    sur_ele = msh.LINES3;
end

% Find surface IEN
IEN_s = [];
for ee = 1 : num_total_sur_ele
    if sur_ele(ee, sur_nLocBas + 1) == PG_id    % Catch surface element
        localien = sur_ele(ee, 1 : sur_nLocBas);
        IEN_s = [IEN_s, localien];  % Insert local IEN
    end
end

DiriNode = unique(IEN_s);

end

