function [K_ele, F_ele] = LocAssem(ee, DataArrays, f)
% Main local assembly routine for Laplacian equation

% Node coordinate
node_ele = zeros(3, DataArrays.nLocBas);
% [x1, x2, ... ;
%  y1, y2, ... ;
%  z1, z2, ... ]

for aa = 1 : DataArrays.nLocBas
    node_id = DataArrays.IEN_v(aa, ee);
    node_ele(:, aa) = DataArrays.Nodes(:, node_id);
end

% Initilization
local_dof = DataArrays.nLocBas * DataArrays.dof;
K_ele = zeros(local_dof, local_dof);
F_ele = zeros(local_dof, 1);

for qua = 1 : DataArrays.triquad.nqp
    triquaddata = TriBasis(DataArrays.ele_order, node_ele, DataArrays.triquad.qp(:, qua));

    % Quadrature point's xyz-coordinate
    qua_coor = node_ele * triquaddata.Basis;
    
    % weight * J
    gwts = DataArrays.triquad.wq(qua) * norm(triquaddata.n_out);

    for A = 1 : DataArrays.nLocBas
        NA = triquaddata.Basis(A);
        NA_x = triquaddata.dN_dx(A);
        NA_y = triquaddata.dN_dy(A);

        F_ele(A) = F_ele(A) + gwts * NA * f(qua_coor(1), qua_coor(2));

        for B = 1 : DataArrays.nLocBas
            NB_x = triquaddata.dN_dx(B);
            NB_y = triquaddata.dN_dy(B);
            K_ele(A, B) = K_ele(A, B) + gwts * (NA_x * NB_x + NA_y * NB_y);
        end
    end
end

end

