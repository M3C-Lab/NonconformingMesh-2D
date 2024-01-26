function [K, F] = WeakBC(K, F, DataArrays, BCs, Paras, u_exact)
% Weakly enforce Dirichlet BC by Nitsche's method

fprintf("  Global assembling for weaky enforced Dirichlet BC ...\n");

% Number of boundary element
num_Diri_elem = size(BCs.Dirichlet.DiriFace, 2);

for ee = 1 : num_Diri_elem
    vol_ele_id = BCs.Dirichlet.DiriFace(1, ee);
    line_id = BCs.Dirichlet.DiriFace(2, ee);

    node_ele = zeros(3, DataArrays.nLocBas);
    for aa = 1 : DataArrays.nLocBas
        node_id = DataArrays.IEN_v(aa, vol_ele_id);
        node_ele(:, aa) = DataArrays.Nodes(:, node_id);
    end

    [n_out, J, TriQuadPts, node_line] = BuildLine(DataArrays.ele_order, node_ele, line_id, DataArrays.linequad.qp);
    
    local_dof = DataArrays.nLocBas * DataArrays.dof;
    K_NitscheBC = zeros(local_dof, local_dof);
    F_NitscheBC = zeros(local_dof, 1);

    for qua = 1 : DataArrays.linequad.nqp
        local_basis = TriBasis(DataArrays.ele_order, node_ele, TriQuadPts(:, qua));

        qua_coor = node_ele * local_basis.Basis;
        g = u_exact(qua_coor(1), qua_coor(2));
        
        n_2d = [n_out(1, qua); n_out(2, qua)];

        hh = 1 / sqrt(dot(n_2d, inv(local_basis.FFF) * n_2d));

        penalty = Paras.penalty_coef / hh;

        gwts = J(qua) * DataArrays.linequad.wq(qua);

        for A = 1 : DataArrays.nLocBas
            NA   = local_basis.Basis(A);
            NA_x = local_basis.dN_dx(A);
            NA_y = local_basis.dN_dy(A);

            F_NitscheBC(A) = F_NitscheBC(A) + gwts * ...
                (- Paras.adjoint_coef * (NA_x * n_2d(1) + NA_y * n_2d(2)) * g ...
                 + penalty * NA * g);

            for B = 1 : DataArrays.nLocBas
                NB   = local_basis.Basis(B);
                NB_x = local_basis.dN_dx(B);
                NB_y = local_basis.dN_dy(B);

                K_NitscheBC(A, B) = K_NitscheBC(A, B) + gwts * ...
                    (- NA * (NB_x * n_2d(1) + NB_y * n_2d(2)) ...
                     - Paras.adjoint_coef * (NA_x * n_2d(1) + NA_y * n_2d(2)) * NB ...
                     + penalty * NA * NB);
            end % B
        end % A

        stop = 1;
    end % qua
    
    % Assembly
    for aa = 1 : local_dof
        LM_a = abs(DataArrays.LM(aa, vol_ele_id));
        F(LM_a) = F(LM_a) + F_NitscheBC(aa);

        for bb = 1 : local_dof
            LM_b = abs(DataArrays.LM(bb, vol_ele_id));
            K(LM_a, LM_b) = K(LM_a, LM_b) + K_NitscheBC(aa, bb);
        end
    end
end % ee

end

