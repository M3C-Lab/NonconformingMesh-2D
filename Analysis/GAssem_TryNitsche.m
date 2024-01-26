function [K, F] = GAssem_TryNitsche(K, F, DataArrays, BCs, Paras, u_exact)
% Global assembly for intergral on the interface

fprintf("  Global assembling on the interface ...\n");

% Number of interface element
num_inter_elem1 = size(BCs.traction.interface1, 2);
num_inter_elem2 = size(BCs.traction.interface2, 2);

% ratio1 + ratio2 = 1.0
ratio1 = 0.0; ratio2 = 0.0;

% Intergral on interface1
for ee = 1 : num_inter_elem1
    vol_ele_id = BCs.traction.interface1(1, ee);
    line_id = BCs.traction.interface1(2, ee);

    node_ele = zeros(3, DataArrays.nLocBas);
    for aa = 1 : DataArrays.nLocBas
        node_id = DataArrays.IEN_v(aa, vol_ele_id);
        node_ele(:, aa) = DataArrays.Nodes(:, node_id);
    end

    [n_out, J, TriQuadPts, node_line] = BuildLine(DataArrays.ele_order, node_ele, line_id, DataArrays.linequad.qp);

    for qua = 1 : DataArrays.linequad.nqp
        local_basis = TriBasis(DataArrays.ele_order, node_ele, TriQuadPts(:, qua));

        qua_coor = node_ele * local_basis.Basis;

        % Find opposite point on interface2
        [oppo_line_ee, oppo_xi] = find_opposite(qua_coor, DataArrays, BCs.traction.interface2, BCs.seeking_priority1(:, ee));
        oppo_vol_ele_id = BCs.traction.interface2(1, oppo_line_ee);
        oppo_line_id = BCs.traction.interface2(2, oppo_line_ee);

        oppo_node_ele = zeros(3, DataArrays.nLocBas);
        for aa = 1 : DataArrays.nLocBas
            node_id = DataArrays.IEN_v(aa, oppo_vol_ele_id);
            oppo_node_ele(:, aa) = DataArrays.Nodes(:, node_id);
        end

        [oppo_n, oppo_J, oppo_TriQuadPts, oppo_node_line] = BuildLine(DataArrays.ele_order, oppo_node_ele, oppo_line_id, oppo_xi);

        oppo_basis = TriBasis(DataArrays.ele_order, oppo_node_ele, oppo_TriQuadPts);

        [K_loc_loc, K_oppo_oppo, K_loc_oppo, K_oppo_loc] = LocAssem_NitscheInterface(local_basis, n_out(:, qua), oppo_basis, oppo_n, DataArrays, Paras);

        % * gwts
        K_loc_loc = ratio1 * K_loc_loc * DataArrays.linequad.wq(qua) * J(qua);
        K_oppo_oppo = ratio1 * K_oppo_oppo * DataArrays.linequad.wq(qua) * J(qua);
        K_loc_oppo = ratio1 * K_loc_oppo * DataArrays.linequad.wq(qua) * J(qua);
        K_oppo_loc = ratio1 * K_oppo_loc * DataArrays.linequad.wq(qua) * J(qua);

        local_dof = size(DataArrays.LM, 1);

        for aa = 1 : local_dof
            LM_a1 = DataArrays.LM(aa, vol_ele_id);
            if LM_a1 > 0
                for bb = 1 : local_dof
                    LM_b1 = DataArrays.LM(bb, vol_ele_id);
                    if LM_b1 > 0
                        K(LM_a1, LM_b1) = K(LM_a1, LM_b1) + K_loc_loc(aa, bb);
                    else
                        node_xyz = DataArrays.Nodes(:, -LM_b1);
                        g = u_exact(node_xyz(1), node_xyz(2));
                        F(LM_a1) = F(LM_a1) - K_loc_loc(aa, bb) * g;
                    end

                    LM_b2 = DataArrays.LM(bb, oppo_vol_ele_id);
                    if LM_b2 > 0
                        K(LM_a1, LM_b2) = K(LM_a1, LM_b2) + K_loc_oppo(aa, bb);
                    else
                        node_xyz = DataArrays.Nodes(:, -LM_b2);
                        g = u_exact(node_xyz(1), node_xyz(2));
                        F(LM_a1) = F(LM_a1) - K_loc_oppo(aa, bb) * g;
                    end
                end
            end

            LM_a2 = DataArrays.LM(aa, oppo_vol_ele_id);
            if LM_a2 > 0
                for bb = 1 : local_dof
                    LM_b2 = DataArrays.LM(bb, oppo_vol_ele_id);
                    if LM_b2 > 0
                        K(LM_a2, LM_b2) = K(LM_a2, LM_b2) + K_oppo_oppo(aa, bb);
                    else
                        node_xyz = DataArrays.Nodes(:, -LM_b2);
                        g = u_exact(node_xyz(1), node_xyz(2));
                        F(LM_a2) = F(LM_a2) - K_oppo_oppo(aa, bb) * g;
                    end

                    LM_b1 = DataArrays.LM(bb, vol_ele_id);
                    if LM_b1 > 0
                        K(LM_a2, LM_b1) = K(LM_a2, LM_b1) + K_oppo_loc(aa, bb);
                    else
                        node_xyz = DataArrays.Nodes(:, -LM_b1);
                        g = u_exact(node_xyz(1), node_xyz(2));
                        F(LM_a2) = F(LM_a2) - K_oppo_loc(aa, bb) * g;
                    end
                end
            end
        end % aa
    end % qua
end % ee

% Intergral on interface2
for ee = 1 : num_inter_elem2
    vol_ele_id = BCs.traction.interface2(1, ee);
    line_id = BCs.traction.interface2(2, ee);

    node_ele = zeros(3, DataArrays.nLocBas);
    for aa = 1 : DataArrays.nLocBas
        node_id = DataArrays.IEN_v(aa, vol_ele_id);
        node_ele(:, aa) = DataArrays.Nodes(:, node_id);
    end

    [n_out, J, TriQuadPts, node_line] = BuildLine(DataArrays.ele_order, node_ele, line_id, DataArrays.linequad.qp);

    for qua = 1 : DataArrays.linequad.nqp
        local_basis = TriBasis(DataArrays.ele_order, node_ele, TriQuadPts(:, qua));

        qua_coor = node_ele * local_basis.Basis;

        % Find opposite point on interface1
        [oppo_line_ee, oppo_xi] = find_opposite(qua_coor, DataArrays, BCs.traction.interface1, BCs.seeking_priority2(:, ee));
        oppo_vol_ele_id = BCs.traction.interface1(1, oppo_line_ee);
        oppo_line_id = BCs.traction.interface1(2, oppo_line_ee);

        oppo_node_ele = zeros(3, DataArrays.nLocBas);
        for aa = 1 : DataArrays.nLocBas
            node_id = DataArrays.IEN_v(aa, oppo_vol_ele_id);
            oppo_node_ele(:, aa) = DataArrays.Nodes(:, node_id);
        end

        [oppo_n, oppo_J, oppo_TriQuadPts, oppo_node_line] = BuildLine(DataArrays.ele_order, oppo_node_ele, oppo_line_id, oppo_xi);

        oppo_basis = TriBasis(DataArrays.ele_order, oppo_node_ele, oppo_TriQuadPts);

        [K_loc_loc, K_oppo_oppo, K_loc_oppo, K_oppo_loc] = LocAssem_NitscheInterface(local_basis, n_out(:, qua), oppo_basis, oppo_n, DataArrays, Paras);

        % * gwts
        K_loc_loc = ratio2 * K_loc_loc * DataArrays.linequad.wq(qua) * J(qua);
        K_oppo_oppo = ratio2 * K_oppo_oppo * DataArrays.linequad.wq(qua) * J(qua);
        K_loc_oppo = ratio2 * K_loc_oppo * DataArrays.linequad.wq(qua) * J(qua);
        K_oppo_loc = ratio2 * K_oppo_loc * DataArrays.linequad.wq(qua) * J(qua);

        local_dof = size(DataArrays.LM, 1);
        
        for aa = 1 : local_dof
            LM_a1 = DataArrays.LM(aa, vol_ele_id);
            if LM_a1 > 0
                for bb = 1 : local_dof
                    LM_b1 = DataArrays.LM(bb, vol_ele_id);
                    if LM_b1 > 0
                        K(LM_a1, LM_b1) = K(LM_a1, LM_b1) + K_loc_loc(aa, bb);
                    else
                        node_xyz = DataArrays.Nodes(:, -LM_b1);
                        g = u_exact(node_xyz(1), node_xyz(2));
                        F(LM_a1) = F(LM_a1) - K_loc_loc(aa, bb) * g;
                    end

                    LM_b2 = DataArrays.LM(bb, oppo_vol_ele_id);
                    if LM_b2 > 0
                        K(LM_a1, LM_b2) = K(LM_a1, LM_b2) + K_loc_oppo(aa, bb);
                    else
                        node_xyz = DataArrays.Nodes(:, -LM_b2);
                        g = u_exact(node_xyz(1), node_xyz(2));
                        F(LM_a1) = F(LM_a1) - K_loc_oppo(aa, bb) * g;
                    end
                end
            end

            LM_a2 = DataArrays.LM(aa, oppo_vol_ele_id);
            if LM_a2 > 0
                for bb = 1 : local_dof
                    LM_b2 = DataArrays.LM(bb, oppo_vol_ele_id);
                    if LM_b2 > 0
                        K(LM_a2, LM_b2) = K(LM_a2, LM_b2) + K_oppo_oppo(aa, bb);
                    else
                        node_xyz = DataArrays.Nodes(:, -LM_b2);
                        g = u_exact(node_xyz(1), node_xyz(2));
                        F(LM_a2) = F(LM_a2) - K_oppo_oppo(aa, bb) * g;
                    end

                    LM_b1 = DataArrays.LM(bb, vol_ele_id);
                    if LM_b1 > 0
                        K(LM_a2, LM_b1) = K(LM_a2, LM_b1) + K_oppo_loc(aa, bb);
                    else
                        node_xyz = DataArrays.Nodes(:, -LM_b1);
                        g = u_exact(node_xyz(1), node_xyz(2));
                        F(LM_a2) = F(LM_a2) - K_oppo_loc(aa, bb) * g;
                    end
                end
            end
        end % aa
    end % qua
end % ee

end

