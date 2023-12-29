function [oppo_line_ee, xi] = find_opposite(anchor_xyz, DataArrays, target_opposite_face, seeking_priority)
% Find the opposite line element id and parametric xi with given anchor point

nElem_seeking = length(seeking_priority);

for ee = 1 : nElem_seeking
    target_elem = seeking_priority(ee);
    vol_ele_id = target_opposite_face(1, target_elem);
    line_id = target_opposite_face(2, target_elem);

    node_ele = zeros(3, DataArrays.nLocBas);
    for aa = 1 : DataArrays.nLocBas
        node_id = DataArrays.IEN_v(aa, vol_ele_id);
        node_ele(:, aa) = DataArrays.Nodes(:, node_id);
    end
    
    % Initial xi = 0.5
    xi = 0.5; start = 0.5;
    [n_out, J, TriQuadPts, node_line] = BuildLine(DataArrays.ele_order, node_ele, line_id, xi);
    
    linequaddata = LineBasis(DataArrays.ele_order, node_line, xi);

    opposite_point = node_line * linequaddata.Basis;

    % Initial residual
    Res0 = norm(anchor_xyz - opposite_point)^2;
    if Res0 < 1e-7
        % Lucky enough
        oppo_line_ee = target_elem;
        return;
    end

    anchor_x = anchor_xyz(1);
    anchor_y = anchor_xyz(2);

    x1 = node_line(1, 1); x2 = node_line(1, 2);
    y1 = node_line(2, 1); y2 = node_line(2, 2);
    N1 = linequaddata.Basis(1); N2 = linequaddata.Basis(2);
    N1_xi = linequaddata.der_xi(1); N2_xi = linequaddata.der_xi(2);

    if DataArrays.ele_order == 1
        % Resx = (N1 * x1 + N2 * x2 - anchor_x)^2
        dResx_dxi = 2 * (x1^2 * N1 * N1_xi + x2^2 * N2 * N2_xi + x1 * x2 * (N1 * N2_xi + N1_xi * N2) ...
                        - anchor_x * (x1 * N1_xi + x2 * N2_xi));

        dResy_dxi = 2 * (y1^2 * N1 * N1_xi + y2^2 * N2 * N2_xi + y1 * y2 * (N1 * N2_xi + N1_xi * N2) ...
                        - anchor_y * (y1 * N1_xi + y2 * N2_xi));

    elseif DataArrays.ele_order == 2
        x3 = node_line(1, 3); y3 = node_line(2, 3);
        N3 = linequaddata.Basis(3); N3_xi = linequaddata.der_xi(3);

        % Resx = (N1 * x1 + N2 * x2 + N3 * x3 - anchor_x)^2
        dResx_dxi = 2 * (x1^2 * N1 * N1_xi + x2^2 * N2 * N2_xi + x3^2 * N3 * N3_xi ...
                         + x1 * x2 * (N1 * N2_xi + N1_xi * N2) + x1 * x3 * (N1 * N3_xi + N1_xi * N3) + x2 * x3 * (N2 * N3_xi + N2_xi * N3) ...
                         - anchor_x * (x1 * N1_xi + x2 * N2_xi + x3 * N3_xi));

        dResy_dxi = 2 * (y1^2 * N1 * N1_xi + y2^2 * N2 * N2_xi + y3^2 * N3 * N3_xi ...
                         + y1 * y2 * (N1 * N2_xi + N1_xi * N2) + y1 * y3 * (N1 * N3_xi + N1_xi * N3) + y2 * y3 * (N2 * N3_xi + N2_xi * N3) ...
                         - anchor_y * (y1 * N1_xi + y2 * N2_xi + y3 * N3_xi));
    end

    dRes_dxi = dResx_dxi + dResy_dxi;

    delta_xi = -Res0 / dRes_dxi;

    s = Res0 / (xi * dRes_dxi);

    xi = xi + s * delta_xi;

    linequaddata = LineBasis(DataArrays.ele_order, node_line, xi);

    opposite_point = node_line * linequaddata.Basis;

    Res1 = norm(anchor_xyz - opposite_point)^2;

    if Res1 < 1e-7
        % Lucky enough
        oppo_line_ee = target_elem;
        return;
    end

    abs_error = abs(Res1 - Res0);
    rel_error =  abs_error / Res0;

    counter_ele = 1;
    
    while (rel_error > 1e-5 && abs_error > 1e-9 && Res1 > 1e-9 && counter_ele < 100)
        N1 = linequaddata.Basis(1); N2 = linequaddata.Basis(2);
        N1_xi = linequaddata.der_xi(1); N2_xi = linequaddata.der_xi(2);

        if DataArrays.ele_order == 1
            % Resx = (N1 * x1 + N2 * x2 - anchor_x)^2
            dResx_dxi = 2 * (x1^2 * N1 * N1_xi + x2^2 * N2 * N2_xi + x1 * x2 * (N1 * N2_xi + N1_xi * N2) ...
                            - anchor_x * (x1 * N1_xi + x2 * N2_xi));
    
            dResy_dxi = 2 * (y1^2 * N1 * N1_xi + y2^2 * N2 * N2_xi + y1 * y2 * (N1 * N2_xi + N1_xi * N2) ...
                            - anchor_y * (y1 * N1_xi + y2 * N2_xi));
    
        elseif DataArrays.ele_order == 2
            N3 = linequaddata.Basis(3); N3_xi = linequaddata.der_xi(3);
    
            % Resx = (N1 * x1 + N2 * x2 + N3 * x3 - anchor_x)^2
            dResx_dxi = 2 * (x1^2 * N1 * N1_xi + x2^2 * N2 * N2_xi + x3^2 * N3 * N3_xi ...
                             + x1 * x2 * (N1 * N2_xi + N1_xi * N2) + x1 * x3 * (N1 * N3_xi + N1_xi * N3) + x2 * x3 * (N2 * N3_xi + N2_xi * N3) ...
                             - anchor_x * (x1 * N1_xi + x2 * N2_xi + x3 * N3_xi));
    
            dResy_dxi = 2 * (y1^2 * N1 * N1_xi + y2^2 * N2 * N2_xi + y3^2 * N3 * N3_xi ...
                             + y1 * y2 * (N1 * N2_xi + N1_xi * N2) + y1 * y3 * (N1 * N3_xi + N1_xi * N3) + y2 * y3 * (N2 * N3_xi + N2_xi * N3) ...
                             - anchor_y * (y1 * N1_xi + y2 * N2_xi + y3 * N3_xi));
        end
    
        dRes_dxi = dResx_dxi + dResy_dxi;
    
        delta_xi = -Res1 / dRes_dxi;
        
        s = Res1 / (xi * dRes_dxi);
        xi = xi + delta_xi;
    
        linequaddata = LineBasis(DataArrays.ele_order, node_line, xi);
    
        opposite_point = node_line * linequaddata.Basis;
    
        Res2 = norm(anchor_xyz - opposite_point)^2;
        
        abs_error = abs(Res2 - Res1);
        rel_error = abs_error / Res0;

        if (Res1 < Res2 && Res1 < 5e-7)
            % Use previous result
            xi = xi - s * delta_xi;
            oppo_line_ee = target_elem;
            return;
        end

        if ((Res1 < Res2 || xi > 1.0 || xi < 0.0) && start == 0.5)
            % Second chance, select another initial value xi = 1.0
            xi = 1.0;
            linequaddata = LineBasis(DataArrays.ele_order, node_line, xi);

            opposite_point = node_line * linequaddata.Basis;

            Res1 = norm(anchor_xyz - opposite_point)^2;
            start = 1.0;
            continue;
        elseif ((Res1 < Res2 || xi > 1.0 || xi < 0.0) && start == 1.0)
            % Third chance, select another initial value xi = 0.0
            xi = 1e-6;
            linequaddata = LineBasis(DataArrays.ele_order, node_line, xi);

            opposite_point = node_line * linequaddata.Basis;

            Res1 = norm(anchor_xyz - opposite_point)^2;
            start = 0.0;
            continue;
        end

        if (xi > 1.0 || xi < 0.0)
            % Fall off from this opposite element
            break;
        end

        Res1 = Res2;
        counter_ele = counter_ele + 1;
    end
    
    % Accept xi in this element
    if xi <= 1.0 && xi >= 0.0
        oppo_line_ee = target_elem;
        return;
    end
    % Else fall off, move on to seek in next opposite element
end

end

