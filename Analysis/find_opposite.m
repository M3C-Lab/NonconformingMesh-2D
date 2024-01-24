function [oppo_line_ee, xi] = find_opposite(anchor_xyz, DataArrays, target_opposite_face, seeking_priority)
% Find the opposite line element id and parametric xi with given anchor point

nElem_seeking = length(seeking_priority);

for ee = 1 : nElem_seeking
    fprintf("ee = %d\n", ee)
    target_elem = seeking_priority(ee);
    vol_ele_id = target_opposite_face(1, target_elem);
    line_id = target_opposite_face(2, target_elem);

    node_ele = zeros(3, DataArrays.nLocBas);
    for aa = 1 : DataArrays.nLocBas
        node_id = DataArrays.IEN_v(aa, vol_ele_id);
        node_ele(:, aa) = DataArrays.Nodes(:, node_id);
    end
    
    % Initial xi = 0.5
    xi = 0.5;
    [n_out, J, TriQuadPts, node_line] = BuildLine(DataArrays.ele_order, node_ele, line_id, xi);
    
    linequaddata = LineBasis(DataArrays.ele_order, node_line, xi);

    opposite_point = node_line * linequaddata.Basis;

    % Initial distance function
    Dist = norm(anchor_xyz - opposite_point);
    % Dist^2 = (x - anchor_x)^2 + (y - anchor_y)^2
    % Res = dDist^2_dxi = 2 * x * dx_dxi - 2 * dx_dxi * anchor_x + 2 * y * dy_dxi - 2 * dy_dxi * anchor_y
    % dRes_dxi = 2 * (x * d2x_dxixi + dx_dxi^2) - 2 * d2x_dxixi * anchor_x + 2 * (y * d2y_dxixi + dy_dxi^2) - 2 * d2x_dxixi * anchor_y

    % But we only reserve the terms of dx_dxi and dy_dxi for tangent.

    anchor_x = anchor_xyz(1);
    anchor_y = anchor_xyz(2);

    x1 = node_line(1, 1); x2 = node_line(1, 2);
    y1 = node_line(2, 1); y2 = node_line(2, 2);
    N1_xi = linequaddata.der_xi(1); N2_xi = linequaddata.der_xi(2);
    % N1_xixi = linequaddata.der_xixi(1); N2_xixi = linequaddata.der_xixi(2);

    x = opposite_point(1);
    y = opposite_point(2);
    
    if DataArrays.ele_order == 1
        dx_dxi = x1 * N1_xi + x2 * N2_xi;
        dy_dxi = y1 * N1_xi + y2 * N2_xi;

        % d2x_dxixi = x1 * N1_xixi + x2 * N2_xixi;
        % d2y_dxixi = y1 * N1_xixi + y2 * N2_xixi;
    elseif DataArrays.ele_order == 2
        x3 = node_line(1, 3); y3 = node_line(2, 3);
        N3_xi = linequaddata.der_xi(3);
        % N3_xixi = linequaddata.der_xixi(3);

        dx_dxi = x1 * N1_xi + x2 * N2_xi + x3 * N3_xi;
        dy_dxi = y1 * N1_xi + y2 * N2_xi + y3 * N3_xi;

        % d2x_dxixi = x1 * N1_xixi + x2 * N2_xixi + x3 * N3_xixi;
        % d2y_dxixi = y1 * N1_xixi + y2 * N2_xixi + y3 * N3_xixi;
    end

    Res0 = 2 * (x * dx_dxi - dx_dxi * anchor_x + y * dy_dxi - dy_dxi * anchor_y);
    fprintf("ee = %d, xi = %e, initial distance = %e, Res0 = %e\n", ee, xi, Dist, Res0);

    if abs(Res0) < 1e-9
        % Lucky enough
        oppo_line_ee = target_elem;
        fprintf("\n");
        return;
    end

    % We only reserve the terms of dx_dxi and dy_dxi for tangent.
    dRes_dxi = 2 * (dx_dxi^2 + dy_dxi^2);

    delta_xi = -Res0 / dRes_dxi;

    xi = xi + delta_xi;

    linequaddata = LineBasis(DataArrays.ele_order, node_line, xi);

    opposite_point = node_line * linequaddata.Basis;

    Dist = norm(anchor_xyz - opposite_point);

    N1_xi = linequaddata.der_xi(1); N2_xi = linequaddata.der_xi(2);
    N1_xixi = linequaddata.der_xixi(1); N2_xixi = linequaddata.der_xixi(2);

    x = opposite_point(1);
    y = opposite_point(2);
    
    if DataArrays.ele_order == 1
        dx_dxi = x1 * N1_xi + x2 * N2_xi;
        dy_dxi = y1 * N1_xi + y2 * N2_xi;

        % d2x_dxixi = x1 * N1_xixi + x2 * N2_xixi;
        % d2y_dxixi = y1 * N1_xixi + y2 * N2_xixi;
    elseif DataArrays.ele_order == 2
        N3_xi = linequaddata.der_xi(3);
        % N3_xixi = linequaddata.der_xixi(3);

        dx_dxi = x1 * N1_xi + x2 * N2_xi + x3 * N3_xi;
        dy_dxi = y1 * N1_xi + y2 * N2_xi + y3 * N3_xi;

        % d2x_dxixi = x1 * N1_xixi + x2 * N2_xixi + x3 * N3_xixi;
        % d2y_dxixi = y1 * N1_xixi + y2 * N2_xixi + y3 * N3_xixi;
    end

    Res1 = 2 * (x * dx_dxi - dx_dxi * anchor_x + y * dy_dxi - dy_dxi * anchor_y);
    fprintf("ee = %d, xi = %e, distance = %e, Res = %e\n", ee, xi, Dist, Res1);

    if abs(Res1) < 1e-9 && xi <= 1.0 && xi >= 0.0 
        % Lucky enough
        oppo_line_ee = target_elem;
        fprintf("\n");
        return;
    end

    abs_error = abs(Res1 - Res0);
    rel_error =  abs_error / abs(Res0);

    counter_ele = 1;
    
    while (rel_error > 1e-5 && abs_error > 1e-9 && abs(Res1) > 1e-9 && counter_ele < 100)
        % We only reserve the terms of dx_dxi and dy_dxi for tangent.
        dRes_dxi = 2 * (dx_dxi^2 + dy_dxi^2);

        delta_xi = -Res1 / dRes_dxi;
    
        xi = xi + delta_xi;
    
        linequaddata = LineBasis(DataArrays.ele_order, node_line, xi);
    
        opposite_point = node_line * linequaddata.Basis;

        Dist = norm(anchor_xyz - opposite_point)^2;

        N1_xi = linequaddata.der_xi(1); N2_xi = linequaddata.der_xi(2);
        % N1_xixi = linequaddata.der_xixi(1); N2_xixi = linequaddata.der_xixi(2);

        x = opposite_point(1);
        y = opposite_point(2);
        
        if DataArrays.ele_order == 1
            dx_dxi = x1 * N1_xi + x2 * N2_xi;
            dy_dxi = y1 * N1_xi + y2 * N2_xi;

            % d2x_dxixi = x1 * N1_xixi + x2 * N2_xixi;
            % d2y_dxixi = y1 * N1_xixi + y2 * N2_xixi;
        elseif DataArrays.ele_order == 2
            N3_xi = linequaddata.der_xi(3);
            % N3_xixi = linequaddata.der_xixi(3);

            dx_dxi = x1 * N1_xi + x2 * N2_xi + x3 * N3_xi;
            dy_dxi = y1 * N1_xi + y2 * N2_xi + y3 * N3_xi;

            % d2x_dxixi = x1 * N1_xixi + x2 * N2_xixi + x3 * N3_xixi;
            % d2y_dxixi = y1 * N1_xixi + y2 * N2_xixi + y3 * N3_xixi;
        end

        Res2 = 2 * (x * dx_dxi - dx_dxi * anchor_x + y * dy_dxi - dy_dxi * anchor_y);
        fprintf("ee = %d, xi = %e, distance = %e, Res = %e", ee, xi, sqrt(Dist), Res2);
        
        abs_error = abs(Res2 - Res1);
        rel_error = abs_error / abs(Res0);

        if (abs(Res1) < abs(Res2) && abs(Res1) < 1e-7)
            % Use previous result
            xi = xi - delta_xi;
            oppo_line_ee = target_elem;
            return;
        end

        Res1 = Res2;
        counter_ele = counter_ele + 1;
        fprintf(", iter = %d\n", counter_ele);
    end
    
    % Accept xi in this element
    if xi <= 1.0 && xi >= 0.0
        oppo_line_ee = target_elem;
        fprintf("\n");
        return;
    else
        % Else fall off, move on to seek in next opposite element
        fprintf("\nFall off \n");
        continue;
    end
end

end

