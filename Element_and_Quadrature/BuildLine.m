function [n_out, J, TriQuadPts, node_line] = BuildLine(ele_order, node_ele, line_id, LineQuadPts)
% Get normals and triangular quadrature rules with boundaryline info

% n_out: Normalized outward normal vector
% J: Line Jacobian
% TriQuadPts: Quadrature point on line with area coordinate in triangle

nqp = length(LineQuadPts);
TriQuadPts = zeros(3, nqp);

if line_id == 1
    if ele_order == 1
        node_line = [node_ele(:, 2), node_ele(:, 3)];
    elseif ele_order == 2
        node_line = [node_ele(:, 2), node_ele(:, 3), node_ele(:, 5)];
    end

    for qua = 1 : nqp
        TriQuadPts(1, qua) = 1.0 - LineQuadPts(qua);
        TriQuadPts(2, qua) = LineQuadPts(qua);
    end

elseif line_id == 2
    if ele_order == 1
        node_line = [node_ele(:, 3), node_ele(:, 1)];
    elseif ele_order == 2
        node_line = [node_ele(:, 3), node_ele(:, 1), node_ele(:, 6)];
    end

    for qua = 1 : nqp
        TriQuadPts(2, qua) = 1.0 - LineQuadPts(qua);
        TriQuadPts(3, qua) = LineQuadPts(qua);
    end

elseif line_id == 3
    if ele_order == 1
        node_line = [node_ele(:, 1), node_ele(:, 2)];
    elseif ele_order == 2
        node_line = [node_ele(:, 1), node_ele(:, 2), node_ele(:, 4)];
    end

    for qua = 1 : nqp
        TriQuadPts(3, qua) = 1.0 - LineQuadPts(qua);
        TriQuadPts(1, qua) = LineQuadPts(qua);
    end

end

n_out = zeros(3, nqp);
J = zeros(1, nqp);
for qua = 1 : nqp
    line = LineBasis(ele_order, node_line, LineQuadPts(qua));
    tan = line.tangent_vector;
    J(qua) = norm(tan);

    triangle = TriBasis(ele_order, node_ele, TriQuadPts(:, qua));
    face_n_out = triangle.n_out;

    line_n_out = cross(tan, face_n_out);
    line_n_out = line_n_out / norm(line_n_out);

    qua_coor = node_line * line.Basis;
    opposite = node_ele(:, line_id);

    inward_vector = opposite - qua_coor;
    if dot(inward_vector, line_n_out) > 0
        line_n_out = -1.0 * line_n_out;
    end
    n_out(:, qua) = line_n_out;
end

end

