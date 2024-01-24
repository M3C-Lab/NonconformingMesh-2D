function Line = LineBasis(ele_order, node_ele, qua)
% Basis functions and derivatives with Gmsh node layout
% node_ele = [x1, x2, ... ;
%             y1, y2, ... ;
%             z1, z2, ... ]
% qua = xi in [0, 1] at quadrature point

if ele_order == 1
    Line.Basis = [1 - qua; qua];
    
    % dN/dxi
    Line.der_xi = [-1.0; 1.0];

    % d2N/dxi2
    Line.der_xixi = [0.0; 0.0];

elseif ele_order == 2
    Line.Basis = [2 * qua^2 - 3 * qua + 1;
                  2 * qua^2 - 1 * qua;
                 -4 * qua^2 + 4 * qua];
    
    % dN/dxi
    Line.der_xi = [4 * qua - 3;
                   4 * qua - 1;
                  -8 * qua + 4];

    % d2N/dxi2
    Line.der_xixi = [4; 4; -8];
end

% [dx/dxi; dy/dxi; dz/dxi]
Line.tangent_vector = node_ele * Line.der_xi;
% norm(tangent_vector) is J

% [dxi/dx; dxi/dy; dxi/dz]
inv_Jac = Line.tangent_vector / (norm(Line.tangent_vector)^2);
Line.dN_dx = Line.der_xi * inv_Jac(1);
Line.dN_dy = Line.der_xi * inv_Jac(2);
Line.dN_dz = Line.der_xi * inv_Jac(3);

end

