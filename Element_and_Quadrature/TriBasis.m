function Tri = TriBasis(ele_order, node_ele, qua)
% Basis functions and derivatives with Gmsh node layout
% node_ele = [x1, x2, x3, ... ;
%             y1, y2, y3, ... ;
%             z1, z2, z3, ... ]
% qua = [r; s; t] at quadrature point

if ele_order == 1
    Tri.Basis = [qua(3); qua(1); qua(2)];
    
    Tri.der_r = [-1; 1; 0];
    Tri.der_s = [-1; 0; 1];

elseif ele_order == 2
    r = qua(1); s = qua(2); t = qua(3);

    Tri.Basis = [t * (2 * t - 1); % 2(r^2 + s^2) + 4rs - 3(r + s) + 1
                 r * (2 * r - 1);
                 s * (2 * s - 1);
                 4 * r * t; % 4r(1 - r - s)
                 4 * r * s;
                 4 * s * t]; % 4s(1 - r - s)
    
    % dN/dr
    Tri.der_r = [4 * r + 4 * s - 3;
             4 * r - 1;
             0;
             -8 * r - 4 * s + 4;
             4 * s;
             -4 * s];

    % dN/ds
    Tri.der_s = [4 * r + 4 * s - 3;
             0;
             4 * s - 1;
             -4 * r;
             4 * r;
             -8 * s - 4 * r + 4];
    
end

Der_1 = [Tri.der_r, Tri.der_s];
% [dN1/dr, dN1/ds;
%  dN2/dr, dN2/ds;
%  ...

nonsymm_Jac = node_ele * Der_1;
% [dx/dr, dx/ds;
%  dy/dr, dy/ds;
%  dz/dr, dz/ds]

% un-normalized outward normal vector
Tri.n_out = cross(nonsymm_Jac(:, 1), nonsymm_Jac(:, 2));

% First fundamental form
FFF = nonsymm_Jac' * nonsymm_Jac;

Grad = zeros(size(node_ele, 2), 3);
% [dN1/dx, dN1/dy, dN1/dz;
%  dN2/dx, dN2/dy, dN2/dz;
%  ...
for aa = 1 : size(node_ele, 2)
    grad_xi = Der_1(aa, :)';
    % [dN/dr; dN/ds]

    grad_x = (nonsymm_Jac * inv(FFF)) * grad_xi;
    % [dN/dx; dN/dy; dN/dz]

    Grad(aa, :) = grad_x';
end
Tri.dN_dx = Grad(:, 1);
Tri.dN_dy = Grad(:, 2);
Tri.dN_dz = Grad(:, 3);

end

