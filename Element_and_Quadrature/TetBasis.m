function Tet = TetBasis(ele_order, node_ele, qua)
% Basis functions and derivatives with Gmsh node layout
% node_ele = [x1, x2, x3, ... ;
%             y1, y2, y3, ... ;
%             z1, z2, z3, ... ]
% qua = [r; s; t; u] at quadrature point

if ele_order == 1
    Tet.Basis = [qua(4); qua(1); qua(2); qua(3)];

    Tet.der_r = [-1; 1; 0; 0];
    Tet.der_s = [-1; 0; 1; 0];
    Tet.der_t = [-1; 0; 0; 1];

elseif ele_order == 2
    r = qua(1); s = qua(2); t = qua(3); u = qua(4);
    Tet.Basis = [u * (2 * u - 1); % 2(r^2 + s^2 + t^2) + 4(rs + rt + st) - 3(r + s + t) + 1
                 r * (2 * r - 1);
                 s * (2 * s - 1);
                 t * (2 * t - 1);
                 4 * u * r; % 4r(1 - r - s - t)
                 4 * r * s;
                 4 * u * s; % 4s(1 - r - s - t)
                 4 * u * t; % 4t(1 - r - s - t)
                 4 * s * t;
                 4 * r * t];
    
    % dN/dr
    Tet.der_r = [4 * r + 4 * (s + t) - 3;
             4 * r - 1;
             0;
             0;
             -8 * r - 4 * (s + t) + 4;
             4 * s;
             -4 * s;
             -4 * t;
             0;
             4 * t];
    
    % dN/ds
    Tet.der_s = [4 * s + 4 * (r + t) - 3;
             0;
             4 * s - 1;
             0;
             -4 * r;
             4 * r;
             -8 * s - 4 * (r + t) + 4;
             -4 * t;
             4 * t;
             0];
    
    % dN/dt
    Tet.der_t = [4 * t + 4 * (r + s) - 3;
             0;
             0;
             4 * t - 1;
             -4 * r;
             0;
             -4 * s;
             -8 * t - 4 * (r + s) + 4;
             4 * s;
             4 * r];

    % d2N/drr
    der_rr = [4; 4; 0; 0; -8; 0; 0; 0; 0; 0];

    % d2N/dss
    der_ss = [4; 0; 4; 0; 0; 0; -8; 0; 0; 0];

    % d2N/dtt
    der_tt = [4; 0; 0; 4; 0; 0; 0; -8; 0; 0];

    % d2N/drs
    der_rs = [4; 0; 0; 0; -4; 4; -4; 0; 0; 0];

    % d2N/drt
    der_rt = [4; 0; 0; 0; -4; 0; 0; -4; 0; 4];

    % d2N/dst
    der_st = [4; 0; 0; 0; 0; 0; -4; -4; 4; 0];

end

Der_1 = [Tet.der_r, Tet.der_s, Tet.der_t];
% [dN1/dr, dN1/ds, dN1/dt;
%  dN2/dr, dN2/ds, dN2/dt;
%  ...

Tet.Jacobian = node_ele * Der_1;
% [dx/dr, dx/ds, dx/dt;
%  dy/dr, dy/ds, dy/dt;
%  dz/dr, dz/ds, dz/dt]
% (matrix)

Grad = Der_1 / (Tet.Jacobian);
% [dN1/dx, dN1/dy, dN1/dz;
%  dN2/dx, dN2/dy, dN2/dz;
%  ...

Tet.dN_dx = Grad(:, 1);
Tet.dN_dy = Grad(:, 2);
Tet.dN_dz = Grad(:, 3);

if ele_order == 1
    Tet.d2N_dxx = zeros(4, 1);
    Tet.d2N_dyy = zeros(4, 1);
    Tet.d2N_dzz = zeros(4, 1);
    Tet.d2N_dxy = zeros(4, 1);
    Tet.d2N_dyz = zeros(4, 1);
    Tet.d2N_dxz = zeros(4, 1);
elseif ele_order == 2
    Tet.d2N_dxx = zeros(10, 1);
    Tet.d2N_dyy = zeros(10, 1);
    Tet.d2N_dzz = zeros(10, 1);
    Tet.d2N_dxy = zeros(10, 1);
    Tet.d2N_dxz = zeros(10, 1);
    Tet.d2N_dyz = zeros(10, 1);

    x_r = Tet.Jacobian(1, 1); x_s = Tet.Jacobian(1, 2); x_t = Tet.Jacobian(1, 3);
    y_r = Tet.Jacobian(2, 1); y_s = Tet.Jacobian(2, 2); y_t = Tet.Jacobian(2, 3);
    z_r = Tet.Jacobian(3, 1); z_s = Tet.Jacobian(3, 2); z_t = Tet.Jacobian(3, 3);

    x_rr = dot(node_ele(1, :), der_rr);
    x_ss = dot(node_ele(1, :), der_ss);
    x_tt = dot(node_ele(1, :), der_tt);
    x_rs = dot(node_ele(1, :), der_rs);
    x_rt = dot(node_ele(1, :), der_rt);
    x_st = dot(node_ele(1, :), der_st);

    y_rr = dot(node_ele(2, :), der_rr);
    y_ss = dot(node_ele(2, :), der_ss);
    y_tt = dot(node_ele(2, :), der_tt);
    y_rs = dot(node_ele(2, :), der_rs);
    y_rt = dot(node_ele(2, :), der_rt);
    y_st = dot(node_ele(2, :), der_st);

    z_rr = dot(node_ele(3, :), der_rr);
    z_ss = dot(node_ele(3, :), der_ss);
    z_tt = dot(node_ele(3, :), der_tt);
    z_rs = dot(node_ele(3, :), der_rs);
    z_rt = dot(node_ele(3, :), der_rt);
    z_st = dot(node_ele(3, :), der_st);

    LHS = [x_r^2, y_r^2, z_r^2, 2 * x_r * y_r, 2 * x_r * z_r, 2 * y_r * z_r;
           x_s^2, y_s^2, z_s^2, 2 * x_s * y_s, 2 * x_s * z_s, 2 * y_s * z_s;
           x_t^2, y_t^2, z_t^2, 2 * x_t * y_t, 2 * x_t * z_t, 2 * y_t * z_t;
           x_r * x_s, y_r * y_s, z_r * z_s, (x_r * y_s + x_s * y_r), (x_r * z_s + x_s * z_r), (y_r * z_s + z_r * y_s);
           x_r * x_t, y_r * y_t, z_r * z_t, (x_r * y_t + x_t * y_r), (x_r * z_t + x_t * z_r), (y_r * z_t + z_r * y_t);
           x_t * x_s, y_t * y_s, z_t * z_s, (x_t * y_s + x_s * y_t), (x_t * z_s + x_s * z_t), (y_t * z_s + z_t * y_s)];

    for aa = 1 : 10
        RHS = [der_rr(aa) - Tet.dN_dx(aa) * x_rr - Tet.dN_dy(aa) * y_rr - Tet.dN_dz(aa) * z_rr;
               der_ss(aa) - Tet.dN_dx(aa) * x_ss - Tet.dN_dy(aa) * y_ss - Tet.dN_dz(aa) * z_ss;
               der_tt(aa) - Tet.dN_dx(aa) * x_tt - Tet.dN_dy(aa) * y_tt - Tet.dN_dz(aa) * z_tt;
               der_rs(aa) - Tet.dN_dx(aa) * x_rs - Tet.dN_dy(aa) * y_rs - Tet.dN_dz(aa) * z_rs;
               der_rt(aa) - Tet.dN_dx(aa) * x_rt - Tet.dN_dy(aa) * y_rt - Tet.dN_dz(aa) * z_rt;
               der_st(aa) - Tet.dN_dx(aa) * x_st - Tet.dN_dy(aa) * y_st - Tet.dN_dz(aa) * z_st];

        GradGrad = LHS \ RHS;
        % [d2N/dxx; d2N/dyy; d2N/dzz; d2N/dxy; d2N/dxz; d2N/dyz]
        Tet.d2N_dxx(aa) = GradGrad(1);
        Tet.d2N_dyy(aa) = GradGrad(2);
        Tet.d2N_dzz(aa) = GradGrad(3);
        Tet.d2N_dxy(aa) = GradGrad(4);
        Tet.d2N_dxz(aa) = GradGrad(5);
        Tet.d2N_dyz(aa) = GradGrad(6);
    end
end

