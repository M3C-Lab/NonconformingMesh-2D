function [qp, wq, nqp] = TriQuad(degree)
% To generate Hammer's quadrature info on Page 173 for the default 
% parent triangle on Page 121 with parent coordinates (xi, eta).
% Input:
%   degree: The degree of precision.
% Output:
%   qp: The area coodinates of quadrature points.
%       qp(*, e): The 'e'th quadrature points.
%       qp(1, *): The 'r' coordinate, N1.
%       qp(2, *): The 's' coordinate, N2.
%       qp(3, *): The 't' coordinate, N3.
%   wq: The quadrature weight.
%       wq(e): The weight of the 'e'th quadrature point.
%   nqp: The number of quadrature points.

% !!!!!! Note: We will make the number of quadrature points (nqp) as little 
%   as possible. The available choices are:
%   (nqp, degree) = 
%   (3, 2)  (4, 3)  (6, 4)  (7, 5)  (12, 6) (13, 7)

% !!!!!! Note: The quadrature rule with the triangular coordinates (r, s, t)
%   should be use in the (r, s) space of Figure 3.I.4 on Page 166.
%   In order words, the mapping rule from (r, s) space to the physical space
%   is necessary.

if degree == 2
    nqp = 3;
    qp = zeros(3, nqp); % qp with triangular coordinates
    wq = zeros(1, nqp);
    for ii = 1 : 3
        for jj = 1 : 3
            if ii == jj
                qp(ii, jj) = 0.666666666666667;
            else
                qp(ii, jj) = 0.166666666666667;
            end
            wq(1, jj) = 0.333333333333333;
        end
    end
    
elseif degree == 3
    nqp = 4;
    qp = zeros(3, nqp);
    wq = zeros(1, nqp);
    for ii = 1 : 3
        qp(ii, 1) = 0.333333333333333;
        for jj = 2 : 4
            if ii == jj - 1
                qp(ii, jj) = 0.6;
            else
                qp(ii, jj) = 0.2;
            end
            wq(1, jj) = 0.520833333333333;
        end
    end
    wq(1, 1) = -0.56250;

elseif degree == 4
    nqp = 6;
    qp = zeros(3, nqp);
    wq = zeros(1, nqp);
    for ii = 1 : 3
        for jj = 1 : 3
            wq(1, jj) = 0.109951743655322;
            if ii == jj
                qp(ii, jj) = 0.816847572980459;
            else
                qp(ii, jj) = 0.091576213509771;
            end
        end
        for jj = 4 : 6
            wq(1, jj) = 0.223381589678011;
            if ii == jj - 3
                qp(ii, jj) = 0.108103018168070;
            else
                qp(ii, jj) = 0.445948490915965;
            end
        end
    end
    
elseif degree == 5
    nqp = 7;
    qp = zeros(3, nqp);
    wq = zeros(1, nqp);
    wq(1, 1) = 0.22500;
    for ii = 1 : 3
        qp(ii, 1) = 0.333333333333333;
        for jj = 2 : 4
            if ii == jj - 1
                qp(ii, jj) = 0.797426985353087;
            else
                qp(ii, jj) = 0.101286507323456;
            end
            wq(1, jj) = 0.125939180544827;
        end
        for jj = 5 : 7
            if ii == jj - 4
                qp(ii, jj) = 0.059715871789770;
            else
                qp(ii, jj) = 0.470142064105115;
            end
            wq(1, jj) = 0.132394152788506;
        end
    end
    
elseif degree == 6
    nqp = 12;
    qp = zeros(3, nqp);
    wq = zeros(1, nqp);
    for ii = 1 : 3
        for jj = 1 : 3
            wq(1, jj) = 0.050844906370207;
            if ii == jj
                qp(ii, jj) = 0.873821971016996;
            else
                qp(ii, jj) = 0.063089014491502;
            end
        end
        for jj = 4 : 6
            wq(1, jj) = 0.116786275726379;
            if ii == jj - 3
                qp(ii, jj) = 0.501426509658179;
            else
                qp(ii, jj) = 0.249286745170910;
            end
        end
        for jj = 7 : 12
            wq(1, jj) = 0.082851075618374;
        end
        r = 0.636502499121399;
        s = 0.310352451033785;
        t = 0.053145049844816;
        qp(:, 7) = [r, s, t]';
        qp(:, 8) = [s, t, r]';
        qp(:, 9) = [t, r, s]';
        qp(:, 10) = [r, t, s]';
        qp(:, 11) = [t, s, r]';
        qp(:, 12) = [s, r, t]';
    end
    
elseif degree == 7
    nqp = 13;
    qp = zeros(3, nqp);
    wq = zeros(1, nqp);
    wq(1, 1) = -0.149570044467670;
    for ii = 1 : 3
        qp(ii, 1) = 0.333333333333333;
        for jj = 2 : 4
            wq(1, jj) = 0.175615257433204;
            if ii == jj - 1
                qp(ii, jj) = 0.479308067841923;
            else
                qp(ii, jj) = 0.260345966079038;
            end
        end
        for jj = 5 : 7
            wq(1, jj) = 0.053347235608839;
            if ii == jj - 4
                qp(ii, jj) = 0.869739794195568;
            else
                qp(ii, jj) = 0.065130102902216;
            end
        end
        for jj = 8 : 13
            wq(1, jj) = 0.077113760890257;
        end
        r = 0.638444188569809;
        s = 0.312865496004875;
        t = 0.086903154253160;
        qp(:, 7) = [r, s, t]';
        qp(:, 8) = [s, t, r]';
        qp(:, 9) = [t, r, s]';
        qp(:, 10) = [r, t, s]';
        qp(:, 11) = [t, s, r]';
        qp(:, 12) = [s, r, t]';
    end
    
else
    disp('TriangularQuad: Please check the data input.');
    disp("The degree input should be a interger from 2 to 7.");
    return;
end

% For triangles in a plane.
wq = 0.5 * wq;
end

% EOF

