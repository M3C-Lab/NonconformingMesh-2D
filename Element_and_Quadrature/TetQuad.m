function [qp, wq, nqp] = TetQuad(degree)
% Quadrature rules of tetrahedral element

if degree == 1
    nqp = 1;
    qp = [0.25; 0.25; 0.25; 0.25];
    wq = [1.0];

elseif degree == 2
    nqp = 4;
    qp = zeros(4, nqp);
    wq = zeros(1, nqp);
    a = 0.5854101966249685;
    b = (1 - a) / 3.0;
    qp(:, 1) = [a; b; b; b];
    qp(:, 2) = [b; a; b; b];
    qp(:, 3) = [b; b; a; b];
    qp(:, 4) = [b; b; b; a];
    for qua = 1 : nqp
        wq(qua) = 0.25;
    end

elseif degree == 3
    nqp = 5;
    qp = zeros(4, nqp);
    wq = zeros(1, nqp);
    qp(:, 1) = [0.25; 0.25; 0.25; 0.25];
    a = 0.5;
    b = (1 - a) / 3.0;
    qp(:, 2) = [a; b; b; b];
    qp(:, 3) = [b; a; b; b];
    qp(:, 4) = [b; b; a; b];
    qp(:, 5) = [b; b; b; a];
    wq(1) = -0.8;
    for qua = 2 : nqp
        wq(qua) = 0.45;
    end

elseif degree == 4
    nqp = 16;
    qp = zeros(4, nqp);
    wq = zeros(1, nqp);
    a = 0.7716429020672371;
    b = (1 - a) / 3.0;
    qp(:, 1) = [a; b; b; b];
    qp(:, 2) = [b; a; b; b];
    qp(:, 3) = [b; b; a; b];
    qp(:, 4) = [b; b; b; a];
    c = 0.1197005277978019;
    d = 0.07183164526766925;
    e = (1 - c - d) / 2.0;
    qp(:, 5) = [c; d; e; e];
    qp(:, 6) = [e; c; d; e];
    qp(:, 7) = [e; e; c; d];
    qp(:, 8) = [d; e; e; c];
    qp(:, 9)  = [c; e; d; e];
    qp(:, 10) = [e; c; e; d];
    qp(:, 11) = [d; e; c; e];
    qp(:, 12) = [e; d; e; c];
    qp(:, 13) = [c; e; e; d];
    qp(:, 14) = [d; c; e; e];
    qp(:, 15) = [e; d; c; e];
    qp(:, 16) = [e; e; d; c];
    for qua = 1 : 4
        wq(qua) = 0.05037379410012282;
    end
    for qua = 5 : 16
        wq(qua) = 0.06654206863329239;
    end

elseif degree == 5
    nqp = 17;
    qp = zeros(4, nqp);
    wq = zeros(1, nqp);
    qp(:, 1) = [0.25; 0.25; 0.25; 0.25];
    a = 0.7316369079576180;
    b = (1 - a) / 3.0;
    qp(:, 2) = [a; b; b; b];
    qp(:, 3) = [b; a; b; b];
    qp(:, 4) = [b; b; a; b];
    qp(:, 5) = [b; b; b; a];
    c = 0.1325810999384657;
    d = 0.0245400379290300;
    e = (1 - c - d) / 2.0;
    qp(:, 6) = [c; d; e; e];
    qp(:, 7) = [e; c; d; e];
    qp(:, 8) = [e; e; c; d];
    qp(:, 9) = [d; e; e; c];
    qp(:, 10) = [c; e; d; e];
    qp(:, 11) = [e; c; e; d];
    qp(:, 12) = [d; e; c; e];
    qp(:, 13) = [e; d; e; c];
    qp(:, 14) = [c; e; e; d];
    qp(:, 15) = [d; c; e; e];
    qp(:, 16) = [e; d; c; e];
    qp(:, 17) = [e; e; d; c];
    wq(1) = 0.1884185567365411;
    for qua = 2 : 5
        wq(qua) = 0.06703858372604275;
    end
    for qua = 6 : 17
        wq(qua) = 0.04528559236327399;
    end

elseif degree == 6
    nqp = 29;
    qp = zeros(4, nqp);
    wq = zeros(1, nqp);
    qp(:, 1) = [0.25; 0.25; 0.25; 0.25];
    a = 0.8277192480479295;
    b = (1 - a) / 3.0;
    qp(:, 2) = [a; b; b; b];
    qp(:, 3) = [b; a; b; b];
    qp(:, 4) = [b; b; a; b];
    qp(:, 5) = [b; b; b; a];
    c = 0.05135188412556341;
    d = 0.4860510285706072;
    e = (1 - c - d) / 2.0;
    qp(:, 6) = [c; d; e; e];
    qp(:, 7) = [e; c; d; e];
    qp(:, 8) = [e; e; c; d];
    qp(:, 9) = [d; e; e; c];
    qp(:, 10) = [c; e; d; e];
    qp(:, 11) = [e; c; e; d];
    qp(:, 12) = [d; e; c; e];
    qp(:, 13) = [e; d; e; c];
    qp(:, 14) = [c; e; e; d];
    qp(:, 15) = [d; c; e; e];
    qp(:, 16) = [e; d; c; e];
    qp(:, 17) = [e; e; d; c];
    f = 0.2967538129690260;
    g = 0.6081079894015281;
    h = (1 - f - g) / 2.0;
    qp(:, 18) = [f, g, h, h];
    qp(:, 19) = [h, f, g, h];
    qp(:, 20) = [h, h, f, g];
    qp(:, 21) = [g, h, h, f];
    qp(:, 22) = [f, h, g, h];
    qp(:, 23) = [h, f, h, g];
    qp(:, 24) = [g, h, f, h];
    qp(:, 25) = [h, g, h, f];
    qp(:, 26) = [f, h, h, g];
    qp(:, 27) = [g, f, h, h];
    qp(:, 28) = [h, g, f, h];
    qp(:, 29) = [h, h, g, f];
    wq(1) = 0.09040129046014750;
    for qua = 2 : 5
        wq(qua) = 0.01911983427899124;
    end
    for qua = 6 : 17
        wq(qua) = 0.04361493840666568;
    end
    for qua = 18 : 29
        wq(qua) = 0.02581167596199161;
    end
    
else
    return;
end

end

