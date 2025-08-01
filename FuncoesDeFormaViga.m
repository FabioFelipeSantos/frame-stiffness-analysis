function [N, dNdx, B, dBdx] = FuncoesDeFormaViga(x, L)
% Funcões de Forma
N(:, 1) = (1 / L^3) * (2*x.^3 - 3*x.^2*L + L^3);
N(:, 2) = (1 / L^3) * (x.^3*L - 2*x.^2*L^2 + x*L^3);
N(:, 3) = (1 / L^3) * (-2*x.^3 + 3*x.^2*L);
N(:, 4) = (1 / L^3) * (x.^3*L - x.^2*L^2);

% Derivada das funções de Forma
dNdx(:, 1) = (6*x.^2 - 6*L*x) / L^3;
dNdx(:, 2) = (L^2 - 4*L*x + 3*x.^2) / L^2;
dNdx(:, 3) = (-6*x.^2 + 6*L*x) / L^3;
dNdx(:, 4) = (3*x.^2 - 2*L*x)/L^2;

% Derivada Segunda das funções de Forma
B(:, 1) = (12*x - 6*L) / L^3;
B(:, 2) = (6*x - 4*L) / L^2;
B(:, 3) = (6*L - 12*x) / L^3;
B(:, 4) = (6*x - 2*L) / L^2;

% Derivada Terceira das funções de Forma
dBdx(:, 1) = 12 / L^3;
dBdx(:, 2) = 6/L^2;
dBdx(:, 3) = -12 / L^3;
dBdx(:, 4) = 6 / L^2;
end