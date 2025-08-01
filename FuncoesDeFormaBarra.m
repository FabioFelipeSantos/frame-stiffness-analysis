function [N, dNdx] = FuncoesDeFormaBarra(x, L)
% Funcões de Forma
N(:, 1) = (1 / L) * (L - x);
N(:, 2) = (1 / L) * x;

% Derivada das funções de Forma
dNdx(:, 1) = -1/L * ones(length(x), 1);
dNdx(:, 2) = 1/L * ones(length(x), 1);
end