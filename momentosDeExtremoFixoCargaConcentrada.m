%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função para o cálculo dos momentos de extremo fixo devido à cargas
% concentradas ao longo de um único elemento.
%
% Entrada: carga - Vetor contendo as informações da carga aplicada
%              L - comprimento da viga em que a carga está aplicada
%        noIViga - Nó inicial da viga em que a carga está aplicada
%
% Saída: Qf - vetor 6x1 contendo os esforços de extremo fixo do elemento
%
% Autor: Fábio Felipe dos Santos
% Data: 29/04/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Qf = momentosDeExtremoFixoCargaConcentrada(carga, L, noIViga)
% Inicializa o vetor de esforços de extremo fixo
Qf = zeros(6, 1);

% Determinação dos valores "a" e "b" para a carga concentrada
coorI = carga([4, 5]);
a = sqrt(sum((coorI - noIViga).^2));
b = L - a;

if carga(2) == 1    % A carga é concentrada
    % Recolhe os valores da carga para a força normal
    P = carga(8) * cos(deg2rad(carga(10)));

    % Colocação das reações horizontais (forças normais)
    Qf(1) = - (P*b) / L;
    Qf(4) = - (P*a) / L;

    % Recolhe os valores da carga para a força cortante (força vertical)
    P = carga(8) * sin(deg2rad(carga(10)));

    % Colocação das cargas de momento
    Qf(3) = ((-P) * a * b^2) / (L^2);
    Qf(6) = (P * a^2 * b) / (L^2);

    % Colocação das cargas Cortante
    Qf(2) = -(P*(- L^3 + L^2*a + a^2*b - a*b^2)) / L^3;
    Qf(2) = -Qf(2);
    Qf(5) = (P*a*(L^2 - b^2 + a*b)) / L^3;
    Qf(5) = -Qf(5);
else        % A carga é um momento concentrado
    % Recolhe os valores do momento
    M = carga(8);

    % Colocação das cargas de momento
    Qf(3) = (M * b * (2*a - b)) / L^2;
    Qf(6) = (M * a * (2*b - a)) / L^2;

    % Colocação das cargas Cortante
    Qf(2) = (6 * M * a * b) / L^3;
    Qf(5) = -(6 * M * a * b) / L^3;
end
end