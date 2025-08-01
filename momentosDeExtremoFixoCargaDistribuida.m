%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função para o cálculo dos momentos de extremo fixo devido à cargas
% distribuídas ao longo de um único elemento.
%
% Entrada: carga - Vetor contendo as informações da carga aplicada
%              L - comprimento da viga em que a carga está aplicada
%        noIViga - Nó inicial da viga em que a carga está aplicada
%
% Saída: Qf - vetor 6x1 contendo os esforços de extremo fixo do elemento
%
% Autor: Fábio Felipe dos Santos
% Data: 28/04/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Qf = momentosDeExtremoFixoCargaDistribuida(carga, L, noIViga)
% Inicializa o vetor de esforços de extremo fixo
Qf = zeros(6, 1);

% Determinação dos valores "a" e "b" para a carga linearmente distribuída
coorI = carga([4, 5]); coorF = carga([6, 7]);
a = sqrt(sum((coorI - noIViga).^2));
b = sqrt(sum((coorF - noIViga).^2));

% Recolhe os valores da carga inicial e final para a força normal
w0 = carga(8) * cos(deg2rad(carga(10)));
if isnan(carga(9))
    w1 = w0;
else
    w1 = carga(9) * cos(deg2rad(carga(10)));
end

% Colocação das reações horizontais (forças normais)
Qf(1) = ((a - b)*(2*a*w0 - 3*L*w1 - 3*L*w0 + a*w1 + b*w0 + 2*b*w1))/(6*L);
Qf(1) = -Qf(1);
Qf(4) = -((a - b)*(2*a*w0 + a*w1 + b*w0 + 2*b*w1))/(6*L);
Qf(4) = -Qf(4);

% Recolhe os valores da carga inicial e final para a força cortante (força
% vertical)
w0 = carga(8) * sin(deg2rad(carga(10)));
if isnan(carga(9))
    w1 = w0;
else
    w1 = carga(9) * sin(deg2rad(carga(10)));
end

% Colocação das cargas de momento
Qf(3) = -((a - b)*(12*a^3*w0 + 3*a^3*w1 + 3*b^3*w0 + 12*b^3*w1 - 30*L*a^2*w0...
    + 20*L^2*a*w0 - 10*L*a^2*w1 + 10*L^2*a*w1 - 10*L*b^2*w0 + 10*L^2*b*w0...
    - 30*L*b^2*w1 + 20*L^2*b*w1 + 6*a*b^2*w0 + 9*a^2*b*w0 + 9*a*b^2*w1...
    + 6*a^2*b*w1 - 20*L*a*b*w0 - 20*L*a*b*w1))/(60*L^2);
Qf(3) = -Qf(3);
Qf(6) = ((a - b)*(12*a^3*w0 + 3*a^3*w1 + 3*b^3*w0 + 12*b^3*w1 - 15*L*a^2*w0...
    - 5*L*a^2*w1 - 5*L*b^2*w0 - 15*L*b^2*w1 + 6*a*b^2*w0 + 9*a^2*b*w0 +...
    9*a*b^2*w1 + 6*a^2*b*w1 - 10*L*a*b*w0 - 10*L*a*b*w1))/(60*L^2);

% Colocação das cargas Cortante
Qf(2) = -((a - b)*(10*L^3*w0 + 10*L^3*w1 + 8*a^3*w0 + 2*a^3*w1 + 2*b^3*w0...
    + 8*b^3*w1 - 15*L*a^2*w0 - 5*L*a^2*w1 - 5*L*b^2*w0 - 15*L*b^2*w1 +...
    4*a*b^2*w0 + 6*a^2*b*w0 + 6*a*b^2*w1 + 4*a^2*b*w1 - 10*L*a*b*w0 -...
    10*L*a*b*w1))/(20*L^3);
Qf(2) = -Qf(2);
Qf(5) = ((a - b)*(8*a^3*w0 + 2*a^3*w1 + 2*b^3*w0 + 8*b^3*w1 - 15*L*a^2*w0...
    - 5*L*a^2*w1 - 5*L*b^2*w0 - 15*L*b^2*w1 + 4*a*b^2*w0 + 6*a^2*b*w0 +...
    6*a*b^2*w1 + 4*a^2*b*w1 - 10*L*a*b*w0 - 10*L*a*b*w1))/(20*L^3);
Qf(5) = -Qf(5);
end