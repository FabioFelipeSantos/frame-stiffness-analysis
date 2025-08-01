%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função que calcula o esforço cortante nos elementos
%
% Entrada: x - Coordenada x da viga que se deseja o valor esforço normal;
%          cargas - vetor com as cargas aplicadas sobre o elemento;
%          coordI - vetor com as coordenadas x e y do ponto inicial do
%                   elemento;
%          f1y - Valor do esforço cortante local no primeiro nó do elemento;
%
% Saída: V - vetor com os valores dos esforços cortantes.
%
% Autor: Fábio Felipe dos Santos
% Data: 05/05/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = EsforcoCortante(x, cargas, coordI, f1y, L)
% Inicialização do vetor de esforços cortantes
V = f1y * ones(length(x), 1);

% Inicializa os valores de a e b para a carga distribuída
a = []; b = [];

% Valores da carga distribuída
auxCarga = cargas(cargas(:, 2) == 2, :);
if ~isempty(auxCarga)
    % Coordenada inicial da aplicação da carga distribuída
    coordCargaDistI = auxCarga([4, 5]);
    % Se houver carga distribuída, calcula o valor a
    a = sqrt(sum((coordCargaDistI - coordI).^2));
    % Cálculo do ponto final de aplicação da carga distribuída
    coordCargaDistF = auxCarga([6, 7]);
    if ~isempty(coordCargaDistF)
        % Se houver um ponto final para a aplicação da carga distribuída
        % calcula o valor de b
        b = sqrt(sum((coordCargaDistF - coordI).^2));
    else
        % Caso contrário b = L
        b = L;
    end
    % Recebe o valor da carga inicial distribuída
    w1 = auxCarga(8) * sin(deg2rad(auxCarga(10)));
    if isnan(auxCarga(9))
        % Se o valor da carga final for vazio então trata-se de uma carga
        % uniforme
        w2 = w1;
    else
        % Caso contrário, a carga é linear
        w2 = auxCarga(9) * sin(deg2rad(auxCarga(10)));
    end
end

% Verifica se ambos a é vazio, se for, não há aplicação de carga distribuída
if isempty(a)
    % Logo, o valor do esforço cortante será igual ao valor da força
    % vertical local do nó 1 da viga
    return;
else
    % Nesse caso existe uma carga distribuída vertical ao longo do elemento
    auxx = x(sum([x>=a, x<=b], 2) == 2);
    % Valores da carga distribuída
    auxInt = ((a - auxx).*(a*(w1 + w2) - 2*b*w1 + (w1 - w2)*auxx))/(2*(a - b));
    
    % Valor da carga cortante de a até b
    V(sum([x>=a, x<=b], 2) == 2) = V(sum([x>=a, x<=b], 2) == 2) - auxInt;
    
    % Se caso b for menor que L então o restante do esforço cortante será
    % igual ao valor de V em b
    if b < L
        % Armazena o valor de V em b
        auxV = V(x<=b);
        auxV = auxV(end);
        % Valor do esforço cortante
        V(x>b) = auxV;
    end
end
end





















