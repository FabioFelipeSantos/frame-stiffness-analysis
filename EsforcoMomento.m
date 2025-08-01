%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função que calcula o esforço momento fletor nos elementos
%
% Entrada: x - Coordenada x da viga que se deseja o valor esforço normal;
%          cargas - vetor com as cargas aplicadas sobre o elemento;
%          coordI - vetor com as coordenadas x e y do ponto inicial do
%                   elemento;
%          m1 - Valor do esforço normal local no primeiro nó do elemento;
%          f1y - Valor do esforço cortante no local do primeiro nó do
%                elemento;
%          L - Comprimento do Elemento.
%
% Saída: M - vetor com os valores dos esforços momentos fletores.
%
% Autor: Fábio Felipe dos Santos
% Data: 05/05/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = EsforcoMomento(x, cargas, coordI, m1, f1y, L)
% Inicialização do vetor de esforços momentos
M = m1 * ones(length(x), 1);

% Inicializa os valores de a e b para a carga distribuída e c para a carga
% concentrada
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
    if a == L
        a = 0;
        if b == 0
            b = L;
        end
    elseif b < a
        temp = b;
        b = a;
        a = temp;
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
    % Logo, o valor do esforço momento fletor será igual ao valor do
    % momento do nó 1 do elemento mais o momento causado pela carga
    % vertical do nó 1.
    M = M + f1y * x;
else
    M(x<a) = M(x<a) + f1y * x(x<a);
    auxM = M(x<=a);
    auxM = auxM(end);

    % Nesse caso existe uma carga distribuída vertical ao longo do elemento
    auxx = x(sum([x>=a, x<=b], 2) == 2);
    % Valores da carga distribuída
    auxInt = -(((a - auxx).^2).*(2*a*w1 + a*w2 - 3*b*w1 + (w1 - w2)*auxx))/(6*(a - b));
    
    % Valor da carga normal de a até b
    M(sum([x>=a, x<=b], 2)==2) = auxM + f1y * x - auxInt;
    
    % Se caso b for menor que L então o restante do esforço normal será
    % igual ao valor de N em b
    if b < L
        % Armazena o valor de M em b
        auxM = M(x<=b);
        auxM = auxM(end);
        % Valor do esforço normal
        M(x>b) = auxM;
    end
end
end