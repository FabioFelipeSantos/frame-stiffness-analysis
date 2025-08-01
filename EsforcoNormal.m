%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Função que calcula o esforço normal nos elementos
%
% Entrada: x - Coordenada x da viga que se deseja o valor esforço normal;
%          cargas - vetor com as cargas aplicadas sobre o elemento;
%          coordI - vetor com as coordenadas x e y do ponto inicial do
%                   elemento;
%          f1x - Valor do esforço normal local no primeiro nó do elemento;
%
% Saída: N - vetor com os valores dos esforços normais.
%
% Autor: Fábio Felipe dos Santos
% Data: 03/05/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = EsforcoNormal(x, cargas, coordI, f1x, L)
% Inicialização do vetor de esforços normais
N = f1x * ones(length(x), 1);

% Verifica se há cargas concentradas no elemento
%auxCarga = cargas(cargas(:, 2) == 1, :);

% Inicializa os valores de a e b para a carga distribuída e c para a carga
% concentrada
a = []; b = []; c = [];

% Determina o valor c desde o ponto inicial de onde a carga concentrada
% está aplicada
% if ~isempty(auxCarga)
%     coordCargaConcentrada = auxCarga([4, 5]);
%     % Coordenada c da carga concentrada
%     c = sqrt(sum((coordCargaConcentrada - coordI).^2));
%     % Valor da carga concentrada
%     P = auxCarga(8) * cos(deg2rad(auxCarga(10)));
% end

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
    w1 = auxCarga(8) * cos(deg2rad(auxCarga(10)));
    if isnan(auxCarga(9))
        % Se o valor da carga final for vazio então trata-se de uma carga
        % uniforme
        w2 = w1;
    else
        % Caso contrário, a carga é linear
        w2 = auxCarga(9) * cos(deg2rad(auxCarga(10)));
    end
end

% Verifica se ambos a e c são vazios, se for, não há aplicação de carga
% concentrada ou distribuída
if isempty(a) && isempty(c)
    % Logo, o valor do esforço normal será igual ao valor da força
    % horizontal local do nó 1 da viga
    return;
elseif isempty(a)
    % Caso somente a seja vazio, então há carga concentrada, dessa forma o
    % esforço normal irá sofrer uma variação no ponto de aplicação da carga
    N(x >= c) = N(x >= c) - P;
elseif isempty(c)
    % Nesse caso existe uma carga distribuída normal ao longo do elemento
    auxx = x(sum([x>=a, x<=b], 2) == 2);
    % Valores da carga distribuída
    auxInt = ((a - auxx).*(a*(w1 + w2) - 2*b*w1 + (w1 - w2)*auxx))/(2*(a - b));
    
    % Valor da carga normal de a até b
    N(sum([x>=a, x<=b], 2) == 2) = N(sum([x>=a, x<=b], 2) == 2) + auxInt;
    
    % Se caso b for menor que L então o restante do esforço normal será
    % igual ao valor de N em b
    if b < L
        % Armazena o valor de N em b
        auxN = N(x<=b);
        auxN = auxN(end);
        % Valor do esforço normal
        N(x>b) = auxN;
    end
else
    % Nesse caso ambos a e c são não vazios, logo há os dois tipos de
    % cargas, distribuídas e concentradas
    if c <= a
        % A carga concentrada está antes da carga distribuída, logo da
        % aplicação da carga concentrada até o inicío da carga distribuída
        % aplicamos a carga concentrada
        N(sum([x>=c,x<a], 2)==2) = N(sum([x>=c,x<a], 2)==2) - P;

        % Armazena os valores de N para pegar o último valor calculado
        auxN = N(x<a);
        auxN = auxN(end);

        % Valores de x que serão utilizados no cálculo da integral da carga
        % distribuída
        auxx = x(sum([x>=a, x<=b], 2) == 2);
        % Valores da carga distribuída
        auxInt = ((a - auxx).*(a*(w1 + w2) - 2*b*w1 + (w1 - w2)*auxx))/(2*(a - b));
        
        % Armazena os valores da normal de a até b
        N(sum([x>=a, x<=b], 2) == 2) = auxN + auxInt;

        % Verifica se o valor de b é igual a L
        if b < L
            % Se não for, b é menor que L, dessa forma a carga distribuída
            % não vai até o fim do elemento. Logo, coloca-se a carga
            % somente até o valor de b
            % Para os valores de N onde x é maior que b e menor que L
            % consideramos o valor de N em b
            auxN = N(x <= b);
            N(x>b) = auxN(end);
        end
    elseif c <= b
        % Nesse caso a carga concentrada está aplicada dentro da carga
        % distribuída
        % Calcula a integral para valores de a até c
        auxx = x(sum([x>=a, x<c], 2)==2);
        auxInt = ((a - auxx).*(a*(w1 + w2) - 2*b*w1 + (w1 - w2)*auxx))/(2*(a - b));
        N(sum([x>=a, x<c], 2)==2) = N(sum([x>=a, x<c], 2)==2) + auxInt;
        
        % Recolhe o último valor da Normal (no caso em c)
        auxN = N(x < c);
        auxN = auxN(end);
        
        % Calcula a integral para valores de c até b
        auxx = x(sum([x>=c, x<=b], 2)==2);
        auxInt = ((a - auxx).*(a*(w1 + w2) - 2*b*w1 + (w1 - w2)*auxx))/(2*(a - b));
        
        % Calcula o valor da normal de c até b
        N(sum([x>=c, x<=b], 2)==2) = auxN + auxInt - P;

        if b < L
            % Se não for, b é menor que L, dessa forma a carga distribuída
            % não vai até o fim do elemento. Logo, coloca-se a carga
            % somente até o valor de b
            % Para os valores de N onde x é maior que b e menor que L
            % consideramos o valor de N em b
            auxN = N(x <= b);
            N(x>b) = auxN(end);
        end
    else
        % No caso c é maior que a e b, ou seja, a carga concentrada está
        % aplicada depois da carga distribuída
        % Calcula a integral para valores de a até b
        auxx = x(sum([x>=a, x<=b], 2)==2);
        auxInt = ((a - auxx).*(a*(w1 + w2) - 2*b*w1 + (w1 - w2)*auxx))/(2*(a - b));
        N(sum([x>=a, x<=b], 2)==2) = N(sum([x>=a, x<=b], 2)==2) + auxInt;

        % De b até c mantém o valor de N em b
        auxN = N(x<=b);
        auxN = auxN(end);
        N(sum([x>=b, x<c], 2)==2) = auxN;

        % Após c aplica a carga concentrada
        N(x>=c) = auxN - P;
    end
end
end





















