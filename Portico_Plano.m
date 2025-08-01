%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An�lise de p�rticos via M�todo da Rigidez
%
% No arquivo P�rticos as cargas s�o dos seguintes tipos
%   1 - Concentrada;
%   2 - Distribu�da;
%   3 - Momento concentrado;
% Para a dire��o das cargas use: 1 - Dire��o X, 2 - Dire��o Y
%
% Autor: F�bio Felipe dos Santos
% Data: 17/04/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

file = fopen('Dados Estrutura.txt','w');

%% Leitura dos dados da estrutura
nome = 'P�rticos4.xlsx';
nos = readmatrix(nome, 'Sheet', 'N�s');
vigas = readmatrix(nome, 'Sheet', 'Elementos (Vigas)');
cargas = readmatrix(nome, 'Sheet', 'For�as');

%% Propriedades geom�tricas dos elementos finitos da estrutura
% N�mero de vigas
numEl = size(vigas, 1);

% N�mero de n�s da estrutura
numNos = size(nos, 1);

% N�mero de graus de liberdade da estrutura
numGL = numNos * 3;

% Recolhe as coordenadas dos n�s da estrutura
coord = nos(:, [2, 3]);

% C�lculo das coordenadas da malha de elemento finitos
compViga = [];      % Comprimento das vigas
cosSenViga = [];    % Cossenos e Senos da viga [c, s]
iconec = [];        % Conectivadade dos elementos (n�s de liga��o)

for elem = 1:numEl
    % Recolhe os n�s inicial e final da viga
    noI = vigas(elem, 2); noF = vigas(elem, 3);
    
    % Recolhe as coordenadas do ponto inicial e final da viga
    coorI = nos(noI, [2,3]); coorF = nos(noF, [2,3]);
    
    % Calcula o comprimento da viga
    compViga = [compViga; sqrt(sum((coorF - coorI).^2))];
    
    % Calcula os cossenos e senos do �ngulo da viga
    c = (coorF(1) - coorI(1)) / compViga(elem);
    s = (coorF(2) - coorI(2)) / compViga(elem);
    
    cosSenViga = [cosSenViga; [c, s]];

    % Monta a conectividade de cada elemento (viga) dos finitos
    % Armazena as conectividades do elemento (viga) em tr�s n�s cada
    iconec = [iconec; [vigas(elem, 2), vigas(elem, 3)]];
end

%% Processamento: Cria��o e c�lculo da matriz do sistema

% Para a estrutura
f = zeros(numGL, 1);       % Vetor de for�as nodais
K = zeros(numGL);          % Matriz de Rigidez
Qf = f;                    % Vetor de momentos de extremo fixo para a extrura
Kelem = zeros(6, 6, numEl);
T = Kelem;
KelemT = Kelem;
% Armazena somente os n�meros dos GL do esfor�o axial
% numGlElemento_Normal = (1:3:9)';
% 
% Calcula os n�meros dos GL da transla��o vertical e rota��o
% numGlElemento_Viga = setdiff((1:9)', numGlElemento_Normal);

% Montagem da Matriz de Rigidez da Estrutura
for elem = 1:numEl
    % Recolhe as caracter�sticas de material da viga
    E = vigas(elem, 4); A = vigas(elem, 5); I = vigas(elem, 6);

    % C�lculo dos graus de liberdade do elemento (no caso 3)
    glElem = iconec(elem, :);      % Numera��o dos graus de liberdade do elemento

    % Monta o vetor com a numera��o dos graus de liberdade do elemento
    numGlElemento(:,3) = glElem' * 3;
    numGlElemento(:,2) = numGlElemento(:,3) - 1;
    numGlElemento(:,1) = numGlElemento(:,2) - 1;

    % Armazena todos os GL do elemento
    numGlElemento = reshape(numGlElemento', [numel(numGlElemento), 1]);
    %numGlElemento = sort(numGlElemento);

    nn = 3*length(glElem);        % N�mero de graus de liberdade por elemento

    % Matriz de Rigidez para o elemento de p�rtico plano
    Kelem(:, :, elem) = matRigVigaPortico(E, A, I, compViga(elem));
    
    fprintf(file,'Matrix de Rigidez do elemento %d\n', elem);
    fprintf(file,'-------------------------------------------------------------------------------------------\n');
    for i = 1:6
        fprintf(file, '\t');
        fprintf(file, '%15.7e ', Kelem(i, :, elem));
        fprintf(file, '\n');
    end
    fprintf(file,'-------------------------------------------------------------------------------------------\n');
    % Recebe a matriz e rota��o da viga
    T(:,:,elem) = matrizRotacao(cosSenViga(elem, :));
    fprintf(file,'Matrix de Rota��o do elemento %d\n', elem);
    fprintf(file,'-------------------------------------------------------------------------------------------\n');
    for i = 1:6
        fprintf(file, '\t');
        fprintf(file, '%15.7e ', T(i,:,elem));
        fprintf(file, '\n');
    end
    KelemT(:,:,elem) = T(:,:,elem)' * Kelem(:, :, elem) * T(:,:,elem);
    fprintf(file,'Matrix de Rigidez Rotacionada do elemento %d\n', elem);
    fprintf(file,'-------------------------------------------------------------------------------------------\n');
    for i = 1:6
        fprintf(file, '\t');
        fprintf(file, '%15.7e ', KelemT(i, :, elem));
        fprintf(file, '\n');
    end
    fprintf(file,'-------------------------------------------------------------------------------------------\n');
    % Composi��o da matriz de rigidez global
    K(numGlElemento, numGlElemento) = K(numGlElemento, numGlElemento) + KelemT(:,:,elem);

    clear numGlElemento;
end

fprintf(file,'Matrix de Rigidez da Estrutura\n');
fprintf(file,'-------------------------------------------------------------------------------------------\n');
for i = 1:size(K,1)
    fprintf(file, '\t');
    fprintf(file, '%15.7e ', K(i, :));
    fprintf(file, '\n');
end
fprintf(file,'-------------------------------------------------------------------------------------------\n');

%% Aplica��o das for�as de extremo fixo dos elementos da estrutura
qf = zeros(6, 1, numEl);

% Percorre cada uma das cargas
for i = 1:size(cargas, 1)
    if cargas(i, 2) == 2
        % Existe carregamentos distribu�dos!!!
        
        % Recolhe a viga que possui a carga distribu�da
        elem = cargas(i, 3);
        
        % Numera��o dos graus de liberdade do elemento
        glElem = iconec(elem, :);
        
        % Monta o vetor com a numera��o dos graus de liberdade do elemento
        numGlElemento(:,3) = glElem' * 3;
        numGlElemento(:,2) = numGlElemento(:,3) - 1;
        numGlElemento(:,1) = numGlElemento(:,2) - 1;
        
        % Armazena todos os GL do elemento
        numGlElemento = reshape(numGlElemento', [numel(numGlElemento), 1]);
        %numGlElemento = sort(numGlElemento);

        % N�mero de graus de liberdade por elemento
        nn = 3*length(glElem);

        qf(:, :, elem) = momentosDeExtremoFixoCargaDistribuida(cargas(i, :),...
            compViga(elem), coord(vigas(elem,2), [1,2]));
        fprintf(file,'Vetor de For�as equivalentes de carregamento distribu�do (sinal contr�rio) do elemento %d\n', elem);
        fprintf(file,'-------------------------------------------------------------------------------------------\n');
        fprintf(file, '\t%15.7e\n', qf(:,:,elem));
        fprintf(file,'-------------------------------------------------------------------------------------------\n');
        % Recebe a matriz de rota��o da viga
        TT = matrizRotacao(cosSenViga(elem, :));
        
        % Armazena as for�as nos respectivos graus de liberdade do vetor
        % global
        Qf(numGlElemento) = Qf(numGlElemento) + TT'*qf(:, :, elem);
        fprintf(file,'Vetor de For�as equivalentes ROTACIONADO de carregamento distribu�do (sinal contr�rio) do elemento %d\n', elem);
        fprintf(file,'-------------------------------------------------------------------------------------------\n');
        fprintf(file, '\t%15.7e\n', TT'*qf(:,:,elem));
        fprintf(file,'-------------------------------------------------------------------------------------------\n');
        clear numGlElemento;
    end
end

%% Aplica��o de cargas concentradas ou momentos concentrados
% Verifica se h� cargas concentradas ou momentos concentrados
if any(cargas(:, 2) == 1) || any(cargas(:, 2) == 3)
    % Determina em quais linhas da planilha est�o as cargas concentradas ou
    % momentos
    aux = sort([find(cargas(:, 2) == 1); find(cargas(:, 2) == 3)]);

    for i = 1:size(aux, 1) % Percorre todas as cargas concentradas ou momentos
        if ~isnan(cargas(aux(i), 3))
            % A carga est� no meio de algum elemento da estrutura
            % Recolhe a viga que possui a carga distribu�da
            elem = cargas(aux(i), 3);

            % Numera��o dos graus de liberdade do elemento
            glElem = iconec(elem, :);

            % Monta o vetor com a numera��o dos graus de liberdade do elemento
            numGlElemento(:,3) = glElem' * 3;
            numGlElemento(:,2) = numGlElemento(:,3) - 1;
            numGlElemento(:,1) = numGlElemento(:,2) - 1;

            % Armazena todos os GL do elemento
            numGlElemento = reshape(numGlElemento', [numel(numGlElemento), 1]);
            %numGlElemento = sort(numGlElemento);
            
            auxqf = momentosDeExtremoFixoCargaConcentrada(cargas(aux(i), :),...
            compViga(elem), coord(vigas(elem,2), [1,2]));

            % Recebe a matriz de rota��o da viga
            TT = matrizRotacao(cosSenViga(elem, :));
            
            % Armazena as for�as nos respectivos graus de liberdade do vetor
            % global
            Qf(numGlElemento) = Qf(numGlElemento) + TT'*auxqf;
            qf(:, :, elem) = qf(:, :, elem) + auxqf;
            clear numGlElemento;
        else
            % A carga est� aplicada em um n� da estrutura
            % Determina o n� em que est� a carga
            no = cargas(aux(i), [4, 5]);
            
            if cargas(aux(i), 2) == 1
                % A carga � concentrada no ponto
                % Recebe o m�dulo da carga
                P = cargas(aux(i), 8);

                % Recebe a dire��o da carga aplicada
                theta = deg2rad(cargas(aux(i), 10));

                % Determina o grau de liberdade para carga na horizontal
                glCargaConc = 3*find(sum(coord == no, 2) == 2) - 2;
                % Armazena a carga na horizontal
                f(glCargaConc) = f(glCargaConc) + P * cos(theta);
                % Determina o grau de liberdade para carga na vertical
                glCargaConc = 3*find(sum(coord == no, 2) == 2) - 1;
                % Armazena a carga na horizontal
                f(glCargaConc) = f(glCargaConc) + P * sin(theta);
            else
                % Trata-se de um momento concentrado
                % Recebe o valor do momento
                M = cargas(aux(i), 8);

                % Determina o grau de liberdade para carga na horizontal
                glCargaConc = 3*find(sum(coord == no, 2) == 2);
                % Armazena o momento concentrado
                f(glCargaConc) = f(glCargaConc) + M;
            end
        end
    end
end
%% Aplica��o das condi��es de contorno nos extremos e resolu��o do sistema Ku = f

% Restri��es na dire��o x aplicadas aos n�s da estrutura 
coordRestringidas = nos(nos(:,4) == 1, [2, 3]);
% Recolhe os n�s da malha com essas coordenadas
nosMalhaRestringidos = zeros(size(coordRestringidas, 1), 1);
for i = 1:size(coordRestringidas, 1)
    nosMalhaRestringidos(i) = find(sum(coord == coordRestringidas(i, :), 2) == 2);
end

glRestringidos = 3*nosMalhaRestringidos - 2;

% Restri��es na dire��o y aplicadas aos n�s da estrutura 
coordRestringidas = nos(nos(:,5) == 1, [2, 3]);
% Recolhe os n�s da malha com essas coordenadas
nosMalhaRestringidos = zeros(size(coordRestringidas, 1), 1);
for i = 1:size(coordRestringidas, 1)
    nosMalhaRestringidos(i) = find(sum(coord == coordRestringidas(i, :), 2) == 2);
end

glRestringidos = [glRestringidos; 3*nosMalhaRestringidos - 1];

% Restri��es de momento aplicados aos n�s da estrutura 
coordRestringidas = nos(nos(:,6) == 1, [2, 3]);
if ~isempty(coordRestringidas)
    % Recolhe os n�s da malha com essas coordenadas
    nosMalhaRestringidos = zeros(size(coordRestringidas, 1), 1);
    for i = 1:size(coordRestringidas, 1)
        nosMalhaRestringidos(i) = find(sum(coord == coordRestringidas(i, :), 2) == 2);
    end
    
    glRestringidos = [glRestringidos; 3*nosMalhaRestringidos];
end

glRestringidos = sort(glRestringidos);

% Armazena os graus de liberdade sem restri��o
glAtivos = setdiff((1:numGL)', glRestringidos);

% C�lculo dos deslocamentos
u = zeros(numGL, 1);
%det(K(glAtivos, glAtivos))
fprintf(file,'Vetor de For�as\n');
fprintf(file,'-------------------------------------------------------------------------------------------\n');
fprintf(file, '\t%15.7e\n', f - Qf);
fprintf(file,'-------------------------------------------------------------------------------------------\n');
u(glAtivos) = K(glAtivos, glAtivos) \ (f(glAtivos) - Qf(glAtivos));

fprintf(file,'Vetor de Deslocamentos Nodais\n');
fprintf(file,'-------------------------------------------------------------------------------------------\n');
fprintf(file, '\t%15.7e\n', u);
fprintf(file,'-------------------------------------------------------------------------------------------\n');

freacao = zeros(numGL, 1);
freacao(glRestringidos) = K(glRestringidos, :)*u + Qf(glRestringidos);

fprintf(file,'Vetor de For�as de Rea��o\n');
fprintf(file,'-------------------------------------------------------------------------------------------\n');
fprintf(file, '\t%15.7e\n', freacao);
fprintf(file,'-------------------------------------------------------------------------------------------\n');

fclose(file);
%% Mostra as rea��es no console
nosComReacao = nos(sum(nos(:, [4, 5, 6]) == 1, 2) ~= 0, 1);
for i = 1:length(nosComReacao)
    fprintf('As for�as de rea��o do n� %d s�o:\n', nosComReacao(i))
    auxReacao = nos(nosComReacao(i), [4, 5, 6]);
    if auxReacao(1) == 1
        Rx = freacao(nosComReacao(i) * 3 - 2);
        fprintf('Rx = %.8f\n', Rx);
    end
    if auxReacao(2) == 1
        Ry = freacao(nosComReacao(i) * 3 - 1);
        fprintf('Ry = %.8f\n', Ry);
    end 
    if auxReacao(3) == 1
        M = freacao(nosComReacao(i) * 3);
        fprintf('M = %.8f\n', M);
    end
    fprintf('\n')
end

%% Desenho das curvas de deforma��o da estrutura
% % Desenho da curva exata com a aproxima��o
numNosInterp = 50;  % N�mero de pontos para o desenho

% Fator de escala para o desenho
escala = 0.01;

figure('units','normalized','outerposition',[0 0 1 1])
title('Curva de deslocamento da Estrutura','Interpreter','Latex','FontSize',24)
ax1 = axes;
hold on
box on
grid on
ax1.FontSize = 18;
ax1.TickLabelInterpreter = 'Latex';

% Vetor de for�as locais dos elementos
fLocal = zeros(6, 1, size(vigas, 1));
for elem = 1:size(vigas, 1)
    % Recolhe os n�s do elemento
    glElem = iconec(elem, :);

    % Monta o vetor com a numera��o dos graus de liberdade do elemento
    numGlElemento(:,3) = glElem' * 3;
    numGlElemento(:,2) = numGlElemento(:,3) - 1;
    numGlElemento(:,1) = numGlElemento(:,2) - 1;

    % Armazena todos os GL do elemento
    numGlElemento = reshape(numGlElemento', [numel(numGlElemento), 1]);

    % Interpola os pontos para desenho para xi entre -1 e 1
    xi = (0:(compViga(elem)/(numNosInterp-1)):compViga(elem))';
    %xi = linspace(0, compViga(elem),numNosInterp)';
    
    % Calcula os pontos ao longo da viga que foram interpolados
    vetor = (1/compViga(elem)) * (coord(vigas(elem, 3), :) - coord(vigas(elem, 2), :));
    pontos = coord(vigas(elem, 2), :) + vetor .* xi;
    
    % Calcula o vetor de deslocamentos nas coordenadas locais
    T = matrizRotacao(cosSenViga(elem, :)); 
    uLocal = T * u(numGlElemento);

    % Calcula as fun��es de forma do elemento
    [N, dNdxi] = FuncoesDeFormaViga(xi, compViga(elem));    

    % Calcula a deflex�o da viga
    v = N * uLocal([2;3;5;6]);

    % Para a barra, calcula o deslocamento
    Nbarra = FuncoesDeFormaBarra(xi, compViga(elem));
    uu = Nbarra * uLocal([1;4]);
    
    % Pontos deslocados
    Taux = [cosSenViga(elem, 1), cosSenViga(elem, 2); -cosSenViga(elem, 2) cosSenViga(elem, 1)];

    pontos = pontos + escala * (Taux'*[uu, v]')';

    % Faz o desenho da viga
    xElem = [coord(vigas(elem, 2), 1) coord(vigas(elem, 3), 1)];
    yElem = [coord(vigas(elem, 2), 2) coord(vigas(elem, 3), 2)];
    plot(ax1, xElem, yElem,'--b', 'LineWidth', 2.5)
    plot(ax1, pontos(:, 1), pontos(:, 2),'-r', 'LineWidth', 2)

    % Desenho do esfor�o normal na estrutura
    aux = cargas(cargas(:, 3) == elem, :);

    fLocal(:, 1, elem) = Kelem(:, :, elem) * uLocal + qf(:, 1, elem);
    
    clear numGlElemento
end
%% Desenho do esfor�o normal da estrutura
fprintf('--------------------------------------------------\n')
% Determina os valores normais m�ximos e m�nimos para se fazer uma escala
% para o desenho do esfor�o normal
auxNormalMin = min(min(abs(fLocal([1;4], 1, :))));
auxNormalMax = max(max(abs(fLocal([1;4], 1, :))));
auxNormalEscala = max(auxNormalMin, auxNormalMax);

% Comprimento dos pontos para escala do eixo no desenho do esfor�o normal
escalaEsforcos = 1;

numNosInterp = 10;  % N�mero de pontos para o desenho

figure('units','normalized','outerposition',[0 0 1 1])
title('Esfor\c{c}o normal','Interpreter','Latex','FontSize',28)
ax1 = axes;
hold on
box on
grid on
ax1.FontSize = 18;
ax1.TickLabelInterpreter = 'Latex';
for elem = 1:size(vigas, 1)
    % Recolhe os n�s do elemento
    glElem = iconec(elem, :);

    % Interpola os pontos para desenho para xi entre -1 e 1
    xi = (0:(compViga(elem)/(numNosInterp-1)):compViga(elem))';
    %xi = linspace(0, compViga(elem), numNosInterp)';
    
    % Calcula os pontos ao longo da viga que foram interpolados
    vetor = (1/compViga(elem)) * (coord(vigas(elem, 3), :) - coord(vigas(elem, 2), :));
    pontos = coord(vigas(elem, 2), :) + vetor .* xi;
    
    % Calcula os esfor�os normais no elemento
    N = EsforcoNormal(xi, cargas(cargas(:, 3)==elem, :),...
        coord(vigas(elem, 2), :), -fLocal(1, 1, elem), compViga(elem));
    
    % Faz a mudan�a de escala para o gr�fico do esfor�o normal
    Nescala = (escalaEsforcos / auxNormalEscala) * N;

    % Determina o �ngulo que o elemento faz com o eixo x global
    if cosSenViga(elem, 1) ~= 0
        theta = atan(cosSenViga(elem, 2) / cosSenViga(elem, 1));
    else
        theta = pi/2;
    end

    % Determina o vetor unit�rio que dar� a dire��o do ponto normal para o
    % desenho do esfor�o normal
    alpha = theta + pi/2;
    vetorNormal = [cos(alpha), sin(alpha)];

    % Pontos do elemento como local para o esfor�o normal
    pontosNormais = pontos + Nescala.*vetorNormal;

    % Faz o desenho da viga
    xElem = [coord(vigas(elem, 2), 1) coord(vigas(elem, 3), 1)];
    yElem = [coord(vigas(elem, 2), 2) coord(vigas(elem, 3), 2)];
    plot(xElem, yElem,'-b', 'LineWidth', 1.5)
    % Faz o desenho da curva normal
    plot(pontosNormais(:, 1), pontosNormais(:, 2),'-r', 'LineWidth', 2)
    % Faz o desenho dos segmentos que ligam os pontos iniciais e finais do
    % elemento ao valor do esfor�o normal
    plot([xElem(1); pontosNormais(1, 1)], [yElem(1); pontosNormais(1, 2)],'-r', 'LineWidth', 2)
    plot([xElem(2); pontosNormais(end, 1)], [yElem(2); pontosNormais(end, 2)],'-r', 'LineWidth', 2)
    
    % Faz a sa�da dos valores no console    
    fprintf('O valor da normal para o elemento %d �:\n', elem)
    fprintf('N[%d] = %.8f\n', glElem(1), N(1))
    fprintf('N[%d] = %.8f\n', glElem(2), N(end))
end

%% Desenho do esfor�o cortante da estrutura
fprintf('\n--------------------------------------------------\n')
% Determina os valores normais m�ximos e m�nimos para se fazer uma escala
% para o desenho do esfor�o normal
auxCortanteMin = min(min(abs(fLocal([2;5], 1, :))));
auxCortanteMax = max(max(abs(fLocal([2;5], 1, :))));
auxCortanteEscala = max(auxCortanteMin, auxCortanteMax);

% Comprimento dos pontos para escala do eixo no desenho do esfor�o normal
escalaEsforcos = 1;

numNosInterp = 10;  % N�mero de pontos para o desenho

figure('units','normalized','outerposition',[0 0 1 1])
title('Esfor\c{c}o Cortante','Interpreter','Latex','FontSize',28)
ax1 = axes;
hold on
box on
grid on
ax1.FontSize = 18;
ax1.TickLabelInterpreter = 'Latex';
for elem = 1:size(vigas, 1)
    % Recolhe os n�s do elemento
    glElem = iconec(elem, :);

    % Interpola os pontos para desenho
    xi = (0:(compViga(elem)/(numNosInterp-1)):compViga(elem))';
    
    % Calcula os pontos ao longo da viga que foram interpolados
    vetor = (1/compViga(elem)) * (coord(vigas(elem, 3), :) - coord(vigas(elem, 2), :));
    pontos = coord(vigas(elem, 2), :) + vetor .* xi;
    
    % Calcula os esfor�os normais no elemento
    V = EsforcoCortante(xi, cargas(cargas(:, 3)==elem, :),...
        coord(vigas(elem, 2), :), fLocal(2, 1, elem), compViga(elem));
    
    % Faz a mudan�a de escala para o gr�fico do esfor�o normal
    Vescala = (escalaEsforcos / auxCortanteEscala) * V;

    % Determina o �ngulo que o elemento faz com o eixo x global
    if cosSenViga(elem, 1) ~= 0
        theta = atan(cosSenViga(elem, 2) / cosSenViga(elem, 1));
    else
        theta = pi/2;
    end

    % Determina o vetor unit�rio que dar� a dire��o do ponto normal para o
    % desenho do esfor�o normal
    alpha = theta + pi/2;
    vetorNormal = [cos(alpha), sin(alpha)];

    % Pontos do elemento como local para o esfor�o normal
    pontosCortante = pontos + Vescala.*vetorNormal;

    % Faz o desenho da viga
    xElem = [coord(vigas(elem, 2), 1) coord(vigas(elem, 3), 1)];
    yElem = [coord(vigas(elem, 2), 2) coord(vigas(elem, 3), 2)];
    plot(xElem, yElem,'-b', 'LineWidth', 1.5)
    % Faz o desenho da curva normal
    plot(pontosCortante(:, 1), pontosCortante(:, 2),'-r', 'LineWidth', 2)
    % Faz o desenho dos segmentos que ligam os pontos iniciais e finais do
    % elemento ao valor do esfor�o normal
    plot([xElem(1); pontosCortante(1, 1)], [yElem(1); pontosCortante(1, 2)],'-r', 'LineWidth', 2)
    plot([xElem(2); pontosCortante(end, 1)], [yElem(2); pontosCortante(end, 2)],'-r', 'LineWidth', 2)
    
    % Faz a sa�da dos valores no console
    fprintf('O valor do esfor�o cortante para o elemento %d �:\n', elem)
    fprintf('V[%d] = %.8f\n', glElem(1), V(1))
    fprintf('V[%d] = %.8f\n', glElem(2), V(end))
end

%% Desenho do esfor�o momento fletor da estrutura
fprintf('\n--------------------------------------------------\n')
% Determina os valores normais m�ximos e m�nimos para se fazer uma escala
% para o desenho do esfor�o normal
auxMomentoMin = min(min(abs(fLocal([3;6], 1, :))));
auxMomentoMax = max(max(abs(fLocal([3;6], 1, :))));
auxMomentoEscala = max(auxMomentoMin, auxMomentoMax);

% Comprimento dos pontos para escala do eixo no desenho do esfor�o normal
escalaEsforcos = 1;

numNosInterp = 20;  % N�mero de pontos para o desenho

figure('units','normalized','outerposition',[0 0 1 1])
title('Esfor\c{c}o Momento Fletor','Interpreter','Latex','FontSize',28)
ax1 = axes;
hold on
box on
grid on
ax1.FontSize = 18;
ax1.TickLabelInterpreter = 'Latex';
for elem = 1:size(vigas, 1)
    % Recolhe os n�s do elemento
    glElem = iconec(elem, :);

    % Interpola os pontos para desenho
    xi = (0:(compViga(elem)/(numNosInterp-1)):compViga(elem))';
    
    % Calcula os pontos ao longo da viga que foram interpolados
    vetor = (1/compViga(elem)) * (coord(vigas(elem, 3), :) - coord(vigas(elem, 2), :));
    pontos = coord(vigas(elem, 2), :) + vetor .* xi;
    
    % Calcula os esfor�os normais no elemento
    M = EsforcoMomento(xi, cargas(cargas(:, 3)==elem, :),...
        coord(vigas(elem, 2), :), -fLocal(3, 1, elem), fLocal(2, 1, elem),...
        compViga(elem));
    
    % Faz a mudan�a de escala para o gr�fico do esfor�o normal
    Mescala = (escalaEsforcos / auxMomentoEscala) * M;

    % Determina o �ngulo que o elemento faz com o eixo x global
    if cosSenViga(elem, 1) ~= 0
        theta = atan(cosSenViga(elem, 2) / cosSenViga(elem, 1));
    else
        theta = pi/2;
    end

    % Determina o vetor unit�rio que dar� a dire��o do ponto normal para o
    % desenho do esfor�o normal
    alpha = theta + pi/2;
    vetorNormal = [cos(alpha), sin(alpha)];

    % Pontos do elemento como local para o esfor�o normal
    pontosMomento = pontos - Mescala.*vetorNormal;

    % Faz o desenho da viga
    xElem = [coord(vigas(elem, 2), 1) coord(vigas(elem, 3), 1)];
    yElem = [coord(vigas(elem, 2), 2) coord(vigas(elem, 3), 2)];
    plot(xElem, yElem,'-b', 'LineWidth', 1.5)
    % Faz o desenho da curva normal
    plot(pontosMomento(:, 1), pontosMomento(:, 2),'-r', 'LineWidth', 2)
    % Faz o desenho dos segmentos que ligam os pontos iniciais e finais do
    % elemento ao valor do esfor�o normal
    plot([xElem(1); pontosMomento(1, 1)], [yElem(1); pontosMomento(1, 2)],'-r', 'LineWidth', 2)
    plot([xElem(2); pontosMomento(end, 1)], [yElem(2); pontosMomento(end, 2)],'-r', 'LineWidth', 2)

    % Faz a sa�da dos valores no console
    fprintf('O valor do momento fletor para o elemento %d �:\n', elem)
    fprintf('M[%d] = %.8f\n', glElem(1), M(1))
    fprintf('M[%d] = %.8f\n', glElem(2), M(end))
end

































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Faz o desenho do esfor�o normal
%     plot(ax5, xx, Nx, '-c', 'LineWidth', 3)
%     title(ax5, 'Esfor�o Normal', 'FontSize', 26)
%     
%     % Calcula as fun��es de forma e derivada do elemento
%     [N, dNdxi, B, dBdxi] = funcoesForma(xi,compElemento);
%     B = (invJ)^2*B;
%     dBdxi = (invJ)^3*dBdxi;
%     
%     % Calcula os deslocamentos interpolados ao longo dos elementos
%     uu = N * u(numGlElemento);
%     
%     % Faz o desenho do deslocamento colorido com a escala da tens�o
%     plot(ax1, xx, uu, '-b', 'LineWidth', 3)
%     title(ax1, 'Deslocamento u(x)', 'FontSize', 26)
%     
%     % Calcula a rota��o na barra interpolada ao longo dos elementos
%     duudx = (1 / detJ) * dNdxi * u(numGlElemento);
%     plot(ax2, xx, duudx,'-b','LineWidth',2)
%     title(ax2, 'Rota��o dos Eixos', 'FontSize', 26)
%     
%     % Calcula o momento fletor aplicado em cada ponto da estrutura
%     M = E*I*B*u(numGlElemento);
%     plot(ax3, xx, M,'-r','LineWidth',2)
%     title(ax3, 'Momento Fletor', 'FontSize', 26)
%     
%     % Calcula o cortante aplicado sobre a estrutura
%     V = E*I*dBdxi*u(numGlElemento);
%     plot(ax4, xx, V,'-m','LineWidth',2)
%     title(ax4, 'Cortante', 'FontSize', 26)
%     
%     % Calcula o momento fletor aplicado em cada ponto da estrutura
%     sigma_max =  (ymax / I) * M;
%     patch(ax6, [xx; NaN], [sigma_max; NaN], [sigma_max; NaN],'LineWidth',2,'EdgeColor','interp')
%     title(ax6, 'Tens�o m�xima na se��o transversal', 'FontSize', 26)