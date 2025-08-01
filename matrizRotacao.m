function T = matrizRotacao(cs)
% Recolhe o seno e o cosseno do elemento
c = cs(1); s = cs(2);

% Monta a matriz de rotação
T = [c s 0 0 0 0;
    -s c 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 c s 0;
     0 0 0 -s c 0;
     0 0 0 0 0 1;
     ];
end