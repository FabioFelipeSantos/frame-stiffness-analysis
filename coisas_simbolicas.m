%% Funções de forma para elemento de barra com 2 nós
clear
clc
syms x L u1 u2 u(x) a1 a2 E A

u(x) = a1*x + a2;

eq1 = u(0) == u1;
eq2 = u(L) == u2;

aS = solve([eq1, eq2], [a1, a2]);

a1 = aS.a1;
a2 = aS.a2;

u(x) = subs(u)

N1(x, L) = 1 - x/L;
N2(x, L) = x/L;

N = [N1, N2]
dN = diff(N, x, 1)

K = collect(int(E * A * transpose(dN) * dN, x, 0, L), E*A/L)

%% Momento de extremo fixo para carga distribuída em uma parte da estrutura
clear
clc

syms w a L

exp1 = (w * a^2 * (3*a^2 - 8*a*L + 6*L^2))/(12*L^2);
exp2 = (w*a^3*(4*L - 3*a)) / (12*L^2);
exp3 = w * a * (a/2);

Rb = (-exp1 + exp2 + exp3)/L;
pretty(simplify(Rb))

Ra = -Rb + w*a;
pretty(simplify(Ra))

%% Momento de extremo fixo para carga distribuída linear em uma parte da estrutura
clear
clc

syms a b w0 w1 x L

w(x) = (b*w0 - a*w1 + (w1 - w0)*x)/(b-a);

Ma = -(1/L^2) * int(w * x * (L-x)^2, x, a, b);
disp("Ma")
pretty(simplify(Ma))

Mb = (1/L^2) * int(w * x^2 * (L-x), x, a, b);
disp("Mb")
pretty(simplify(Mb))

% Ma - Mb + Rb*L - w0*(b-a)*(a + (b-a)/2) - ((w1-w0)*(b-a)/2)*(a + (2/3)*(b-a))=0
Rb = -(1/L) * (Ma + Mb + w0*(b-a)*(a + (b-a)/2) + ((w1-w0)*(b-a)/2)*(a + (2/3)*(b-a)));
disp("Rb")
pretty(simplify(Rb))

Ra = - Rb - ((w1 + w0)*(b-a)/2);
disp("Ra")
pretty(simplify(Ra))

Maa = subs(Ma, [a,b,L,w0,w1], [1,4,4,-20,0])
Mbb = subs(Mb, [a,b,L,w0,w1], [1,4,4,-20,0])
Raa = subs(Ra, [a,b,L,w0,w1], [1,4,4,-20,0])
Rbb = subs(Rb, [a,b,L,w0,w1], [1,4,4,-20,0])

%% Momento de extremo fixo para carga distribuída linear normal em uma parte da estrutura
clear
clc

syms a b w0 w1 x L

w(x) = (b*w0 - a*w1 + (w1 - w0)*x)/(b-a);

Ra = (1/L) * int(w*(L-x), a, b);
disp("Ra")
pretty(simplify(Ra))

Rb = - Ra + int(w, a, b);
disp("Rb")
pretty(simplify(Rb))

double(subs(Ra, [a,b,L,w0,w1], [1,4,6,10,20]))
double(subs(Rb, [a,b,L,w0,w1], [1,4,6,10,20]))

%% Momento de extremo fixo para carga concentrada em uma parte da estrutura
clear
clc

syms a b P L

Ma = (P*a*b^2) / L^2;
Mb = (P*a^2*b) / L^2;

Rb = (1/L) * (-Ma + Mb + P*a);
pretty(simplify(Rb))
pretty(simplify(subs(Rb, L, a+b)))


Ra = P - Rb;
pretty(simplify(Ra))
pretty(simplify(subs(Ra, L, a+b)))

double(subs(Ra, [a,b,L,P], [1,5,6,20]))
double(subs(Rb, [a,b,L,P], [1,5,6,20]))
double(subs(Ma, [a,b,L,P], [1,5,6,20]))
double(subs(Mb, [a,b,L,P], [1,5,6,20]))

%% Momento de extremo fixo para momento
clear
clc

syms a b M L

Ma = (M*b*(2*a-b)) / L^2;
Mb = (M*a*(2*b-a)) / L^2;

Rb = (1/L) * (-Ma - Mb - M);
pretty(simplify(subs(Rb, L, a+b)))

Ra = - Rb;
pretty(simplify(subs(Ra, L, a+b)))

double(subs(Ra, [a,b,L,M], [1,5,6,20]))
double(subs(Rb, [a,b,L,M], [1,5,6,20]))
double(subs(Ma, [a,b,L,M], [1,5,6,20]))
double(subs(Mb, [a,b,L,M], [1,5,6,20]))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Funções de forma para elemento finito de viga
clear
clc

syms v1 v2 o1 o2 x L a1 a2 a3 a4 E I P M(x) V(x) q u a K v3 o3 t

u = [v1;o1;v2;o2];

v(x) = a1*x^3 + a2*x^2 + a3*x + a4;
eq1 = subs(v, x, 0) == v1;
eq2 = subs(diff(v, x, 1), x, 0) == o1;
eq3 = subs(v, x, L) == v2;
eq4 = subs(diff(v, x, 1), x, L) == o2;

S = solve([eq1, eq2, eq3, eq4], [a1, a2, a3, a4]);

v(x) = S.a1*x^3 + S.a2*x^2 + S.a3*x + S.a4;

vSimp = collect(v, [v1, v2, o1, o2]);

N1(x) = (2*x^3)/L^3 - (3*x^2)/L^2 + 1;
N2(x) = x - (2*x^2)/L + x^3/L^2;
N3(x) = (3*x^2)/L^2 - (2*x^3)/L^3;
N4(x) = x^3/L^2 - x^2/L;

N(x) = [N1 N2 N3 N4];

B(x) = diff(N, x, 2);

Naux1 = subs(N, L, a);
Baux1 = diff(Naux1, x, 2);
Naux2 = subs(N, L, L - a);
Baux2 = diff(Naux2, x, 2);

K = sym(zeros(6,6));

u = [v1;o1;v2;o2;v3;o3];

K(1:4, 1:4) = int(E*I*transpose(Baux1)*Baux1, x, 0, a);
K(3:6, 3:6) = K(3:6, 3:6) + int(E*I*transpose(Baux2)*Baux2, x, 0, L-a);

f = sym(zeros(6, 1));

%f(1:4) = int(q*transpose(Naux1),x,0,a);
%f(3:6) = f(3:6) + int(q*transpose(Naux2),x,0,L-a);
f(3) = f(3) + P;

uSol = sym(zeros(6, 1));

uSol(3:4) = K(3:4, 3:4)\f(3:4);

freac(a, L) = simplify(K*uSol - f)

M(x) = E*I*Baux1*uSol(1:4)
M(x) = subs(M, [a, P, L], [1, -50, 4])
M2(x) = E*I*Baux2*uSol(3:6)
M2(x) = subs(M2, [a, P, L], [1, -50, 4])

% 
% uS1 = Naux1 * uSol(1:4);
% uS1aux(x) = simplify(subs(uS1,[a,L,P,q,E,I],[3,6,-1,-1,205e6,(0.2*0.4^3)/12]));
% uS2 = Naux2 * uSol(3:6);
% uS2aux(x) = simplify(subs(uS2,[a,L,P,q,E,I],[3,6,-1,-1,205e6,(0.2*0.4^3)/12]));
% figure(1)
% fplot(uS1aux,[0,3])
% hold on
% fplot(subs(uS2aux,x,x-3),[3,6])
% hold off




%% Equações para o esforço normal nas vigas do pórtico
clear
clc
syms N x q L f1x a b w1 w2 w t

m = (w2 - w1) / (b - a);
q(x) = w1 + m*(x - a);

N(x) = simplify(int(-q(x), x, a, t))

simplify(subs(N, [w1, w2, a, b, x], [800, 100, 0, 8, 8]))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
syms a x y N

A = [x, y, 1; -N -a 1;N a 1]

pretty(collect(simplify(solve(det(A) == 0, y)), x))

%% Equações para o esforço cortante nas vigas do pórtico
clear
clc
syms M x q L f1y a b w1 w2 w t

m = (w2 - w1) / (b - a);
q(x) = w1 + m*(x - a);

M(x) = simplify(int(int(-q(x), x, a, x),x,a,x))

simplify(subs(M, [w1, w2, a, b, x], [800, 100, 0, 8, 8]))

%% Funções de Forma para elemento de viga com 3 nós
clear
clc

syms a1 a2 a3 a4 a5 a6 x L v1 v2 v3 m1 m2 m3 E I w0 w1 a b w P

% Definição da função de deflexão em função de seis ceficientes
v(x) = a1*x^5 + a2*x^4 + a3*x^3 + a4*x^2 + a5*x + a6;
df = diff(v, x, 1);     % Derivada primeira da deflexão v
% Equações para que v e dvdx satisfaçam os deslocamentos nodais
eq1 = v(0) == v1;       % Deflexão nó 1
eq2 = df(0) == m1;      % Rotação nó 1
eq3 = v(L/2) == v2;     % Deflexão nó 2
eq4 = df(L/2) == m2;    % Rotação nó 2
eq5 = v(L) == v3;       % Deflexão nó 3
eq6 = df(L) == m3;      % Rotação nó 3
% Resolve o sistema formado pelas equações anteriores
S = solve([eq1, eq2, eq3, eq4, eq5, eq6], [a1, a2, a3, a4, a5, a6]);
% Substitui os coeficientes a1, a2, a3, a4, a5 e a6 na função v
v(x) = subs(v, [a1,a2,a3,a4,a5,a6], [S.a1,S.a2,S.a3,S.a4,S.a5,S.a6]);
% Retira as funções coeficientes dos deslocamentos v1, m1, v2, m2, v3 e m3
collect(v,[v1,m1,v2,m2,v3,m3]);

% Recolhe as funções de forma para o elemento de viga com 3 nós
N1 = ((L + 6*x)*(L^2 - 3*L*x + 2*x^2)^2)/L^5;
N2 = (x*(L^2 - 3*L*x + 2*x^2)^2)/L^4;
N3 = (16*x^2*(L - x)^2)/L^4;
N4 = -(8*x^2*(L - x)^2*(L - 2*x))/L^4;
N5 = (x^2*(L - 2*x)^2*(7*L - 6*x))/L^5;
N6 = -(x^2*(L - x)*(L - 2*x)^2)/L^4;

% Monta a matriz de funções de forma N
N = [N1, N2, N3, N4, N5, N6];
dN = diff(N,x,1);       % Derivada das funções de forma

% Reescreve v e dvdx em função das funções de forma
v = N * [v1;m1;v2;m2;v3;m3];
dv = dN * [v1;m1;v2;m2;v3;m3];

% Determina a matriz B, derivada segunda de N
B = simplify(diff(N, x, 2));

% Calcula o produto B' * B para aferição
simplify(transpose(B) * B);

% Calcula a matriz de Rigidez K
K = E*I * simplify(int(transpose(B) * B, x, 0, L));

% Reescreve cada um dos 36 componentes da matriz de rigidez do elemento
K11 = (5092*E*I)/(35*L^3);
K12 = (1138*E*I)/(35*L^2);
K13 = -(512*E*I)/(5*L^3);
K14 = (384*E*I)/(7*L^2);
K15 = -(1508*E*I)/(35*L^3);
K16 = (242*E*I)/(35*L^2);
K21 = (1138*E*I)/(35*L^2);
K22 = (332*E*I)/(35*L);
K23 = -(128*E*I)/(5*L^2);
K24 = (64*E*I)/(7*L);
K25 = -(242*E*I)/(35*L^2);
K26 = (38*E*I)/(35*L);
K31 = -(512*E*I)/(5*L^3);
K32 = -(128*E*I)/(5*L^2);
K33 = (1024*E*I)/(5*L^3);
K34 = 0;
K35 = -(512*E*I)/(5*L^3);
K36 = (128*E*I)/(5*L^2);
K41 = (384*E*I)/(7*L^2);
K42 = (64*E*I)/(7*L);
K43 = 0;
K44 = (256*E*I)/(7*L);
K45 = -(384*E*I)/(7*L^2);
K46 = (64*E*I)/(7*L);
K51 = -(1508*E*I)/(35*L^3);
K52 = -(242*E*I)/(35*L^2);
K53 = -(512*E*I)/(5*L^3);
K54 = -(384*E*I)/(7*L^2);
K55 = (5092*E*I)/(35*L^3);
K56 = -(1138*E*I)/(35*L^2);
K61 = (242*E*I)/(35*L^2);
K62 = (38*E*I)/(35*L);
K63 = (128*E*I)/(5*L^2);
K64 = (64*E*I)/(7*L);
K65 = -(1138*E*I)/(35*L^2);
K66 = (332*E*I)/(35*L);

%% Momentos de extremo fixo de carga concentrada em elementos de 3 nós para viga
clc

fint = sym(zeros(6, 1));

fint = simplify(P * transpose(subs(N, x, a)));

u = sym(zeros(6, 1));

u(3:4) = K(3:4, 3:4) \ fint(3:4)

FR = simplify(K * u - fint);

double(subs(FR, [a,L,P],[1,4,-50]))

double(subs(u, [a,L,P, E, I],[1,4,-50, 205e6, (0.2 * 0.4^3)/12]))

v(x) = simplify(N * u);
dv(x) = simplify(dN * u);
M(x) = simplify(E*I*B*u);

M(x) = subs(M, [a,L,P], [1,4,-50])

double(M(1))

simplify(subs(Mx(0), a, L/2))

vaux(x) = subs(vaux, [a, L, P, E, I], [2,4,-50000,205e9,2/1875])
dvdx(x) = subs(dvdx, [a, L, P, E, I], [2,4,-50000,205e9,2/1875])
Mx(x) = subs(Mx, [a, L, P, E, I], [2,4,-50000,205e9,2/1875])
double(vaux(2))
double(dvdx(2))
double(Mx(0))
fplot(Mx, [0, 4])


subs(FR, [a, L, P], [2, 4, -50])

%% Momentos de extremo fixo de momento concentrado em elementos de 3 nós para viga
clc
syms M

fint = simplify(M * transpose(subs(dN, x, a)))

u = sym(zeros(6, 1));

u(3:4) = K(3:4, 3:4) \ fint(3:4);

FR = simplify(K * u - fint);

subs(FR, [a, L, M], [1, 4, -100])

%% Calcula a função da carga distribuída que começa em a com magnitude w0 e termina em b com magnitude w1
clc
syms h u
A = [x w 1;a w0 1;(L-b) w1 1];
wsol = (collect(solve(det(A) == 0, w), x));

fint = simplify(int(-wsol * transpose(N), x, a, L-b));

u = sym(zeros(6, 1));

u(3:4) = K(3:4, 3:4) \ fint(3:4);

FR = simplify(K * u - fint);

f = simplify(subs(FR, a, a));
Ryt = simplify(collect(f(1) + f(3) + f(5), [w0, w1]));
Mr = simplify(f(2) + f(4) + f(6));

pretty(Ryt)
pretty(Mr)

%pretty(collect(simplify(subs(wsol,x,L-b)), x));
%% Funções de forma para elemento finito de viga em elementos isoparamétricos
clear
clc

syms v1 v2 o1 o2 x L a1 a2 a3 a4 E I xi q P
u = [v1;o1;v2;o2];

A = [x, xi 1;0 -1 1;L 1 1];

Sol = solve(det(A) == 0, x);
Solxi = solve(det(A) == 0, xi);

N1(x) = (2*x^3)/L^3 - (3*x^2)/L^2 + 1;
N2(x) = x - (2*x^2)/L + x^3/L^2;
N3(x) = (3*x^2)/L^2 - (2*x^3)/L^3;
N4(x) = x^3/L^2 - x^2/L;

N1(xi) = simplify(subs(N1, x, Sol));
N2(xi) = simplify(subs(N2, x, Sol));
N3(xi) = simplify(subs(N3, x, Sol));
N4(xi) = simplify(subs(N4, x, Sol));

dN1 = diff(N1, xi, 1);
dN2 = diff(N2, xi, 1);
dN3 = diff(N3, xi, 1);
dN4 = diff(N4, xi, 1);

N = [N1, N2, N3, N4];
dN = [dN1, dN2, dN3, dN4];

B = diff(N, xi, 2);

K = ((8*E*I)/L^3) * int(transpose(B)*B, xi, -1, 1)
f = (L/2)* int(-q*transpose(N), xi, -1, 1)

Kg = sym(zeros(6,6));
Kg(1:4, 1:4) = subs(K, L, L/2);
Kg(3:6, 3:6) = Kg(3:6, 3:6) + subs(K, L, L/2);

fg = sym(zeros(6,1));
fg(1:4) = subs(f, L, L/2);
fg(3:6) = fg(3:6) + subs(f, L, L/2);
fg(3) = fg(3) - P;

ug = sym(zeros(6, 1));
ug(3:4) = Kg(3:4, 3:4)\fg(3:4);
ug

freac = Kg*ug - fg

u = N*ug(1:4)
u = subs(u, [L, P, q, E, I], [1,1,1,200e6,0.004])
figure(1)
fplot(u, [-1, 1])
hold on
u = N*ug(3:6)
u = subs(u, [L, P, q, E, I], [1,1,1,200e6,0.004])
fplot(u, [-1, 1])
hold off

theta1 = (2/L) * dN * ug(1:4);
theta1 = subs(theta1, [L, P, q, E, I], [1,1,1,200e6,0.004]);

figure(2)
fplot(theta1, [-1, 1])
hold on

theta1 = (2/L) * dN * ug(3:6);
theta1 = subs(theta1, [L, P, q, E, I], [1,1,1,200e6,0.004]);
fplot(theta1, [-1, 1])

M1 = ((4*E*I)/(L^2)) * B * ug(1:4)
M1 = subs(M1, [L, P, q, E, I], [1,1,1,200e6,0.004]);
figure(3)
fplot(M1, [-1,1])
hold on

M1 = ((4*E*I)/(L^2)) * B * ug(3:6);
M1 = subs(M1, [L, P, q, E, I], [1,1,1,200e6,0.004]);
fplot(M1, [-1,1])

%% Solução de uma viga com 2 elementos de 3 nós
clear
clc

syms L E I x v1 m1 v2 m2 v3 m3 P a

% Recolhe as funções de forma para o elemento de viga com 3 nós
N1 = ((L + 6*x)*(L^2 - 3*L*x + 2*x^2)^2)/L^5;
N2 = (x*(L^2 - 3*L*x + 2*x^2)^2)/L^4;
N3 = (16*x^2*(L - x)^2)/L^4;
N4 = -(8*x^2*(L - x)^2*(L - 2*x))/L^4;
N5 = (x^2*(L - 2*x)^2*(7*L - 6*x))/L^5;
N6 = -(x^2*(L - x)*(L - 2*x)^2)/L^4;

% Monta a matriz de funções de forma N
N(x, L) = [N1, N2, N3, N4, N5, N6];
dN(x, L) = diff(N,x,1);       % Derivada das funções de forma

% Determina a matriz B, derivada segunda de N
B(x, L) = simplify(diff(N, x, 2));

% Calcula o produto B' * B para aferição
simplify(transpose(B) * B);

% Calcula a matriz de Rigidez K
Kg = sym(zeros(10,10));
K(L) = E*I * simplify(int(transpose(B(x, L)) * B(x, L), x, 0, L));
Kg(1:6, 1:6) = K(L/2);
Kg(5:10, 5:10) = Kg(5:10, 5:10) + K(L/2);

% Vetor de Forças
f = sym(zeros(10,1))
f(1:6) = P * transpose(N(L/4, L/2))

u = sym(zeros(10,1));
u(3:8) = Kg(3:8, 3:8) \ f(3:8)

double(subs(u, [L, P, E, I], [4, -50, 200e6, (0.2*0.4^3)/12]))

FR(a, L) = simplify(Kg * u - f)

% Reescreve v e dvdx em função das funções de forma
v1(x) = N(x, L/2) * u(1:6);
dv1(x) = dN(x, L/2) * u(1:6);
Mx(x, a) = E*I*B(x, L/2)*u(1:6);

%% Exercício do Canal do Clayton Pettit
clear
clc

syms L E A I pt(x) y

N1v(x, L) = (2*x^3)/L^3 - (3*x^2)/L^2 + 1;
N2v(x, L) = x - (2*x^2)/L + x^3/L^2;
N3v(x, L) = (3*x^2)/L^2 - (2*x^3)/L^3;
N4v(x, L) = x^3/L^2 - x^2/L;

Nv = [N1v, N2v, N3v, N4v];
dNv = diff(Nv, x, 1);
Bv = diff(Nv, x, 2);

N1b(x, L) = 1 - x/L;
N2b(x, L) = x/L;

Nb = [N1b, N2b];
dNb = diff(Nb, x, 1);

pt(x) = 80/3 - (20*x)/3;

Kg = sym(zeros(11));
fg = sym(zeros(11, 1));
fe = sym(zeros(4, 1, 2));
u = fg;
Kev = sym(zeros(4, 4, 2));
Keb = sym(zeros(2, 2, 4));
T = sym(zeros(2, 4, 4));

% Elemento 1
Keb(:, :, 1) = E*A*int(transpose(dNb(x, 5))*dNb(x, 5), x, 0, 5);
T(:,:,1) = [1 0 0 0;0 0 1 0];
Kg([1,2,4,5],[1,2,4,5]) = Kg([1,2,4,5],[1,2,4,5]) + transpose(T(:,:,1)) * Keb(:, :, 1) * T(:,:,1);
Kev(:, :, 1) = E*I*int(transpose(Bv(x, 5))*Bv(x, 5), x, 0, 5);
Kg([2,3,5,6],[2,3,5,6]) = Kg([2,3,5,6],[2,3,5,6]) + Kev(:, :, 1);

fe(:, 1, 1) = int(-20*transpose(Nv(x, 5)), x, 1, 5);
fg([2,3,5,6], 1) = fg([2,3,5,6], 1) + fe(:, 1, 1);

% Elemento 2
Keb(:, :, 2) = E*A*int(transpose(dNb(x, 4))*dNb(x, 4), x, 0, 4);
T(:,:,2) = [1 0 0 0;0 0 1 0];
Kg([4,5,7,8],[4,5,7,8]) = Kg([4,5,7,8],[4,5,7,8]) + transpose(T(:,:,2)) * Keb(:, :, 2) * T(:,:,2);
Kev(:, :, 2) = E*I*int(transpose(Bv(x, 4))*Bv(x, 4), x, 0, 4);
Kg([5,6,8,9],[5,6,8,9]) = Kg([5,6,8,9],[5,6,8,9]) + Kev(:, :, 2);

fe(:, 1, 2) = int(-pt(x)*transpose(Nv(x, 4)), x, 1, 4);
fg([5,6,8,9], 1) = fg([5,6,8,9], 1) + fe(:, 1, 2);

% Elemento 3
Keb(:, :, 3) = E*A*int(transpose(dNb(x, 5))*dNb(x, 5), x, 0, 5);
T(:,:,3) = [4/5 -3/5 0 0;0 0 4/5 -3/5];
Kg([10,11,7,8],[10,11,7,8]) = Kg([10,11,7,8],[10,11,7,8]) + transpose(T(:,:,3)) * Keb(:, :, 3) * T(:,:,3);

% Elemento 4
Keb(:, :, 4) = E*A*int(transpose(dNb(x, 4))*dNb(x, 4), x, 0, 4);
T(:,:,4) = [0 1 0 0;0 0 0 1];
Kg([4,5,10,11],[4,5,10,11]) = Kg([4,5,10,11],[4,5,10,11]) + transpose(T(:,:,4))*Keb(:, :, 4)*T(:,:,4);

% Forças Nodais
fg(3) = fg(3) - 20;
fg(10) = fg(10) - 40;

% Deslocamentos totais
gLiberdade = (1:11)';
gRestringidos = [1;2;5];
gSemRest = setdiff(gLiberdade, gRestringidos);

% Cálculo dos deslocamentos nodais
u(gSemRest) = Kg(gSemRest, gSemRest) \ fg(gSemRest)

% Cálculo das reações de apoio
Fr = simplify(Kg(gRestringidos, :) * u - fg(gRestringidos))

% Cálculo das forças em cada elemento da estrutura
Forcas_Barra_Elem_1 = simplify(Keb(:,:,1) * T(:,:,1) * u([1,2,4,5]))
Forcas_Barra_Elem_2 = simplify(Keb(:,:,2) * T(:,:,2) * u([4,5,7,8]))
Forcas_Barra_Elem_3 = simplify(Keb(:,:,3) * T(:,:,3) * u([10,11,7,8]))
Forcas_Barra_Elem_4 = simplify(Keb(:,:,4) * T(:,:,4) * u([4,5,10,11]))

Forcas_Viga_Elem_1 = simplify(Kev(:,:,1) * u([2,3,5,6]) - fe(:,:,1))
Forcas_Viga_Elem_2 = simplify(Kev(:,:,2) * u([5,6,8,9]) - fe(:,:,2))

%%
clear
clc

syms X1 X2 X3 x(X1,X2,X3) u(X1, X2, X3)

x(X1,X2,X3) = [0.8*X1 + 0.2*X2 + 0.01*X3;0.01*X1 + 0.9*X2;1.4*X3]

F = [diff(x,X1, 1), diff(x,X2, 1), diff(x,X3, 1)]

u(X1, X2, X3) = x - [X1; X2; X3]

nablau = [diff(u,X1, 1), diff(u,X2, 1), diff(u,X3, 1)]

epsilon = 0.5 * (nablau + transpose(nablau))

