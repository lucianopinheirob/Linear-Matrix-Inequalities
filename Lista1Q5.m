function [relacao_entrada_saida, norma_hinf_LMI, norma_h2, norma_h2_LMI] = Luciano_cinco()

clear all
clc

%Sistema
A = [0 1; -16 -1.2];
B = [1; 1];
C = [1 2];
D = 0;
sis = ss(A, B, C, D);

%5 b)
bode(sis);

T = 0.01;
t = 0:T:10;

%5 c)
%LMI para computo da norma Hinf
P = sdpvar(2);
mi = sdpvar(1);
F = [A'*P+P*A+C'*C P*B+C'*D; B'*P+D'*C D'*D-mi*eye(1)];
des1 = P;
des2 = F;
prob = [des1 >= 0, des2 <= 0];
obj = mi;
S = sdpsettings('solver', 'lmilab');
optimize(prob, obj, S);
mi = double(mi);
norma_hinf_LMI = sqrt(mi);


%Relação energia saída/entrada
w = 3.8:0.05:4.2;
res = [];
for i = 1:length(w)
    [r] = rel(w(i), sis, t, T);
    res = [res r];
end

relacao_entrada_saida = res

%5 d)
%Cálculo da norma H2
h = impulse(sis, t)
norma_h2 = sqrt(trapz(h.^2)*T)

%5 e)
%LMI para computo da norma H2
P = sdpvar(2);
des1 = P;
des2 = A'*P+P*A+C'*C;
prob = [des1 >= 0, des2 <= 0];
obj = trace(B'*P*B);
S = sdpsettings('solver', 'lmilab');
optimize(prob, obj, S);
P = double(P);
norma_h2_LMI = sqrt(B'*P*B);




end

%5 a)
function [r] = rel(w, sis, t, T)
%Sinal de Entrada
entrada = sin(w*t).*exp(-0.1*t);
energia_entrada = sqrt(trapz(entrada.^2)*T);

%Sinal de saída
saida = lsim(sis,entrada,t);
energia_saida = sqrt(trapz(saida.^2)*T);

%figure(1)
%lsim(sis,entrada,t)

%Relação eneria saída/entrada
r = energia_saida/energia_entrada;
end



