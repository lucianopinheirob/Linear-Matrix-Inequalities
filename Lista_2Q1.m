clc
clear all

n = 3;
N = 2;

P = rolmipvar(n, n, 'P', 'symmetric', N, 1);
F = rolmipvar(n, n, 'F', 'full', N, 0);
G = rolmipvar(n, n, 'G', 'full', N, 1);

A_{1} = [0 1 0; 0 0 1; -1 -3 -2];
A_{2} = [0 1 0; 0 0 1; -1 -22 -2];

A = rolmipvar(A_, 'A', 2, 1);

%testar vários valores de gama
gamma = 2.076;
while (1)
  dotP = diff(P, 'dot_P', [-gamma gamma; -gamma gamma]);
  des1 = [dotP+F*A+A'*F' P-F+A'*G'; P+G*A-F' -G-G'];
  des2 = P;
  obj = [];
  prob = [des1<=0; des2>=0];
  S = sdpsettings('solver', 'sedumi');
  solvesdp(prob,obj,S);
  k = min(checkset(prob))
  if (k>0)
    gamma = gamma+0.0001
  else
    break
  end
end
