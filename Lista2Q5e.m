clc
clear all

A = [1 -2 6; 7 2 -5; -9 7 4];
B = [2 0; -3 -1; 2 1];
C = eye(3);

eps = [0.001 0.01 0.1 1 10 100];

for i = 1:length(eps)
    
    Z = sdpvar(2,3);
    W = sdpvar(3,3,'symmetric');
    X = sdpvar(3,3,'full');
    
    des1 = [A*X+X'*A'+B*Z+Z'*B' W-X'+eps(i)*A*X+eps(i)*B*Z; W-X+eps(i)*X'*A'+eps(i)*Z'*B' -eps(i)*X-eps(i)*X'];
    des2 = W;
    prob = [des1<=0 des2>=0];
    
    obj = [];
    
    S = sdpsettings('solver', 'lmilab');
    diagnostics = optimize(prob,obj,S);

    Z = double(Z);
    W = double(W);
    X = double(X);

    K{i} = Z*inv(X);
  
end

for i  = 1:length(eps)
    
    P = sdpvar(3,3,'symmetric');
    F = sdpvar(3,3,'full');
    G = sdpvar(3,3,'full');
    
    h11 = sdpvar(1);
    h22 = sdpvar(1);
    j11 = sdpvar(1);
    j22 = sdpvar(1);
    j13 = sdpvar(1);
    H = [h11 0; 0 h22];
    J = [j11 0 j13; 0 j22 0];  
    
    M1 = [A'*F'+F*A+K{i}'*B'*F'+F*B*K{i} P-F+A'*G'+K{i}'*B'*G' F*B+C'*J'-K{i}'*H'];
    M2 = [(P-F+A'*G'+K{i}'*B'*G')' -G-G' G*B];
    M3 = [(F*B+C'*J'-K{i}'*H')' (G*B)' -H-H'];   
    des1 = [M1; M2; M3];
    
    des2 = P;   
    
    obj = [];
    prob = [des1<=0; des2>=0];
    S = sdpsettings('solver', 'lmilab');
    diagnostics = optimize(prob,obj,S);
    
    H = double(H);
    J = double(J);
    
    L{i} = inv(H)*J
end

disp ('Ls (do menor para o maior epsilon):')
for i = 1:length(eps)
    L{i}
end



