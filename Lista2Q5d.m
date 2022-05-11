clc
clear all

A = [1 -2 6; 7 2 -5; -9 7 4];
B = [2 0; -3 -1; 2 1];

eps = [0.001 0.01 0.1 1 10 100];

for i =1:length(eps)
    
    z11 = sdpvar(1);
    z22 = sdpvar(1);
    z13 = sdpvar(1);

    x11 = sdpvar(1);
    x13 = sdpvar(1);
    x22 = sdpvar(1);
    x31 = sdpvar(1);
    x33 = sdpvar(1);

    Z = [z11 0 z13; 0 z22 0];
    X = [x11 0 x13; 0 x22 0; x31 0 x33];
    
    W = sdpvar(3,3,'symmetric');

    des1 = [A*X+X'*A'+B*Z+Z'*B' W-X'+eps(i)*A*X+eps(i)*B*Z; W-X+eps(i)*X'*A'+eps(i)*Z'*B' -eps(i)*X-eps(i)*X'];
    des2 = W;
    
    prob = [des1<=0 des2>=0];
    obj = [];
    S = sdpsettings('solver', 'lmilab');
    diagnostics = optimize(prob,obj,S);

    Z = double(Z);
    W = double(W);
    X = double(X);

    K{i} = Z*inv(X)
end

disp ('Ks (do menor para o maior epsilon):')
for i = 1:length(eps)
    K{i}
end


