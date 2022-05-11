function [resp] = Luciano_dois(A, N, n)

[maximo, ponto_max] = Luciano_um(A, 1);
maior = maximo;

res{1} = [n, N, maior];
[res{2}] = Teorema1(A,N,n);
[res{3}] = Lema1(A,N,n);
[res{4}] = Lema2(A,N,n);
[res{5}] = Lema3(A,N,n);
[res{6}] = Lema4(A,N,n);

resp = [];
for i = 1:6
    resp = [resp res{i}];
end
        


end

function [res] = Teorema1(A,N,n)

P = sdpvar(n);
obj = trace(P);
%obj = [];
des1 = P;
prob = [des1>=0];
for i = 1:N
    des2{i} = A{i}'*P + P*A{i};
    prob = [prob, des2{i}<=0];
end
S = sdpsettings('solver', 'lmilab');
diagnostics = optimize(prob,obj,S);
E = ~diagnostics.problem;
V = length(depends(prob));
L = n*(N+1);
P = double(P);

res = [E V L];

end

function [res] = Lema1(A,N,n)

S1 = sdpvar(n,n,'full');
S2 = sdpvar(n,n,'full');
%obj = [];
prob = [];
for i = 1:N
    P{i} = sdpvar(n);
    F{i} = [S1*A{i}+A{i}'*S1' P{i}-S1+A{i}'*S2'; P{i}+S2*A{i}-S1' -S2-S2']; 
    des1{i} = F{i};
    des2{i} = P{i};
    prob = [prob, des1{i}<=0, des2{i}>=0];
end
obj = trace(P{1});
S = sdpsettings('solver', 'lmilab');
diagnostics = optimize(prob,obj,S);
E = ~diagnostics.problem;
V = length(depends(prob));
L = 3*n*N;

res = [E V L];

end

function [res] = Lema2(A,N,n)

V1 = sdpvar(n,n,'full');
%obj = [];
prob = [];
for i = 1:N
    P{i} = sdpvar(n);
    F{i} = [-V1-V1' V1'*A{i}+P{i} V1'; A{i}'*V1+P{i} -P{i} zeros(n); V1 zeros(n) -P{i}]; 
    des1{i} = P{i};
    des2{i} = F{i};
    prob = [prob, des1{i}>=0, des2{i}<=0];
end
obj = trace(P{1});
S = sdpsettings('solver', 'lmilab');
diagnostics = optimize(prob,obj,S);
E = ~diagnostics.problem;
V = length(depends(prob));
L = 4*n*N;

res = [E V L];

end

function [res] = Lema3(A,N,n)

%obj = [];
prob1 = [];
prob2 = [];
for i = 1:N
    P{i} = sdpvar(n);
end
for i = 1:N
    des1{i} = P{i};
    des2{i} = A{i}'*P{i} + P{i}*A{i}+eye(n);
    prob1 = [prob1, des1{i}>=0, des2{i}<=0];
    for j = i+1:N
        des3{i} = A{i}'*P{j}+P{j}*A{i}+A{j}'*P{i}+P{i}*A{j}+(2*eye(n))/(1-N);
        prob2 = [prob2, des3{i}<=0];
    end
end
obj = trace(P{1});
prob = [prob1, prob2]
S = sdpsettings('solver', 'lmilab');
diagnostics = optimize(prob,obj,S);
E = ~diagnostics.problem;
V = length(depends(prob));
L = 2*n*N + nchoosek(N,2)*n;

res = [E V L];

end

function [res] = Lema4(A,N,n)

%obj = [];
prob1 = [];
prob2 = [];
for i = 1:N
    P{i} = sdpvar(n);
    S1{i} = sdpvar(n,n,'full');
    S2{i} = sdpvar(n,n,'full');
end
obj = trace(P{2});
for i = 1:N
    F{i} = [S1{i}*A{i}+A{i}'*S1{i}' P{i}-S1{i}+A{i}'*S2{i}'; P{i}+S2{i}*A{i}-S1{i}' -S2{i}-S2{i}'];
    des1{i} = P{i};
    des2{i} = F{i};
    prob1 = [prob1, des1{i}>=0, des2{i}<=0];
    for j = i+1:N
        G{i} = [A{i}'*S1{j}'+S1{j}*A{i}+A{j}'*S1{i}'+S1{i}*A{j} ...
            P{i}+P{j}-S1{i}-S1{j}+A{i}'*S2{j}'+A{j}'*S2{i}';...
            P{i}+P{j}-S1{i}'-S1{j}'+S2{j}*A{i}+S2{i}*A{j} ...
            -S2{i}-S2{i}'-S2{j}-S2{j}'];
        des3{i} = G{i};
        prob2 = [prob2, des3{i}<=0];
    end
end
prob = [prob1, prob2]
S = sdpsettings('solver', 'lmilab');
diagnostics = optimize(prob,obj,S);
E = ~diagnostics.problem;
V = length(depends(prob));
L = 3*n*N + nchoosek(N,2)*2*n;

res = [E V L];

end















