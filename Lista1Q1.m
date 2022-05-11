function [maximo, ponto_max] = Luciano_um(A,flag)

N = length(A); %Número de vértices

%Pontos uniformemente distribuídos
for i = 1:100    
    alfa{i}(1)=1-rand^(1/(N-1));
    for k = 2:N-1
        alfa{i}(k)=(1-sum(alfa{i}(1:k-1)))*(1-rand^(1/(N-k)));
    end
    alfa{i}(N)=1-sum(alfa{i}(1:N-1));
    A_de_alfa{i} = 0;
    for k = 1:N
        A_de_alfa{i} = A_de_alfa{i} + alfa{i}(k)*A{k};
    end
    autoa{i} = eig(A_de_alfa{i});
end

%Pontos igualmente espaçados para cada segmento de reta     
L = 0;
for j = 1:N
    for k = j:N-1
        for i = 1:100
            beta{L*100+i}(j) = 1-0.01*i;
            beta{L*100+i}(k+1) = 1-beta{L*100+i}(j);
            A_de_beta{L*100+i} = beta{L*100+i}(j)*A{j} + beta{L*100+i}(k+1)*A{k+1};
            autob{L*100+i} = eig(A_de_beta{i});
            if k ~= N-1
                for count = k+2:N
                    beta{L*100+i}(count) = 0;
                end
            end
        end
        L = L+1;
    end
end

%Pontos uniformemente distribuídos no subpolitopo
T = 0;
if (N>3)
    for j = 1:N
        for k = j:N-1
            for l = k:N-2
                for i = 1:150
                    gama{T*150+i}(j)=1-rand^(1/2);
                    gama{T*150+i}(k+1)=(1-gama{T*150+i}(j))*(1-rand^(1));
                    gama{T*150+i}(l+2)=1-sum(gama{T*150+i}(1:k+1));
                    A_de_gama{T*150+i} = gama{T*150+i}(j)*A{j} + gama{T*150+i}(k+1)*A{k+1} + gama{T*150+i}(l+2)*A{l+2};
                    autog{T*150+i} = eig(A_de_gama{T*150+i});
                    if l ~= N-2
                        for count = l+3:N
                            gama{T*150+i}(count) = 0;
                        end
                    end
                end
                T = T+1;
            end
        end
    end
end

%Concatena autovalores
if (N>3)
    auto = [autoa autob autog];
    ponto = [alfa beta gama];
else
    auto = [autoa autob];
    ponto = [alfa beta];
end

%Plota autovalores
%figure(1)
%hold on
for i = 1:length(auto)
    for j = 1:length(A{1})
        %plot(real(auto{i}(j)),imag(auto{i}(j)),'bx')
    end
end

%plota círculo unitário
if flag == 1
    for i = 1:length(auto)
        maximo_auto(i) = max(real(auto{i}));
    end
    [maximo, ind_max] = max(maximo_auto(i));
    ponto_max = ponto{ind_max};
elseif flag == 0
    for i = 1:length(auto)
        maximo_auto(i) = max(abs(auto{i}));
    end
    [maximo, ind_max] = max(maximo_auto(i));
    ponto_max = ponto{ind_max};
%     figure(1)
%     hold on
%     ang=0:0.01:2*pi; 
%     xp=1*cos(ang);
%     yp=1*sin(ang);
%     plot(xp,yp,'r');
end

end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        