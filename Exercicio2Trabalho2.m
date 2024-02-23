clc; clear;

%Definir Parâmetros
M = rand() * 10^-4;
R = rand() * 9 + 1;
R = round(R, 3);
T = [250, 300, 350];
n_vals = [200 400 600];

for i = 1:length(T)

    fprintf('Valor de T = %d K:\n', T(i));
    fprintf('Valor de R = %d K:\n', R);
    fprintf('Valor de M = %d K:\n', M);

    for j = 1:length(n_vals)
        
    f = @(v) ((4.*sqrt(pi)).*(M./(2*R*T(i))).^(3./2)) .* (v.^3 .* exp(-M*v.^2 ./ (2*R*T(i)))); %define a função

    %define a malha
    a = 0; %Limite inferior
    b = 1500;
    n = n_vals(j);
    h(j) = (b-a)/n;

    x = linspace(a,b,n+1); % M+1 pontos
    pm = (x(1:n)+x(2:n+1))/2; % pontos médios
    
    %Ponto Médio
    I_PM = h(j)*sum(f(pm));
    %Trapézio
    I_T = h(j)/2*(f(x(1))+2*sum(f(x(2:n)))+f(x(n+1)));
    %Simpson
    I_S = h(j)/6*(f(x(1))+4*sum(f(pm))+2*sum(f(x(2:n)))+f(x(n+1)));
    %Valor real do integral 
    V_R = integral(f, a, b);

    E_PM(j) = abs(V_R - I_PM);
    E_T(j) = abs(V_R - I_T);
    E_S(j) = abs(V_R - I_S);

    %fprintf dos resultados
    fprintf('n = %d:\n', n);
    fprintf('Regra do Ponto Médio: %.6f m/s\n', I_PM);
    fprintf('Regra do Trapézio: %.6f m/s\n', I_T);
    fprintf('Regra de Simpson: %.6f m/s\n', I_S);
    fprintf('Valor real do integral: %.6f m/s\n', V_R);
    fprintf('Erro do PM: %.6f\n', E_PM(j));
    fprintf('Erro do T: %.6f\n', E_T(j));
    fprintf('Erro de S: %.6f\n', E_S(j))
    fprintf('\n');
    

    end

    %Calcular a ordem de convergência para cada método
    for k = 1:length(n_vals)-1
        
        p_PM = (log(E_PM(k+1)/E_PM(k)))/(log(h(k+1)/h(k)));
        p_T = (log(E_T(k+1)/E_T(k)))/(log(h(k+1)/h(k)));
        p_S = (log(E_S(k+1)/E_S(k)))/(log(h(k+1)/h(k)));

    end

    fprintf('Ordem de convergência da regra do Ponto Médio: %.2f\n', p_PM);
    fprintf('Ordem de convergência da regra do Trapézio: %.2f\n', p_T);
    fprintf('Ordem de convergência da regra de Simpson: %.2f\n', p_S);

end

