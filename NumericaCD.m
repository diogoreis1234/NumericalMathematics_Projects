%Dados e partição
N =1000;
h = 1/(N+1);
x = (0:h:1)';
f = ones(N, 1);
e = ones(N,1);

num_eps = 3; % Número de valores de epsilon
uh_values = cell(num_eps, 1); % Célula para armazenar as soluções
u_exata_values = cell(num_eps, 1); % Célula para armazenar as soluções exatas
for k = 1:num_eps
    eps = 10^(-k);
    
    % Criar a matriz de rigidez esparsa
    diagonal = (2*eps/(h^2)) - (1/h);
    diagonalinf = -eps/(h^2);
    diagonalsup = -eps/(h^2) + (1/h);

    A = spdiags([diagonalinf*e, diagonal*e, diagonalsup*e], -1:1, N, N);

    % Vetor de carga b
    b = f;

    % Resolução do sistema de equações lineares A*uh = b
    uh = A\b;
    uh = [0; uh; 0];
    uh_values{k} = uh; % Armazenar a solução para o valor atual de epsilon

    % Cálculo da solução exata
    u_exata = @(x) x - (exp((x-1)/eps) - exp(-1/eps)) / (1 - exp(-1/eps));
    u_exata_vals = u_exata(x); % Avaliar a função nos pontos da rede
    u_exata_values{k} = u_exata_vals; % Armazenar a solução exata para o valor atual de epsilon
    
    % Cálculo do erro
    erro = norm(uh-u_exata_vals, inf);
    erro1 = h*norm(uh-u_exata_vals, 1);
    erro2 = (h^(1/2))*norm(uh-u_exata_vals, 2);

    fprintf('Erro para epsilon = 1e-%d:\n', k);
    fprintf('Erro infinito: %.6f\n', erro);
    fprintf('Erro 1: %.6f\n', erro1);
    fprintf('Erro 2: %.6f\n', erro2);
    fprintf('\n');

end

% Plot dos gráficos individuais
colors = {'b', 'b', 'b'};
for k = 1:num_eps
    figure;
plot(x, uh_values{k}, colors{k}, 'LineWidth', 2);
hold on
plot(x, u_exata_values{k}, 'r--', 'LineWidth', 2);
hold off
legend('Diferenças Finitas', 'Solução Exata');
xlabel('x');
ylabel('u(x)');
title(sprintf('Solução por Diferenças Finitas vs. Solução Exata (Epsilon = 1e-%d)', k));
end
