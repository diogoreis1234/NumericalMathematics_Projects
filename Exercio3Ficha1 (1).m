clc; close all; clear;

n       = 4;                    % Ordem da matriz
A       = randi([0 10],n);      % Matriz arbitrária de ordem n com elemntos a variar entre 0 e 10
I       = eye(n);               % Matriz identidade de ordem n
x       = rand(n, 1);           % Vector arbitrário
x       = x / norm(x);          % Normalizaçao do vector arbitrário anterior, sendo norm(x) a norma do vector
y       = x'*A*x;               % Valor próprio de A, aproximacao inicial
tol     = 1e-10;                % Tolerância
kmax    = 100;                  % #máximo de iteracões  
k       = 0;                    % Iteração inicial
erro    = tol + 1;              % Iniciar o ciclo
aux     = 1;                    % Variável auxiliar para o cálculo do vetor próprio 

fprintf(' \n  Matriz A: \n');
    disp(A);

disp('---------------------------------');
disp(' Iteração         Valor prório   ');
disp('---------------------------------');
fprintf('    %2.0f   %20.16f \n',k,y);

%Ciclo de Newton
while (erro > tol) && (k < kmax)
    F       = Funcao(A,x,y);
    J       = Jacobiana(A,x,y,I);
    pk      = - J \ F;
    v       = [x ; y];              % Vector normalizado anteriormente com y na ultima linha         
    b       = v + pk;               % Soma do vector v com o vector pk, ficando b com n+1 elems          
    x       = b(1:n);               % x toma os valores do n primeiro elementos de b
    y       = b(n+1);               % y toma o último valor de b       
    erro    = norm(F);                      
    k       = k + 1;                        
    fprintf('    %2.0f   %20.16f \n',k,y);  % print do valor próprio na iteração k
 
end

    fprintf(' \n  Vector próprio: \n');

        % Impressão do vector próprio
while (aux<n+1)                 
    fprintf('   %20.16f \n',x(aux));
    aux=aux+1;
end
     
function F = Funcao(A, x, y) % Funcao F conforme o enunciado
F = [A*x - y*x; x'*x - 1];
end

function J = Jacobiana(A, x, y, I) % Jacobiana conforme o enunciado
J = [A - y*I, -x; 2*x', 0];
end