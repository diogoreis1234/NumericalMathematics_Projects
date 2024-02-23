clc; clear;

for k = 1:3
    eps = 10^(-k);
    
    u = @(x) x - (exp((x-1)/eps) - exp(-1/eps)) / (1 - exp(-1/eps));
  
    fplot(u, [0, 1], 'LineWidth', 1.5)
    hold on
end

ylim([0, 1])
xlim([0, 1])
legend('\epsilon = 10^{-1}', '\epsilon = 10^{-2}', '\epsilon = 10^{-3}');
title('Soluções numéricas para diferentes \epsilon')
xlabel('x')
ylabel('u(x)')

