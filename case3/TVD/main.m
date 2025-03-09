clc;
clear;
% �?��计算的限制器
Limiters = {'Superbee', 'Minmod', 'Woodward', 'Van Leer', 'UMIST', 'WACEB', 'Albada', 'OSPRE', 'TCDF', 'Koren'};
% �?��计算的格�?
Methods = {'MUSCL', 'MTVDLF'};
for i = 1:length(Limiters)
    for j = 1:length(Methods)
        limiter = Limiters{i};
        method = Methods{j};
        disp("limiter: " +  limiter  + "         method:  " + method);
        case1_TVD_f(limiter, method);
    end
end

