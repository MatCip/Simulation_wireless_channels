clear;
clc;
%%
t_sim = 15;
delta_t = 1e-3;
N_step = t_sim/delta_t;
fd = 10;
lags = (0:N_step-1)*(fd*delta_t); lags = lags';
N = 8;
%%
f_c = zeros(N_step,1);
f_s = zeros(N_step,1);
I_k_cos = zeros(N_step,1);
I_k_sin = zeros(N_step,1);
for i=1:N_step
    if mod(i,100)==0
        disp(i)
    end
    for k=1:N
        fun_cos = @(x) cos(2*pi*lags(i)*cos(x));
        fun_sin = @(x) sin(2*pi*lags(i)*cos(x));
        I_k_cos(i) = integral(fun_cos,(2*pi*k-pi)/N,(2*pi*k+pi)/N);
        I_k_sin(i) = integral(fun_sin,(2*pi*k-pi)/N,(2*pi*k+pi)/N);
        f_c(i) = f_c(i) + (1/(2*pi)*I_k_cos(i))^2;
        f_s(i) = f_s(i) + (1/(2*pi)*I_k_sin(i))^2;
    end
end

%%
save('data/fc_fs.mat','f_c','f_s');