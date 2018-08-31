clear;
load('data/var_01.mat');
load('data/var_05.mat');
load('data/var_002.mat');
%% smooth
LCR_4_01_sm = smooth(LCR_4_01,0.1);
LCR_4_05_sm = smooth(LCR_4_05,0.1);
LCR_4_002_sm = smooth(LCR_4_002,0.05);
l_db_4_05 = 20*log10(lambda_4_05);
l_db_4_01 = 20*log10(lambda_4_01);
l_db_4_002 = 20*log10(lambda_4_002);
sp = 1000;
l = linspace(-40,8,length(LCR_4_01(1:sp:end)));
%% ind
ind_05 = zeros(length(l),1);
for i=1:length(l)
    zz = find(abs(l(i)-l_db_4_05) < 0.1);
    ind_05(i) = zz(1);
end

ind_01 = zeros(length(l),1);
for i=1:length(l)
    zz = find(abs(l(i)-l_db_4_01) < 0.1);
    ind_01(i) = zz(1);
end

ind_002 = zeros(length(l),1);
for i=1:length(l)
    zz = find(abs(l(i)-l_db_4_002) < 0.1);
    ind_002(i) = zz(1);
end
%%
ln_wdt = 2; %line width for plot
pl = figure('Name','Level Crossing Rate','NumberTitle','off');
set(pl, 'Position', [200, 100, 800, 600]);
semilogy(l,LCR_4_05_sm(ind_05),'^b','LineWidth',ln_wdt/2,'MarkerSize',10);grid on; hold on;
semilogy(l,LCR_4_01_sm(ind_01),'sb','LineWidth',ln_wdt/2,'MarkerSize',10); 
semilogy(l,LCR_4_002_sm(ind_002),'ob','LineWidth',ln_wdt/2,'MarkerSize',10);
semilogy(20*log10(lambda_4_05),LCR_id_05,'r','LineWidth',ln_wdt);
semilogy(20*log10(lambda_4_01),LCR_id_01,'r','LineWidth',ln_wdt);
semilogy(20*log10(lambda_4_002),LCR_id_002,'r','LineWidth',ln_wdt);
xlabel({'$\lambda$ [dB]'},'Interpreter','latex','FontSize',20);
ylabel({'Normalized LCR'},'Interpreter','latex','FontSize',20);
xlim([-40,10]);ylim([1e-5,1e-1]);
legend({strcat('Sim. $f_d \tau$ = ',num2str(fd_05*delta_t_05)),strcat('Sim. $f_d \tau$ = ',num2str(fd_01*delta_t_01)),strcat('Sim. $f_d \tau$ = ',num2str(fd_002*delta_t_002)),'Theory'},'Interpreter','latex','FontSize',20,'Location','northwest');
%% save
saveas(pl,'graphs/komninakis_var_rho','epsc');
