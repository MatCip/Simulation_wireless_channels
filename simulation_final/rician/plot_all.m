%% this function plots all the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot all 
function plot_all()
    close all;
    load('data/all_data.mat')
    clc; disp('plotting...')
    global delta_t;
    lags = (0:N_step-1)*(fd*delta_t); lags = lags';
    ln_wdt = 2; %line width for plot
    f_size = 20;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rician
    pl1 = figure('Name','Rician - Second Order Statistic','NumberTitle','off');
    set(pl1, 'Position', [200, 100, 800, 600]);
    grid on; hold on; 
    plot(lags,real(R_sim_0{5}),':','LineWidth',ln_wdt);
    plot(lags,real(R_sim_1{5}),':','LineWidth',ln_wdt);
    plot(lags,real(R_sim_3{5}),':','LineWidth',ln_wdt);
    plot(lags,real(R_id_0{5}),'k','LineWidth',ln_wdt/2);
    plot(lags,real(R_id_1{5}),'k','LineWidth',ln_wdt/2);
    plot(lags,real(R_id_3{5}),'k','LineWidth',ln_wdt/2);
    xlabel({'Normalized Time: $f_d \tau$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'$Re[R_{Z,Z}(\tau)]$'},'Interpreter','latex','FontSize',f_size); 
    xlim([0,10]); ylim([-1,1.2]);
    legend({'$K=0$','$K=1$','$K=3$','$Theory$'},'Interpreter','latex','FontSize',f_size);
    
    pl2 = figure('Name','Rician - Second Order Statistic','NumberTitle','off');
    set(pl2, 'Position', [200, 100, 800, 600]);
    grid on; hold on; 
    plot(lags,R_sim_0{6},':','LineWidth',ln_wdt);
    plot(lags,R_sim_1{6},':','LineWidth',ln_wdt);
    plot(lags,R_sim_3{6},':','LineWidth',ln_wdt);
    plot(lags,R_id_0{6}/max(R_id_0{6}),'k','LineWidth',ln_wdt/2);
    plot(lags,R_id_1{6}/max(R_id_1{6}),'k','LineWidth',ln_wdt/2);
    plot(lags,R_id_3{6}/max(R_id_3{6}),'k','LineWidth',ln_wdt/2);
    xlabel({'Normalized Time: $f_d \tau$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'Normalized $R_{|Z|^2,|Z|^2}(\tau)$'},'Interpreter','latex','FontSize',f_size); 
    xlim([0,10]); ylim([0.4,1.2]);
    legend({'$K=0$','$K=1$','$K=3$','$Theory$'},'Interpreter','latex','FontSize',f_size);
    
    pl3 = figure('Name','Rician - Second Order Statistic','NumberTitle','off');
    set(pl3, 'Position', [200, 100, 800, 600]);
    grid on; hold on; 
    plot(lags,imag(R_sim_0{5}),':','LineWidth',ln_wdt);
    plot(lags,imag(R_sim_1{5}),':','LineWidth',ln_wdt);
    plot(lags,imag(R_sim_3{5}),':','LineWidth',ln_wdt);
    plot(lags,imag(R_id_0{5}),'k','LineWidth',ln_wdt/2);
    plot(lags,imag(R_id_1{5}),'k','LineWidth',ln_wdt/2);
    plot(lags,imag(R_id_3{5}),'k','LineWidth',ln_wdt/2);
    xlabel({'Normalized Time: $f_d \tau$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'$Im[R_{Z,Z}(\tau)]$'},'Interpreter','latex','FontSize',f_size); 
    xlim([0,10]); ylim([-1,1.2]);
    legend({'$K=0$','$K=1$','$K=3$','$Theory$'},'Interpreter','latex','FontSize',f_size);
    
    pl4 = figure('Name','Rician - Fading Envelope','NumberTitle','off');
    set(pl4, 'Position', [200, 100, 800, 600]);
    grid on; hold on;
    x = linspace(0,4,500);
    plot(x,pdf_sim_0{1},':','LineWidth',ln_wdt);
    plot(x,pdf_sim_1{1},':','LineWidth',ln_wdt);
    plot(x,pdf_sim_3{1},':','LineWidth',ln_wdt);
    plot(x,pdf_sim_5{1},':','LineWidth',ln_wdt);
    plot(x,pdf_sim_10{1},':','LineWidth',ln_wdt);
    plot(x,pdf_id_0{1},'k','LineWidth',ln_wdt/2);
    plot(x,pdf_id_1{1},'k','LineWidth',ln_wdt/2);
    plot(x,pdf_id_3{1},'k','LineWidth',ln_wdt/2);
    plot(x,pdf_id_5{1},'k','LineWidth',ln_wdt/2);
    plot(x,pdf_id_10{1},'k','LineWidth',ln_wdt/2);
    xlabel({'$z$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'$f_{|Z|}(z)$'},'Interpreter','latex','FontSize',f_size); 
    xlim([0,4]); ylim([-0.2,2]);
    legend({'$K=0 (N=8)$','$K=1 (N=8)$','$K=3 (N=8)$','$K=5 (N=8)$','$K=10 (N=8)$','$Theory (N=\infty)$'},'Interpreter','latex','FontSize',f_size);
    
    pl5 = figure('Name','Rician - Fading Envelope','NumberTitle','off');
    set(pl5, 'Position', [200, 100, 800, 600]);
    grid on; hold on;
    x1 = linspace(-1,1,50);
    x = linspace(-1,1,500);
    plot(x1,smooth(pdf_sim_0{2},0.9),':','LineWidth',ln_wdt);
    plot(x1,smooth(pdf_sim_1{2},0.5),':','LineWidth',ln_wdt);
    plot(x1,smooth(pdf_sim_3{2},0.5),':','LineWidth',ln_wdt);
    plot(x1,smooth(pdf_sim_5{2},0.5),':','LineWidth',ln_wdt);
    plot(x1,smooth(pdf_sim_10{2},0.5),':','LineWidth',ln_wdt);
    plot(x,pdf_id_0{2},'k','LineWidth',ln_wdt/2);
    plot(x,pdf_id_1{2},'k','LineWidth',ln_wdt/2);
    plot(x,pdf_id_3{2},'k','LineWidth',ln_wdt/2);
    plot(x,pdf_id_5{2},'k','LineWidth',ln_wdt/2);
    plot(x,pdf_id_10{2},'k','LineWidth',ln_wdt/2);
    xlabel({'$\psi (\times \pi)$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'$f_{\Psi}(\psi)$'},'Interpreter','latex','FontSize',f_size); 
    xlim([-1,1]); ylim([0.14,0.18]);
    legend({'$K=0 (N=8)$','$K=1 (N=8)$','$K=3 (N=8)$','$K=5 (N=8)$','$K=10 (N=8)$','$Theory (N=\infty)$'},'Interpreter','latex','FontSize',f_size);
    
    pl6 = figure('Name','Rician - Level Crossing Rate','NumberTitle','off');
    rho_db = 20*log10(rho);
    set(pl6, 'Position', [200, 100, 800, 600]);
    semilogy(rho_db, LCR_sim_0,':','LineWidth',ln_wdt); grid on; hold on;
    semilogy(rho_db, LCR_sim_1,':','LineWidth',ln_wdt); grid on; hold on;
    semilogy(rho_db, LCR_sim_3,':','LineWidth',ln_wdt); grid on; hold on;
    semilogy(rho_db, LCR_sim_5,':','LineWidth',ln_wdt); grid on; hold on;
    semilogy(rho_db, LCR_sim_10,':','LineWidth',ln_wdt); grid on; hold on;
    semilogy(rho_db, LCR_id_0,'k','LineWidth',ln_wdt/2);
    semilogy(rho_db, LCR_id_1,'k','LineWidth',ln_wdt/2);
    semilogy(rho_db, LCR_id_3,'k','LineWidth',ln_wdt/2);
    semilogy(rho_db, LCR_id_5,'k','LineWidth',ln_wdt/2);
    semilogy(rho_db, LCR_id_10,'k','LineWidth',ln_wdt/2);
    xlabel({'Normalized fading envelope level $\rho$ [dB]'},'Interpreter','latex','FontSize',f_size);
    ylabel({'Normalized LCR'},'Interpreter','latex','FontSize',f_size);
    xlim([-25,10]); ylim([1e-3,10^0.5]);
    legend({'$K=0 (N=8)$','$K=1 (N=8)$','$K=3 (N=8)$','$K=5 (N=8)$','$K=10 (N=8)$','$Theory (N=\infty)$'},'Interpreter','latex','FontSize',f_size,'Location','northwest');
   
    pl7 = figure('Name','Rician - Average Fading Duration','NumberTitle','off');
    rho_db = 20*log10(rho);
    set(pl7, 'Position', [200, 100, 800, 600]);
    semilogy(rho_db, AFD_sim_0,   ':', 'LineWidth', ln_wdt); grid on; hold on;
    semilogy(rho_db, AFD_sim_1,   ':', 'LineWidth', ln_wdt);
    semilogy(rho_db, AFD_sim_3,   ':', 'LineWidth', ln_wdt);
    semilogy(rho_db, AFD_sim_5,   ':', 'LineWidth', ln_wdt);
    semilogy(rho_db, AFD_sim_10,  ':', 'LineWidth', ln_wdt);
    semilogy(rho_db, AFD_id_0,    'k', 'LineWidth', ln_wdt/2); 
    semilogy(rho_db, AFD_id_1,    'k', 'LineWidth', ln_wdt/2);
    semilogy(rho_db, AFD_id_3,    'k', 'LineWidth', ln_wdt/2);
    semilogy(rho_db, AFD_id_5,    'k', 'LineWidth', ln_wdt/2);
    semilogy(rho_db, AFD_id_10,   'k', 'LineWidth', ln_wdt/2);
    xlabel({'Normalized fading envelope level $\rho$ [dB]'},'Interpreter','latex','FontSize',f_size);
    ylabel({'Normalized AFD'},'Interpreter','latex','FontSize',f_size);
    xlim([-20,5]); ylim([10^-1.3,10^1]);
    legend({'$K=0 (N=8)$','$K=1 (N=8)$','$K=3 (N=8)$','$K=5 (N=8)$','$K=10 (N=8)$','$Theory (N=\infty)$'},'Interpreter','latex','FontSize',f_size,'Location','northwest');
    
%% save plots
    plots = cell(1,7);
    plots{1} = pl1;
    plots{2} = pl2;
    plots{3} = pl3;
    plots{4} = pl4;
    plots{5} = pl5;
    plots{6} = pl6;
    plots{7} = pl7;

    %save_plots(plots); %save plots
    clc;disp('Simulation is finished!')
end
%% save plots
function save_plots(pl)
    clc;disp('Saving Plots...')
    saveas(pl{1}, 'graphs/rician_re_rgrg','epsc');
    saveas(pl{2}, 'graphs/rician_r2r2','epsc');
    saveas(pl{3}, 'graphs/rician_im_rgrg','epsc');
    saveas(pl{4}, 'graphs/rician_pdf_abs','epsc');
    saveas(pl{5}, 'graphs/rician_pdf_ph','epsc');
    saveas(pl{6}, 'graphs/rician_lcr','epsc');
    saveas(pl{7}, 'graphs/rician_afd','epsc');
end