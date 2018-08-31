%% this function plots all the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot all 
function plot_all()
    close all;
    load('data/all_data.mat')
    clc; disp('plotting...')
    global delta_t;
    ln_wdt = 2; %line width for plot
    f_size = 20;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Pop - Beaulieu
    pl1 = figure('Name','Pop/Beaulieu - Fading Envelope','NumberTitle','off');
    set(pl1, 'Position', [200, 100, 800, 600]);
    grid on; hold on;
    x = linspace(0,4,500);
    plot(x,pdf_2_id{1},'k','LineWidth',ln_wdt);
    plot(x,pdf_2_3{1},':','LineWidth',ln_wdt);
    plot(x,pdf_2_4{1},'-.','LineWidth',ln_wdt);
    plot(x,pdf_2_8{1},'--','LineWidth',ln_wdt);
    xlabel({'$r$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'$f_{|g|}(r)$'},'Interpreter','latex','FontSize',f_size); 
    xlim([0,4]);
    legend({'Reference','$M=3$','$M=4$','$M=8$'},'Interpreter','latex','FontSize',f_size);
    
    %% Xiao - Zheng
    pl2 = figure('Name','Xiao/Zheng - Second Order Statistic','NumberTitle','off');
    set(pl2, 'Position', [200, 100, 800, 600]);
    lags_3 = (0:N_step-1)*(fd_3*delta_t); lags_3 = lags_3';
    grid on; hold on;
    plot(lags_3,R_3_id{1},'k','LineWidth',ln_wdt);
    plot(lags_3,R_3_10{1},':','LineWidth',ln_wdt);
    plot(lags_3,R_3_50{1},'-.','LineWidth',ln_wdt);
    plot(lags_3,R_3_100{1},'--','LineWidth',ln_wdt);
    xlabel({'Normalized Time: $f_d \tau$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'$R_{X_c,X_c}(\tau)$'},'Interpreter','latex','FontSize',f_size);
    xlim([0,t_sim]);
    legend({'Reference','$N_{stat}=10$','$N_{stat}=50$','$N_{stat}=100$'},'Interpreter','latex','FontSize',f_size);

    pl3 = figure('Name','Xiao/Zheng - Second Order Statistic','NumberTitle','off');
    set(pl3, 'Position', [200, 100, 800, 600]);
    lags_3 = (0:N_step-1)*(fd_3*delta_t); lags_3 = lags_3';
    grid on; hold on; 
    plot(lags_3,R_3_id{3},'k','LineWidth',ln_wdt);
    plot(lags_3,R_3_10{3},':','LineWidth',ln_wdt);
    plot(lags_3,R_3_50{3},'-.','LineWidth',ln_wdt);
    plot(lags_3,R_3_100{3},'--','LineWidth',ln_wdt);
    xlabel({'Normalized Time: $f_d \tau$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'$R_{X_c,X_s}(\tau)$'},'Interpreter','latex','FontSize',f_size);
    ylim([-0.5,0.5]); xlim([0,t_sim]);
    legend({'Reference','$N_{stat}=10$','$N_{stat}=50$','$N_{stat}=100$'},'Interpreter','latex','FontSize',f_size);

    pl4 = figure('Name','Xiao/Zheng - Second Order Statistic','NumberTitle','off');
    set(pl4, 'Position', [200, 100, 800, 600]);
    lags_3 = (0:N_step-1)*(fd_3*delta_t); lags_3 = lags_3';
    grid on; hold on; 
    plot(lags_3,R_3_id{2},'k','LineWidth',ln_wdt);
    plot(lags_3,R_3_10{2},':','LineWidth',ln_wdt);
    plot(lags_3,R_3_50{2},'-.','LineWidth',ln_wdt);
    plot(lags_3,R_3_100{2},'--','LineWidth',ln_wdt);
    xlabel({'Normalized Time: $f_d \tau$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'$R_{X_s,X_s}(\tau)$'},'Interpreter','latex','FontSize',f_size); 
    xlim([0,t_sim]);
    legend({'Reference','$N_{stat}=10$','$N_{stat}=50$','$N_{stat}=100$'},'Interpreter','latex','FontSize',f_size);

    pl5 = figure('Name','Xiao/Zheng - Second Order Statistic','NumberTitle','off');
    set(pl5, 'Position', [200, 100, 800, 600]);
    lags_3 = (0:N_step-1)*(fd_3*delta_t); lags_3 = lags_3';
    grid on; hold on; 
    plot(lags_3,real(R_3_id{5}),'k','LineWidth',ln_wdt);
    plot(lags_3,real(R_3_10{5}),':','LineWidth',ln_wdt);
    plot(lags_3,real(R_3_50{5}),'-.','LineWidth',ln_wdt);
    plot(lags_3,real(R_3_100{5}),'--','LineWidth',ln_wdt);
    xlabel({'Normalized Time: $f_d \tau$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'Re[$R_{X,X}(\tau)$]'},'Interpreter','latex','FontSize',f_size); 
    xlim([0,t_sim]);ylim([-1,2]);
    legend({'Reference','$N_{stat}=10$','$N_{stat}=50$','$N_{stat}=100$'},'Interpreter','latex','FontSize',f_size);

    pl6 = figure('Name','Xiao/Zheng - Second Order Statistic','NumberTitle','off');
    set(pl6, 'Position', [200, 100, 800, 600]);
    lags_3 = (0:N_step-1)*(fd_3*delta_t); lags_3 = lags_3';
    grid on; hold on; 
    plot(lags_3,imag(R_3_id{5}),'k','LineWidth',ln_wdt);
    plot(lags_3,imag(R_3_10{5}),':','LineWidth',ln_wdt);
    plot(lags_3,imag(R_3_50{5}),'-.','LineWidth',ln_wdt);
    plot(lags_3,imag(R_3_100{5}),'--','LineWidth',ln_wdt);
    xlabel({'Normalized Time: $f_d \tau$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'Im[$R_{X,X}(\tau)$]'},'Interpreter','latex','FontSize',f_size); 
    xlim([0,t_sim]);ylim([-0.5,0.5]);
    legend({'Reference','$N_{stat}=10$','$N_{stat}=50$','$N_{stat}=100$'},'Interpreter','latex','FontSize',f_size);

    pl7 = figure('Name','Xiao/Zheng - Second Order Statistic','NumberTitle','off');
    set(pl7, 'Position', [200, 100, 800, 600]);
    lags_3 = (0:N_step-1)*(fd_3*delta_t); lags_3 = lags_3';
    grid on; hold on; 
    plot(lags_3,R_3_id{6},'k','LineWidth',ln_wdt);
    plot(lags_3,R_3_10{6},':','LineWidth',ln_wdt);
    plot(lags_3,R_3_50{6},'-.','LineWidth',ln_wdt);
    plot(lags_3,R_3_100{6},'--','LineWidth',ln_wdt);
    xlabel({'Normalized time: $f_d \tau$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'$R_{|X|^2,|X|^2}(\tau)$'},'Interpreter','latex','FontSize',f_size); 
    xlim([0,t_sim]);
    legend({'Reference','$N_{stat}=10$','$N_{stat}=50$','$N_{stat}=100$'},'Interpreter','latex','FontSize',f_size);

    pl8 = figure('Name','Xiao/Zheng - Fading Envelope','NumberTitle','off');
    set(pl8, 'Position', [200, 100, 800, 600]);
    grid on; hold on; 
    x = linspace(0,4,500);
    plot(x,pdf_3_id{1},'k','LineWidth',ln_wdt);
    plot(x,pdf_3_10{1},':','LineWidth',ln_wdt);
    plot(x,pdf_3_50{1},'-.','LineWidth',ln_wdt);
    plot(x,pdf_3_100{1},'--','LineWidth',ln_wdt);
    xlabel({'$x$'},'Interpreter','latex','FontSize',f_size);
    ylabel({'$f_{|X|}(x)$'},'Interpreter','latex','FontSize',f_size); 
    xlim([0,4]);
    legend({'Reference','$N_{stat}=10$','$N_{stat}=50$','$N_{stat}=100$'},'Interpreter','latex','FontSize',f_size);

    pl9 = figure('Name','Xiao/Zheng - Fading Envelope','NumberTitle','off');
    set(pl9, 'Position', [200, 100, 800, 600]);
    grid on; hold on; xlabel({'$\theta_X$ ($\times \pi$)'},'Interpreter','latex','FontSize',f_size);
    x = linspace(-1,1,30);
    x1 = linspace(-1,1,500);
    plot(x1,pdf_3_id{2},'k','LineWidth',ln_wdt);
    plot(x,smooth(pdf_3_10{2}',1-1e-12),':','LineWidth',ln_wdt);
    plot(x,smooth(pdf_3_50{2}',0.9999999),'-.','LineWidth',ln_wdt);
    plot(x,smooth(pdf_3_100{2}',0.8),'--','LineWidth',ln_wdt);
    ylabel({'$f_{\theta_X}(\theta_X)$'},'Interpreter','latex','FontSize',f_size); 
    xlim([-1,1]); ylim([0.13,0.19]);
    legend({'Reference','$N_{stat}=10$','$N_{stat}=50$','$N_{stat}=100$'},'Interpreter','latex','FontSize',f_size);
    
    pl10 = figure('Name','Xiao/Zheng - Level Crossing Rate','NumberTitle','off');
    set(pl10, 'Position', [200, 100, 800, 600]);
    semilogy(20*log10(lambda_3_10),LCR_3_id,'k','LineWidth',ln_wdt);  grid on; hold on;
    semilogy(20*log10(lambda_3_10),LCR_3_10,':','LineWidth',ln_wdt);
    semilogy(20*log10(lambda_3_50),LCR_3_50,'-.','LineWidth',ln_wdt); 
    semilogy(20*log10(lambda_3_100),LCR_3_100,'--','LineWidth',ln_wdt); 
    xlabel({'Normalized fading envelope level $\rho$ [dB]'},'Interpreter','latex','FontSize',f_size);
    ylabel({'Normalized LCR'},'Interpreter','latex','FontSize',f_size);
    xlim([-30,10]); ylim([10^-3,10^1]);
    legend({'Reference','$N_{stat}=10$','$N_{stat}=50$','$N_{stat}=100$'},'Interpreter','latex','FontSize',f_size,'Location','northwest');
    
    pl11 = figure('Name','Xiao/Zheng - Average Fading Duration','NumberTitle','off');
    set(pl11, 'Position', [200, 100, 800, 600]);
    semilogy(20*log10(lambda_3_10),AFD_3_id,'k','LineWidth',ln_wdt); grid on; hold on;
    semilogy(20*log10(lambda_3_10),AFD_3_10,':','LineWidth',ln_wdt);
    semilogy(20*log10(lambda_3_50),AFD_3_50,'-.','LineWidth',ln_wdt);
    semilogy(20*log10(lambda_3_100),AFD_3_100,'--','LineWidth',ln_wdt);
    xlabel({'Normalized fading envelope level $\rho$ [dB]'},'Interpreter','latex','FontSize',f_size);
    ylabel({'Normalized AFD'},'Interpreter','latex','FontSize',f_size);
    xlim([-30,5]); ylim([10^-2,10^1.5]);
    legend({'Reference','$N_{stat}=10$','$N_{stat}=50$','$N_{stat}=100$'},'Interpreter','latex','FontSize',f_size,'Location','northwest');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Komninakis paper
    pl12 = figure('Name','Komninakis - Magnitude Response of the Filter','NumberTitle','off');
    set(pl12, 'Position', [200, 100, 800, 600]);
    M = 500; W = linspace(0,1,M+1);
    hold on; grid on;
    plot(W,20*log10(Y),'k','LineWidth',ln_wdt);
    plot(W,20*log10(abs(H_abs)),'r','LineWidth',ln_wdt/2);
    ylim([-80,20]); xlim([0,1]);
    xlabel({'Normalized Frequency $f_d\tau$ (1 = Fc/2)'},'Interpreter','latex','FontSize',f_size);
    ylabel({'$|$H$|$ [dB]'},'Interpreter','latex','FontSize',f_size);
    legend({'Ideal Deideres Response','Designed Response'},'Interpreter','latex','FontSize',f_size);
    
    pl13 = figure('Name','Komninakis - Second Order Statistic','NumberTitle','off');
    set(pl13, 'Position', [200, 100, 800, 600]);
    grid on; hold on; 
    R_id_shift = besselj(0,1/pi*linspace(-100,100,N_step));
    plot(-N_step+1:N_step-1,real(R_4{5}),'--','LineWidth',ln_wdt);
    plot(-N_step+1:N_step-1,imag(R_4{5}),'--','LineWidth',ln_wdt);
    plot(linspace(-100,100,N_step),real(R_id_shift),'k','LineWidth',ln_wdt/2);
    plot(linspace(-100,100,N_step),imag(R_id_shift),'k','LineWidth',ln_wdt/2);
    xlabel({'Correlation lag in samples'},'Interpreter','latex','FontSize',f_size);
    xlim([-100,100]);
    legend({'Simulated Re[$R_{g,g}$]','Simulated Im[$R_{g,g}$]','Theory'},'Interpreter','latex','FontSize',f_size);
    
%%
    pl = cell(13,1);
    pl{1 } = pl1; 
    pl{2 } = pl2; 
    pl{3 } = pl3; 
    pl{4 } = pl4; 
    pl{5 } = pl5;
    pl{6 } = pl6; 
    pl{7 } = pl7; 
    pl{8 } = pl8; 
    pl{9 } = pl9; 
    pl{10} = pl10;
    pl{11} = pl11; 
    pl{12} = pl12; 
    pl{13} = pl13; 

    %save_plots(pl); %save plots
    clc;disp('Simulation is finished!')
end
%% save plots
function save_plots(pl)
    clc;disp('Saving Plots...')
    saveas(pl{1 }, 'graphs/pop_beaulieu_fading_abs','epsc');
    saveas(pl{2 }, 'graphs/xiao_zheng_rgcgc','epsc');
    saveas(pl{3 }, 'graphs/xiao_zheng_rgcgs','epsc');
    saveas(pl{4 }, 'graphs/xiao_zheng_rgsgs','epsc');
    saveas(pl{5 }, 'graphs/xiao_zheng_re_rgg','epsc');
    saveas(pl{6 }, 'graphs/xiao_zheng_im_rgg','epsc');
    saveas(pl{7 }, 'graphs/xiao_zheng_rg2g2','epsc');
    saveas(pl{8 }, 'graphs/xiao_zheng_fading_abs','epsc');
    saveas(pl{9 }, 'graphs/xiao_zheng_fading_ang','epsc');
    saveas(pl{10}, 'graphs/xiao_zheng_lcr','epsc');
    saveas(pl{11}, 'graphs/xiao_zheng_afd','epsc');
    saveas(pl{12}, 'graphs/komninakis_filter','epsc');
    saveas(pl{13}, 'graphs/komninakis_rgg','epsc');
end