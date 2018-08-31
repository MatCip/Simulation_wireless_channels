%% this function implement the channel based on the work of Jakes,  Pop-Beaulieu and Xiao-Zheng.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: 
% str_mode : string for selecting which channel has to be computed.
% N_stat   : number of statistical trials.
% K        : ratio between power of LoS and power of non-LoS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT: 
% g        : fading,
% pdf      : a cell that contains the pdf of the absolute value of g and
%            the pdf of the phase of g.
% R        : a cell that contains the autocorrelations and the
%            crosscorrelations of the fading g.
% rho      : normalized levels for lcr and afd.
% LCR      : level crossing rate.
% AFD      : average fading duration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate channels
function [g,pdf,R,rho,LCR,AFD] = gen_channel(str_mode,N_stat,K)
    global N_step;
    global delta_t;
    global theta_0;
    global N;
    global fd;
    % R{1} = Rgcgc, R{2} = Rgsgs, R{3} = Rgcgs, R{4} = Rgsgc, R{5} = Rgg,
    % R{6} = Rg2g2.
    g = zeros(N_step,1);
    pdf = cell(1,2); 
    pdf{1} = zeros(1,500); 
    pdf{2} = zeros(1,50);
    R = cell(1,6);
    for i=1:6
        R{i} = zeros(N_step,1);
    end
    N_lev = 1000;
    normalize = true;
    rho = zeros(N_lev,1);
    LCR = zeros(N_lev,1);
    AFD = zeros(N_lev,1);

    switch str_mode
        case 'sim_rician'
            for i=1:N_stat
                clc; disp(strcat('Generating Simulated Rician Channel... (N_stat_max=',num2str(N_stat),')'))
                disp(strcat('N_stat=',num2str(i)))
                disp(strcat('K=',num2str(K)))
                %generate the temporary channel
                [g_temp,pdf_temp,R_temp] = gen_si_rician_ch(N,fd,K);
                [rho, LCR_temp, AFD_temp] = get_LCR_AFD(g_temp,fd,'simulated',0,normalize,N_lev);
                %update
                g = g + g_temp/N_stat;
                for j=1:2
                    pdf{j} = pdf{j} + pdf_temp{j}/N_stat;
                end
                for j=1:6
                    R{j} = R{j} + R_temp{j}/N_stat;
                end
                LCR = LCR + LCR_temp/N_stat;
                AFD = AFD + AFD_temp/N_stat;
            end
            
        case 'id_rician'
            clc; disp('Generating Ideal Rician Channel...')
            disp(strcat('K=',num2str(K)))
            r = linspace(0,4,500);
            pdf = cell(1,2);
            pdf{1} = 2*(1+K)*r.*exp(-K-(1+K)*r.^2).*besseli(0,2*r*sqrt(K*(K+1)));
            pdf{2} = 1/(2*pi)*ones(1,size(r,2));
            g=0;
            R = cell(1,6);
            lags = (0:N_step-1)*(fd*delta_t); lags = lags';
            wd_tau = 2*pi*lags;
            bes = besselj(0,wd_tau);
            R{1} = (bes + K*cos(wd_tau*cos(theta_0)))/(2+2*K);
            R{2} = R{1};
            R{3} = K*sin(wd_tau*cos(theta_0))/(2+2*K);
            R{4} = -R{3};
            R{5} = (bes + K*cos(wd_tau*cos(theta_0)) + 1i*K*sin(wd_tau*cos(theta_0)))/(1+K);
            load('data/fc_fs.mat');
            R{6} = (1 + bes.^2 + K^2 - f_c - f_s + 2*K*(1 + bes.*cos(wd_tau*cos(theta_0))))/(1+K)^2;
            [lev, LCR, AFD] = get_LCR_AFD(g,fd,'id_rician',K,normalize,N_lev);
            rho = lev;
        otherwise disp('Error in string mode.');
    end
end

%% Simulated Rician
function [z,pdf,R] = gen_si_rician_ch(N,fd,K)
    global time; global theta_0;
    %------------------------

    w_0 = 2*pi*fd;
    y_c = zeros(size(time,1),1);
    y_s = zeros(size(time,1),1);
    for n = 1:N
        theta_n = rand(1,1)*2*pi-pi;
        phi_n = rand(1,1)*2*pi-pi;
        alpha_n = (2*pi*n + theta_n)/N;
        y_c = y_c + 1/sqrt(N)*cos(w_0*time*cos(alpha_n) + phi_n);
        y_s = y_s + 1/sqrt(N)*sin(w_0*time*cos(alpha_n) + phi_n);
    end
    phi_0 = rand(1,1)*2*pi-pi;
    z_c = (y_c + sqrt(K)*cos(w_0*time*cos(theta_0)+phi_0))/sqrt(1+K);
    z_s = (y_s + sqrt(K)*sin(w_0*time*cos(theta_0)+phi_0))/sqrt(1+K);
    z = z_c + 1i*z_s;
    
    F = abs(z); %fading power
    [pdf_abs,~] = ksdensity(F,linspace(0,4,500));% pdf of fading power
    TH = angle(z);
    pdf_ph = get_pdf_ph(TH,50);
    pdf = cell(1,2);
    pdf{1} = pdf_abs; pdf{2} = pdf_ph;
    R = get_corr(z,z_c,z_s);
    R{5} = R{5}/2;
    R{6} = R{6}/8;
    
end

%% get cross correlation
function R = get_corr(g,gr,gi)
    global N_step;
    R = cell(6,1);
    Rgcgc = xcorr(gr,gr,'coeff'); Rgcgc = Rgcgc(N_step:end);
    Rgsgs = xcorr(gi,gi,'coeff'); Rgsgs= Rgsgs(N_step:end);
    Rgcgs = xcorr(gr,gi,'coeff'); Rgcgs = Rgcgs(N_step:end);
    Rgsgc = xcorr(gi,gr,'coeff'); Rgsgc = Rgsgc(N_step:end);
    Rgg = 2*xcorr(g,g,'coeff'); Rgg = Rgg(N_step:end);
    Rg2g2 = 8*xcorr(abs(g).^2,abs(g).^2,'coeff'); Rg2g2 = Rg2g2(N_step:end);
    R{1} = Rgcgc;
    R{2} = Rgsgs;
    R{3} = Rgcgs;
    R{4} = Rgsgc;
    R{5} = Rgg;
    R{6} = Rg2g2;
end

%% get pdf phase
function pdf_ph = get_pdf_ph(TH,n_points)
    edges = linspace(-pi,pi,n_points+1);
    d_theta = edges(2)-edges(1);
    [counts,edges] = histcounts(TH,edges);
    pdf_ph = counts/(sum(counts)*d_theta);
end