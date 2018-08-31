%% this function implement the channel based on the work of Jakes,  Pop-Beaulieu and Xiao-Zheng.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: 
% M        : N = 4*M+2 where N is the number of oscillators in the Jakes
%            simulator,
% fd       : the doppler frequency,
% str_mode : string for selecting which channel has to be computed,
% N_stat   : number of statistical trials.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT: 
% g        : fading,
% pdf      : a cell that contains the pdf of the absolute value of g and
%            the pdf of the phase of g,
% R        : a cell that contains the autocorrelations and the
%            crosscorrelations of the fading g.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate channels
function [g,pdf,R,rho,LCR_id,LCR,AFD_id,AFD] = gen_channel(M,fd,str_mode,N_stat)
    global N_step;
    global delta_t;
    % R{1} = Rgcgc, R{2} = Rgsgs, R{3} = Rgcgs, R{4} = Rgsgc, R{5} = Rgg,
    % R{6} = Rg2g2.
    g = zeros(N_step,1);
    pdf = cell(1,2); 
    pdf{1} = zeros(1,500); 
    pdf{2} = zeros(1,30);
    R = cell(1,6);
    for i=1:6
        R{i} = zeros(N_step,1);
    end
    N_lev = 1000;
    normalize = true;
    rho = zeros(N_lev,1);
    LCR = zeros(N_lev,1);
    AFD = zeros(N_lev,1);
    LCR_id = zeros(N_lev,1);
    AFD_id = zeros(N_lev,1);

    switch str_mode
        case 'ideal'
            clc; disp('Generating Ideal Channel...')
            R = cell(1,6);
            g=0;
            lags = (0:N_step-1)*(fd*delta_t); lags = lags';
            R{1} = besselj(0,2*pi*lags);
            R{2} = besselj(0,2*pi*lags);
            R{3} = zeros(N_step,1);
            R{4} = zeros(N_step,1);
            R{5} = 2*besselj(0,2*pi*lags);
            R{6} = 4 + 4*besselj(0,2*pi*lags).^2;
            
            r = linspace(0,4,500);
            theta = ones(1,500);
            pdf_abs = r.*exp(-r.^2/2);
            pdf_ph = 1/(2*pi)*theta;
            pdf = cell(1,2);
            pdf{1} = pdf_abs;
            pdf{2} = pdf_ph;
            
        case 'jakes'
            clc; disp('Generating Jakes Channel...')
            [g,pdf,R] = gen_rayleigh_ch_1(M,fd); 
            
        case 'pop_beaulieu'
            for i=1:N_stat
                clc; disp(strcat('Generating Pop-Beaulieu Channel... (N_stat_max=',num2str(N_stat),')',' M=',num2str(M)))
                disp(strcat('N_stat=',num2str(i)))
                [g_temp,pdf_temp,R_temp] = gen_rayleigh_ch_2(M,fd); %g= fading
                g = g + g_temp/N_stat;
                for j=1:2
                    pdf{j} = pdf{j} + pdf_temp{j}/N_stat;
                end
                for j=1:6
                    R{j} = R{j} + R_temp{j}/N_stat;
                end
            end
            
        case 'zheng_xiao'
            for i=1:N_stat
                clc; disp(strcat('Generating Zheng-Xiao Channel... (N_stat_max=',num2str(N_stat),')'))
                disp(strcat('N_stat=',num2str(i)))
                %generate the temporary channel
                [g_temp,pdf_temp,R_temp] = gen_rayleigh_ch_3(M,fd);
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
            [rho,LCR_id,AFD_id] = get_LCR_AFD(g,fd,'id_rayleigh',0,normalize,N_lev);
        
        otherwise disp('Error in string mode.');
    end
end

%% ch_1 Jakes
function [g,pdf,R] = gen_rayleigh_ch_1(M,fd)
    global time;
    N=4*M+2;
    %------------------------
    beta_0 = pi/4;
    a_0 = sqrt(2)*cos(beta_0);
    b_0 = sqrt(2)*sin(beta_0);
    w_0 = 2*pi*fd;
    phi = 0;
    %------------------------
    gr = a_0*cos(w_0*time + phi);
    gi = b_0*cos(w_0*time + phi);
    %------------------------
    for n = 1:M
        %------------------------
        beta_n = n*pi/M;
        a_n = 2*cos(beta_n);
        b_n = 2*sin(beta_n);
        w_n = w_0*cos((2*pi*n)/N);
        %------------------------
        gr = gr + a_n*cos(w_n*time + phi);
        gi = gi + b_n*cos(w_n*time + phi);
    end
    g = 2/sqrt(N)*(gr + 1i*gi);
    F = abs(g); %fading power
    [pdf_abs,~] = ksdensity(F,linspace(0,4,100));% pdf of fading power
    TH = angle(g);
    [pdf_ph,~] = ksdensity(TH,linspace(-pi,pi,100));
    pdf = cell(1,2);
    pdf{1} = pdf_abs; pdf{2} = pdf_ph;
    R = get_corr(g,gr,gi);
    
end

%% ch_2 Pop-Beaulieau
function [g,pdf,R] = gen_rayleigh_ch_2(M,fd)
global time;
N=4*M+2;
phi= 2*pi*rand(4*M+2,1);

X_c=zeros(length(time),1);
 for n = 1:M
     A_n=2*pi*n/N;
     w_n=2*pi*fd*cos(A_n);
       Term1=(cos(phi(n))+cos(phi(2*M+1-n))+cos(phi(2*M+1+n))+cos(phi(4*M+2-n)))*cos(w_n*time);
       Term2=(sin(phi(n))-sin(phi(2*M+1-n))-sin(phi(2*M+1+n))+sin(phi(4*M+2-n)))*sin(w_n*time);
       X_c=X_c+Term1-Term2;
 end
 Term3=(cos(phi(2*M+1))+cos(phi(4*M+2)))*cos(2*pi*fd*time);
 Term4=(sin(phi(2*M+1))-sin(phi(4*M+2)))*sin(2*pi*fd*time);
 X_c=X_c+Term3+Term4;
 X_c=X_c*sqrt(2/(N));
 
 X_s=zeros(length(time),1);
 for n = 1:M
     A_n=2*pi*n/N;
     w_n=2*pi*fd*cos(A_n);
       Term1=(cos(phi(n))-cos(phi(2*M+1-n))-cos(phi(2*M+1+n))+cos(phi(4*M+2-n)))*sin(w_n*time);
       Term2=(sin(phi(n))+sin(phi(2*M+1-n))+sin(phi(2*M+1+n))+sin(phi(4*M+2-n)))*cos(w_n*time);
       X_s=X_s+Term1+Term2;
 end
 Term3=(cos(phi(2*M+1))-cos(phi(4*M+2)))*sin(2*pi*fd*time);
 Term4=(sin(phi(2*M+1))+sin(phi(4*M+2)))*cos(2*pi*fd*time);
 X_s=X_s-Term3+Term4;
 X_s=X_s*sqrt(2/(N));
 
 g=X_c+1i*X_s;
 gc=X_c;
 gs=X_s;
 
  F = abs(g); %fading power
    [pdf_abs,~] = ksdensity(F,linspace(0,4,500));% pdf of fading power
    TH = angle(g);
    [pdf_ph,~] = ksdensity(TH,linspace(-pi,pi,30),'support',[-pi,pi]);
    pdf=cell(1,2); pdf{1} = pdf_abs; pdf{2} = pdf_ph;
    R = get_corr(g,gc,gs);
 
end

%% ch_3 Zheng-Xiao
function [g,pdf,R] = gen_rayleigh_ch_3(M,fd)
    global time;
    gc = zeros(length(time),1);
    gs = zeros(length(time),1);
    wd = 2*pi*fd;
    theta = 2*pi*rand - pi; 
    phi = 2*pi*rand - pi;

    for n = 1:M
        %------------------------
        psi_n = 2*pi*rand - pi;
        alpha_n = (2*pi*n - pi + theta)/(4*M);
        wn = wd*cos(alpha_n);
        %------------------------
        gc = gc + cos(psi_n)*cos(wn*time + phi);
        gs = gs + sin(psi_n)*cos(wn*time + phi);
    end;
    g = 2/sqrt(M)*(gc + 1i*gs);
    F = abs(g); %fading power
    [pdf_abs,~] = ksdensity(F,linspace(0,4,500));% pdf of fading power
    TH = angle(g);
    pdf_ph = get_pdf_ph(TH,30);
    pdf=cell(1,2); pdf{1} = pdf_abs; pdf{2} = pdf_ph;
    R = get_corr(g,gc,gs);
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
    [counts, edges] = histcounts(TH,edges);
    pdf_ph = counts/(sum(counts)*d_theta);
end