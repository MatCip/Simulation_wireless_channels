%% this function compute the channel based on the work of Komninakis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT: 
% g       : fading,
% pdf_abs : pdf of the absolute value of g,
% pdf_ph  : pdf of the phase of g,
% R       : a cell that contains the autocorrelations and crosscorrelations,
% H_abs   : amplitude of the computed filter,
% Y       : ideal amplitude of the filter,
% zer     : zeros of the computed filter,
% pol     : poles of the computed filter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate channel with an IIR filter
function [g,pdf_abs,pdf_ph,R ,H_abs,Y,zer,pol] = gen_channel_IIR(fd)
    global delta_t;
    %filter parameters
    K = 7; % filter order is 2K 
    M = 500;
    rho_0 = 0.2;
    N= floor(rho_0/(fd*delta_t));
    rho=N*(fd*delta_t);
    Y = get_Y(M,rho);
    W = linspace(0,1,M+1);
    z = exp((1i*pi).*W);
    
    %initialization of the paramaters for the ellip. alg.
    x = zeros(4*K,1);
    B = 100*eye(4*K); %usata  per aggiornare x ad ogni iterazione
    eps = 0.01; %accuracy for the convergence
    F = get_F(M,K,z,x);
    A0 = ((abs(F'))*(Y))/((abs(F'))*(abs(F)));
    R_gd = get_R_gd(M,K,F,z,x,A0,Y);
    beta = sqrt(R_gd'*B*R_gd);
    
    iter = 0;
    while(beta>eps)
        g=R_gd/beta;
        x=x-(B*g)/(4*K+1);
        B=(4*K)^2/((4*K)^2-1)*(B-2/(4*K+1)*B*g*(g'*B));  
        F = get_F(M,K,z,x);
        A0=((abs(F'))*(Y))/((abs(F'))*(abs(F)));
        R_gd = get_R_gd(M,K,F,z,x,A0,Y);
        beta=abs(sqrt(R_gd'*B*R_gd));
        
        if mod(iter,50) == 0
            clc;disp(strcat('Filter Design Convergence... (threshold=',num2str(eps),')'))
            disp(strcat('beta=',num2str(beta)))
            disp(strcat('iteration=',num2str(iter)))
        end
        iter = iter+1;
    end
    H_abs = A0*F;
    g = A0*get_fading(x,K,delta_t,rho,fd); %fading
    [pdf_abs,~ ] = ksdensity(abs(g));
    [pdf_ph,~ ] = ksdensity(angle(g));
    [zer,pol] = get_zer_pol(x);
    R = get_corr(g,real(g),imag(g));
end

%% compute the gradient of R
function R_gd = get_R_gd(M,K,F,z,x,A0,Y)
    G=ones(4*K,M+1);
    R_gd=ones(4*K,1);
    for m=1:M+1
        for k=1:K-1
            a_k = x(4*k+1);
            b_k = x(4*k+2);
            c_k = x(4*k+3);
            d_k = x(4*k+4);
            G(4*k+1,m)= abs(F(m))*real(z(m)^(-1)/(1+a_k*z(m)^(-1)+b_k*z(m)^(-2)));
            G(4*k+2,m)= abs(F(m))*real(z(m)^(-2)/(1+a_k*z(m)^(-1)+b_k*z(m)^(-2)));
            G(4*k+3,m)=-abs(F(m))*real(z(m)^(-1)/(1+c_k*z(m)^(-1)+d_k*z(m)^(-2)));
            G(4*k+4,m)=-abs(F(m))*real(z(m)^(-2)/(1+c_k*z(m)^(-1)+d_k*z(m)^(-2)));
        end
    end

    for k=1:4*K
        for m=1:M+1
            R_gd(k)=R_gd(k)+(A0*abs(F(m))-Y(m))*G(k,m);
        end
        R_gd(k)=R_gd(k)*2*A0;
    end
end

%% compute the function F
function F = get_F(M,K,z,x)
    F=ones(M+1,1);
    for m=1:(M+1) %for each zi
        for k=1:K-1
            a_k = x(4*k+1);
            b_k = x(4*k+2);
            c_k = x(4*k+3);
            d_k = x(4*k+4);
            F(m)=F(m)*((1+a_k*z(m)^-1+b_k*z(m)^-2)/(1+c_k*z(m)^-1+d_k*z(m)^-2));
        end
    end
end

%% compute the target filter amplitude Y
function Y = get_Y(M,rho)
    
    L = floor(2*rho*M);
    Y = zeros(M+1,1);
    for i=0:M
        if i<=L-1
            Y(i+1) = sqrt(1/(sqrt(1-(i/L)^2)));
        end
        if i==L
            Y(i+1) = sqrt(L*(pi/2-asin((L-1)/L)));
        end
        if i>=L+1 
            Y(i+1) = 0;
        end
    end 
end
%% compute the fading
function g = get_fading(x,K,T,rho,fd)
    global N_step;
    
    I=rho/(fd*T);
    % generate gaussian process
    input_dim=ceil(N_step/I);
    w_real=randn(input_dim,1);
    w_im=randn(input_dim,1);
    w=w_real+1i*w_im;

%     h = ifft(H_abs);
%     g = conv(h,w);
    g=w;
    for i=1:K
        a=[1,x(3+4*(i-1)),x(4+4*(i-1))];
        b=[1,x(1+4*(i-1)),x(2+4*(i-1))];
        a = polystab(a);
        b = polystab(b);

        g = filter(b,a,g);
    end
    g=interp(g,I);
    g = g(1:N_step);
end

%% get the zeros and the poles of the designed filter
function [zer,pol] = get_zer_pol(x)
    K = length(x)/4;
    zer = zeros(2*K,1);
    pol = zeros(2*K,1);
    for i=1:K
        a=[1,x(3+4*(i-1)),x(4+4*(i-1))];
        b=[1,x(1+4*(i-1)),x(2+4*(i-1))];
        a = polystab(a);
        b = polystab(b);
        zer(2*(i-1)+1:2*(i-1)+2) = roots(b);
        pol(2*(i-1)+1:2*(i-1)+2) = roots(a);
    end
end

%% get cross correlation
function R = get_corr(g,gr,gi)
    R = cell(1,6);
    Rgcgc = xcorr(gr,gr,'coeff'); %Rgcgc = Rgcgc(N_step:end);
    Rgsgs = xcorr(gi,gi,'coeff'); %Rgsgs = Rgsgs(N_step:end);
    Rgcgs = xcorr(gr,gi,'coeff'); %Rgcgs = Rgcgs(N_step:end);
    Rgsgc = xcorr(gi,gr,'coeff'); %Rgsgc = Rgsgc(N_step:end);
    Rgg = xcorr(g,g,'coeff'); %Rgg = Rgg(N_step:end);
    Rg2g2 = 8*xcorr(abs(g).^2,abs(g).^2,'coeff'); %Rg2g2 = Rg2g2(N_step:end);
    R{1} = Rgcgc;
    R{2} = Rgsgs;
    R{3} = Rgcgs;
    R{4} = Rgsgc;
    R{5} = Rgg;
    R{6} = Rg2g2;
end