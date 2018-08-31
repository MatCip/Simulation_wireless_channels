%% this function compute the LCR and the AFD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: 
% g      : fading,
% fd     : doppler frequency.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT: 
% lambda : lambda = R/Rrms where R is the level,
% LCR_id : ideal level crossing rate,
% LCR    : computed level crossing rate,
% AFD_id : ideal average fade duration,
% AFD    : computed average fade duration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [rho, LCR, AFD] = get_LCR_AFD(g,fd,str_mode,K,normalize,N_dim)
    global N_step; global delta_t; global theta_0;
    F = abs(g); %fading
    rho_db = linspace(-30,10,N_dim);
    rho = 10.^(rho_db/20);
    R = rho*rms(F);
    
    %% ideal LCR and AFD
    LCR = zeros(N_dim,1);
    AFD = zeros(N_dim,1);
    switch str_mode
        case 'id_rayleigh'
            LCR = sqrt(2*pi)*fd*rho.*exp(-rho.^2); 
            AFD = (exp(rho.^2)-1)./(rho*fd*sqrt(2*pi));
            if normalize
                AFD = AFD*fd;
                LCR = LCR/fd;
            end
            
        case 'id_rician'
            %LCR_id rician
            alpha=0:0.001:pi;
            delta_area=zeros(length(rho),length(alpha)-1);
            for j=1:length(rho)
                for i=1:1:length(alpha)-1
                    delta_area(j,i)=(alpha(i+1)-alpha(i))*(1+(2/rho(j))*sqrt(K/(1+K))*(cos(theta_0))^2*cos(alpha(i+1)))*exp(2*rho(j)*sqrt(K*(1+K))*cos(alpha(i+1))-2*K*(cos(theta_0))^2*(sin(alpha(i+1)))^2);
                end
            end
            sum_area=sum(delta_area,2);
            LCR = (sqrt(2*(1+K)/pi)*rho*fd).*exp(-K-(1+K)*(rho.^2)).*(sum_area)';
            
            %AFD_id rician
            for i=1:1:length(AFD)
                AFD(i)= (1 - marcumq(sqrt(2*K),sqrt(2*(1+K)*(rho(i)^2)),1))/LCR(i);
            end
            %%%%%
            %normalization
            AFD = AFD*fd;
            LCR = LCR/fd;
            
        case 'simulated'
            %% simulated LCR and AFD
            for i=1:N_dim
                g_trasl = F-R(i);
                g_trasl_abs = abs(F-R(i));
                g_diff = g_trasl-g_trasl_abs;
                for j=1:(N_step-1)
                    if g_diff(j)==0 && g_diff(j+1)~=0;
                        LCR(i) = LCR(i)+1;
                    end
                end
                LCR(i) = LCR(i)/(N_step*delta_t);
                AFD(i) = mean(g_trasl<0)/LCR(i);
            end
            if normalize
                AFD = AFD*fd;
                LCR = LCR/fd;
            end
    end
end