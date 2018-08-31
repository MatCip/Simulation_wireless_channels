clear;
close all;
clc;

%% simulation parameters:
global t_sim; t_sim = 15;  %duration of simulation (s)
global delta_t; delta_t = 1e-3; %time step (s)
global N_step; N_step = t_sim/delta_t; %number of step in time in the simulation  
global time; time = linspace(0,t_sim,N_step)'; %time vector
N_stat = 100; %number of statistical trials

%% jakes basic simulator
M_1 = 50;  %N = 4*M+2 where N is the number of oscillators
fd_1 = 50; %doppler frequency (Hz)
[g_1    , pdf_1    , R_1    ] = gen_channel(M_1,fd_1,'jakes',N_stat);

%% from Pop-Beaulieu paper
fd_2 = 50; %doppler frequency (Hz)
[~      , pdf_2_id , R_2_id ] = gen_channel(0, fd_2, 'ideal',        1     );
[g_2_3  , pdf_2_3  , R_2_3  ] = gen_channel(3, fd_2, 'pop_beaulieu', N_stat); %M=3
[g_2_4  , pdf_2_4  , R_2_4  ] = gen_channel(4, fd_2, 'pop_beaulieu', N_stat); %M=4
[g_2_8  , pdf_2_8  , R_2_8  ] = gen_channel(8, fd_2, 'pop_beaulieu', N_stat); %M=8

%% from Xiao-Zheng paper
M_3 = 8; %N = 4*M+2 where N is the number of oscillators
fd_3 = 25; %doppler frequency (Hz)
[~      , pdf_3_id , R_3_id   ] = gen_channel(0,fd_3,'ideal',1);
[g_3_10 , pdf_3_10 , R_3_10 , lambda_3_10,  LCR_3_id, LCR_3_10,  AFD_3_id, AFD_3_10 ] = gen_channel(M_3, fd_3, 'zheng_xiao', 10 );
[g_3_50 , pdf_3_50 , R_3_50 , lambda_3_50,  ~,        LCR_3_50,  ~,        AFD_3_50 ] = gen_channel(M_3, fd_3, 'zheng_xiao', 50 );
[g_3_100, pdf_3_100, R_3_100, lambda_3_100, ~,        LCR_3_100, ~,        AFD_3_100] = gen_channel(M_3, fd_3, 'zheng_xiao', 100);

%% from Komninakis paper
fd_4 = 50; %doppler frequency (Hz)
[g_4,pdf_abs_4,pdf_ph_4,R_4,H_abs,Y,zer,pol] = gen_channel_IIR(fd_4); 

%% plots
save('data/all_data.mat');
plot_all();
