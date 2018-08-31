clear;
close all;
clc;

%% simulation parameters:
global t_sim; t_sim = 15;  %duration of simulation (s)
global delta_t; delta_t = 1e-3; %time step (s)
global N_step; N_step = t_sim/delta_t; %number of step in time in the simulation  
global time; time = linspace(0,t_sim,N_step)'; %time vector
global theta_0; theta_0 = pi/4; %angle of arrival
global N; N = 8; %number of the sinusoids
global fd; fd = 10; %doppler frequency (Hz)
N_stat = 100; %number of statistical trials


%% computation of the channels
[~,        pdf_id_0,   R_id_0,   rho, LCR_id_0,      AFD_id_0  ] = gen_channel('id_rician',  1,      0 );
[~,        pdf_id_1,   R_id_1,   ~,   LCR_id_1,      AFD_id_1  ] = gen_channel('id_rician',  1,      1 );
[~,        pdf_id_3,   R_id_3,   ~,   LCR_id_3,      AFD_id_3  ] = gen_channel('id_rician',  1,      3 ); 
[~,        pdf_id_5,   R_id_5,   ~,   LCR_id_5,      AFD_id_5  ] = gen_channel('id_rician',  1,      5 );
[~,        pdf_id_10,  R_id_10,  ~,   LCR_id_10,     AFD_id_10 ] = gen_channel('id_rician',  1,      10);
[g_sim_0,  pdf_sim_0,  R_sim_0,  ~,   LCR_sim_0,     AFD_sim_0 ] = gen_channel('sim_rician', N_stat, 0 );
[g_sim_1,  pdf_sim_1,  R_sim_1,  ~,   LCR_sim_1,     AFD_sim_1 ] = gen_channel('sim_rician', N_stat, 1 );
[g_sim_3,  pdf_sim_3,  R_sim_3,  ~,   LCR_sim_3,     AFD_sim_3 ] = gen_channel('sim_rician', N_stat, 3 );
[g_sim_5,  pdf_sim_5,  R_sim_5,  ~,   LCR_sim_5,     AFD_sim_5 ] = gen_channel('sim_rician', N_stat, 5 );
[g_sim_10, pdf_sim_10, R_sim_10, ~,   LCR_sim_10,    AFD_sim_10] = gen_channel('sim_rician', N_stat, 10);

%% plots
save('data/all_data.mat');
plot_all();
