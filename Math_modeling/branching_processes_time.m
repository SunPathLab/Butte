function [time, type, event, z_0, z_1]=branching_processes_time(r_0, d_0, r_1, d_1, z_0, z_1, v_0, v_1)

% edited 02/24/2021
% This function is used to get the next event

%% r_i: birth rate of type i
%% d_i: death rate of type i
%% v_0: neutral mutation rate
%% v_1: advantageous mutation rate
%% type 0 cells only have neutral mutation
%% type 1 cells have advantageous mutation

%% time: lenght of time for the next event to happen





rate_0_b=r_0*z_0; % total birth rate of type 0 cells
rate_0_d=d_0*z_0; % total death rate of type 0 cells

rate_1_b=r_1*z_1; % total birth rate of type 1 cells
rate_1_d=d_1*z_1; % total death rate of type 1 cells



rate_mutate_advantage=z_0*v_1; % total driver mutation rate
rate_mutate_neutral_0=z_0*v_0; % total passenger mutation rate for type 0 cells
rate_mutate_neutral_1=z_1*v_0; % total passenger mutation rate for type 1 cells

total_rate=rate_0_b+rate_0_d+rate_1_b+rate_1_d+rate_mutate_advantage+rate_mutate_neutral_0+rate_mutate_neutral_1;

time=exprnd(1/total_rate);

u=rand;

if u<rate_0_b/total_rate
    % birth of neutral cells
    type=0;
    event=0;
    z_0=z_0+1;
elseif u<(rate_0_b+rate_0_d)/total_rate
    % death of neutral cells
    type=0;
    event=1;
    z_0=z_0-1;
elseif u<(rate_0_b+rate_0_d+rate_1_b)/total_rate
    % birth of advantageous cells
    type=1;
    event=0;
    z_1=z_1+1;
elseif u<(rate_0_b+rate_0_d+rate_1_b+rate_1_d)/total_rate
    % death of advantageous cells
    type=1;
    event=1;
    z_1=z_1-1;
elseif u<(rate_0_b+rate_0_d+rate_1_b+rate_1_d+rate_mutate_advantage)/total_rate
    % advantageous mutation
    type=2;
    event=1;
    z_0=z_0-1; % in countinuous mutation model, a type 0 switches to a type 1
    z_1=z_1+1;
elseif u<(rate_0_b+rate_0_d+rate_1_b+rate_1_d+rate_mutate_advantage+rate_mutate_neutral_0)/total_rate
    % neutral mutation 0
    type=3;
    event=1;
else
    % neutral mutation 1
    type=4;
    event=1; 
end
    