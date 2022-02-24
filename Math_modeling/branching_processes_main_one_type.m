function [population_0,population, T, mutation_num, num_list_0, num_list, z_0, z]=branching_processes_main_one_type(r, d, v, n)

% edited 02/24/2021

% r: birth rate
% d: death rate 
% v: neutral mutation rate
% n: targeted population size


population_0={};
population={};
population{end+1} = [1;-1;0];
z=1; % initial population size
T=0;
mutation_num=0;
num_list=[1];
num_list_0=[];

r_0=0;
d_0=0;
v_1=0;
z_0=0;

while z<n % iterate until the population size reaches n
    [time, type, event, z_0, z]=branching_processes_time(r_0, d_0, r, d, z_0, z, v, v_1);
    T=T+time;
    [population_0, population, num_list_0, num_list, mutation_num]=branching_processes_evolution(type, event, z_0, z, population_0, population, T, num_list_0, num_list, mutation_num);
    if z_0==0 && z==0
        z_0=0;
        z=1;
        mutation_num=0;
        T=0;
        population_0={};
        population={};
        population{end+1} = [1;-1;0];
        num_list_0=[];
        num_list=[1];
    end        
end