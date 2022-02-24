function [population_0,population_1, T, mutation_num, num_list_0, num_list_1, z_0, z_1]=branching_processes_main(r_0, d_0, r_1, d_1, v_0, v_1)

% edited 02/24/2021
% this is the main function

% r_i: birth rate of type i
% d_i: death rate of type i
% type 0 cells only have neutral mutation
% type 1 cells have another advantageous mutation
% v_0: passenger mutation rate
% v_1: driver mutation rate






population_0={};
population_0{end+1} = [1;-1;0];
population_1={};
%% (1,1) stores the number of cells in this clone
%% the second row stores the index of mutations
%% the third row stores the acquisition time of each mutation

z_0=1;
z_1=0;
% start from a single type zero cell

T=0;

mutation_num=0;
% used for recording index of mutations

num_list_0=[1];
% used to record the number of cells in each subclone (type 0 cells)
num_list_1=[];
% used to record the number of cells in each subclone (type 1 cells)

first_type_two=0;

while first_type_two==0 % iterate until the first non-extinct type 2 
    [time, type, event, z_0, z_1]=branching_processes_time(r_0, d_0, r_1, d_1, z_0, z_1, v_0, v_1);
    T=T+time;
    [population_0, population_1, num_list_0, num_list_1, mutation_num]=branching_processes_evolution(type, event, z_0, z_1, population_0, population_1, T, num_list_0, num_list_1, mutation_num);
    if z_0==0 && z_1==0
        % if the pupulatin dies out, start over from the initial state
        z_0=1;
        z_1=0;
        mutation_num=0;
        T=0;
        population_0={};
        population_1={};
        population_0{end+1} = [1;-1;0];
        num_list_0=[1];
        num_list_1=[];
    elseif z_1==1 && z_0==0
        temp=rand;
        if temp< (r_1-d_1)/(r_1)
            % check if this type 1 cell has infinite lineage
            first_type_two=1;
        else
            z_0=1;
            z_1=0;
            mutation_num=0;
            T=0;
            population_0={};
            population_1={};
            population_0{end+1} = [1;-1;0];
            num_list_0=[1];
            num_list_1=[];
        end
    elseif z_1==1 && z_0>0
        temp=rand;
        if temp< (r_1-d_1)/(r_1)
            first_type_two=1;
        else
            z_1=0;
            population_1={};
            num_list_1=[];
        end
    end
        
end