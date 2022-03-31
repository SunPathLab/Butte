function [time_to_first_type_1_cell, num_of_passenger_mutation]=main(a_0,b_0,a_1,b_1,u_0,u_1,threshold)


%% a_i: birth rate of type i
%% b_i: death rate of type i
%% u_0: neutral mutation rate
%% u_1: advantageous mutation rate




[~,population_1, T, mutation_num, ~, num_list_1, ~, ~]=branching_processes_main(a_0, b_0, a_1, b_1, u_0, u_1);
time_to_first_type_1_cell=T;

num_of_mutation_in_type_1_before=0;
[~, SFS_neutral_threshold]=identify_SFS_in_second_type(population_1, num_list_1, mutation_num, threshold);
if ~isempty(SFS_neutral_threshold)
    [temp,~]=size(SFS_neutral_threshold);
    num_of_mutation_in_type_1_before=temp; % number of neutral mutations accumulated before acquiring the driver mutation
end
        
n=10000; % number of target size of type 1 cell to ensure that it does not go extinct
[~,population, ~, mutation_num, ~, num_list, ~, ~]=branching_processes_main_one_type(a_1, b_1, u_0, n);
[~, SFS_threshold]=identify_SFS_time(population, num_list, mutation_num, n, threshold);
num_of_mutation_in_type_1_after=0;
if ~isempty(SFS_threshold)
    [temp,~]=size(SFS_threshold);
    num_of_mutation_in_type_1_after=temp; % number of neutral mutations accumulated after acquiring the driver mutation
end

num_of_passenger_mutation=num_of_mutation_in_type_1_before+num_of_mutation_in_type_1_after;
