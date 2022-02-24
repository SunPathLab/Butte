function [population_0, population_1, num_list_0, num_list_1, mutation_num]=branching_processes_evolution(type, event, z_0, z_1, population_0, population_1, T, num_list_0, num_list_1, mutation_num)

% edited 02/24/2021
% This function is used to simulate the evolution of the process given the
% next event

%% population_0: stores information of advantageous cells.
%% it is a list, each element stores information of a clone (distinguished by different mutations)
%% (1,1) stores the number of cells in this clone
%% the second row stores the index of mutations
%% the third row stores the adquisition time of each mutation


%% num_list_0:  vector which stores the number of cells in each clone of neutral cells
%% num_list_1:  vector which stores the number of cells in each clone of advantageous cells

%% mutation_num: stores the current index of mutation

%% type and event: information about what type of event just happened

%% z_0: total number of neutral cells
%% z_1: total number of advantageous cells

%% T: current time



if type==4 && event==1 % neutral mutation for type 1 cells happens
    mutation_num=mutation_num+1;
    k=find_position(num_list_1, z_1);
    % find which clone does the cell with the new mutation belong to
    
    temp=population_1{k};
    temp(1,1)=1; % set the clone size to be 1
    temp=[temp,[0;mutation_num;T]]; % add the new mutation
    
    num_list_1(k)=num_list_1(k)-1;
    if num_list_1(k)>0
        population_1{k}(1,1) = num_list_1(k);
    else
        num_list_1(k)=[];
        population_1(:,k)=[];
    end
    
    num_list_1=[num_list_1,1];
    population_1{end+1} = temp;

elseif type==3 && event==1 % neutral mutation for type 0 cells happens
    mutation_num=mutation_num+1;
    k=find_position(num_list_0, z_0);
    
    temp=population_0{k};
    temp(1,1)=1;
    temp=[temp,[0;mutation_num;T]];
    
    num_list_0(k)=num_list_0(k)-1;
    if num_list_0(k)>0
        population_0{k}(1,1) = num_list_0(k);
    else
        num_list_0(k)=[];
        population_0(:,k)=[];
    end     
    
    num_list_0=[num_list_0,1];
    population_0{end+1} = temp;
    
elseif type==2 && event==1 % advantageous mutation
    k=find_position(num_list_0, z_0+1); % we have already updated z_0 in branching_processes_time, thus need to plus one here
    num_list_1=[num_list_1,1];
     
    temp=population_0{k};
    temp(1,1)=1;
    temp=[temp,[0;0;T]];  % 0 represents the advantageous mutation
    
    num_list_0(k)=num_list_0(k)-1;
    if num_list_0(k)>0
        population_0{k}(1,1) = num_list_0(k);
    else
        num_list_0(k)=[];
        population_0(:,k)=[];
    end     
    
    population_1{end+1} = temp;
    
elseif type==1 && event==1 % death of advantageous cells
    k=find_position(num_list_1, z_1+1);
    num_list_1(k)=num_list_1(k)-1;
    if num_list_1(k)>0
        population_1{k}(1,1) = num_list_1(k);
    else
        num_list_1(k)=[];
        population_1(:,k)=[];
    end
    
elseif type==1 && event==0 % birth of advantageous cells
    k=find_position(num_list_1, z_1-1);
    num_list_1(k)=num_list_1(k)+1;
    population_1{k}(1,1)=num_list_1(k);

    
elseif type==0 && event==1 % death of neutral cells
    k=find_position(num_list_0, z_0+1);
    num_list_0(k)=num_list_0(k)-1;
    if num_list_0(k)>0
        population_0{k}(1,1) = num_list_0(k);
    else
        num_list_0(k)=[];
        population_0(:,k)=[];
    end     
    
else  % birth of neutral cells
    k=find_position(num_list_0, z_0-1);
    num_list_0(k)=num_list_0(k)+1;
    population_0{k}(1,1)=num_list_0(k);
    
end
                    