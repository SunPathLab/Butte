% In this code, we change the neutral mutation rate 
num=1000; % number of samples for each set of parameters
time_to_first_2=zeros(num*10,1);
num_of_mutation_2_1=zeros(num*10,1);
num_of_mutation_2_2=zeros(num*10,1);
neutral_mutation_rate_2=zeros(num*10,1);
count=0;

n=10000; % number of target size of type 1 cell to ensure that it does not go extinct

for k=1:10

    fitness=1;
    r_0=1.2;
    d_0=1;
    r_1=2;
    d_1=1;
    v_0=0.05+k*0.05;
    v_1=0.0001;
    threshold=0.9;

    for j=1:num
        
        count=count+1;
        [population_0,population_1, T, mutation_num, num_list_0, num_list_1, z_0, z_1]=branching_processes_main(r_0, d_0, r_1, d_1, v_0, v_1);
        time_to_first_2(count)=T;

        [SFS_neutral, SFS_neutral_threshold]=identify_SFS_in_second_type(population_1,  num_list_1, mutation_num, threshold);
        if ~isempty(SFS_neutral_threshold)
            [temp,~]=size(SFS_neutral_threshold);
            num_of_mutation_2_1(count)=temp; % neutral mutation before
        end
        
        [~,population, ~, mutation_num, ~, num_list, ~, ~]=branching_processes_main_one_type(r_1, d_1, v_0, n);
        [SFS, SFS_threshold]=identify_SFS_time(population, num_list, mutation_num, n, threshold);
        if ~isempty(SFS_threshold)
            [temp,~]=size(SFS_threshold);
            num_of_mutation_2_2(count)=temp; % neutral mutation after
        end
        
        neutral_mutation_rate_2(count)=v_0;
    end

end

csvwrite('time_to_first_2_test22.txt',time_to_first_2);
csvwrite('num_of_mutation_2_1_test22.txt',num_of_mutation_2_1);
csvwrite('num_of_mutation_2_2_test22.txt',num_of_mutation_2_2);
csvwrite('neutral_mutation_rate_2_test22.txt',neutral_mutation_rate_2);




