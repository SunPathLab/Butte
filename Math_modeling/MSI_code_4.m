% In this code, we change the advantageous mutation rate
num=10000; % number of samples for each set of parameters
time_to_first_4=zeros(num*10,1);
num_of_mutation_4_1=zeros(num*10,1);
num_of_mutation_4_2=zeros(num*10,1);
advantageous_mutation_rate_4=zeros(num*10,1);
count=0;

n=10000; % number of target size of type 1 cell to ensure that it does not go extinct

for k=1:10

    fitness=1;
    r_0=1.2;
    d_0=1;
    r_1=2;
    d_1=1;
    v_0=0.1;
    v_1=0.0001+(k-1)*0.0002;
    threshold=0.9;

    for j=1:num
        
        count=count+1;
        [population_0,population_1, T, mutation_num, num_list_0, num_list_1, z_0, z_1]=branching_processes_main(r_0, d_0, r_1, d_1, v_0, v_1);
        time_to_first_4(count)=T;

        [SFS_neutral, SFS_neutral_threshold]=identify_SFS_in_second_type(population_1,  num_list_1, mutation_num, threshold);
        if ~isempty(SFS_neutral_threshold)
            [temp,~]=size(SFS_neutral_threshold);
            num_of_mutation_4_1(count)=temp; % neutral mutation before
        end
        
        [~,population, ~, mutation_num, ~, num_list, ~, ~]=branching_processes_main_one_type(r_1, d_1, v_0, n);
        [SFS, SFS_threshold]=identify_SFS_time(population, num_list, mutation_num, n, threshold);
        if ~isempty(SFS_threshold)
            [temp,~]=size(SFS_threshold);
            num_of_mutation_4_2(count)=temp; % neutral mutation after
        end
        
        advantageous_mutation_rate_4(count)=v_1;
    end

end

csvwrite('time_to_first_4.txt',time_to_first_4);
csvwrite('num_of_mutation_4_1.txt',num_of_mutation_4_1);
csvwrite('num_of_mutation_4_2.txt',num_of_mutation_4_2);
csvwrite('advantageous_mutation_rate_4.txt',advantageous_mutation_rate_4);


