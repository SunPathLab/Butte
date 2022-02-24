function [SFS_neutral, SFS_neutral_threshold]=identify_SFS_in_second_type(population_1,  num_list_1, mutation_num, threshold)

% edited 02/24/2021
% this function is used to generate SFS result for type 1 cells


total_type_two=sum(num_list_1);

[~,x1]=size(population_1);

data=[];
count=0;

for i=1:mutation_num
    total_num_1=0;
    if x1>0
        for k=1:x1
            temp = find(population_1{k}(2,:)==i);
            if ~isempty(temp)
                total_num_1=total_num_1+num_list_1(k);
                t=population_1{k}(3,temp);
            end
        end
    end
    if total_num_1>0
        count=count+1;
        data(count,1)=i; % store mutation index
        data(count,2)=t; % store mutation time
        data(count,3)=total_num_1/total_type_two; % store site frequency
    end
end
SFS_neutral=data;

SFS_neutral_threshold=[];
for i=1:count
    if SFS_neutral(i,3)>=threshold
        SFS_neutral_threshold=[SFS_neutral_threshold; SFS_neutral(i,:)];
    end
end