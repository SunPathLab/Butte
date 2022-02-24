function [SFS, SFS_threshold]=identify_SFS_time(population, num_list, mutation_num, n, threshold)

% edited 02/24/2021
% this function is used to generate SFS result given total population size
% n

[~,x]=size(population);


data=[];
count=0;

for i=1:mutation_num
    total_num=0;
    if x>0
        for k=1:x
            temp = find(population{k}(2,:)==i);
            if ~isempty(temp)
                total_num=total_num+num_list(k);
                t=population{k}(3,temp);
            end
        end
    end
    if total_num>0
        count=count+1;
        data(count,1)=i;
        data(count,2)=t;
        data(count,3)=total_num/n;
    end
end
SFS=data;

SFS_threshold=[];
for i=1:count
    if SFS(i,3)>=threshold
        SFS_threshold=[SFS_threshold; SFS(i,:)];
    end
end
    



                
            
    