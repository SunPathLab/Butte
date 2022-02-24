function k=find_position(num_list, total_num)

% num_list: n by 1 matrix



u=rand;

i=1;
current_num=num_list(1);
while current_num/total_num < u
    i=i+1;
    current_num=current_num+num_list(i);
end

k=i;
    
    
