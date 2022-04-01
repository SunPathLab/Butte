## Model Description:
The tumor starts with a single tumor-initiating cell. \
The tumor-initiating cell and its descendants with only passenger mutations form the **Type 0** population. \
**Type 0** cells give birth at a rate of `a_0` and die at a rate of `b_0`. \
**Type 0** cells mutate to **Type 1** cells at a rate of `u_1`. \
**Type 1** cells give birth at a rate of `a_1` and die at a rate of `b_1`. \
Both **Type 0** and **Type 1** cells accumulate passenger mutations at a rate of `u_0`. \
This package simulates the number of passenger mutations shared by more than a certain percentage (denoted by a threshold) of descendents of the first **Type 1** cell when the size of its living descendents reaches 10,000.


## Run (in Matlab environment):
```sh
#1. assign values to parameters. A code example is as follows:
    a_0=1;
    b_0=1;
    a_1=2;
    b_1=1;
    u_0=0.1;
    u_1=0.0001;
    threshold=0.9;

#2. run
[time_to_first_type_1_cell, num_of_passenger_mutation]=main(a_0,b_0,a_1,b_1,u_0,u_1,threshold)
```
*output:* \
`time_to_first_type_1_cell`: the birth time of the first type 1 cell with at least 10,000 living descendents at some time point.\
`num_of_passenger_mutation`: the number of passenger mutations shared by more than [threshold] of descendents of the first type 1 cell when the size of its living descendents reaches 10,000.


