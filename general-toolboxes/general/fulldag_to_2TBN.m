function [real_TBN,drifts] = fulldag_to_2TBN(real_dag,O,L,E,C)

%convert a flatten model into 2TBN
num_structures=E*(C+1)-1;
num_var=2*(O+L);
num_obs=2*O;
real_TBN=cell(1,num_structures);
for i=1:num_structures
    start=(num_structures-1)*O+(i-1)*L+1;
    real_TBN{i}=zeros(num_var,num_var);
    real_TBN{i}(1:num_var,1:num_obs)=real_dag(start:start+num_var-1,(i-1)*O+1:i*O+O);
    real_TBN{i}(num_obs+1:num_var,num_obs+1:num_var)=real_dag(start+num_obs:start+num_obs+2*L-1,start+num_obs:start+num_obs+2*L-1);
end
%find drift points
drifts=find_drifts(real_TBN,E,C);