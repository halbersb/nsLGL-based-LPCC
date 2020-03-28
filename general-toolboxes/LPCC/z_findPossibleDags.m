function dags_array = z_findPossibleDags(base_dag,SOC,Observed,Latent)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find all possible dags for a baseline dag that will be used for search and score
%
% *assume SOC can have only 2 parents
% *assume the 2 parents are not children of the same latent parent
% *assume the 2 parents are not both latent
%
% input:
% [base_dag]   - (matrix) the dag for which we search neighbors
% [SOC]        - (structure array) suspected observed colliders
% [Observed]   - (array) array of structures
% [Latent]     - (array) array of structures
%
% output:
% [dags_array] - (array) array of matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
v=[]; %vector that hold all possible parents of SOC

for i=1:length(SOC.PPS_index)
    if SOC.PPS_index(i)<=length(Observed) 
        %observed parent
        v=[v SOC.PPS_index(i)];
    else
        %latent parent
        v=[v Latent(SOC.PPS_index(i)-length(Observed)).PCC_index]; %add the latent itself
        v=[v Latent(SOC.PPS_index(i)-length(Observed)).OCH]; %add its children
    end
end

combos=combntns(v,2); %all possible pairs of parents
for i=1:size(combos,1)
    dag=base_dag;
    %if one of the combinations (but not both) is latent or both are observed but not of the same latent parent!!!
    if (combos(i,1)>length(Observed) || combos(i,2)>length(Observed) || Observed(combos(i,1)).PA~=Observed(combos(i,2)).PA) && (combos(i,1)<=length(Observed) || combos(i,2)<=length(Observed))
        %add parents
        dag(combos(i,1),SOC.PCC_index)=1;
        dag(combos(i,2),SOC.PCC_index)=1;
        dags_array{i}=dag; 
    end
end