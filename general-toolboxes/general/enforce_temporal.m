function [dag] = enforce_temporal(dag,num_obs_slice)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enforce latent-latent edges from t to t+1
% for 2-TBN (mainly for undirected edges)
% assumptions:
% 1) dag has latents and their indices starts after num_obs_slice*2
% 2) latent variable cannot have children from both slices
%
% input:
% [dag]           - (matrix) dircected acyclic graph
% [num_obs_slice] - (scalar) number of observed in slice
%
% output:
% [dag]           - (matrix) the new dag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    num_obs=2*num_obs_slice;
    latent_mat=dag(num_obs+1:end,num_obs+1:end);
    %for each edge in that matrix
    for i=1:size(latent_mat,2)
         for j=1:size(latent_mat,2)
             if i~=j
                %for i->j edge that is not self oriented
                L1_min_obs_indx = min(find(dag(num_obs+i,1:num_obs))); %from (min)
                L1_max_obs_indx = max(find(dag(num_obs+i,1:num_obs))); %from (max)
                L2_min_obs_indx = min(find(dag(num_obs+j,1:num_obs))); %to (min)
                L2_max_obs_indx = max(find(dag(num_obs+j,1:num_obs))); %to (max)
                if ~isempty(L1_min_obs_indx) && ~isempty(L1_max_obs_indx) && ~isempty(L2_min_obs_indx) && ~isempty(L2_max_obs_indx)
                    if L1_max_obs_indx<=num_obs_slice && L2_min_obs_indx>num_obs_slice %from t to t+1
                        %edge must be i->j
                        latent_mat(j,i)=0;
                    end
                end
             end
         end
    end
    dag(num_obs+1:end,num_obs+1:end)=latent_mat;
end