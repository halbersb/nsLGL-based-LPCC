function [clust_new,clusters_new,clusters_size_new] = findNewClusters(clusters,clusters_size,clusters_unrounded,Observed,Latent,l_ch,ch_v,l1,pcc_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find subset from clusters to be the new clusters relative to specific latent
% a cluster is selected if the columns that equal to the latent children are
% equal to one of the latent children value configuration (as a whole)
%
% input:
% [clusters]           - (matrix) the clusters centers
% [clusters_size]      - (vector) the clusters sizes
% [clusters_unrounded] - (matrix) real: number of clusters -X- number of observed
% [Observed]           - (structure array) observed variables
% [Latent]             - (structure array) latents variables
% [l_ch]               - (array) array of vectors with the latent children
% [ch_v]               - (array) array of matrices with the latent children value configurations
% [l1]                 - (scalar) index of current latent (first latent is 1..)
% [pcc_type]           - [scalar] indicator for how to run pairwise cluster comparison
%
% output:
% [clust_new]          - (vector) the indeces of the clusters to be used in clusters_new
% [clusters_new]       - (matrix) clusters_new is a subset of clusters
% [clusters_size_new]  - (vector) izes of the new clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initializations
clust_new=[];
clusters_new=[];
clusters_size_new=[];
w=0;
if pcc_type>1, clusters=clusters_unrounded; end %if pcc_type is not regular use the unrounded clusters

%loop over all the clusters
for c=1:size(clusters,1)
    flag1=1;
    %loop over all the latents
    for l2=1:length(Latent)
        if Latent(l2).PCC_index~=Latent(l1).PCC_index
            latent_indx=Latent(l2).PCC_index-length(Observed); %the latent position in the latent array
            cl=[]; %vector
            t=0;
            %loop over all children of the curr latent
            for ch=1:length(l_ch{latent_indx})
                if l_ch{latent_indx}(ch)<=length(Observed) %if it is an observed child
                    t=t+1;
                    cl(t)=clusters(c,l_ch{latent_indx}(ch)); %extract from cluster 'c' (one row) only the values of children of the current latent (columns) 
                end
            end
            flag2=0;
            for i=1:size(ch_v{latent_indx},1)
                %if cl is equal to at least one of the latent's children value configurations
                %then 'c' is not eliminated from the new clusters
                %however if it does not equal to any of them, this cluster 'c' is eliminated
                %Dan changed condition here...
                if isequal(createPcc([cl;ch_v{latent_indx}(i,:)],pcc_type),zeros(1,length(l_ch{latent_indx})))
                    flag2=1;
                    break;
                end
            end
            %if cl never equals to any of ch_v{latent_indx} - for one of the latetns!!!
            %then cluster 'c' will not go into clusters_new
            if flag2==0
                flag1=0;
                break;
            end
        end
    end
    %add 'c' to the new clusters
    if flag1==1
        clust_new=[clust_new c];
        w=w+1;
        clusters_new(w,:)=clusters(c,:);
        clusters_size_new=[clusters_size_new clusters_size(c)];
    end
end