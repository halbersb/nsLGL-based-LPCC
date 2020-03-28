function [cardinality,l_ch,ch_v] = findLatentCardinality(dag,latents,clusters,clusters_unrounded,major_clusters,pcc_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find latents cardinality and find the corresponding joint  
% value of the latents children for each if its possible states 
%
% input:
% [dag]                - (matrix) the current directed acyclic graph
% [latents]            - (vector) indices of latents
% [clusters]           - (matrix) the clusters centers
% [clusters_unrounded] - (matrix) real: number of clusters -X- number of observed
% [major_clusters]     - (vector) array of major clusters indexes
% [pcc_type]           - [scalar] indicator for how to run cluster comparison
%
% output:
% [cardinality]        - (vector) a list with the cardinality of each latent variable
% [l_ch]               - (array) array of vectors with the latent children
% [ch_v]               - (array) array of matrices with the latent children value configurations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
l_ch=[]; %latents children - array of cells, where each cell in the array is a latent and the content of the cell is a vector of its children
ch_v=[]; %children values
if pcc_type>1, clusters=clusters_unrounded; end %if pcc_type is not regular use the unrounded clusters

%loop to find the observed children of each latent based on the graph (dag)
for k=1:length(latents)
    r=0;
    e=children(dag,latents(k));
    for t=1:length(children(dag,latents(k)))
        %make sure you don't insert other latents as children of current latent (only observed)
        if e(t)<=length(clusters(1,:))
            r=r+1;
            l_ch{k}(r)=e(t);
        end
    end
end

%loop to find the cardinality of each latent
cardinality=zeros(1,length(latents)); %initialization
for f=1:length(latents)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find all ch_v confugurations based on all clusters %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %mc_ch are f's children (subset of the observed) in the first major cluster only
    mc_ch=clusters(major_clusters(1),l_ch{f}); %major_clusters(1) is the first major cluster - it is used as a baseline cluster
    %loop over all clusters and compare each cluster to mc_ch based on the same subset of observed
    %add '+1' to cardinality if you find a cluster which is completely different (all its observed values are different from the mc_ch values)
    %that is, the cardinality is equal to the number of times that when we move from the baseline cluster (mc_ch) to another cluster, all the observed of the current latent change their value
    for k=1:size(clusters,1)
        cc_ch=clusters(k,l_ch{f}); %subset of the current cluster: only variables which are children of current latent
        %for each latent - first index in the for loop (each f...) increases the cardinality from 0 to 1
        if k==1
            ch_v{f}(cardinality(f)+1,:)=mc_ch;
            cardinality(f)=cardinality(f)+1; %increase the cardinality by 1
        end
        check_pcc=[cc_ch;mc_ch];
        [cc1,~]=createPcc(check_pcc,pcc_type);
        if cc1==ones(1,length(l_ch{f})) %if the clusters are completely different
            ch_v{f}(cardinality(f)+1,:)=cc_ch;
            cardinality(f)=cardinality(f)+1; %increase the cardinality by 1
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % remove ch_v confugurations based on major clusters and calculate cardinality %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ch_v{f}=union(ch_v{f},ch_v{f},'rows'); %equal to sql join on row i = row j. remove duplications for discrete case.
    to_remove{f}=[];
    t=0;
    for v=1:size(ch_v{f},1)
        cc_ch=ch_v{f}(v,:);
        for i=1:length(major_clusters)
            mc_ch=clusters(major_clusters(i),l_ch{f});
            %Dan changed condition here...
            check_pcc=[cc_ch;mc_ch];
            [cc1,~]=createPcc(check_pcc,pcc_type);
            if cc1==zeros(1,length(l_ch{f}))
                if isequal(cc_ch,mc_ch)
                    %if curr children configuration is from major cluster then don't remove the children configuration from ch_v{f}
                    break;
                else
                    %remove duplication for the continious case
                    to_remove{f}(t+1)=v;
                    t=t+1;
                    break;
                end
            end
            %remove only if cc_ch~=mc_ch and regarding one of the observed they are equal (0)
            if length(find(cc1==0))>0 %the pcc contains at least one zero 
                to_remove{f}(t+1)=v;
                t=t+1;
            end
        end
    end
    
    to_remove_2{f}=unique(to_remove{f});
    for k=1:length(to_remove_2{f})
        ch_v{f}(to_remove_2{f}(k),:)=[]; %here we delete from ch_v
        for j=1:length(to_remove_2{f})
            if to_remove_2{f}(j)>to_remove_2{f}(k)
                to_remove_2{f}(j)=to_remove_2{f}(j)-1;
            end
        end
    end
    cardinality(f)=size(ch_v{f},1); %calculate final cardinality
end