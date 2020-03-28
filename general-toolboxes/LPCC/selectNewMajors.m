function [new_majors] = selectNewMajors(EXO_L,Observed,clusters,clusters_size,clusters_unrounded,l_ch,ch_v,ns,pcc_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select a new set of major clusters
% the entrance rule to the new set of major clusters is as follow:
% if the cluster center vector is equal to exogenous children values combination
%
% input:
% [EXO_L]              - (vector) latents index
% [Observed]           - (structure array) observed variables
% [clusters]           - (matrix) clusters centers
% [clusters_size]      - (vector) clusters sizes
% [clusters_unrounded] - (matrix) real: number of clusters -X- number of observed
% [l_ch]               - (cell array of matrices) latents children
% [ch_v]               - (cell array of matrices) children values 
% [ns]                 - (vector) nodes size - both observed and latents
% [pcc_type]           - [scalar] indicator for how to run cluster comparison
%
% output:
% [new_majors]         - (vector) indeces if the new major clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
EXV=calcEXV(ns(:,EXO_L)); %create all combinations of the latents size
EXO_L_index=EXO_L-length(Observed);
t=0; %index for new majors
new_majors=[];
if pcc_type>1, clusters=clusters_unrounded; end %if pcc_type is not regular use the unrounded clusters

%loop over all children configurations (of all exogenous latents together)
for i=1:size(EXV,1)
    v1=zeros(1,length(Observed));
    for j=1:length(EXO_L)
        for h=1:length(l_ch{EXO_L_index(j)}) %loop over all of the exogenous children
            v1(l_ch{EXO_L_index(j)}(h))=ch_v{EXO_L_index(j)}(EXV(i,j),h);
        end
    end
    v2=zeros(1,length(Observed));
    max=0; %initialize the size of the biggest cluster per children configuration
    n=0;
    for m=1:size(clusters,1)
        for j=1:length(EXO_L)
            for h=1:length(l_ch{EXO_L_index(j)})
                v2(l_ch{EXO_L_index(j)}(h))=clusters(m,l_ch{EXO_L_index(j)}(h));
            end
        end
        %if the exogenous latents configuration is equal to the current cluster center then add it to the new major set
        if isequal(createPcc([v2;v1],pcc_type),zeros(1,length(Observed))) %dan changed the condition
            if clusters_size(m)>max %if the equal current cluster center is physically larger than the previous found then take it and update n
                max=clusters_size(m); %update new biggest cluster size (for specific combination)
                n=m;
                %a major cluster from previous iteration can be dropped out
                %here if another major cluster (which is bigger) has the same children
                %configuration (that could happen if not all observed variables were assigned under a latent variable)
            end
        end
    end
    if (n~=0)
        t=t+1;
        new_majors(t)=n;
    end
end