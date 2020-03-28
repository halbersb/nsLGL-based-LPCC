function [distinct_clusters,distinct_clusters_size,distinct_clusters_unrounded] = z_findDistinctClusters(clusters,clusters_size,clusters_unrounded,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% join duplicate clusters due to the round() and sum up their sizes
%
% input:
% [clusters]                    - (matrix) cluster centers
% [clusters_size]               - (vector) clusters size
% [clusters_unrounded]          - (matrix) unrounded cluster centers
% [flag]
%
% output:
% [distinct_clusters]           - (matrix) new cluster centers
% [distinct_clusters_size]      - (vector) new clusters size
% [distinct_clusters_unrounded] - (matrix) new unrounded cluster centers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
ind=1;
to_ignore=[];
distinct_clusters=[];
distinct_clusters_size=[];
distinct_clusters_unrounded=[];
if ~isempty(varargin), flag=varargin{1}; else flag=0; end
if ~isnumeric(flag), error('flag should be 0 or 1'); end %'0' for sum and '1' for max

%validation
if ~isequal(size(clusters),size(clusters_unrounded)), error('clusters and clusters_unrounded are not the same size'); end

%loop over all clusters
for i=1:size(clusters,1)
    %if 'i'th center was already used then ignore it
    if ismember(i,to_ignore)
        continue;
    end
    joint_size=0;
    distinct_clusters_size(1,ind)=clusters_size(1,i);
    %loop over all other clusters
    for j=(i+1):size(clusters,1)
        if clusters(i,:)==clusters(j,:)
            if flag
                distinct_clusters_size(1,ind)=max(distinct_clusters_size(1,ind),clusters_size(1,j)); %max between their sizes
            else
                distinct_clusters_size(1,ind)=distinct_clusters_size(1,ind)+clusters_size(1,j); %sum their sizes
            end
            to_ignore=[to_ignore j]; %do not use the 'j'th center anymore
        end
    end
    distinct_clusters(ind,:)=clusters(i,:);
    distinct_clusters_unrounded(ind,:)=clusters_unrounded(i,:);
    ind=ind+1;
end