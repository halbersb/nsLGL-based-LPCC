function [new_clusters,new_clusters_size,new_clusters_unrounded,new_label] = z_reCalcCentroids(data,clusters_unrounded)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recalculate clusters directly from data based on SOM centers
%
% input:
% [data]                   - (matrix) the original data set
% [clusters_unrounded]     - (matrix) real: number of clusters -X- number of observed
%
% output:
% [new_clusters]           - (matrix) the new found clusters centers
% [new_clusters_size]      - (vector) new sizes of the clusters
% [new_clusters_unrounded] - (matrix) the new found clusters centers without rounding
% [new_label]              - (vector) the new cluster label for each sample in the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initializations
new_label=zeros(size(data,1),1);
new_clusters_unrounded=zeros(size(clusters_unrounded));
new_clusters_size=zeros(1,size(clusters_unrounded,1));

%loop over data
for i=1:size(data,1)
    sample=data(i,:);
    dist=[];
    %loop over clusters
    for j=1:size(clusters_unrounded,1)
        center=clusters_unrounded(j,:);
        dist=[dist norm(sample - center)]; %calculate dist to current cluster
    end
    new_label(i)=find(dist==min(dist),1,'first'); %assign closest cluster
end
data=[data new_label]; %add new column to data with labels

%caluculate centroids and their respective size
for i=1:size(clusters_unrounded,1)
    avg=mean(data(data(:,end)==i,:),1);
    new_clusters_unrounded(i,:)=avg(1:end-1); %ignore the label column
    new_clusters_size(i)=length(data(data(:,end)==i,end)); %calculate new sizes
end

%round final new_clusters_unrounded to nearest interger
new_clusters=round(new_clusters_unrounded);