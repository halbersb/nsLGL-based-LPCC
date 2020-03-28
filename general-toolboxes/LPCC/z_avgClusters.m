function [clusters,clusters_size_new,clusters_unrounded_new] = z_avgClusters(clusters_size,clusters_unrounded)
%this function implements smoothing for the clusters

%workaround if on matrix has NAN
for i=length(clusters_unrounded):-1:1 
    if sum(isnan(clusters_unrounded{i}(:))) %there is at least one zero
        clusters_unrounded(i)=[];
        clusters_size(i)=[];
    end
end

base_centers=clusters_unrounded{1};
num_centroids=length(clusters_unrounded);
clusters_unrounded_new=zeros(size(base_centers,1),size(base_centers,2));
clusters_size_new=zeros(size(base_centers,1),1);

for i=1:size(base_centers,1) %for each cluster of base_centroids
    sample=base_centers(i,:);
    closest_cluster=[];
    closest_cluster_size=[];
    for j=2:num_centroids %for each other centroid matrix
        dist=[];
        for k=1:size(clusters_unrounded{j},1) %for each cluster of current centroid metrix
            center=clusters_unrounded{j}(k,:);
            dist=[dist norm(sample - center)]; %calculate dist to current cluster
        end
        closest_cluster(j-1,:)=clusters_unrounded{j}(find(dist==min(dist),1,'first'),:); %assign closest cluster
        closest_cluster_size(j-1)=clusters_size{j}(find(dist==min(dist),1,'first')); %assign closest cluster
    end
    clusters_unrounded_new(i,:)=mean(closest_cluster,1);
    clusters_size_new(i)=round(mean(closest_cluster_size));
end

clusters=round(clusters_unrounded_new);
clusters_size_new=clusters_size_new.';

