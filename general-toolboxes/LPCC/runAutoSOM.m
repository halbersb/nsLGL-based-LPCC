function [clusters,clusters_size,clusters_unrounded,init_clusters,SM_250,clusters_250] = runAutoSOM(data,max_k,min_k,server_ind,units)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find clusters with automatic som and optimal k with
% Davis Bouldin index
%
% input:
% [data]               - (matrix) the data to be clustered
% [max_k]              - (scalar) upper bound for my_kmeans_clusters
% [max_k]              - (scalar) lower bound for my_kmeans_clusters
% [server_ind]         - (scalar) indicator if running on PC or Server
%
% output:
% [clusters[           - (matrix) the found clusters centers
% [clusters_size]      - (vector) the size of each cluster
% [clusters_unrounded] - (matrix) the found clusters centers without rounding
% [init_clusters]      - (matrix) the initial data points that were randomly selected
% [SM_250]             - (structure) the som model (250 is historical for 250 units)
% [clusters_250]       - (vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read input
if nargin>4
    if ~isnumeric(units), error('units is not an integer'); end
else
    units=250; %default
end
if ~isnumeric(max_k), error('max_k is not an integer'); end
if ~isnumeric(min_k), error('min_k is not an integer'); end
if ~isnumeric(server_ind), error('server_ind must be 0 or 1'); end

SD_250=som_data_struct(data);
SM_250=som_make(SD_250,'munits',units);
bmus_250=som_bmus(SM_250,SD_250);
SD_250.labels=num2cell(bmus_250);
SM_250=som_autolabel(SM_250,SD_250,'vote');

%plot maps (on PC only if 1)
if(server_ind)
    K=som_show(SM_250,'umat','all','norm','d');
    h=som_show_add('label',SM_250,'subplot',1);
end

%run k-means
[c,ic,lbl,~,ind]=my_kmeans_clusters(SM_250,max_k,min_k); %find clusterings. ic='init_clusters' is dan addition. 'lbl' stands for labels
[~,i]=min(ind); %select the one with smallest davies-bouldin index

%update plot (on PC only if 1)
if(server_ind)
    som_show(SM_250,'color',{lbl{i},sprintf('%d clusters',i)}); %visualize
    colormap(jet(i)), som_recolorbar %change colormap
end

%round final clusters to nearest interger
clusters_unrounded=c{i};
clusters=round(c{i});
init_clusters=ic{i};
%find to which cluster each vector was mapped
for j=1:length(bmus_250)
    clusters_250(j)=lbl{i}(bmus_250(j));
end
clusters_size=zeros(1,i);
for j=1:i %i represent also the number of clusters (K)
    for k=1:length(clusters_250)
        if clusters_250(k)==j
            clusters_size(j)=clusters_size(j)+1;
        end
    end
end
%option 2 for sizes - take directtly the SOM sizes (not recommended!)
%clusters_size=zeros(1,i);
%for j=1:length(labels)
%    clusters_size(lbl{i}(j))=clusters_size(lbl{i}(j))+1;
%end