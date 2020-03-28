function [uuc,uuc_s] = unUsedClusters(mc,c,c_s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the clusters that are not in the major clusters (part of 1-order minor clusters search)
%
% input:
% [mc]     - (vector) list of major clusters
% [c]      - (vector) list of cluster
% [c_s]    - (vector) list of clusters size
%
% output:
% [uuc]    - (vector) cluster minus major clusters (uuc = un-used clusters)
% [uuc_s]  - (vector) sizes of uuc (uuc_s = un-used clusters size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
uuc=[]; 
uuc_s=[];
t=0;

for i=1:length(c)
    flag=0;
    for j=1:length(mc)
        if mc(j)==c(i)
            flag=1;
            break;
        end
    end
    if flag==0 %if c(i) is not part of mc set
        t=t+1;
        uuc_s(t)=c_s(i);
        uuc(t)=c(i);
    end
end