function [dag] = z_covertToDag(pdag,un_directed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert pdag to dag by orientation rules if possible
% otherwise orient randomly
%
% input:
% [pdag]        - (matrix) pdag
% [un_directed] - (matrix) nodes indices of undirected edges
%
% output:
% [dag]         - (matrix) dag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% loop over all undirected edges
for i=1:size(un_directed,1)
    if sum(find(pdag(:,un_directed(i,2))==1))==0 %the new edge will not create a new collider
        pdag(un_directed(i,1),un_directed(i,2))=1;
        pdag(un_directed(i,2),un_directed(i,1))=0;
    else
        if sum(find(pdag(:,un_directed(i,1))==1))==0 %try the opposite way
            pdag(un_directed(i,2),un_directed(i,1))=1;
            pdag(un_directed(i,1),un_directed(i,2))=0;
        else
            warning('edge %d-%d cannot be oriented without introducing a new V-structure',un_directed(i,1),un_directed(i,2))
            pdag(un_directed(i,2),un_directed(i,1))=0;
            pdag(un_directed(i,1),un_directed(i,2))=0;
        end
    end
end
dag=pdag;