function dag = createDag(Observed,Latent)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine the observed and latents into a dag matrix
%
% input:
% [Observed] - (structure) observed variables
% [Latent]   - (structure) latents variables
%
% output:
% [dag]      - (matrix) the new directed acyclic graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
dag=zeros(length(Observed)+length(Latent),length(Observed)+length(Latent));

for i=1:length(Latent)
    for j=1:length(Latent(i).CH)
        dag(Latent(i).PCC_index,Latent(i).CH(j))=1; %edges between latent and his children
    end
end
for i=1:length(Observed)
    for j=1:length(Observed(i).CH)
        dag(Observed(i).PCC_index,Observed(i).CH(j))=1; %edges between observed and his children
    end
end
for i=1:size(dag,1)
    for j=1:size(dag,2)
        %bi-directed edges are removed
        if dag(i,j)==1 && dag(j,i)==1
            dag(i,j)=0;
            dag(j,i)=0;
        end
    end
end