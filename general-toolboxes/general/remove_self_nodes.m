function [dag] = remove_self_nodes(dag)

%remove self countaining nodes
for i=1:size(dag,2)
    if dag(i,i)
        warning('A diagonal elenemt is not zero');
        dag(i,i)=0;
    end
end