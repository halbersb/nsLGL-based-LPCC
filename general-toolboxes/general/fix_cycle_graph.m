function new_graph = fix_cycle_graph(graph)

%replace all undirected edges with a directed one
for a=1:size(graph,2)
    for b=a+1:size(graph,2)
        if graph(a,b) && graph(b,a) %undirected edge
            graph(b,a)=0;
        end
    end
end
%validate acyclic graph
if acyclic(graph)
    new_graph=graph;
else
    for i=size(graph,2):-1:1 %remove only latent-latent edges...
        temp_graph=graph;
        prt=find(graph(:,i));
        if ~isempty(prt)
            temp_graph(prt(1),i)=0;%try to remove the first parent
            if acyclic(temp_graph)
                new_graph=temp_graph;
                break;
            end
        end
    end
end