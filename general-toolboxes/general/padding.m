function [dag] = padding(dag,req_size,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill dag with zeros (add unconnected latent)
% for 2-TBN
%
% input:
% [dag]      - (matrix) dircected acyclic graph
% [req_size] - (scalar) the required new size
%
% output:
% [dag]      - (matrix) the new dag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dag_size=size(dag,1);
    if (dag_size>=req_size)
        warning('no padding is needed');
    else
        if nargin==2
            if size(dag,1)~=req_size
                for i=dag_size+1:req_size
                    dag(:,i)=0; %add column
                    dag(i,:)=0; %add row
                end
            end
        else
            for i=1:2:length(varargin)
                switch varargin{i},
                    case 'num_O', num_obs_slice=varargin{i+1};
                    otherwise, error(['invalid argument name: ' varargin{i}]);       
                end
            end
            L1=count_latents(dag(num_obs_slice*2+1:end,1:num_obs_slice));
            L2=count_latents(dag(num_obs_slice*2+1:end,num_obs_slice+1:num_obs_slice*2));
            N=L1-L2;
            switch N
                case num2cell(0), % L1==L2, add (req_size-dag_size)/2 to each slice
                    for i=1:(req_size-dag_size) %req_size is always odd
                        dag=insert_vertex(dag,dag_size+i);%add at the end
                    end
                case num2cell(1:100), % L1 > L2, add to slice 2
                    dag=insert_vertex(dag,dag_size+1);%add at the end to L2
                    dag_size=size(dag,1); %update dag size
                    for i=1:(req_size-dag_size) %req_size is always odd
                        dag=insert_vertex(dag,dag_size+i);%add at the end
                    end
                case num2cell(-100:-1), % L1 < L2, add to slice 1
                    dag=insert_vertex(dag,num_obs_slice*2+1);%add at the middle
                    dag_size=size(dag,1); %update dag size
                    for i=1:(req_size-dag_size) %req_size is always odd
                        dag=insert_vertex(dag,dag_size+i);%add at the end
                    end
            end 
        end
    end
end

%%
function count = count_latents(mat)
%this function counts the number of row with at least one '1'
    
    count=0;
    for i=1:size(mat,1)
        if sum(find(mat(i,:)))>0, count=count+1; end
    end
end

%%
function new_mat = insert_vertex(mat,ind)
%this function insert new vertex in the middle of a dag

    new_mat=zeros(size(mat)+1);
    %four part to fill
    for i=1:ind-1
        for j=1:ind-1
            new_mat(i,j)=mat(i,j);
        end
        for j=ind:size(mat,1)
            new_mat(i,j+1)=mat(i,j);
        end
    end
    for i=ind:size(mat,1)
        for j=1:ind-1
            new_mat(i+1,j)=mat(i,j);
        end
        for j=ind:size(mat,1)
            new_mat(i+1,j+1)=mat(i,j);
        end
    end
end