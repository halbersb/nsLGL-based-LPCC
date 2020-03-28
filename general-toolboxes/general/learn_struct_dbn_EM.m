function [curr_2TBN,init_2TBN,bestscore,bestengine] = learn_struct_dbn_EM(init_2TBN,data,ns,num_latents,max_iter,PMM_ind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% learn Dynamic Bayesian Network based on S&S - hill climbing procedure
% assumptions: 1. number of latent variables is known (as well as their cardinalities)
%              2. observed variable cannot be a parent of latent
%              3. variables cannot be directed from 't+1' to 't'
%
% input:
% [init_2TBN]     - (matrix) initial 2TBN (could be a random graph)
% [data]          - (matrix) rows are samples and columns are variables
% [ns]            - (vector) the cardinality of each of the observed and latent variables
% [num_latents]   - (scalar) number of latent variable in slice t
% [max_iter]      - (scalar) maximum nuber of iteration for EM procedure
% [PMM_ind]       - (scalar) indicator to search over pure structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
convergence=1;
bestscore=-inf;
bestengine=[];
iter=1;
nodes_per_slice=size(init_2TBN,1)/2; %observed and latents per slice
onodes=(num_latents+1):nodes_per_slice; %observed variables per slice
hnodes=1:num_latents; %latent variables per slice
dnodes=1:nodes_per_slice; %latent and observed variables per slice
curr_2TBN=init_2TBN;
%reshape data
for i=1:size(data,1)
    sample=reshape(data(i,:),length(onodes),length(data(i,:))/length(onodes));
    sample_w_h=cell(1,size(sample,2)); %init record
    sample=num2cell(sample); %convert to cell matrix
    sample(cellfun(@isnan,sample))={''}; %replace NaNs with zeros
    sample_w_h(onodes,:)=sample;
    data2{i}=sample_w_h;
end
data=data2;clear sample sample_w_h data2;

while convergence
    
    %find all neighbours in a single slice
    intra_2TBN=curr_2TBN(1:nodes_per_slice,1:nodes_per_slice); %represents single slice
    inter_2TBN=curr_2TBN(1:nodes_per_slice,nodes_per_slice+1:end); %all inter edges
    paths=reachability_graph(intra_2TBN+transpose(intra_2TBN));
    [Gs,op,nodes]=mk_nbrs_of_dag(intra_2TBN);
    
    %remove neighbours according to constraints
    to_remove=[]; 
    for i=1:length(Gs)
        %don't try to add observed->latent
        if strcmp(op{i},'add') && ismember(nodes(i,1),onodes) && ismember(nodes(i,2),hnodes)
            to_remove=[to_remove i];
        end
        %both are observed and none of them have a path to any latent
        if ismember(nodes(i,1),onodes) && ismember(nodes(i,2),onodes) && sum(sum(paths(hnodes,nodes(i,:))))==0 && sum(sum(paths(nodes(i,:),hnodes)))==0
            to_remove=[to_remove i];
        end
        %don't try to reverse latent->observed
        if strcmp(op{i},'rev') && ismember(nodes(i,1),hnodes) && ismember(nodes(i,2),onodes)
            to_remove=[to_remove i];
        end
        %don't try to delete latent->observed if it has only 2 (or less) observed children
        if strcmp(op{i},'del') && length(parents(intra_2TBN,nodes(i,2)))<3
            to_remove=[to_remove i];
        end
        %don't try to add latent->latent if one of the latent has no children
        if strcmp(op{i},'add') && isempty(children(intra_2TBN,nodes(i,1))) && isempty(children(intra_2TBN,nodes(i,2)))
            to_remove=[to_remove i];
        end
        if PMM_ind %serach over PMM space...
            %observed->latent/observed
            if strcmp(op{i},'add') && ismember(nodes(i,1),onodes)
                to_remove=[to_remove i];
            end
            %latent->observed that already has a parent (avoid observed collider)
            if strcmp(op{i},'add') && ismember(nodes(i,2),onodes) && sum(Gs{i}(:,nodes(i,2)))>1
                to_remove=[to_remove i];
            end
        end
        %for even iterations, only latent changes are allowed (to reduce complexity)
        %if mod(iter,2) && ismember(nodes(i,1),onodes) && ismember(nodes(i,2),onodes)
        %    to_remove=[to_remove i];
        %end
    end
    Gs(unique(to_remove))=[]; %delete the invalid neighbours
    
    %preserve 2-TBN
    for i=1:length(Gs)
        Gs{i}(nodes_per_slice+1:2*nodes_per_slice,nodes_per_slice+1:2*nodes_per_slice)=Gs{i};
        %preserve inter edges in necessary
        if sum(sum(inter_2TBN)), Gs{i}(1:nodes_per_slice,nodes_per_slice+1:end)=inter_2TBN; end
    end
    
    %score all neighbours
    flag=0; %indicator if we found better 2TBN
    for i=1:length(Gs)
        try
            intra=Gs{i}(1:nodes_per_slice,1:nodes_per_slice);
            inter=Gs{i}(1:nodes_per_slice,nodes_per_slice+1:end);
            bnet=mk_dbn(intra,inter,ns,'discrete',dnodes,'observed',onodes);
            for j=1:max(bnet.eclass2), bnet.CPD{j}=tabular_CPD(bnet,j); end %here we add random parameters
            engine=smoother_engine(jtree_2TBN_inf_engine(bnet));
            [~,LLtrace,engine]=learn_params_dbn_em(engine,data,'max_iter',max_iter); %here we update the parameters using EM
            %update the curr best found 2TBN
            if LLtrace(end)>bestscore
                flag=1;
                bestengine=engine;
                bestscore=LLtrace(end);
                curr_2TBN=Gs{i};
                warning(['better neighbor: ' int2str(i)])
                break %super hill climbing
            end
        catch
            warning(['error for neighbor: ' int2str(i)])
            continue
        end
    end
    
    %two stop conditions
    if (iter>=max_iter || flag==0), convergence=0; end 
    iter=iter+1 %increase number of iterations by one and print to screen
    warning(['num neighbors: ' int2str(length(Gs))]) %disp neighbors 
end

%move latent variables to the end
temp=curr_2TBN(:,[1:num_latents,(nodes_per_slice+1):(nodes_per_slice+num_latents)]);
curr_2TBN(:,[1:num_latents,(nodes_per_slice+1):(nodes_per_slice+num_latents)])=[];
curr_2TBN(:,end+1:end+num_latents*2)=temp;
temp=curr_2TBN([1:num_latents,(nodes_per_slice+1):(nodes_per_slice+num_latents)],:);
curr_2TBN([1:num_latents (nodes_per_slice+1):(nodes_per_slice+num_latents)],:)=[];
curr_2TBN(end+1:end+num_latents*2,:)=temp;