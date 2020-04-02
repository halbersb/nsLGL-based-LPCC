function [allTBN,init_2TBN,bestscore,bestengine] = learn_struct_nsdbn_EM(init_2TBN,data,ns,L,C,max_iter,PMM_ind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% learn Dynamic Bayesian Network based on S&S - hill climbing procedure (non-stationary)
% assumptions: 1. number of latent variables is known (as well as their cardinalities)
%              2. observed variable cannot be a parent of latent
%              3. variables cannot be directed from 't+1' to 't'
%
% input:
% [init_2TBN]     - (matrix) initial 2TBN (could be a random graph)
% [data]          - (matrix) rows are samples and columns are variables
% [ns]            - (vector) the cardinality of each of the observed and latent variables
% [L]             - (scalar) number of latent variable in slice t
% [C]             - (scalar) number of changes in the DBN
% [max_iter]      - (scalar) maximum nuber of iteration for EM procedure
% [PMM_ind]       - (scalar) indicator to search over pure structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% initialization
    convergence=1;
    count=0;
    bestscore=-inf;
    bestengine=[];
    iter=1;
    N=size(init_2TBN,1)/2; %number of observed and latents per slice
    O=N-L;
    S=size(data,2)/O; %number of slices
    T=S/(C+1); %epoch size
    allTBN=cell(1); for i=1:S-1, allTBN{i}=init_2TBN; end
    %reshape data
    for i=size(data,2)+1:S*N, data(:,i)=nan; end
    data=num2cell(data');
    data(cellfun(@isnan,data))={''};  
        
    while convergence
        
        curr_allTBN=allTBN;
        
        for c=1:C+1  
            
            %find valid neighbours for the c's change point
            my_clock=clock;
            Gs{c}=find_valid_neighbors(curr_allTBN{(c-1)*T+1},O,L,PMM_ind,round(my_clock(end),0)*C*sum(sum(init_2TBN)));
                        
            %score all neighbours
            flag=0; %indicator if we found better 2TBN
            for i=1:length(Gs{c}) %for each neighbor
                try
                    for j=1:T, curr_allTBN{(c-1)*T+j}=Gs{c}{i}; end
                    curr_dag=flatten(curr_allTBN,S,N,O);
                    bnet=mk_bnet(curr_dag,ns);
                    engine = jtree_inf_engine(bnet);
                    for j=1:length(bnet.dnodes), bnet.CPD{j}=tabular_CPD(bnet,j); end %here we add random parameters
                    engine = jtree_inf_engine(bnet);
                    [~,LLtrace,engine]=learn_params_em(engine,data,max_iter); %here we update the parameters using EM
                    %update the curr best found 2TBN
                    count=count+1;
                    if LLtrace(end)>bestscore
                        flag=1;
                        bestengine=engine;
                        bestscore=LLtrace(end);
                        allTBN=curr_allTBN;
                        display(['better neighbor: ' int2str(i)])
                        break %super hill climbing
                    end
                catch
                    warning(['error for neighbor: ' int2str(i)])
                    continue
                end
            end
        end
        
        %two stop conditions
        if (iter>=max_iter || flag==0), convergence=0; end 
        iter=iter+1 %increase number of iterations by one and print to screen
        display(['total (comulative) num neighbors: ' int2str(count)]) %disp neighbors 
    end

end

%%
function [dag] = flatten(TBN,num_slices,num_var_slice,num_obs_slice)
% this function creates a flattened dag from 2-TBN
    
    %check final_2TBN size vs. num_var_slice
    for i=1:length(TBN)
        if size(TBN{i},1)<2*num_var_slice
            TBN{i}=padding(TBN{i},2*num_var_slice,'num_O',num_obs_slice); %padding
        end
    end
    %start falttening
    dag=zeros(num_slices*num_var_slice,num_slices*num_var_slice);
    num_lats_slice=num_var_slice-num_obs_slice;
    total_obs=num_obs_slice*num_slices;
    %add latent => observed edges
    for i=1:num_slices-1
        dag(total_obs+num_lats_slice*(i-1)+1:total_obs+num_lats_slice*i+num_lats_slice,(num_obs_slice*i-num_obs_slice+1):num_obs_slice*(i+1))=TBN{i}(2*num_obs_slice+1:end,1:2*num_obs_slice);
    end
    %add latent => latent edges
    for i=1:num_slices-1
        dag(total_obs+num_lats_slice*(i-1)+1:total_obs+num_lats_slice*i+num_lats_slice,total_obs+num_lats_slice*(i-1)+1:total_obs+num_lats_slice*(i-1)+2*num_lats_slice)=TBN{i}(2*num_obs_slice+1:end,2*num_obs_slice+1:end);
    end
end

%%
function [Gs] = find_valid_neighbors(graph,O,L,PMM_ind,my_seed)

% this function remove invalid neighbours

    O1=1:O;
    O2=O+1:2*O;
    L1=(2*O+1:2*O+L);
    L2=(2*O+L+1:2*(O+L));
    O=[O1 O2];
    L=[L1 L2];
    paths=reachability_graph(graph+transpose(graph));
    [Gs,op,nodes]=mk_nbrs_of_dag(graph);
    rng(my_seed,'twister');
    vec=randperm(length(Gs));
    Gs=Gs(vec);op=op(vec);nodes=nodes(vec,:);
    to_remove=[]; 
    for i=1:length(Gs)
        %don't try to add observed->latent
        if strcmp(op{i},'add') && ismember(nodes(i,1),O) && ismember(nodes(i,2),L)
            to_remove=[to_remove i];
        end
        %both are observed and none of them have a path to any latent
        if ismember(nodes(i,1),O) && ismember(nodes(i,2),O) && sum(sum(paths(L,nodes(i,:))))==0 && sum(sum(paths(nodes(i,:),L)))==0
            to_remove=[to_remove i];
        end
        %don't try to reverse latent->observed
        if strcmp(op{i},'rev') && ismember(nodes(i,1),L) && ismember(nodes(i,2),O)
            to_remove=[to_remove i];
        end
        %don't try to delete latent->observed if it has only 2 (or less) observed children
        if strcmp(op{i},'del') && length(parents(graph,nodes(i,2)))<3
            to_remove=[to_remove i];
        end
        %don't try to add latent->latent if one of the latent has no children
        if strcmp(op{i},'add') && isempty(children(graph,nodes(i,1))) && isempty(children(graph,nodes(i,2)))
            to_remove=[to_remove i];
        end
        %temporal: avoid latent->observed from different slice
        if strcmp(op{i},'add') && ((ismember(nodes(i,1),L1) && ismember(nodes(i,2),O2)) || (ismember(nodes(i,1),L2) && ismember(nodes(i,2),O1)))
            to_remove=[to_remove i];
        end
        if PMM_ind %serach over PMM space...
            %observed->latent/observed
            if strcmp(op{i},'add') && ismember(nodes(i,1),O)
                to_remove=[to_remove i];
            end
            %latent->observed that already has a parent (avoid observed collider)
            if strcmp(op{i},'add') && ismember(nodes(i,2),O) && sum(Gs{i}(:,nodes(i,2)))>1
                to_remove=[to_remove i];
            end
        end
    end
    Gs(unique(to_remove))=[];
end