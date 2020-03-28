function [final_2TBN] = combine_local_struct(pdags,num_obs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combile local structures into a single 2TBN
% assumptions: 1. each observed variable has a single latent patents (MIM)
%              2. latent cannot be a parent of observed variables from both slices
%              3. all pdags are ordered: obs slc1, obs slc2, lat slc1, lat slc2
%
% input:
% [pdags]   - (cell array) each cell is a matrix represnting a single pdag
% [num_obs] - (scalar) number of observed variables in both slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stage 1: find observed groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbrs=zeros(num_obs,num_obs); %init neighbour counter matrix. S in the paper
for i=1:length(pdags)
    if ~isempty(pdags{i})
        for j=1:num_obs
            %here we assume that observed can have only one parent
            if length(parents(pdags{i},j))>1, error('observed can have only one parent'); end
            cld=children(pdags{i},parents(pdags{i},j));
            for k=1:length(cld)
                if(cld(k)<=num_obs && cld(k)~=j) %the child is not a latent
                    nbrs(j,cld(k))=nbrs(j,cld(k))+1;
                end
            end
        end
    end
end

%search for the most frequent children of the same latent parent
closest_nbrs=[]; %LP in the paper
for i=1:num_obs
    if max(nbrs(i,:))>0 %nonzero
        closest_nbrs{i}=find(nbrs(i,:)==max(nbrs(i,:)));
    else
        closest_nbrs{i}=[];
    end
end

%create groups of MSO based on "closest_nbrs"
groups=cell(1,0); %GP in the paper
closest=[];
for i=1:num_obs
    flag=0;
    if ~isempty(closest_nbrs), closest=closest_nbrs{i}; end
    for k=1:length(groups)
        if ismember(i,groups{k})
            groups{k}=[groups{k} closest];
            flag=1; %indicator that we added the observed to an existing group
            break;
        end
    end
    if flag==0
        %open a new group (i.e., latent with set of observed children)
        groups{end+1}=[i closest];
    end
end
%remove groups with single observed variable
groups=clean_groups(groups);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stage 2: validate that variables are in their slice (PMM assumption)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%seperate groups with children from both slices 1 and 2
new_group=cell(1,0);
latent_edges=cell(1,0); %we'll keep split as potential for stage 2
num_obs_per_slice=num_obs/2;
for i=1:length(groups)
    %the group has variables from both slices 1 and 2
    if max(groups{i})>num_obs_per_slice && min(groups{i})<=num_obs_per_slice
        to_remove=[];
        %if there are more in slice 2, then disconnect observed variables from slice 1
        if length(find(groups{i}>num_obs_per_slice))>length(find(groups{i}<=num_obs_per_slice))
            for j=1:length(groups{i})
                if groups{i}(j)<=num_obs_per_slice, to_remove=[to_remove j]; end
            end
        %otherwise, there are more in slice 1, then disconnect observed variables from slice 2
        else
            for j=1:length(groups{i})
                if groups{i}(j)>num_obs_per_slice, to_remove=[to_remove j]; end
            end
        end
        if length(to_remove)>1
            %split the group (create a new group)
            new_group{end+1}=groups{i}(to_remove);
            latent_edges{end+1}={groups{i}(~ismember(groups{i},groups{i}(to_remove))),groups{i}(to_remove)};
        end
        %remove from the original group
        groups{i}(to_remove)=[];
    end
end
groups=[groups new_group]; %add the split groups
groups=clean_groups(groups);
[groups1,groups2]=split_groups(groups,num_obs_per_slice); %groups1 for slice 1 and groups2 for slice 2
%try to merge
[groups1]=merge_groups(groups1,unique(cell2mat(groups2)),num_obs_per_slice);
[groups2]=merge_groups(groups2,unique(cell2mat(groups1)),num_obs_per_slice);
groups=[groups1 groups2];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % stage 2.1: find unconnected obvserved variables and try to connect them
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% unconnected=[];
% for i=1:num_obs
%     flag=0;
%     for j=1:length(groups)
%         if ismember(i,groups{j})
%             flag=1;
%         end
%     end
%     if ~flag, unconnected=[unconnected i]; end
% end
% %try to connect unconnected obvserved variables 
% for i=1:length(unconnected)
%     if unconnected(i)<=num_obs_per_slice
%         %unconnected variable is in slice 1
%         to_check=unconnected(i)+num_obs_per_slice;
%         for j=1:length(groups2) %search in slice 2
%             if ismember(to_check,groups2{j})
%                 %serach for the most similar group
%                 ind=find_most_similar_mso(groups2{j}-num_obs_per_slice,groups1);
%                 if ~isempty(ind)
%                     groups1{ind}=sort([groups1{ind} unconnected(i)]); %add to the relevant group
%                     break;
%                 end
%             end
%         end
%     else
%         to_check=unconnected(i)-num_obs_per_slice;
%         for j=1:length(groups1) %search in slice 1
%             if ismember(to_check,groups1{j})
%                 %serach for the most similar group
%                 ind=find_most_similar_mso(groups1{j}+num_obs_per_slice,groups2);
%                 if ~isempty(ind)
%                     groups2{ind}=sort([groups2{ind} unconnected(i)]); %add to the relevant group
%                     break;
%                 end
%             end
%         end
%     end
% end
% groups=[groups1 groups2];
groups=handle_observed_collider(groups,closest_nbrs,num_obs);
[groups1,groups2]=split_groups(groups,num_obs_per_slice); %split again to slice 1 and slice 2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stage 3: build the graph with latents based on groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_lats=length(groups); %find latent set size
final_2TBN=zeros(num_obs+num_lats,num_obs+num_lats);
ind=1;
for i=1:length(groups)
    if length(groups{i})>1
        for j=1:length(groups{i})
            final_2TBN(num_obs+ind,groups{i}(j))=1;
        end
        ind=ind+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stage 4: find latent connections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lats_mat=zeros(num_lats,num_lats); %init latents matrix. CO in the paper for KDD2019
for i=1:length(pdags)
    for j=num_obs+1:size(pdags{i},2)
        for k=num_obs+1:size(pdags{i},2)
            if pdags{i}(j,k)
                cld_from=find(pdags{i}(j,1:num_obs)); %observed children
                from=find_most_similar_mso(cld_from,groups);
                cld_to=find(pdags{i}(k,1:num_obs)); %observed children
                to=find_most_similar_mso(cld_to,groups);
                if ~isempty(from*to)
                    %oreintaion rule t => t+1
                    if from>length(groups1) && to<=length(groups2)
                        temp=from;
                        from=to;
                        to=temp;
                    end
                    if from~=to
                        lats_mat(from,to)=lats_mat(from,to)+1;
                    end
                end
            end
        end
    end
end
%for undirected edges take majority vote
for i=1:size(lats_mat,1)
    for j=i:size(lats_mat,1)
        if lats_mat(i,j)>0 && lats_mat(j,i)>0 && lats_mat(i,j)~=lats_mat(j,i)
            if lats_mat(i,j)>lats_mat(j,i)
                lats_mat(j,i)=0;
            else
                lats_mat(i,j)=0;
            end
        end
    end
end
lats_mat(find(lats_mat>1))=1;
%recover letent connction from observed splits
if ~isempty(latent_edges)
    for i=1:length(latent_edges)
        from=find_most_similar_mso(latent_edges{i}{1},groups);
        to=find_most_similar_mso(latent_edges{i}{2},groups);
        %if to<from, warning('"to" cannot be smaller than "from"'); end
        if from~=to && to>from , lats_mat(from,to)=1; end
    end
end
%oreintaion rule: intra structures must be equal
for i=1:length(groups1)
    for j=1:length(groups1)
        if lats_mat(i,j) %an intra edge
            from=find_most_similar_mso(groups1{i}+num_obs/2,groups2);
            to=find_most_similar_mso(groups1{j}+num_obs/2,groups2);
            if ~isempty(from) && ~isempty(to)
                lats_mat(length(groups1)+from(1),length(groups1)+to(1))=1;
            end
        end
    end
end
final_2TBN(num_obs+1:end,num_obs+1:end)=lats_mat;
%validated that the diagonal is zero
for i=1:size(final_2TBN,2)
    if final_2TBN(i,i)
        warning('A diagonal elenemt is not zero');
        final_2TBN(i,i)=0;
    end
end

%%
function groups = clean_groups(groups)
%remove groups with single observed variable

to_remove=[];
for i=1:length(groups)
    groups{i}=sort(unique(groups{i})); %remove duplications from the set
    if length(groups{i})<2, to_remove=[to_remove i]; end %remove group
end
groups(to_remove)=[];
%%
function index = find_most_similar_mso(set,groups)
%find index in groups (i.e., a latent variable) which is the most similar to set

cnt=0;
index=[];
for i=1:length(groups)
    if sum(ismember(set,groups{i}))>cnt
        cnt=sum(ismember(set,groups{i}));
        index=i;
    end
end
%%
function [groups1,groups2] = split_groups(groups,num_obs_per_slice)
% split the groups of MSO to groups from slice 1 and from slice 2

groups1=cell(1,0);groups2=cell(1,0);
for i=1:length(groups)
    if length(find(groups{i}>num_obs_per_slice))~=length(groups{i}) && length(find(groups{i}<=num_obs_per_slice))~=length(groups{i})
        error(['group ' int2str(i) ' has children from both slices']);
    end
    if groups{i}(1)<=num_obs_per_slice
        groups1{end+1}=groups{i};
    else
        groups2{end+1}=groups{i};
    end
end
%%
function [new_groups] = merge_groups(groups,other,num_obs_per_slice)
% merge two groups if they have more than 50% overlap in terms of shared observed

if length(groups)==1 
    new_groups=groups;
else
    for i=1:length(groups)-1
        cnt=0;with=0;
        for j=i+1:length(groups)
            if sum(ismember(groups{i},groups{j}))>cnt
                cnt=sum(ismember(groups{i},groups{j}));
                with=j;
            end
        end
        %do merge (for i and with)
        if with>0 && cnt>=0.5*min(length(groups{i}),length(groups{with})) %more than 50% overlap
            diff_var=[setdiff(groups{i},groups{with}) setdiff(groups{with},groups{i})]; %the non-common
            if min(cellfun(@min,groups))>num_obs_per_slice
                diff_var_other=diff_var-num_obs_per_slice; %to search in the second slice
            else
                diff_var_other=diff_var+num_obs_per_slice;
            end
            diff_var(ismember(diff_var_other,other))=[]; %do not add them to the merged set (they should be left unconnected)
            temp=unique([groups{i} groups{with}]);
            groups([i with])=[];
            groups=[groups temp(~ismember(temp,diff_var))];
            groups=merge_groups(groups,other,num_obs_per_slice); %recursive call
            break;
        end
    end
    new_groups=groups;
end
%%
function [groups] = handle_observed_collider(groups,closest_nbrs,num_obs)
% assign observed collider to the one most probable latent parent

for i=1:num_obs
    find_in_cell = cellfun(@(x) find(x==i),groups,'Uni',0);
    %observed collider was found
    if numel([find_in_cell{:}])>1 
        latent_parents=find(~cellfun(@isempty,find_in_cell));
        %choose one of the parents
        %try 1
        i_nbrs=closest_nbrs{i};
        for j=1:length(latent_parents)
            cnt(j)=sum(ismember(i_nbrs,groups{latent_parents(j)}));
        end
        %try 2 - change nbrs to the same sline
        if range(cnt)==0
            if i<=num_obs/2
                i_nbrs(i_nbrs>=num_obs/2)=i_nbrs-15;
            else
                i_nbrs(i_nbrs<num_obs/2)=i_nbrs+15;
            end
            for j=1:length(latent_parents)
                cnt(j)=sum(ismember(i_nbrs,groups{latent_parents(j)}));
            end
        end
        %try 3 - select random parent
        if range(cnt)==0
            for j=2:length(latent_parents) %j=1 is the random parent selected
                groups{latent_parents(j)}(groups{latent_parents(j)}==i)=[];
            end
        else
            to_leave=find(cnt==max(cnt));
            for j=1:length(latent_parents) %if length(to_leave)>1 then to_leave(1) is the random parent selected
                if j~=to_leave(1)
                    groups{latent_parents(j)}(groups{latent_parents(j)}==i)=[];
                end
            end
        end
    end
end
%%