function [pdag_final,Latent_final] = z_connectLatents(pdag,Latent,order,data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if unconnected latents exist - connect them using PC algorithm 
%
% input:
% [pdag]       - (matrix) partial directed acyclic graph
% [Latent]     - (structure array) latents variables
% [order]      - (vector) the order of the latents per pcc_index
% [data]       - (cell) matrix data
%
% output:
% [pdag_final] - (matrix) new pdag after edge additions
% [Latent]     - (structure array) updated latents variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
pdag_final=pdag;
Latent_final=Latent;
count=0;
n_obs=size(pdag,1)-length(Latent);

%validation
L_pdag=pdag(size(pdag,1)-length(Latent)+1:end,size(pdag,1)-length(Latent)+1:end); %extract the sub plot of pdag
L_pdag=L_pdag(order-n_obs,order-n_obs);
for i=1:size(L_pdag,1)
    if (sum(L_pdag(:,i))+sum(L_pdag(i,:)))==0, count=count+1; end
end
if count==0, return; end

%learn pc algorithm for latents
data=cell2mat(data');
[nrow,ncol]=size(data);
try
    PC_pdag=learn_struct_pdag_pc('cond_indep_fisher_z',ncol,ncol-1,cov(data),nrow,0.01);
    PC_pdag=PC_pdag(order-n_obs,order-n_obs); %reorder by latents order so pcc_index could be used
catch
    warning('pc algorithm failed');
    return; %if pc failed
end

%find unconected latents in original pdag
for i=1:length(Latent)
    if isempty(Latent(i).PA) && isempty(setdiff(Latent(i).CH,Latent(i).OCH)) && isempty(Latent(i).UD) %unconnected latent: no parents, no latent children and no undirected edges
        %add parents if exists in PC_pdag
        prt=find(PC_pdag(:,i)==-1);
        for j=1:length(prt)
            Latent_final(i).PA=[Latent_final(i).PA Latent_final(prt(j)).PCC_index];
            Latent_final(prt(j)).CH=[Latent_final(prt(j)).CH Latent_final(i).PCC_index]; %add curr latent also as child of the parent
            pdag_final(Latent_final(prt(j)).PCC_index,Latent_final(i).PCC_index)=1; %add an edge to new pdag
        end
        %add children
        chl=find(PC_pdag(i,:)==-1);
        for j=1:length(chl)
            Latent_final(i).CH=[Latent_final(i).CH Latent_final(chl(j)).PCC_index];
            Latent_final(chl(j)).PA=[Latent_final(chl(j)).PA Latent_final(i).PCC_index]; %add curr latent also as parent of the child
            pdag_final(Latent_final(i).PCC_index,Latent_final(chl(j)).PCC_index)=1; %add an edge to new pdag
        end
        %add bi-directed edges
        bi=find(PC_pdag(:,i)==1); %'1' will appear symmetrically
        for j=1:length(bi)
            Latent_final(i).UD=[Latent_final(i).UD Latent_final(bi(j)).PCC_index];
            Latent_final(bi(j)).UD=[Latent_final(bi(j)).UD Latent_final(i).PCC_index];
            Latent_final(i).UD=unique(Latent_final(i).UD','rows')';
            Latent_final(bi(j)).UD=unique(Latent_final(bi(j)).UD','rows')';
            pdag_final(Latent_final(i).PCC_index,Latent_final(bi(j)).PCC_index)=2; %add a bi-directed edge to new pdag
            pdag_final(Latent_final(bi(j)).PCC_index,Latent_final(i).PCC_index)=2; %add symmetrically
        end
    end
end
%print
%names={};
%for i=1:size(PC_pdag,2), names=[names,['L' int2str(order(i))]]; end
%view(biograph(L_pdag,names));
%view(biograph(PC_pdag,names));
