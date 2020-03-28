function [emcc] = extendPCCMatrix(mcc,Latent)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the extended cluster compartion (emcc) matrix
% added latent will have the value of its children for each cc
% note that the reason all its children were grouped under it
% was that they have the same value for each mcc
%
% input:
% [mcc]    - (matrix) major pairwise cluster comparison
% [Latent] - (structure array) latents variables
%
% output:
% [emcc]   - (matrix) the extended pcc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add columns to cc (one coulmn for each latent)
emcc=[mcc zeros(size(mcc,1),length(Latent))];

%loop to update the added columns
for i=1:length(Latent) %num column
    for j=1:length(emcc(:,1)) %num row
        flag=0;
        for k=1:length(Latent(i).CH)
            if emcc(j,Latent(i).CH(k))==0
                flag=1;
                break;
            end
        end
        if flag==0
            %emcc for latent column gets 1 only if all his children have 1 for that specific major cluater
            %that is, this column represent the CC for this group of observed
            emcc(j,Latent(i).PCC_index)=1;
        end
    end
end