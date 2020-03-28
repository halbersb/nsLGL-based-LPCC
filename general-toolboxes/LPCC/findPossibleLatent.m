function [Observed,Latent,OHLP] = findPossibleLatent(cc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find possible latent based on maximum set of observed (MSO)
% this icludes both latent colliders and latent exogenous
%
% input:
% [cc]       - (matrix) clusters comparison
%
% output:
% [Observed] - (structure array) observed variables
% [Latent]   - (structure array) latents variables
% [OHLP]     - (vector) observed that have a latent parent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize - read the observed variables from the data (clusters comparsion)
Observed=struct([]);
for k=1:length(cc(1,:))
    Observed(k).name=['O',num2str(k)];
    Observed(k).PCC_index=k;
    Observed(k).parentAdded=false;
    Observed(k).SLP=[]; %share latent parent (SLP) - vector with other observed that share the same latent
    Observed(k).CH=[];
    Observed(k).PA=[]; %the parent of the observed - only one perant is allowed (hence it is an integer)
    Observed(k).type='Observed';
    Observed(k).parentAdded2=false; %for the latent split
    Observed(k).SLP2=[]; %for the latent split
    Observed(k).CH2=[]; %for the latent split
    Observed(k).PA2=[]; %for the latent split
end

%find possible latent
OHLP=[];
ohlp_index=0;
latent_index=0;
gap=0;
Latent=struct([]);

%first phase - find maximum set of observed (MSO)
for i=1:length(Observed)
    slp_index=0; %at first no other observed is marked to have the same parent
    for j=1:length(Observed)
        if i==j
            continue;
        end
        flag=0;
        for t=1:length(cc(:,1)) %if the two observed are equal for all major clusters in CC then flag will remain 0
            if cc(t,i)~=cc(t,j)
                flag=1;
                break;
            end
        end
        if(flag==0) %if the two observed are equal for all major clusters in CC then they share the same latent
            slp_index=slp_index+1;
            Observed(i).SLP(slp_index)=j; %add to the SLP vector the current observed j (since they are equal for all major clusters in CC)
            found=0;
            for k=1:length(OHLP)
                if(OHLP(k)==i) %if i already exists in OHLP then don't add it (found=1)
                    found=1;
                    break;
                end
            end
            if(found==0)
                ohlp_index=ohlp_index+1;
                OHLP(ohlp_index)=i;
            end
        end
    end
end

%second phase - create the latents
if length(OHLP)~=0
    for k=1:length(OHLP(1,:)) %loop over all observed that have a latent parent
        if Observed(OHLP(k)).parentAdded==false %if the observed does not have latent yet then create new latent
            latent_index=latent_index+1;
            Latent(latent_index).name=['L',num2str(latent_index)]; %assign sequential name
            Latent(latent_index).PCC_index=length(Observed)+latent_index;
            Latent(latent_index).CH=[];
            Latent(latent_index).OCH=[];
            Latent(latent_index).CH_2S=[];
            Latent(latent_index).PA=[]; %parent set
            Latent(latent_index).K={}; %number of changes in PCCS_2S_tag - one for each path
            Latent(latent_index).W={};
            Latent(latent_index).D=[];
            Latent(latent_index).UD=[];
            Latent(latent_index).PATH={};
            Latent(latent_index).PPS={};
            Latent(latent_index).type='Latent';
            Latent(latent_index).CH(1)=Observed(OHLP(k)).PCC_index; %add the current observed as a child of the new latent
            Latent(latent_index).OCH(1)=Observed(OHLP(k)).PCC_index; %add the current observed as an observed child of the new latent
            Observed(OHLP(k)).parentAdded=true; %change the flag to true - the observed has a parent now
            Observed(OHLP(k)).PA(1)=Latent(latent_index).PCC_index; %PA(1) since only one parent is allowed
            %go over the SLP vector and add the other observed (that share the same latent) as the current new latent
            %also, add this new latent as a parent - PA(1) to the observed set in the SLP 
            for j=1:length(Observed(OHLP(k)).SLP)
                Latent(latent_index).CH(j+1)=Observed(OHLP(k)).SLP(j);
                Latent(latent_index).OCH(j+1)=Observed(OHLP(k)).SLP(j);
                Observed(Observed(OHLP(k)).SLP(j)).parentAdded=true;
                Observed(Observed(OHLP(k)).SLP(j)).PA(1)=Latent(latent_index).PCC_index;
            end
        end
    end
 end