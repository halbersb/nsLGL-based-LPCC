function [Observed,Latent_temp,OHLP2]=findPossibleLatentSplit3(cc,l,Observed,l_ch)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find possible latent
%
% input:
% [cc]          - (vector) pairwise cluster comparison
% [l]           - (scalar) current latent
% [Observed]    - (structure array) observed variables
% [l_ch]        - (vector) list of current latent children
%
% output:
% [Observed]    - (structure array) observed variables after update
% [Latent_temp] - (structure array) latents variables after update
% [OHLP2]       - (vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
OHLP2=[];
ohlp_index=0;
latent_index=0;
Latent_temp=struct([]);
for i=1:length(l_ch)
    Observed(l_ch(i)).parentAdded2=false;
    Observed(l_ch(i)).SLP2=[];
    Observed(l_ch(i)).CH=[];
    Observed(l_ch(i)).PA=[];
end

%first step - for each observed fill in its shared observed
for i=1:length(l_ch)
    slp_index=0; %in the beginning no other Observed identified to have the same parent
    for j=1:length(l_ch)
        if i==j
            continue;
        end
        flag=0;
        for t=1:length(cc(:,1))
          if cc(t,i)~=cc(t,j)
            flag=1;
            break
          end
        end
        if(flag==0)
            slp_index=slp_index+1;
            Observed(l_ch(i)).SLP2(slp_index)=l_ch(j);
            found=0;
            for k=1:length(OHLP2)
                if(OHLP2(k)==l_ch(i))
                    found=1;
                    break;
                end
            end
            if(found==0)
                ohlp_index=ohlp_index+1;
                OHLP2(ohlp_index)=l_ch(i);
            end
        end
    end
end

%second step - update observed set and create new latents
if ~isempty(OHLP2)
    for k=1:length(OHLP2(1,:))
        if Observed(OHLP2(k)).parentAdded2==false
            latent_index=latent_index+1;
            Latent_temp(latent_index).name=['L',num2str(l),num2str(latent_index)];
            Latent_temp(latent_index).CH=[];
            Latent_temp(latent_index).OCH=[];
            Latent_temp(latent_index).CH_2S=[];
            Latent_temp(latent_index).PA=[];
            Latent_temp(latent_index).K={}; %number of changes in PCCS_2S_tag - one for each path
            Latent_temp(latent_index).W={};
            Latent_temp(latent_index).D=[];
            Latent_temp(latent_index).UD=[];
            Latent_temp(latent_index).PATH={};
            Latent_temp(latent_index).PPS={};
            Latent_temp(latent_index).type='Latent';
            Latent_temp(latent_index).CH(1)=OHLP2(k);
            Latent_temp(latent_index).OCH(1)=OHLP2(k);
            Latent_temp(latent_index).PCC_index=length(Observed)+latent_index;
            Observed(OHLP2(k)).parentAdded2=true;
            Observed(OHLP2(k)).PA2(1)=Latent_temp(latent_index).PCC_index;
            for j=1:length (Observed(OHLP2(k)).SLP2)
                Latent_temp(latent_index).CH(j+1)=Observed(OHLP2(k)).SLP2(j);
                Latent_temp(latent_index).OCH(j+1)=Observed(OHLP2(k)).SLP2(j);
                Observed(Observed(OHLP2(k)).SLP2(j)).parentAdded2=true;
                Observed(Observed(OHLP2(k)).SLP2(j)).PA2(1)=Latent_temp(latent_index).PCC_index;
            end
        end
    end
 end