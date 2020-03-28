function SOC = z_findSuspectedObservedColliders(Observed,Latent,emcc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find for all observed with no (latent) parents their possible (candidate) parents
%
% input:
% [latents]    - (vector)
% [observeds]  - (vector)
% [emcc]       - (matrix) extended major pairwise clusters comparison
%
% output:
% [SOC]        - (array) each element is observed with no parents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(Latent), SOC=[]; end %return if there are no latents

%% initialization
SOC=struct([]); %SOC stands for suspected observed colliders
to_remove=[];
s=1;

%first main phase - find SOC
for i=1:length(Observed)
    %no parents
    if Observed(i).parentAdded==false 
        SOC(s).index=i; %i is a potential observed collider
        SOC(s).type=false; %false for partial, true for full
        SOC(s).PPS=[]; %init potential parent set (PPS)
        SOC(s).UPPS=[]; %init unique potential parent set (unique PPS)
        for j=1:size(emcc,1)
            for k=1:length(Latent)
                %the observed changes value with the k latent in the j pcc and it doesn't already exist in PPS
                if emcc(j,i)==1 && emcc(j,i)==emcc(j,length(Observed)+k) && sum(ismember(SOC(s).PPS,length(Observed)+k))==0
                    %add the latent index to PPS
                    SOC(s).PPS=[SOC(s).PPS length(Observed)+k];
                end
            end
            if length(SOC(s).PPS)==length(Latent)
                %all latent are in PPS
                break;
            end
        end
        if length(SOC(s).PPS)<2
            %remove the observed from the potential set
            SOC(s)=[];
        else
            %increase the index of potential observed collidrs
            s=s+1; 
        end
    end
end

%second phase - clean SOC
for i=1:length(SOC)
    for j=1:size(emcc,1)
        %the potential observed collider changes alone, thus its not a collider
        if emcc(j,SOC(i).index)==1 && sum(emcc(j,SOC(i).PPS))==0
            %remove from SOC
            to_remove=[to_remove i];
            break;
        end
        for k=1:length(SOC(i).PPS)
            %the observed changes value with the k latent in the j pcc. also all other latents are zero (no change) and it doesn't already exist in UPPS
            if emcc(j,SOC(i).index)==1 && emcc(j,SOC(i).PPS(k))==1 && sum(emcc(j,SOC(i).PPS))==1 && sum(ismember(SOC(i).UPPS,length(Observed)+k))==0
                SOC(i).UPPS=[SOC(i).UPPS length(Observed)+k];
            end
        end
        if length(SOC(i).PPS)==length(SOC(i).UPPS)
            SOC(i).type=true;
            break;
        end
    end
end
SOC(to_remove)=[];