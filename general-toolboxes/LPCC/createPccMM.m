function [PCCS,IJ] = createPccMM(C,M,pcc_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create pcc between two center sets (major major): 
% try_big centers (C) and the original major clusters (M)
%
% input:
% [C]        - (matrix) cluster centers set A
% [M]        - (matrix) cluster centers set B
% [pcc_type] - [scalar] indicator for how to run pairwise cluster comparison
%
% output:
% [PCCS]     - (matrix) pcc between each row of C and M
% [IJ]       - (scalar) the comparsion indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initializtion
t=0;
PCCS=[];
[nrowC,ncolC]=size(C);
[nrowM,~]=size(M);

for i=1:nrowC
    for j=1:nrowM
        t=t+1;
        IJ(t,1)=i;
        IJ(t,2)=j;
        [PCCS(t,:),~]=createPcc([C(i,:);M(j,:)],pcc_type);
        %for k=1:ncolC
        %    if C(i,k)~=M(j,k)
        %        PCCS(t,k)=1;
        %    else
        %        PCCS(t,k)=0;
        %    end    
        %end
    end
end    