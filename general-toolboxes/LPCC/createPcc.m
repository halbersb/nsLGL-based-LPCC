function [CC,IJ] = createPcc(C,pcc_type,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the pairwise cluster compartion matrix
%
% input:
% [C]        - (matrix) major clusters centers
% [pcc_type] - (scalar) indicator: 1 for regulat pcc
%                                  2 for ratio pcc
%                                  3 for interval pcc
%                                  4 for z-test pcc
%                                  5 for KS-test pcc
%
% output:
% [CC]       - (matrix) pairwise cluster comparison (0/1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
t=0;
CC=[]; %'CC' stands for cluster comparison
IJ=[];
thers=0.3; %for case '3'

switch pcc_type,
    %original PCC
    case 1
        for i=1:length(C(:,1))
            for j=i+1:length(C(:,1))
                t=t+1;
                IJ(t,1)=i;
                IJ(t,2)=j;
                %loop over all observed
                for k=1:length(C(i,:))
                    if C(i,k)~=C(j,k) %if two clusters i & j are not equal with respect to observed k then 1
                        CC(t,k)=1;
                    else
                        CC(t,k)=0;
                    end
                end
            end
        end
    %ratio pcc
    case 2
        for i=1:length(C(:,1))
            for j=i+1:length(C(:,1))
                t=t+1;
                IJ(t,1)=i;
                IJ(t,2)=j;
                %loop over all observed
                for k=1:length(C(i,:))
                    if abs(C(i,k)-C(j,k))>=0.5
                        %first condition: if distance is greater than 0.5 then PCC is 1
                        CC(t,k)=1;
                        continue;
                    end
                    if C(i,k)>2 || C(j,k)>2
                        %second condition if one of them is greater than 2
                        if floor(C(i,k))~=floor(C(j,k))
                            %convert both to the same floor
                            C(i,k)=C(i,k)-0.5;
                            C(j,k)=C(j,k)-0.5;
                        end
                        %distributed between 1 and 2
                        C(i,k)=C(i,k)-floor(C(i,k))+1;
                        C(j,k)=C(j,k)-floor(C(j,k))+1;
                    end
                    %the following condition fits to the case where both numbers are between 1 and 2
                    if sqrt(max(C(i,k),C(j,k))/min(C(i,k),C(j,k))-1)>0.5 
                        %if two clusters i & j are not equal with respect to observed k then 1
                        if abs(C(i,k)-C(j,k))>=0.3, CC(t,k)=1; else CC(t,k)=0; end
                    else
                        if abs(C(i,k)-C(j,k))<0.3, CC(t,k)=0; else CC(t,k)=1; end
                    end
                end
            end
        end
    %interval pcc
    case 3
        for i=1:length(C(:,1))
            for j=i+1:length(C(:,1))
                t=t+1;
                IJ(t,1)=i;
                IJ(t,2)=j;
                %loop over all observed
                for k=1:length(C(i,:))
                    if abs(C(i,k)-C(j,k))>thers %if two clusters i & j are not equal with respect to observed k then 1
                        CC(t,k)=1;
                    else
                        CC(t,k)=0;
                    end
                end
            end
        end
    %statistical test - z (for binary variables)
    case 4
        data=varargin{1};
        clusters=varargin{2};
        major_clusters=varargin{3};
        for i=1:length(C(:,1))
            for j=i+1:length(C(:,1))
                t=t+1;
                IJ(t,1)=i;
                IJ(t,2)=j;
                %loop over all observed
                for k=1:length(C(i,:))
                    data1=data(find(clusters==major_clusters(i)),k);n1=length(data1);p1=length(find(data1==1))/n1;
                    data2=data(find(clusters==major_clusters(j)),k);n2=length(data2);p2=length(find(data2==1))/n2;
                    p=(n1*p1+n2*p2)/(n1+n2);
                    if abs(p2-p1)/sqrt(p*(1-p)*(1/n1+1/n2))>norminv(0.95);
                        CC(t,k)=1;
                    else
                        CC(t,k)=0;
                    end
                end
            end
        end
    %statistical test - KS (for ordinal variables)
    case 5
        for i=1:length(C(:,1))
            for j=i+1:length(C(:,1))
                t=t+1;
                IJ(t,1)=i;
                IJ(t,2)=j;
                %loop over all observed
                for k=1:length(C(i,:))        
                    data=varargin{1};
                    clusters=varargin{2};
                    major_clusters=varargin{3};
                    if kstest2(data(find(clusters==major_clusters(i)),k),data(find(clusters==major_clusters(j)),k));
                        CC(t,k)=1;
                    else
                        CC(t,k)=0;
                    end
                end
            end
        end
    otherwise, error('invalid pcc type: '); 
end    