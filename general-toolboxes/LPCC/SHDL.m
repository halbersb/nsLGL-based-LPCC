function [SHD_Score,SHD_Observed,SHD_Latent] = SHDL(Real_Dag,Found_Dag,First_latent,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate SHD scores with possible Latents
% we wish to minimize this score
% missing, revese or unnecessary edge would cause a high score
% the code assumes Real_Dag Found_Dag are pdag.
%
% input:
% [Real_Dag]     - (matrix) is the real graph to compare to
% [Found_Dag]    - (matrix) is the graph to score
% [First_latent] - if populated then measure SHD for latent model - divide SHD to L and O
%
% output:
% [SHD_Score]    - (scalar) the overall perdormance
% [SHD_Observed] - (scalar) the overall perdormance
% [SHD_Latent]   - (scalar) the overall perdormance
% if the model is not latent then SHD_Score=SHD_Observed and SHD_Latent is empty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read input
if size(Real_Dag,1)~=size(Found_Dag,1), 
    error('the two compared dags must be in the same size'); 
end
if nargin>2 %check if model is latent
    if ~isnumeric(First_latent) || First_latent>size(Real_Dag,2), error('latent is out of bound or is not an integer'); end
else
    First_latent=0; %indicator for non latent model
end
for i=1:2:length(varargin)
    switch varargin{i},
        case 'threshold', threshold=varargin{i+1};
        otherwise, error(['invalid argument name: ' varargin{i}]);       
    end
end

%initialization
SHD_Score=0;
SHD_Observed=0;
if First_latent==0 
    SHD_Latent=[]; 
else
    SHD_Latent=0;
end
[x1,y1]=find(Real_Dag);
[x2,y2]=find(Found_Dag);

%go over all edges in Real_Dag
for i=1:length(x1)
    %missing edge
    if Found_Dag(x1(i),y1(i))==0 && Found_Dag(y1(i),x1(i))==0
        if Real_Dag(x1(i),y1(i))*Real_Dag(y1(i),x1(i))==0, flag1=1; else flag1=0; end %the missing edge is directed in Real_Dag
        if Real_Dag(x1(i),y1(i))*Real_Dag(y1(i),x1(i))==1 && intersect(find(x1==y1(i)),find(y1==x1(i)))<i, flag2=1; else flag2=0; end %the missing edge is undirected in Real_Dag so ignore the second direction to avoid double counting
        if flag1 || flag2 
            SHD_Score=SHD_Score+1;
            %handle latent model (split SHD into latents and observed)
            if First_latent && x1(i)>=First_latent && y1(i)>=First_latent
                SHD_Latent=SHD_Latent+1;
            else
                SHD_Observed=SHD_Observed+1;
                if exist('threshold','var') && SHD_Observed>=threshold, break; end
            end
        end
    end
    %undirected edge in Found_Dag that is directed in Real_Dag
    if Found_Dag(x1(i),y1(i))==1 && Found_Dag(y1(i),x1(i))==1 && Real_Dag(x1(i),y1(i))*Real_Dag(y1(i),x1(i))==0 % directed in Real_Dag 1*0=0
        if First_latent && y1(i)>=First_latent %validate that we only fine latent-latent in these cases (for latent graphs. if not that First_latent=0...)
            SHD_Score=SHD_Score+1;
            %handle latent model (split SHD into latents and observed)
            if First_latent && x1(i)>=First_latent && y1(i)>=First_latent
                SHD_Latent=SHD_Latent+1;
            else
                SHD_Observed=SHD_Observed+1;
                if exist('threshold','var') && SHD_Observed>=threshold, 
                    break; 
                end
            end
        end
    end  
    %revese edge: directed edge in Found_Dag that is directed the opposite in Real_Dag
    if Found_Dag(x1(i),y1(i))==0 && Found_Dag(y1(i),x1(i))==1 && Real_Dag(y1(i),x1(i))==0
        SHD_Score=SHD_Score+1;
        %handle latent model (split SHD into latents and observed)
        if First_latent && x1(i)>=First_latent && y1(i)>=First_latent
            SHD_Latent=SHD_Latent+1;
        else
            SHD_Observed=SHD_Observed+1;
            if exist('threshold','var') && SHD_Observed>=threshold, break; end
        end
    end
end

%go over all edges in Found_Dag
for i=1:length(x2)
    if exist('threshold','var') && SHD_Observed>=threshold, break; end
    %unnecessary edge (redundant)
    if Real_Dag(x2(i),y2(i))==0 && Real_Dag(y2(i),x2(i))==0
        if Found_Dag(x2(i),y2(i))*Found_Dag(y2(i),x2(i))==0, flag1=1; else flag1=0; end %the missing edge is directed in Found_Dag
        if Found_Dag(x2(i),y2(i))*Found_Dag(y2(i),x2(i))==1 && intersect(find(x2==y2(i)),find(y2==x2(i)))<i, flag2=1; else flag2=0; end %the missing edge is undirected in Found_Dag so ignore the second direction to avoid double counting
        if flag1 || flag2
            SHD_Score=SHD_Score+1;
            %handle latent model (split SHD into latents and observed)
            if First_latent && x2(i)>=First_latent && y2(i)>=First_latent
                SHD_Latent=SHD_Latent+1;
            else
                SHD_Observed=SHD_Observed+1;
            end
        end
    end
end