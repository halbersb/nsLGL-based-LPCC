function [Latent] = findColliderLatents(Latent,emcc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find latent colliders in latent set
% a latent will be marked as collider if it satisfies the following:
% 1) it changes value (in emcc) together with another latent
%    and no other latent change value (="1") in that row
%    (each row refer to comparison between two major clustrs)
%    it is enough to find one row that satisfy this condition
% 2) this latent never change value (="1") alone, that is its the
%    only "1" in the row - this means that latent collider
%    will always change value in EX comb due to a change in one or
%    more of its parents
% 3) it doesn't create cycle - parent of parent is the latent
%
% input:
% [Latent] - (structure array) latents variables
% [emcc]   - (matrix) extended pairwise cluster comparison
%
% output:
% [Latent] - (structure array) new latents variables with latents PA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
num_Obs=size(emcc,2)-length(Latent); %number of observed variables
CPS={}; %CPS stands for copy set of 'Latent'
a=0; 
for i=1:length(Latent)
    a=a+1;
    CPS{a}=Latent(i);
end

%first phase - mark PPS (stands for possible parents set) for each latent
for i=1:length(Latent)
    f=0;
    for j=1:length(CPS)
        if CPS{j}.PCC_index~=Latent(i).PCC_index %if CPS is not the current latent
            flag1=0;
            flag3=0;
            for k=1:length(emcc(:,1))
                %if the two compared latents both have "1" (change value together) for the 'k'th pairwise cluster comparsion
                if emcc(k,Latent(i).PCC_index)==1 && emcc(k,CPS{j}.PCC_index)==1
                    flag1=1;
                    flag2=0;
                    for t=1:length(CPS)
                        %if there is another latent with "1" (i.e., change value together with i and j latents) 
                        %then do not consider the 'k'th emcc for adding PPS to latent i
                        if CPS{t}.PCC_index~=Latent(i).PCC_index && ...
                            CPS{t}.PCC_index~=CPS{j}.PCC_index && ...
                            emcc(k,CPS{t}.PCC_index)==1
                            flag2=1;
                            break;
                        end
                    end
                    %if for at least one of the 'k'th pairwise cluster comparsion, i and j latents equal "1" 
                    %and all other latents are "0" then it is enough to add j are PPS to i (flag3=1)
                    if flag1==1 && flag2==0
                        flag3=1;
                        break;
                    end
                end
             end
             if flag3==1
                f=f+1;
                Latent(i).PPS{f}=CPS{j};
            end
        end
    end
end

%second phase - set PA (stands for parents set) out of PPS for each latent (second condition)
for i=1:length(Latent)
    %check if there is more than one possible parent (otherwise it cannot be collider)
    if length(Latent(i).PPS)>1
        flag1=0;
        flag2=0;
        for k=1:length(emcc(:,1))
            %if the current latent has "1" in emcc for the 'k'th pairwise cluster comparsion
            if emcc(k,Latent(i).PCC_index)==1
                flag1=0;
                for f=1:length(Latent(i).PPS)
                    %if one of the PPS has also "1" for the 'k'th pairwise cluster comparsion
                    %then do not disqualified it as PA of the 'i'th latent (keep flag2=0)
                    if emcc(k,Latent(i).PPS{f}.PCC_index)==1
                        flag1=1;
                        break;
                    end
                end
                %if at least for one 'k'th in emcc, all the 'i'th latent's PPS are "0" (and the latent gets "1" for that 'k' in emcc)
                %then do not add PA (the 'i' latent is not a collider)
                if flag1==0
                    flag2=1;
                    break;
                end
            end
        end
        %if the latent was qualified to be a colider
        if flag2==0
            for f=1:length(Latent(i).PPS)
                if strcmp(Latent(i).PPS{f}.type,'Latent')
                    Latent(Latent(i).PPS{f}.PCC_index-num_Obs).CH(length(Latent(Latent(i).PPS{f}.PCC_index-num_Obs).CH)+1)=Latent(i).PCC_index;
                    Latent(i).PA(f)=Latent(i).PPS{f}.PCC_index; %add the parent index to PA set
                end          
            end
        end
    end
end

%third phase - clean PA set: remove unnecessary parents
for i=1:length(Latent) 
    to_remove=[];
    r=0;
    for f=1:length(Latent(i).PA)  
        ff=Latent(i).PA(f);
        %if the parent has parents (i.e., the parent itself is a colider)
        if length(Latent(ff-num_Obs).PA)>0 
            %for each of the parent parents
            for j=1:length(Latent(ff-num_Obs).PA)
                %if the parent of the parent is the latent itself
                if Latent((Latent(i).PA(f))-num_Obs).PA(j)==Latent(i).PCC_index
                    r=r+1;
                    to_remove(r,1)=f;
                    to_remove(r,2)=j;
                    break;
                end
            end
        end
    end
    for k=1:size(to_remove,1)
        h2=0;
        %remove from children set (CH)
        for h=1:length(Latent(i).CH)
            h2=h2+1;
            if Latent(i).CH(h2)==Latent(i).PA(to_remove(k,1))
                Latent(i).CH(h2)=[];
                h2=h2-1;
            end
        end
        m2=0;
        for m=1:length(Latent(Latent(i).PA(to_remove(k,1))-num_Obs).CH)
            m2=m2+1;
            if Latent(Latent(i).PA(to_remove(k,1))-num_Obs).CH(m2)==Latent(i).PCC_index
                Latent(Latent(i).PA(to_remove(k,1))-num_Obs).CH(m2)=[];
                m2=m2-1;
            end
        end
        %remove from parent set (PA)
        Latent(Latent(i).PA(to_remove(k,1))-num_Obs).PA(to_remove(k,2))=[];
        Latent(i).PA(to_remove(k,1))=[];
        %remove from remove list
        for l=1:size(to_remove,1)
            if to_remove(l,1)>to_remove(k,1)
                to_remove(l,1)=to_remove(l,1)-1;
            end
        end 
    end
end         