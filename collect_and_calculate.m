function collect_and_calculate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function collects the SHD and F1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alg={'After_Weights','After_Sliding','After_LPCC','After_SEM_DBN','After_TLGL'};
path1='databases\';
path2='after_databases\';

for a=1:length(alg)
    tic
    display(['algorithm ' alg{a}]);
    folders=dir([path2 alg{a}]);folders(1:2)=[];
    for d=1:length(folders)
        files=dir([path2 alg{a} '\' folders(d).name]);files(1:2)=[];
        for f=1:length(files)
            load([path2 alg{a} '\' folders(d).name '\' files(f).name]);
            %collect parameters from file name
            G=str2double(files(f).name(strfind(files(f).name,'G')+1));
            C=str2double(files(f).name(strfind(files(f).name,'_C')+2));
            E=str2double(files(f).name(strfind(files(f).name,'_E')+2));
            O=str2double(files(f).name(strfind(files(f).name,'_O')+2));
            L=str2double(files(f).name(strfind(files(f).name,'_L')+2));
            if a==2
                load([path1 'G' int2str(G) '\' folders(d).name(1:end-6) '\' files(f).name(1:end-12) '_' files(f).name(end-4) '.mat']); %the original data for real_dag                   
            else
                load([path1 'G' int2str(G) '\' folders(d).name(1:end) '\' files(f).name]); %the original data for real_dag
            end
            
            %convert to 2TBNs (they are already DAGs by assumption)
            [real_TBN,real_drifts]=fulldag_to_2TBN(real_dag,O,L,E,C); %real_dag comes from the original data

            %complete first and last learned structures (padding of graphs) to be equal to length(real_TBN)
            if ~exist('TBN','var'), TBN{1}=final_2TBN; end %for stationary LGL-LPCC and SEM-DBN
            if length(real_TBN)-length(TBN)>2 && a~=3 && a~=4
                W=str2double(files(f).name(strfind(files(f).name,'_W')+2));
                copy_TBN=cell(1);
                for i=1:ceil(W/2), copy_TBN{i}=TBN{1}; end %pad the first graph to 1:window/2
                for i=2:length(TBN), copy_TBN{end+1}=TBN{i}; end %add all TBN
                for i=1:length(real_TBN)-length(copy_TBN), copy_TBN{end+1}=TBN{end}; end %padd the last graph to window/2
                TBN=copy_TBN;
            else
                W=2;
                for i=1:length(real_TBN)-length(TBN), TBN{end+1}=TBN{end}; end %pad last TBN
            end
            if length(TBN)>length(real_TBN), TBN(end-(length(TBN)-length(real_TBN))+1:end)=[]; end %remove from TBN if needed

            %process each TBN
            for i=1:length(TBN)
                TBN{i}=enforce_temporal(TBN{i},O);
                temp_TNB_i=pdag_to_dag(TBN{i});
                if isempty(temp_TNB_i), TBN{i}=fix_cycle_graph(TBN{i}); else TBN{i}=temp_TNB_i; end
                TBN{i}=dag_to_cpdag(TBN{i});
                %padding if necessary
                if size(TBN{i},2)<size(real_TBN{i},2), TBN{i}=padding(TBN{i},2*(O+L),'num_O',O); end
                if size(TBN{i},2)>size(real_TBN{i},2), real_TBN{i}=padding(real_TBN{i},size(TBN{i},2),'num_O',O); end
            end
            drifts=find_drifts(TBN,E,C);
            num_drifts_Recall(f) = measure_drifts(real_drifts,drifts,W,1);
            num_drifts_Precision(f) = measure_drifts(real_drifts,drifts,W,0);
            %here SHD starts
            for i=1:length(TBN)
                best_SHD=inf;
                %try all combinations for SHD score by replacing latent indeces and take the min
                comb=perms(O*2+1:(O*2+size(TBN{i},2)-2*O)); %L=size(TBN{i},2)-2*O!!
                pdag=TBN{i};
                for m=1:size(comb,1)
                    curr_pdag=pdag([1:O*2 comb(m,:)],[1:O*2 comb(m,:)]);
                    curr_pdag=enforce_temporal(curr_pdag,O); %direct undirected edges from t to t+1
                    [~,curr_SHD,~]=SHDL(dag_to_cpdag(real_TBN{i}),curr_pdag,O*2+1,'threshold',best_SHD);
                    if curr_SHD<best_SHD
                        best_SHD=curr_SHD;
                        All_curr_SHD(f,i)=curr_SHD;
                        if best_SHD==0, break; end
                    end
                end
            end
            SHD(f)=mean(All_curr_SHD(f,:));
            clear data LL BIC real_dag TBN data2 full_dag pdag pdags O L S P m comb curr_pdag real_drifts real_TBN drifts curr_SHD best_SHD
        end
        if C==3, C=5; end
        FinalTable((G-1)*8+C+E-5,a)=mean(SHD);
        FinalTable((G-1)*8+C+E-5,length(alg)+a)=2*mean(num_drifts_Precision)*mean(num_drifts_Recall)/(mean(num_drifts_Precision)+mean(num_drifts_Recall));
        clear All_curr_SHD SHD
    end
    toc
end
header=cell(1,2*length(alg));
for i=1:length(alg), header{i}=['SHD_' alg{i}]; end
for i=1:length(alg), header{length(alg)+i}=['F1_' alg{i}]; end
fig=figure('Position',[500,500,1050,350]);
uit=uitable(fig,'Data',FinalTable,'ColumnName',header);
uit.Position=[20,20,1010,310];
save('SHD_results.mat','FinalTable');