fast=1;flag=1;
cd general-toolboxes
cd bnt-master
addpath(genpathKPM(pwd))
cd ..
cd somtoolbox
addpath(pwd)
cd ..
cd LPCC
addpath(pwd) 
cd ..
cd general
addpath(pwd)
cd ..
cd ..

if fast
    for c=1:2:3
        for e=5:8
            tic
            wrapper_nsLGL_based_LPCC(1,c,e,5000,1,4,1,flag); %non-stationary sliding window
            times{1,1}(c,e-4)=toc;toc
            save('run_times','times');
            tic
            wrapper_nsLGL_based_LPCC(1,c,e,5000,2,2,1,flag); %non-stationary weights
            times{2,1}(c,e-4)=toc;toc
            save('run_times','times');
            tic
            wrapper_stationary_models(1,c,e,5000,1,flag); %stationary SEM_DBN
            times{3,1}(c,e-4)=toc;toc
            save('run_times','times');
            tic
            wrapper_stationary_models(1,c,e,5000,2,flag); %stationary LPCC
            times{4,1}(c,e-4)=toc;toc
            save('run_times','times');
        end
    end
end
collect_and_calculate
