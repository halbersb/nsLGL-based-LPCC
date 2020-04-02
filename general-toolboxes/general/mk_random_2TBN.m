function dag = mk_random_2TBN(N,num_obs,PMM_ind,flag,my_seed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a random 2TBN (DBN)
%
% input:
% [N]       - (scalar) number of variables per slice
% [num_obs] - (scalar) number of observed variables per slice
% [PMM_ind] - (scalar) indicator to init a pure structure
% [flag]    - [scalar] if 1 put latent at the end of the matrix
% [my_seed] - [scalar] my seed for randomness
%
% output:
% dag     - [matrix] 2TBN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_lat=N-num_obs;
lat_ind=1:num_lat;
if nargin<5, my_seed=1234; end
rng(my_seed,'twister');

if PMM_ind
    %this is a PMM
    slice_dag=zeros(N,N);
    for i=1:num_obs
        prt=round(1+rand*(num_lat-1),0);
        if rand>0.6 %prob of 0.6 for an ege
            slice_dag(prt,num_lat+i)=1;
        end
    end
else
    [slice_dag,~]=mk_rnd_dag(N,1);  %each node can have only one parent!

    %if there is an edge from observed to latent then reverse it
    for i=1:num_lat
        to=find(slice_dag(:,lat_ind(i)));
        if ~ismember(to,lat_ind) %'to' is not a latent
            %reverse the edge
            slice_dag(to,lat_ind(i))=0;
            slice_dag(lat_ind(i),to)=1;
        end
    end
end
%validate at least a single edge
if isempty(find(slice_dag)), slice_dag(1,num_lat+unidrnd(num_obs))=1; end

%convert to 2-TBN
dag=slice_dag;
dag(N+1:2*N,N+1:2*N)=slice_dag;

%add inter adges between latents only
if num_lat>1
    from=lat_ind(sample_discrete(normalise(ones(1,num_lat))));
    lat_ind=(N+1):(N+num_lat);
    to=lat_ind(sample_discrete(normalise(ones(1,num_lat))));
    if rand>0.9 %add inter adge
        dag(from,to)=1;
    end
end

%move latent variables to the end of the matrix
if flag && num_lat>0
    temp=dag(:,[1:num_lat (N+1):(N+num_lat)]);
    dag(:,[1:num_lat (N+1):(N+num_lat)])=[];
    dag(:,end+1:end+num_lat*2)=temp;
    temp=dag([1:num_lat (N+1):(N+num_lat)],:);
    dag([1:num_lat (N+1):(N+num_lat)],:)=[];
    dag(end+1:end+num_lat*2,:)=temp;
end