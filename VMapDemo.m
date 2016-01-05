% This is a demo for generating a variation map using the joint K-means
% algorithm
% Writen by Lei NIE (nieleimail@gmail.com)
% 4 Jan. 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-specified parameters

% The Scan-1 data for Subject-1 
file11 = './HCP/100307/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii';
% The Scan-2 data for Subject-1
file12 = './HCP/100307/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.nii';
% The Scan-1 data for Subject-2 
file21 = './HCP/103414/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii';
% The Scan-2 data for Subject-2
file22 = './HCP/103414/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.nii';

% The number of parcels
numClusters = 150;

% Please set the following variable to 1 for the left hemisphere, 2 for the right hemisphere
hemisphere = 1;
% Load the spatial relationship matrix
if hemisphere == 1
    load('LeftNMap.mat');
elseif hemisphere == 2
    load('RightNMap.mat');
end

% The following variable controls whether estimate lambda using resampling
% It will talk a very long time for resampling
resample = false;
% The parameters for resampling.
if resample == true 
    % The number of resamples
    numResamples = 10;
    % The number of bootstrap variations 
    numVar = 297;
    % The success probability 
    pLeng = 1/(60+1);
else
    % A resampled value for subject 100307 and 103414 from the Human Connectome Project
    lambda = 0.0105;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lambda Estimation
if resample == true    
    % For Subject-1
    LTable1 = zeros(numResamples,1);
    X1 = [Raw2Norm(file11,hemisphere,[]),Raw2Norm(file12,hemisphere,[])];
    timePoints = size(X1,2)/2;
    parfor i = 1:numResamples
        disp(i);
        disp('Preprocessing');
        reIndex1 = BlockResample(timePoints,pLeng);
        reIndex2 = BlockResample(timePoints,pLeng);
        X2 = [Raw2Norm(file11,hemisphere,reIndex1),Raw2Norm(file12,hemisphere,reIndex2)];
        disp('Initilization Generation');
        Z = SCWard([X1,X2],NeiMap);
        init = cluster(Z,'maxclust',numClusters);
        [init,~] = JKmeans(X1,X2,init,inf);
        disp('Lambda Estimation');
        lambda = ULambda(X1,X2,init,numVar);
        LTable1(i) = max(lambda);
    end
    % For Subject-2
    LTable2 = zeros(numResamples,1);
    X1 = [Raw2Norm(file21,hemisphere,[]),Raw2Norm(file22,hemisphere,[])];
    timePoints = size(X1,2)/2;
    parfor i = 1:numResamples
        disp(i);
        disp('Preprocessing');
        reIndex1 = BlockResample(timePoints,pLeng);
        reIndex2 = BlockResample(timePoints,pLeng);
        X2 = [Raw2Norm(file21,hemisphere,reIndex1),Raw2Norm(file22,hemisphere,reIndex2)];
        disp('Initilization Generation');
        Z = SCWard([X1,X2],NeiMap);
        init = cluster(Z,'maxclust',numClusters);
        [init,~] = JKmeans(X1,X2,init,inf);
        disp('Lambda Estimation');
        lambda = ULambda(X1,X2,init,numVar);
        LTable2(i) = max(lambda);
    end
    lambda = max(prctile(LTable1,95),prctile(LTable2,95));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalizing raw data
disp('Preprocessing');
X1 = [Raw2Norm(file11,hemisphere,[]),Raw2Norm(file12,hemisphere,[])];
X2 = [Raw2Norm(file21,hemisphere,[]),Raw2Norm(file22,hemisphere,[])];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating the initialization assignment
disp('Initilization Generation');
Z = SCWard([X1,X2],NeiMap);
init = cluster(Z,'maxclust',numClusters);
[init,~] = JKmeans(X1,X2,init,inf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating functional parcellations
disp('Parcellations Generation');
[pa1,pa2] = JKmeans(X1,X2,init,lambda);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating a file for visualization in Connectome Workbench
RawData = ft_read_cifti('Template.dtseries.nii');
OutData = RawData;
AllVert = 96854;
numF = 3;
DTMatrix = nan(AllVert,numF);
hemisphere = 1;
VCounter = 0;
for i = 1:size(RawData.brainstructure,1)
    if ((RawData.brainstructure(i)==hemisphere)&&(~isnan(RawData.dtseries(i,1))))     
        VCounter = VCounter+1;
        DTMatrix(i,1) = pa1(VCounter);
        DTMatrix(i,2) = pa2(VCounter);
        DTMatrix(i,3) = (pa1(VCounter)~=pa2(VCounter));
    end
end
OutData.time = zeros(1,numF);
for k = 1:numF
     OutData.time(k) = k;
end
OutData.dtseries = DTMatrix;
outfilename = 'VMap.nii';
ft_write_cifti(outfilename,OutData,'parameter','dtseries');