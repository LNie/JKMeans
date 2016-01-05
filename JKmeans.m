function [label1,label2] = JKmeans(X1,X2,init,lambda)

% Joint k-means clustering
% Writen by Lei NIE (nieleimail@gmail.com)
% 23 Nov. 2015

% Input:
% X1 and X2 are matrices of fMRI data; each row is the time series of a vertex
% init is the initialization assignment; each element represents the label
% for a vertex; the label starts from 1 to N; N is the number of clusters
% lambda is the regularization parameter

% Output:
% label1 and label2 are the clustering assignments for X1 and X2 respectively

label1 = [];
label2 = [];
% The threshold for stability
esp = 0.01;
% Maximum number of iterations
maxIter = 300;

[numSamples,numFeatures1] = size(X1);
[numSamples2,numFeatures2] = size(X2);
if numSamples~=numSamples2
    disp('The dimensions of data does not match each other.');
    return;
end
if (size(init,1)~=numSamples)||(size(init,2)~=1)
    disp('The dimensions of initialization does not match the data.');
    return;
end
if min(init) ~= 1
    disp('The labels in initialization do not start from 1.');
    return;
end
numClusters = max(init);
for i = 1:numClusters
    if sum(init==i) == 0
        disp('There are empty clusters in initialization.');
        return;
    end
end
    
changed1 = 1:numClusters;
changed2 = 1:numClusters;
centroids1 = zeros(numClusters,numFeatures1);
centroids2 = zeros(numClusters,numFeatures2);
dists1 = zeros(numSamples,numClusters);
dists2 = zeros(numSamples,numClusters);
label1 = init;
label2 = init;

iter = 0;
while true
    iter = iter + 1;
    
%     Update centroids.
    for k = 1:length(changed1)
        i = changed1(k);
        members = (label1 == i);
        if any(members)
            centroids1(i,:) = sum(X1(members,:),1) / sum(members);
        end
    end
    for k = 1:length(changed2)
        i = changed2(k);
        members = (label2 == i);
        if any(members)
            centroids2(i,:) = sum(X2(members,:),1) / sum(members);
        end
    end
    
%     Compute distances of all samples to each centroid.
    for k = 1:length(changed1)
        i = changed1(k);
        dists1(:,i) = (X1(:,1) - centroids1(i,1)).^2;
        for j = 2:numFeatures1
            dists1(:,i) = dists1(:,i) + (X1(:,j) - centroids1(i,j)).^2;
        end
    end
    for k = 1:length(changed2)
        i = changed2(k);
        dists2(:,i) = (X2(:,1) - centroids2(i,1)).^2;
        for j = 2:numFeatures2
            dists2(:,i) = dists2(:,i) + (X2(:,j) - centroids2(i,j)).^2;
        end
    end
    
%     Terminal conditions.
    if iter > maxIter
        disp('Stop at the maximal iter.');
        return;
    elseif iter >= 2
        cost = sum(dists1((label1-1)*numSamples+(1:numSamples)'));
        cost = cost + sum(dists2((label2-1)*numSamples+(1:numSamples)'));
        if lambda~=inf
            cost = cost + 2*lambda*sum(label1~=label2);
        end
        reduced = preCost-cost;
        if reduced < 0
            label1 = preLabel1;
            label2 = preLabel2;
            return;
        elseif reduced < esp
            return;
        end
    elseif iter == 1
        cost = sum(dists1((label1-1)*numSamples+(1:numSamples)'));
        cost = cost + sum(dists2((label2-1)*numSamples+(1:numSamples)'));
    end
    
%     Assign samples to the clusters.
    preCost = cost;
    preLabel1 = label1;
    preLabel2 = label2;
    [minDists1, label1] = min(dists1, [], 2);
    [minDists2, label2] = min(dists2, [], 2);
    disagreements = (label1~=label2);
    distsSum = dists1(disagreements,:) + dists2(disagreements,:);
    [minDistsSum,candidates] = min(distsSum, [], 2);
    if lambda~=inf
        toAgree = (minDistsSum-minDists1(disagreements)-minDists2(disagreements)<=2*lambda);
        label1(disagreements) = (1-toAgree).*label1(disagreements) + toAgree.*candidates;
        label2(disagreements) = (1-toAgree).*label2(disagreements) + toAgree.*candidates;
    else
        label1(disagreements) = candidates;
        label2(disagreements) = candidates;
    end

    moved1 = find(label1 ~= preLabel1);
    moved2 = find(label2 ~= preLabel2);
    if (~isempty(moved1))&&(~isempty(moved2))
        changed1 = unique([label1(moved1); preLabel1(moved1)])';
        changed2 = unique([label2(moved2); preLabel2(moved2)])';
    else
        return;
    end
    
%     Display information
%     disp(['Iteration: ',num2str(iter),'; Cost: ',num2str(cost),'.' ]);
end