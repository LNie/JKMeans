function lambda = ULambda(X1,X2,init,numDif)

% Estimate a lambda value for joint k-means clustering with a specified number of variations 
% Writen by Lei NIE (nieleimail@gmail.com)
% 24 Nov. 2015

maxIter = 2;
lambda = zeros(maxIter,1);

[numSamples,numFeatures] = size(X1);
if (numSamples~=size(X2,1))||(numFeatures~=size(X2,2))
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
centroids1 = zeros(numClusters,numFeatures);
centroids2 = zeros(numClusters,numFeatures);
dists1 = zeros(numSamples,numClusters);
dists2 = zeros(numSamples,numClusters);
label1 = init;
label2 = init;


iter = 0;
while true
    iter = iter + 1;
    
    % Update centroids.
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
    
    % Compute distances of all samples to each centroid.
    for k = 1:length(changed1)
        i = changed1(k);
        dists1(:,i) = (X1(:,1) - centroids1(i,1)).^2;
        for j = 2:numFeatures
            dists1(:,i) = dists1(:,i) + (X1(:,j) - centroids1(i,j)).^2;
        end
    end
    for k = 1:length(changed2)
        i = changed2(k);
        dists2(:,i) = (X2(:,1) - centroids2(i,1)).^2;
        for j = 2:numFeatures
            dists2(:,i) = dists2(:,i) + (X2(:,j) - centroids2(i,j)).^2;
        end
    end
    
    
    % Assign samples to the clusters.
    preLabel1 = label1;
    preLabel2 = label2;
    [minDists1, label1] = min(dists1, [], 2);
    [minDists2, label2] = min(dists2, [], 2);
    disagreements = (label1~=label2);
    distsSum = dists1(disagreements,:) + dists2(disagreements,:);
    [minDistsSum,candidates] = min(distsSum, [], 2);
    
    difs = (minDistsSum-minDists1(disagreements)-minDists2(disagreements))/2;
    if numDif<=0
        lambda(iter) = max(difs);
    elseif numDif<=length(difs)
        difs = sort(difs,'descend');
        lambda(iter) = difs(numDif);
    else
        lambda(iter) = 0;
    end
    
    if iter>=maxIter
        return;
    end
    
    toAgree = (minDistsSum-minDists1(disagreements)-minDists2(disagreements)<=2*lambda(iter));
    label1(disagreements) = (1-toAgree).*label1(disagreements) + toAgree.*candidates;
    label2(disagreements) = (1-toAgree).*label2(disagreements) + toAgree.*candidates;
    
    moved1 = find(label1 ~= preLabel1);
    moved2 = find(label2 ~= preLabel2);
    changed1 = unique([label1(moved1); preLabel1(moved1)])';
    changed2 = unique([label2(moved2); preLabel2(moved2)])';
end
