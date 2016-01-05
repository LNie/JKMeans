function Z = SCWard(X,connectivity)

% Spatially constrained Ward's method
% Writen by Lei NIE (nieleimail@gmail.com)
% 23 Nov. 2015

if islogical(connectivity)==false
    disp('The data type of connectivity is not logical.')
    return;
end
[numSamples,~] = size(X);
if (numSamples~=size(connectivity,1))||(numSamples~=size(connectivity,1))
    disp('The dimension of connectivity does not match the data.')
    return;
end
connectivity = tril(connectivity,-1);
numEdges = sum(connectivity(:));
Z = zeros(numSamples-1,3);
T = ones(numEdges,3);
W = ones(2*numSamples-1,1);
L = 1:numSamples;

[T(:,1),T(:,2)] = find(connectivity==1);
for i = 1:numEdges
    dif = X(T(i,1),:)-X(T(i,2),:);
    T(i,3) = dif*dif';
end

for step = 1:(numSamples-1)
    [~,pos] = min(T(:,3));
    merge = T(pos,:);
    T(pos,:) = [];
    Z(step,1:2) = sort(merge(1:2));
    Z(step,3) = sqrt(merge(3));
    numEdges = numEdges-1;
    V1 = (T(:,1)==merge(1))|(T(:,1)==merge(2));
    V2 = (T(:,2)==merge(1))|(T(:,2)==merge(2));
    N = [T(V2,2),T(V2,1),T(V2,3)];
    N = [N;T(V1,:)];
    [~,reorder] = sort(N(:,2));
    N = N(reorder,:);
    ID = step+numSamples;
    C = size(N,1);
    M = zeros(C,3);
    counter = 0;
    inc = 0;
    while counter<C
        counter = counter+1;
        inc = inc+1;
        Index1 = N(counter,1);
        w1 = W(Index1);
        D13 = N(counter,3);
        if Index1 == merge(2)
            Index2 = merge(1);
        else
            Index2 = merge(2);
        end
        Index3 = N(counter,2);
        w3 = W(Index3);
        w2 = W(Index2);
        if (counter<C)&&(N(counter+1,2)==N(counter,2))
            counter = counter+1;
            D23 = N(counter,3);
        else
            dif = mean(X(L==Index2,:),1)-mean(X(L==Index3,:),1);
            D23 = (2*(w2*w3)/(w2+w3))*(dif*dif');
        end
        D12 = merge(3);
        DD = ((w1+w3)*D13+(w2+w3)*D23-w3*D12)/(w1+w2+w3);
        M(inc,:) = [ID,Index3,DD];
    end
    M = M(1:inc,:);
    V = V1|V2;
    T(V,:) = [];
    T = [T;M];
    W(ID) = W(merge(1))+W(merge(2));
    L((L==merge(1))|(L==merge(2))) = ID;
end

