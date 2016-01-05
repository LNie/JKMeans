function Data = Raw2Norm(RawFileName,Hemisphere,ReIndex)

% Normalize raw data
% Writen by Lei NIE (nieleimail@gmail.com)
% 08 Nov. 2015

%disp(RawFileName);
Data = ft_read_cifti(RawFileName);
Data = Data.dtseries(Data.brainstructure == Hemisphere,:);
Data(isnan(sum(Data,2)),:) = [];
[NumVert,~] = size(Data);
if ~isempty(ReIndex)
    Data = Data(:,ReIndex);
end
% Normalization
for i = 1:NumVert
    TmpSe = Data(i,:);
    TmpSe = TmpSe-mean(TmpSe);
    TmpSe_norm = norm(TmpSe,2);
    if TmpSe_norm ~= 0
        Data(i,:) = TmpSe/TmpSe_norm;
    else
        disp('Zero Norm');
        disp(i);
    end
end