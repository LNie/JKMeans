function ReIndex = BlockResample(NumTime,PLeng)

% Circular Block Resampling
% Writen by Lei NIE (nieleimail@gmail.com)
% 08 Nov. 2015

ReIndex = zeros(NumTime,1);
TCounter = 0;
while TCounter<NumTime
    Posi = randi(NumTime);
    Leng = geornd(PLeng);
    Leng = min(Leng,NumTime-TCounter);
    ReIndex((TCounter+1):(TCounter+Leng)) = Posi:(Posi+Leng-1);
    TCounter = TCounter + Leng;
end
ReIndex = mod(ReIndex,NumTime);
ReIndex(ReIndex==0) = NumTime;