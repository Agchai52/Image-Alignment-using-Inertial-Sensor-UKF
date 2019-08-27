function [ R,G1,G2,B ] = Separate3Planes(I)
% I 4*4
% I = [11 12 13 14;
%      21 22 23 24;
%      31 32 33 34;
%      41 42 43 44];
 
temp11 = downsample(I,2);
temp1R = downsample(temp11',2);
temp1G = downsample(temp11',2,1);

temp22 = downsample(I,2,1);
temp2G = downsample(temp22',2);
temp2B = downsample(temp22',2,1);

R = temp1R';
G1 = temp1G';
G2 = temp2G';
B = temp2B';

end

