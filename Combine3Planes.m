function [ CombinedImage ] = Combine3Planes(R,G1,G2,B)
% R 4*4
% R = [11 12 13 14;
%      21 22 23 24;
%      31 32 33 34;
%      41 42 43 44];
% G1 = R;
% G2 = R;
% B =R;

CombinedImage = zeros(size(R,1)+size(G2,1),size(R,2)+size(G1,2));
CombinedImage(1:2:end,1:2:end) = R;
CombinedImage(1:2:end,2:2:end) = G1;
CombinedImage(2:2:end,1:2:end) = G2;
CombinedImage(2:2:end,2:2:end) = B;


end

