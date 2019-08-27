function [I1,dy,dx,theta] = PrepareOneOffsetImage(I0,sigma)
%% Noise
I0 = imnoise(I0/255,'gaussian',0,sigma^2)*255;
%% Offset
dy = -20+40*rand(1); % range [-20,20];
dx = -20+40*rand(1); % range [-20,20];

%% Rotation
theta = -15+30*rand(1); % range [-15,15];

%% Create a new image
original = I0;
outputViewo = imref2d(size(original));
T = [cosd(theta) sind(theta) 0;-sind(theta) cosd(theta) 0; 0 0 1]*[1 0 0;0 1 0;dx dy 1];
tr = projective2d(T);
distorted  = imwarp(original,tr,'cubic','OutputView',outputViewo);
I1 = distorted;

%%


end

