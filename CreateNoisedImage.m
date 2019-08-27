clc
clear all
close all
%% Create Raw image

%img0 = imread('peppers.png');
img0 = imread('Peper_reference.png');

[m,n,s] = size(img0);

I0 = zeros(m,n);

for i = 1:m
    for j= 1:n
        if mod(i,2)==1 && mod(i,2)==1
            I0(i,j) = img0(i,j,2);
        elseif mod(i,2)==0 && mod(i,2)==0
            I0(i,j) = img0(i,j,2);
        elseif mod(i,2)==1 && mod(i,2)==0
            I0(i,j) = img0(i,j,1);
        elseif mod(i,2)==0 && mod(i,2)==1
            I0(i,j) = img0(i,j,3);
        end
            
    end
end
I00 = I0;
% figure,imshow(I0/255)
%% Create Image with Noise 
     mu = 0;
  sigma = 0.05;
 
I = I0/255;
I0 = imnoise(I,'gaussian',mu,sigma^2)*255;
I1 = imnoise(I,'gaussian',mu,sigma^2)*255;
I2 = imnoise(I,'gaussian',mu,sigma^2)*255;
I3 = imnoise(I,'gaussian',mu,sigma^2)*255;

% figure,imshow(I0)
save NoisedData I00 I0 I1 I2 I3