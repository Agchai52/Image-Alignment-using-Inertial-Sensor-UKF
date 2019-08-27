function [T] = MatrixExpandRGB(I)
% I 4*4*3
% I = zeros(4,4,3);
% I(:,:,1) = [11 12 13 14;21 22 23 24;31 32 33 34;41 42 43 44];

a = 0.25;

[m0,n0,s] = size(I);
[y0,x0] = meshgrid(1:m0,1:n0);
[y1,x1] = meshgrid(1:a:m0,1:a:n0);

I1 = double(I(:,:,1));
I2 = double(I(:,:,2));
I3 = double(I(:,:,3));

T(:,:,1) = interp2(y0,x0,I1',y1,x1)';
T(:,:,2) = interp2(y0,x0,I2',y1,x1)';
T(:,:,3) = interp2(y0,x0,I3',y1,x1)';

end

