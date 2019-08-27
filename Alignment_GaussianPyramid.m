%% Alignment by Gaussian Pyramid
% Handle Raw Image
disp('Alignment');
du = zeros(1,N);
dv = zeros(1,N);
[m,n]=size(I0);

for i = 1:N-1
    alternative = cell2mat(img_offset(i));
    
    [temp1,du(i),dv(i)] = AlignOneImage(I0,alternative);
    img_aligned(i) = {temp1};
end

umax = max(abs(du));
vmax = max(abs(dv));

I0 = I0(1:m-umax,1:n-vmax);
for i = 1:N-1
    temp2 = cell2mat(img_aligned(i));
    img_aligned(i) = {temp2(1:m-umax,1:n-vmax)};
end

disp('Alignment ends');

% Photo=demosaic(im2uint16(img_aligned1/255),'rggb');
% figure,imshow(Photo,[])