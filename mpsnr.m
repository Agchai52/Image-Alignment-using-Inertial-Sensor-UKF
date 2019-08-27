function PSNR = mpsnr(A,ref )
[m,n]=size(A);
Z = A-ref;
MSE = 1/m*1/n*sum(sum(Z.*Z));
PSNR = 10*log10(255^2/MSE);



end

