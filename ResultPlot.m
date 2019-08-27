

noiselevel = [0.02168687930014;
0.001075788862692;
0.001001009753240;
0.000968381917034;
0.001029331840389];

PSNR = [64.768832984910560;  
77.813533170146243;
78.126420518623846;
78.270336897691962;
78.005249538247512];

ImgNum = [1 5 10 15 18];

figure
plot(ImgNum,PSNR,'LineWidth',2);
xlabel('Merged Image Number')
title('PSNR')

figure
plot(ImgNum,noiselevel,'LineWidth',2);
xlabel('Merged Image Number')
title('Noise Level')

