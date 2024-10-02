function p_val=psnr(img,img1)
[O,P] = size(img);
e=double(img)-double(img1);
mse=sum(e(:).^2)/(O*P);
RMS=sqrt(mse);
PSNR_val=20*log10(255/RMS);
p_val=PSNR_val;
