clear all
clc
%----------------------------------------------------------
f=imread('Lena512.bmp');
f=uint8(f);
f=im2gray(f);
  W=double(imread('LOGO.bmp'));
% W=zeros(W,[32 32]);
%  W = xlsread('test3.csv');
W=rgb2gray(W);
%----------------------------------------------------------
f=imresize(f,[512 512],'bicubic');
f1=double(f);
figure(1);
imshow(mat2gray(f1));title('Original Host Image'); 
%-----------------------------------------------------------
W=round(imresize(W,[32 32],'bicubic'));
figure(2);
imshow(W);title('Original Watermark Image ');
%-----------------------------------------------------------
W1=W;m=10;
W=arnold( W,m);
figure(8);
imshow(W);title('scambled Watermark Image ');
%------------------dither creation----------------------------------------
%   deltaA=15;          deltaB=14;         deltaC=13;          deltaE=12;         deltaD=11;
  deltaA=7;          deltaB=6;         deltaC=5;          deltaE=4;         deltaD=3;
deltaA2=deltaA/2;   deltaB2=deltaB/2;  deltaC2=deltaC/2;  deltaE2=deltaE/2;  deltaD2=deltaD/2;
%------------------dither creation----------------------------------------
 DVecA=DitherCreate(deltaA,deltaA2,16);
% DVecA(1,:)=[4.7209 6.0869 -5.5952 6.2006 1.9854 -6.0369 -3.3225 0.7032 6.8626 6.9733 -5.1358 7.0589 6.8575 -0.2194 4.5042 -5.3717];
% DVecA(2,:)=[-2.7791 -1.4131 1.9048 -1.2994 -5.5146 1.4631 4.1775 -6.7958 -0.6374 -0.5267 2.3642 -0.4411 -0.6425 7.2806 -2.9958 2.1283];
%         %-------------
 DVecB=DitherCreate(deltaB,deltaB2,4);
% DVecB(1,:)=[-1.0953 5.8203 4.0909 6.4329];
% DVecB(2,:)=[5.9047 -1.1797 -2.9091 -0.5671];
%         %--------------
DVecC=DitherCreate(deltaC,deltaC2,1);
% DVecC(1,:)=[2.0246];
% DVecC(2,:)=[-4.4754];
% 
 DVecE=DitherCreate(deltaE,deltaE2,2);
% DVecE(1,:)=[-5.5715 4.1896];
% DVecE(2,:)=[0.4285 -1.8104];

% 
 DVecD=DitherCreate(deltaD,deltaD2,1);
% DVecD(1,:)=[4.7739];
% DVecD(2,:)=[-0.7261];
% %---------------------------------------------------------
DVecM=zeros(2,24);
DVecM(1,:)=[DVecD(1,:) DVecE(1,:)  DVecC(1,:) DVecB(1,:) DVecA(1,:)];
DVecM(2,:)=[DVecD(2,:) DVecE(2,:)  DVecC(2,:) DVecB(2,:) DVecA(2,:)];
Q=[deltaD deltaE deltaE deltaC deltaB deltaB deltaB deltaB deltaA deltaA deltaA deltaA deltaA deltaA deltaA deltaA deltaA deltaA deltaA deltaA deltaA deltaA deltaA deltaA];
%--------------------------------------------------------
% %save 'xyzdwt' DVecM;
% load xyzdwt;
%-------------------------------------------------------------------------
% w=liftwave('haar','Int2Int');
% [CA1,CH1,CV1,CD1] =lwt2(f1,w);
% [ca,chd,cvd,cdd] = swt2(X,2,'db6')
[CA1,CH1,CV1,CD1] =swt2(f1,1,'db1');
[CA2,CH2,CV2,CD2] =swt2(CA1,1,'db1');
[CA3,CH3,CV3,CD3] =swt2(CA2,1,'db1');
%---------------------------------
% Lm=zeros(32,32); L3H=zeros(32,32); H3L=zeros(32,32);D3m=zeros(32,32);
% D2m=zeros(64,64);
% D1m=zeros(128,128);

Lm=CA3; L3H=CH3; H3L=CV3;D3m=CD3;
D2m=CD2;
D1m=CD1;

Counter=1;

for i=1:32
  for j=1:32
      
     % if Counter<=576
      %-------------------------------------
      L=CA3(i,j); LH=CH3(i,j);HL=CV3(i,j);D3=CD3(i,j);
      D2=CD2(2*(i-1)+1:2*i,2*(j-1)+1:2*j);
      D1=CD1(4*(i-1)+1:4*i,4*(j-1)+1:4*j);
      %-------------------------------------
      d2=D2(:)';
      d1=D1(:)';
      Hd=[L LH HL D3 d2 d1];
      %---------------------
      Wbit=W(i,j); 
      dither0=DVecM(1,:);
      dither1=DVecM(2,:);
      %----------------------
      if (Wbit==0)
        s=quantize_vector0(Hd,dither0,Q,24,2);
      else
        s=quantize_vector1(Hd,dither1,Q,24,2);
      end
      Lm(i,j)=s(1); L3H(i,j)=s(2); H3L(i,j)=s(3);
      D3m(i,j)=s(4);
      %-----------------------------------
      a22=s(5:8);
      a22D=reshape(a22,2,2);
      
      a44=s(9:24);
      a44D=reshape(a44,4,4);
      %----------------------------------
      D2m(2*(i-1)+1:2*i,2*(j-1)+1:2*j)=a22D;
      D1m(4*(i-1)+1:4*i,4*(j-1)+1:4*j)=a44D;
      
     %Counter=Counter+1;
      
      %end
     
  end
end
%--------------------------------------------------
% w=liftwave('haar','Int2Int');
x1=iswt2(Lm,L3H,H3L,D3m,'db1');
x2=iswt2(x1,CH2,CV2,D2m,'db1');
x3=iswt2(x2,CH1,CV1,D1m,'db1');
iFf=x3;
%----------------------------------------------------
figure(3);
imshow(mat2gray(x3)),title(' Watermarked Image ');
%------------------------------------------------------
 %MSE=mean(mean((f1-iFf).^2))                 % Doubt about formulae 
% MSE=mse(f1,iFf);
% SNR=10*log10(255^2/MSE)% Doubt about formulae
%--------------- Calculate SSI  -------------------------------
% SSI=ssim_index(f1,iFf)
% snr=psnr(f1,iFf)%PSNR between original image and watermarked image%
% SSI=ssim(iFf,f1)%SSIM between original image and watermarked image%
% disp(SSI);
%----------------------=====================================---------------
 N1=3;
% %------------------------------------------
hostr=uint8(iFf);
% imwrite(hostr,'your.bmp');
% 
% % imwrite(hostr,'hostr100.jpg','quality',100);%//////JPEG
% hostr100=imread('your.bmp');
% figure(10);
% imshow(hostr100);
% hostr100=double(hostr100);

%----------------------------------------------
% hostr=uint8(iFf);
% imwrite(hostr,'hostr100.jpg','quality',70);%//////JPEG
% hostr100=imread('hostr100','jpg');
%------------------------------------------
hostr100= wiener2(hostr,[N1 N1]);
%------------------------------------------
%   hostr=uint8(iFf);
%   hostr100=medfilt2(hostr,[N1 N1]);
% image sharpen
% hostr100=imsharpen(hostr);
% image rotation
% hostr100=imrotate(hostr,45);
% hostr100=imresize(hostr100,[512 512],'bilinear');
%------------------------------------
%  H=fspecial('gaussian',[N1 N1]);   
% hostr100=imfilter(hostr,H);
% % %-----------------------------------------
% highpass=abs(iFf-hostr100);
% edgeh=highpass+iFf*(1.8-1);
% hostr100=edgeh;
% % % H=hpfilter(N1,N1,1.8) ;
% % % hostr100=imfilter(iFf,H);
%--------------------------------------------------
%  IbR=uint8(iFf);
%  hostr100=histeq(hostr);
%-------------------------------------------------
% host=imresize(iFf,.75,'bilinear');%//////SCALED DOWN
% hostr100=imresize(host,[256 256],'bilinear');
%-------------------------------------------------------------
% hostr100 = imtranslate(hostr,[25.3, -10.1],'FillValues',255,'OutputView','full');
% hostr100=imresize(hostr100,[512 512],'bicubic');
%  IbR=uint8(iFf);
%  hostr100=imnoise(hostr,'speckle',0.5); %'salt & pepper' 'speckle'   'gaussian' 
%   hostr100=IbR; %'salt & pepper' 'speckle'   'gaussian' 
% hostr100=imgaussfilt(hostr);%gaussian filter
%---------------------------------------
%------------DYNAMIC RANNGE CHANGE[50 200]-------------------
% Drc=50+(abs(1-iFf)*150)/254;
% hostr100=Drc;
%-------------------------------------------------------------
hostr100=double(hostr100);
iFf=hostr100;
% % %--------------------------------------------------------------------------
% MSE=mean(mean((f-iFf).^2))                 % Doubt about formulae 
% SNR=10*log10(255^2/MSE)% Doubt about formulae
% % %--------------- Calculate SSI  -------------------------------
% SSI1=ssim(f,iFf);
snr=psnr(f1,iFf)%PSNR between original image and attacked watermarked image%
SSI=ssim_index(iFf,f1)%SSIM between original image and attacked watermarked image%
nccV=ncc(f1,iFf)%Attacked Image%
% disp(SSI);
% figure(4);
% imshow(mat2gray(iFf));title('Attacked Image');
%--------------------categorization of block--at destination------------------
% w=liftwave('haar','Int2Int');
[CA1,CH1,CV1,CD1] =swt2(iFf,1,'db1');
[CA2,CH2,CV2,CD2] =swt2(CA1,1,'db1');
[CA3,CH3,CV3,CD3] =swt2(CA2,1,'db1');

%---------------------------------
% Lm=zeros(32,32); L3H=zeros(32,32); H3L=zeros(32,32);
% D3m=zeros(32,32);
% D2m=zeros(64,64);
% D1m=zeros(128,128);

Lm=CA3; L3H=CH3; H3L=CV3;D3m=CD3;
D2m=CD2;
D1m=CD1;


EwmB=zeros(32,32);

Counter=1;

for i=1:32
  for j=1:32
      
     % if Counter<=576
      %-------------------------------------
      L=CA3(i,j); LH=CH3(i,j);HL=CV3(i,j);
      D3=CD3(i,j);
      D2=CD2(2*(i-1)+1:2*i,2*(j-1)+1:2*j);
      D1=CD1(4*(i-1)+1:4*i,4*(j-1)+1:4*j);
      %-------------------------------------
      d2=D2(:)';
      d1=D1(:)';
      Hd=[L LH HL D3 d2 d1];
      %---------------------
      dither0=DVecM(1,:);
      dither1=DVecM(2,:);
      %----------------------
      sum0=mydistance0(Hd,dither0,Q,24);
      sum1=mydistance1(Hd,dither1,Q,24);
      if sum0<sum1
         s=Hd+dither0;EwmB(i,j)=0; 
      else
         s=Hd-dither1;EwmB(i,j)=1; 
      end    
      %----------------------------------
      Lm(i,j)=s(1); L3H(i,j)=s(2); H3L(i,j)=s(3);
      D3m(i,j)=s(4);
      %-----------------------------------
      a22=s(5:8);
      a22D=reshape(a22,2,2);
      
      a44=s(9:24);
      a44D=reshape(a44,4,4);
      %----------------------------------
      D2m(2*(i-1)+1:2*i,2*(j-1)+1:2*j)=a22D;
      D1m(4*(i-1)+1:4*i,4*(j-1)+1:4*j)=a44D;
     
      
     %Counter=Counter+1;
     %end
  end
end
%--------------------------------------------------
% w=liftwave('haar','Int2Int');
x1=iswt2(Lm,L3H,H3L,D3m,'db1');
x2=iswt2(x1,CH2,CV2,D2m,'db1');
x3=iswt2(x2,CH1,CV1,D1m,'db1');
iFf=x3;
% iFf = imresize(x3,1,'bicubic');
iFf = imresize(iFf,1,'bicubic');
%-----------------------------------------------
figure(5);
imshow(EwmB),title(' Scambled Extracted Watermark');
EwmB1=iarnold(EwmB,m);
figure(6);
imshow(EwmB1),title('  Extracted Watermark');
%dif=sum(sum((W-EwmB)));
%  nccV2=ncc(W1,EwmB1)%NCC value bewtween original watermark before scrambling and extracted watermark after inverse scrambling%
%  nccV=ncc(W,EwmB)%NCC value bewtween  watermark after scrambling and extracted watermark before inverse scrambling%
% nccV=normxcorr2(W,EwmB)
% nccV = corr2(f1,iFf)
% disp('NCC');
% disp(nccV);
%-------------------------------------------
figure(7);
imshow(mat2gray(iFf)),title(' Reconstructed ');
%--------------- Calculate PSNR  -------------------------------
% MSE=mean(mean((f1-iFf).^2))                 % Doubt about formulae 
%MSE=mse(iFf,f1)
% SNR=10*log10(255^2/MSE)% Doubt about formulae
  SNR1=psnr(iFf,f1) %PSNR between original image and reconstructed host image%
%--------------- Calculate SSI  -------------------------------
% SSI=ssim_index(f1,iFf)
% disp(SSI);
SSI1=ssim_index(iFf,f1) %SSIM between original image and reconstructed host image%
%--------------- Calculate PSNR  -------------------------------
err = immse(iFf,f1)
st0p=1;

        
        
        









     


