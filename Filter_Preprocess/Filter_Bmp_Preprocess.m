% Anamik Jhunjhunwala
% 10th January 2020
% Preprocessing the BMP Images to increase contrast and have better input
% data

%%
clear all
close all

disp('Start!')


direc = 'C:\Users\W540\Desktop\BMP_US_Images';
Img1 = (imread('C:\Users\W540\Desktop\BMP_US_Images\Sample_Ultrasound_1.bmp')); % Load image
Img2 = (imread('C:\Users\W540\Desktop\BMP_US_Images\Sample_Ultrasound_2.bmp')); % Load image

%%

% figure
% imshow(Img1)
% imshow(Img2)
imtool(Img1)
imtool(Img2)

%%
% X = multibandread(filename,size,precision,offset,interleave,byteorder)
figure
imhist(Img1(:,:,1))
title('Hist of the first image')

% n = 2;  
% Idouble = im2double(Img1); 
% avg = mean2(Idouble);
% sigma = std2(Idouble);
% proImg2 = imadjust(Img1,[avg-n*sigma avg+n*sigma],[]);

Img1_V = Img1(:);
Zeros_Img1 = find(~Img1_V);
for i = 1:length(Zeros_Img1)
    j = Zeros_Img1(i);
    Img1_V(j) = 128;
end

Img1_V_Gray = vec2mat(Img1_V,480);
Img1_V_Gray_true = Img1_V_Gray';
imshow(Img1_V_Gray_true)

proImg1 = imadjust(Img1_V_Gray_true);
proImg2 = imadjust(Img1);
imtool(proImg1)
imtool(proImg2)

imhist(proImg1(:,:,1))

%%
E = entropyfilt(I);
Eim = rescale(E);
figure
imshow(Eim)

BW1 = imbinarize(Eim, .8);
imshow(BW1);

figure
imshow(I)
BWao = bwareaopen(BW1,2000);
imshow(BWao)

nhood = true(9);
closeBWao = imclose(BWao,nhood);
imshow(closeBWao)

roughMask = imfill(closeBWao,'holes');
imshow(roughMask);

figure
imshow(I)

I2 = I;
I2(roughMask) = 0;
imshow(I2)

E2 = entropyfilt(I2);
E2im = rescale(E2);
imshow(E2im)
BW2 = imbinarize(E2im);
imshow(BW2)

mask2 = bwareaopen(BW2,1000);
imshow(mask2)

texture1 = I;
texture1(~mask2) = 0;
texture2 = I;
texture2(mask2) = 0;

imshow(texture1)
figure
imshow(texture2)

boundary = bwperim(mask2);
segmentResults = I;
segmentResults(boundary) = 255;
imshow(segmentResults)

%%

direc = 'C:\Users\W540\Desktop\Matlab_USG\Texture_Recognition';
I = imread('Tendon_Img1.bmp');

S = stdfilt(I,ones(21));
imshowpair(I,S,'montage')

%%

R = rangefilt(I,ones(21));
imshowpair(I,R,'montage')

%%


