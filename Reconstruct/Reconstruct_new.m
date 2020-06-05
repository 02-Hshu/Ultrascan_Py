% Austen Poteet
% 5/26/19
% Reconstruct US Image files into 3D
clear all
close all
clc

disp('Start!')
profile on

crop=1;
adjust=1;


%% Variables

direc = 'C:\Users\W540\Desktop\Matlab_USG\PlantarFasciitis_V2\Scan_3d__bmp';

padSide = 0;
numSlices = 512; % Number of rotational slices
rotTot=255; % Total rotation sweep in degrees
imgW = 640; % slice image width
imgH = 480; % slice image height
thic = 1; % Thicknes of slice in pixels
interval = 1; % Interval between slices, 1 => use all slices
depth=12; % depth of scan in cm
dia=16; % diameter of circle in cm
filtSize1=1;
filtSize2=1;


fileID=fopen(fullfile(direc,'Scan_3d__data.txt'));
A = textscan(fileID,'\');
fclose(fileID);

depth=str2num(A{2}{2})/10;
rotTot=str2num(A{2}{20});
dia=2*str2num(A{2}{21});

b=70;
m=0.5;

d1= round(imgH/(depth*2)); % first 0.5cm of image
d2 = imgH-round((depth-dia/2)*imgH/depth); % center line of circle
d3 = round((dia-depth)*imgH/depth); % pixels to pad array for rotation
if (d3<0)
   d3=abs(d3);
   padSide = 1;
else
   padSide = 0;
end

if crop==1
   disp('Cropping Images: ON');
   imgW=440;
else
   disp('Cropping Images: OFF');    
end
if adjust==1
   disp('Adjusting Images: ON');
else
   disp('Adjusting Images: OFF'); 
end
%% Load Images

%direc = 'Scan_3d_2019_06_02_14_03_59\Scan_3d__bmp\';
%direc = 'Scan_3d_2019_06_02_13_10_03\Scan_3d__bmp\'; % Brush Scan 1
%direc = 'Scan_3d_2019_06_02_11_49_44\Scan_3d__bmp_apple\'; %apple scan
%direc = 'Scan_3d_2019_06_02_13_38_02\Scan_3d__bmp\';
%direc = 'Jun17\Scan_3d_2019_06_17_15_17_18\Scan_3d__bmp\'; % Pot scan
%direc = 'Jun17\Scan_3d_2019_06_17_15_28_20\Scan_3d__bmp\'; % Pot scan
%direc = 'Jun18Scans\Scan_3d_2019_06_18_15_10_06\Scan_3d__bmp\'; % pot with wire
%direc = 'Jun18Scans\Scan_3d_2019_06_18_15_24_29\Scan_3d__bmp\';
%direc = 'Jun24Scans\Scan_3d_2019_06_24_15_14_15\Scan_3d__bmp\'; % Ankle
%direc = 'Jun24Scans\Scan_3d_2019_06_24_15_35_01\Scan_3d__bmp\'; % Elbow
%direc = 'Jun24Scans\Scan_3d_2019_06_24_15_29_05\Scan_3d__bmp\'; % 

%%
nf=dir(fullfile(direc,'Image_*_*_*.bmp'));

[~, reindex] = sort(str2double(regexp({nf.name}, '\d+', 'match', 'once')));
nf = nf(reindex);

numSlices = length(nf);
disp(strcat('Directory: ',direc))
disp(strcat('Number of slices: ',int2str(numSlices)))
disp(strcat('Diameter: ',num2str(dia)))
rotInc = rotTot/numSlices; % degrees rotated between slices

%for i=1:1:numSlices
%    %filename = strcat('Image_',num2str(i-1),'_0_0.bmp'); % name of file to load Image_#.bmp
%    filename = nf(i).name;
%    %rawImg{i}=zeros(imgH,imgW); % fill array with zeros
%    rawImg{i} = imread(filename); % Load image
%end

%% Rotate Images & Construct 3D Matrix

%d3=(imgH/1); % pixels to pad for rotation around center

summed = uint8(zeros(imgH+d3,imgW,imgH+d3)); % Initialize array for summed slices
% maskSummed = int16(zeros(imgH+d3,imgW,imgH+d3)); % Initialize array for summed slices
% summed(:,:,:)=NaN;

[L,W,H]=size(summed);
% tic % output time for processing each slice
disp('Reconstructing...')
for i=1:interval:numSlices
% parfor i=1:1:numSlices % Assemble slices in parallel !Careful of memory! max 2 workers on 16GB

    filename = strcat(direc,nf(i).name);
    rawImg = (imread(filename)); % Load image

    strcat(direc,filename);
    if adjust==1
        rawImg = imguidedfilter(rawImg, 'NeighborhoodSize',[filtSize1 filtSize1], 'DegreeOfSmoothing', 0.01*diff(getrangefromclass(rawImg)).^2);
%         rawImg = imadjust(rawImg,[0.05, 0.8],[0, 1]);
    end
    if crop==1
        rawImg(1:d1,:)=0;
        rawImg(d2:end,:)=0;
        rawImg=rawImg(:,106:545);
    end

    for q=(1:imgH)
        d4=imgW/2-(m*q+b);
        rawImg(q,1:d4)=0;
        rawImg(q,imgW-d4:imgW)=0;
    end
    
    rawImg(1,:)=-1;
    for k=1:size(rawImg')
        for m=1:(size(rawImg)-6)
            if (rawImg(m+6,k)/7+rawImg(m+5,k)/7+rawImg(m+4,k)/7+rawImg(m+3,k)/7+rawImg(m+2,k)/7+rawImg(m+1,k)/7+rawImg(m,k)/7)>90
                break;
            end
            rawImg(m,k)=-1;
        end
    end
    
%     rawImg1=rawImg*2^7;
%     legMask=im2bw(rawImg1,0.04*0.5+0.5);
%     skinMask=im2bw(rawImg1,0.74*0.5+0.5);
%     tendonMask=im2bw(rawImg1,0.6*0.5+0.5)-skinMask;
%     boneMask=-im2bw(rawImg1,0.2*0.5+0.5)+legMask;
%     totalImg=50*boneMask+100*tendonMask+200*skinMask;
    
    temp=uint8(zeros(imgH+d3,imgW,(2*thic)+1)); % initialize 3d array to contain slice
    %temp(:,:,:)=-1;
    %     temp(:,:,:)=NaN;
    for j=-thic:thic
        if padSide==0
        temp(:,:,thic+1+j)=vertcat(rawImg,zeros(d3,imgW)); % Pad slice and create 3d slice with thickness
        else
        temp(:,:,thic+1+j)=vertcat(zeros(d3,imgW),rawImg); % Pad slice and create 3d slice with thickness
        end
    end
    %figure
    %imagesc(rawImg{i});
    %daspect([1 1 1]);
    %figure

    theta=rotInc*(i-1); % Calculate rotation of slice
    temp = permute(imrotate(permute(temp,[1 3 2]),theta,'bilinear','loose'),[1 3 2]); % Permute 3d slice, rotate about vertical axis, permute back
    [L2,~,H2]=size(temp);
    if(L2>L)
        temp=temp(1:L,:,:);
        L2=L;
    end
    if(H2>H)
        temp=temp(:,:,1:H);
        H2=H;
    end
    off1=floor((L-L2)/2);
    off2=floor((H-H2)/2);
%     temp2=summed((1:L2)+off1,:,(1:H2)+off2);
%     zeroInd=(temp>temp2);
%     nonZeroInd=~zeroInd;
% 
%     summed((1:L2)+off1,:,(1:H2)+off2)=temp2(nonZeroInd)+temp(zeroInd);
    summed((1:L2)+off1,:,(1:H2)+off2)=max(summed((1:L2)+off1,:,(1:H2)+off2),temp);
    
%     maskTemp=int16(zeros(imgH+d3,imgW,(2*thic)+1))-1; % initialize 3d array to contain slice
%     %temp(:,:,:)=-1;
%     %     temp(:,:,:)=NaN;
%     for j=-thic/2:thic/2
%         if padSide==0
%         maskTemp(:,:,3+j)=vertcat(totalImg,zeros(d3,imgW)-1); % Pad slice and create 3d slice with thickness
%         else
%         maskTemp(:,:,3+j)=vertcat(zeros(d3,imgW)-1,totalImg); % Pad slice and create 3d slice with thickness
%         end
%     end
%     %figure
%     %imagesc(rawImg{i});
%     %daspect([1 1 1]);
%     %figure
% 
%     theta=rotInc*(i-1); % Calculate rotation of slice
%     maskTemp = permute(imrotate(permute(maskTemp,[1 3 2]),theta,'bilinear','loose'),[1 3 2]); % Permute 3d slice, rotate about vertical axis, permute back
%     [L2,~,H2]=size(maskTemp);
%     if(L2>L)
%         maskTemp=maskTemp(1:L,:,:);
%         L2=L;
%     end
%     if(H2>H)
%         maskTemp=maskTemp(:,:,1:H);
%         H2=H;
%     end
%     off1=floor((L-L2)/2);
%     off2=floor((H-H2)/2);
% %     temp2=summed((1:L2)+off1,:,(1:H2)+off2);
% %     zeroInd=(temp>temp2);
% %     nonZeroInd=~zeroInd;
% % 
% %     summed((1:L2)+off1,:,(1:H2)+off2)=temp2(nonZeroInd)+temp(zeroInd);
%     maskSummed((1:L2)+off1,:,(1:H2)+off2)=max(maskSummed((1:L2)+off1,:,(1:H2)+off2),maskTemp);
% %     summed((1:L2)+off1,:,(1:H2)+off2)=summed((1:L2)+off1,:,(1:H2)+off2)+temp;
% %     if mod(i,10)==0
% %         figure
% %         imagesc(rawImg);
% %         figure
% %         imagesc(totalImg);
% %     end
    disp(i) % output slice number
end
%%

% summed = imguidedfilter(summed, 'NeighborhoodSize', filtSize2); % Filter 3D Data
% Invert Image

% Invert = 255*uint8(ones(L,W,H))-summed;
% toc
%% Export New Slices
% disp('Outputting New Slices...')
% mkdir(strcat(direc,'NewSlices_Fast'));
nanInd=isnan(summed);
summed2=single(summed);
summed2(nanInd)=-1000;
oldDir=cd;
cd(direc);
addpath(oldDir);
summed_small=imresize3(summed,0.5);
Mat2VTK('Construct.vtk',summed,'binary');
Mat2VTK('Construct_small.vtk',summed_small,'binary');
% Mat2VTK('Construct_mask.vtk',maskSummed,'binary');
cd(oldDir);
convert=max(max(max(summed)))/255;


for j=1:H
filename2 = strcat(direc,'NewSlices_Fast\','Image_',num2str(j,'%03d'),'.bmp');
% imwrite(summed(:,:,j),filename2);
end

disp('Done!')

profile viewer


