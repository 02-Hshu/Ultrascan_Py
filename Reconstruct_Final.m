% Anamik Jhunjhunwala & Dr. Ted Selker
% OD - Austen Poteet (5/26/19)
% LM - 01/23/2020
% Takes multiple bmp images in a radial orientation and converts them into
% X,Y and Z slices for 3D reconstruction and analysis using Slicer 4.2
% Linear Probe

clear all %Clears or all variables and assignments
close all %Closes all open images and graphs
clc %Cleans out the command window

disp('Start!')
% profile on %Used to optimize the code and check for computation capability used

crop = 0;
adjust = 1;


%% Importing the different variables from the scan

direc = 'C:\Users\W540\Desktop\Matlab_USG\Scan_3d_2020_01_27_16_15_21\Scan_3d__bmp';
% direc = 'Insert location of the images here - complete address';
    
    %Finding the correct amount to crop and adjust based on depth
    %ArrayImg = dir(fullfile(direc,'Image_*_*_*.bmp')) 
    %originalImage = imread('C:\Users\W540\Desktop\Matlab_USG\Image size\Scan_3d_2020_01_23_13_44_58\Scan_3d__bmp\Image_0_1_0.bmp');
    %[rows, columns, numberOfColorChannels] = size(originalImage);

% transfer the meta data and use that to import values
imgSpecFile = fopen(fullfile(direc,'Scan_3d__data.txt'));
A = textscan(imgSpecFile,'%s %s');
fclose(imgSpecFile);

%Preassigning all the variables for quicker computation
padSide = 0;
numSlices = 0; % Number of rotational slices
rotTot = 0;
probeID = 0; %256 refers to the 10000HZ while 512 refers to the 7500 HZ
imgW = 0; 
imgH = 0; 
thic = 3; % Thicknes of slice in pixels
interval = 1; % Interval between slices, 1 => use all slices
imgDepth = 0;
dia = 0; 
        %filtSize1=1;
        %filtSize2=2;

% Always adjust the number of slices analyzed
% Always adjust the rotation angle of the images

probeID = str2double(A{2}{16}); %Frequency of the transducer used
imgW = str2double(A{2}{17}); %.bmp image width
imgH = str2double(A{2}{18}); %.bmp image height
imgDepth = str2double(A{2}{2})/10;  %depth of scan in cm
rotTot = str2double(A{2}{20}); %Total rotation sweep in degrees
dia = 2*str2double(A{2}{21}); %diameter of circle in cm


%This is to crop the sides of the image in the conical shape for 7.5 Hz
% b is the value that corresponds to guardband and cone edges
% m is the value that corresponds to the diagonal
    gband = 0;
    
if probeID == 512
    if imgDepth == 4
        b = 0; 
        m = 0;
    elseif imgDepth == 8
        b = 70;  %b=70; m=0.5;
        m = 0.5;
    elseif imgDepth == 12  
        b = 0;  %b=70; m=0.5;
        m = 0;
    elseif imgDepth == 16
        b = 0;
        m = 0;
    end
elseif probeID == 256
    if imgDepth == 2
        b = 320; 
    elseif imgDepth == 4
        b = 200+gband;  
    elseif imgDepth == 6  
        b = 132.5+gband; 
    elseif imgDepth == 8
        b = 100+gband;
    elseif imgDepth == 10
        b = 80+gband;
    end
end

    
%Taking the first part of the image to crop the top area where it is convex
cropTop = round(imgH/(imgDepth*4)); %Is it actually taking the first quarter cm of the image?

%Bottom tangent line to the circle of imaging
%Will be used to crop the bottom of the image
cropBottom = imgH-round((imgDepth-dia/2)*imgH/imgDepth); 

%Pixels to pad array for rotation
%So that we do not lose any data while slicing
cropPadSide = round((dia-imgDepth)*imgH/imgDepth);

if (cropPadSide < 0)
   cropPadSide = abs(cropPadSide);
   padSide = 1;
   %Only in the case of 16cm imaging - need to pad sides to prevent loss
else
   padSide = 0;
end

%pixels to take off sides
%Based on the C# code in the software sdk we have settings of 4 cm interval
%We can adjust the C# code and get intervals of 2 cm as well
cropSide_d5 = 0; %d5 which is never used to crop side - d4 is used instead
switch (imgDepth)
    case 4
        imgW = 640;
        startPix=80;
    case 8
        %imgW = 640; %  imgW = 228; %434-207 
        %startPix=200;
    case 12
        imgW = 640;
        startPix=240;
    case 16
        imgW = 640;
        startPix=260;
end

%For our knowledge
if crop==1
   disp('Cropping Images: ON');
   %imgW=440;
else
   disp('Cropping Images: OFF');    
end

%For our knowledge
if adjust==1
   disp('Adjusting Images: ON');
else
   disp('Adjusting Images: OFF'); 
end

disp('Always adjust the number of slices analyzed and Always adjust the rotation angle of the images');
%% Load Images into matlab for slicing

%direc = 'Scan_3d_2019_06_02_14_03_59\Scan_3d__bmp\';
imgArray = dir(fullfile(direc,'Image_*_*_*.bmp'));

[~, reindex] = sort(str2double(regexp({imgArray.name}, '\d+', 'match', 'once')));
imgArray = imgArray(reindex);

numSlices = length(imgArray);
disp(strcat('Directory: ',direc)) %Output = location of images analyzed
disp(strcat('Number of slices: ',int2str(numSlices))) %Output = number of images
disp(strcat('Diameter: ',num2str(dia))) %Output = diameter of measurement
disp(strcat('Rotation Angle: ',num2str(rotTot))) %Output = angle of measurement made

rotInc = rotTot/numSlices; % degrees rotated between slices


%for i=1:1:numSlices
%    %filename = strcat('Image_',num2str(i-1),'_0_0.bmp'); % name of file to load Image_#.bmp
%    filename = nf(i).name;
%    %rawImg{i}=zeros(imgH,imgW); % fill array with zeros
%    rawImg{i} = imread(filename); % Load image
%end

%% Rotate Images & Construct 3D Matrix

%cropPadSide = (imgH/1); % pixels to pad for rotation around center

%Initialize array for summed slices
summed = uint8(zeros(imgH+cropPadSide,imgW,imgH+cropPadSide)); %Old code 
[L,W,H]=size(summed);



# Pick back up from here


% tic % output time for processing each slice
disp('Reconstructing...')
for i=1:interval:numSlices
    filename = strcat(direc,'\',imgArray(i).name);
    rawImg = (imread(filename)); % Load image
    
    %Size of the image before cropping - always 480x640 pixels
    SImg1 = size(rawImg);
    strcat(direc,filename);
    
    %Currently adjust is turned off 
    if adjust==1
%         rawImg = imguidedfilter(rawImg, 'NeighborhoodSize',[filtSize1 filtSize1], 'DegreeOfSmoothing', 0.01*diff([0 255]).^2);
%         rawImg = imadjust(rawImg,[0.05, 0.8],[0, 1]);
    end
    
   if crop==1
       if probeID == 512
            rawImg(1:cropTop,:)=0;
            rawImg(cropBottom:end,:)=0;
            %Here we essentially take a much smaller raw image after cropping
            %the top and the bottom components
            rawImg=rawImg(:,startPix:startPix+imgW-1);
       elseif probeID ==256
           %We do not crop the top or the bottom components at all
           %We only crop from the sides
       end
    end
    
    %Crop the image into a triangle/cone shape while cutting off the sides
    %here we use the topwidth and the bottom width and shave off the sides
    %of the image turning it into a polygon
    if probeID == 512 %7500Hz transducer
        for q=(1:imgH)
            cropTriangle = round(imgW/2-(m*q+b));
            rawImg(q,1:cropTriangle)=0;
            rawImg(q,imgW-cropTriangle:imgW)=0;
        end
    elseif probeID == 256 %10000Hz transducer
        for q = (1:imgH)
            cropSide = round(imgW/2-b);
            rawImg(q,1:cropSide)=0;
            rawImg(q,imgW-cropSide:imgW)=0;
        end
    end
    
    %filtering using an averaging method to make brightness -1 for all bad areas
    rawImg(1,:)=-1;
    for k=1:size(rawImg')
        for m=1:(size(rawImg)-6)
            if (rawImg(m+6,k)/7+rawImg(m+5,k)/7+rawImg(m+4,k)/7+rawImg(m+3,k)/7+rawImg(m+2,k)/7+rawImg(m+1,k)/7+rawImg(m,k)/7)>90
                break;
            end
            rawImg(m,k)=-1;
        end
    end 
    
    %size of the raw image after the cropping
    SImg2 = size(rawImg);
        %     rawImg1=rawImg*2^7;
        %     legMask=im2bw(rawImg1,0.04*0.5+0.5);
        %     skinMask=im2bw(rawImg1,0.74*0.5+0.5);
        %     tendonMask=im2bw(rawImg1,0.6*0.5+0.5)-skinMask;
        %     boneMask=-im2bw(rawImg1,0.2*0.5+0.5)+legMask;
        %     totalImg=50*boneMask+100*tendonMask+200*skinMask;

    %We only PadSide for 16cm arrays and greater    
    temp = uint8(zeros(imgH+cropPadSide,imgW,(2*thic)+1)); % initialize 3d array to contain slice
    % temp(:,:,:)=-1;
    % temp(:,:,:)=NaN;
    
    for j = -thic:thic
        if padSide==0
        temp(:,:,thic+1+j)=vertcat(rawImg,zeros(cropPadSide,imgW)); % Pad slice and create 3d slice with thickness
        else
        temp(:,:,thic+1+j)=vertcat(zeros(cropPadSide,imgW),rawImg); % Pad slice and create 3d slice with thickness
        end
    end
    
    %figure
    %imagesc(rawImg{i});
    %daspect([1 1 1]);
    %figure

    theta = rotInc*(i-1); % Calculate rotation of slice
    temp = permute(imrotate(permute(temp,[1 3 2]),theta,'bilinear','loose'),[1 3 2]); 
    % Permute 3d slice, rotate about vertical axis, permute back
    
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
    summed((1:L2)+off1,:,(1:H2)+off2)= max(summed((1:L2)+off1,:,(1:H2)+off2),temp);
    
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
% %         temp2=summed((1:L2)+off1,:,(1:H2)+off2);
% %         zeroInd=(temp>temp2);
% %         nonZeroInd=~zeroInd;
% % 
% %         summed((1:L2)+off1,:,(1:H2)+off2)=temp2(nonZeroInd)+temp(zeroInd);
%     maskSummed((1:L2)+off1,:,(1:H2)+off2)=max(maskSummed((1:L2)+off1,:,(1:H2)+off2),maskTemp);
% %         summed((1:L2)+off1,:,(1:H2)+off2)=summed((1:L2)+off1,:,(1:H2)+off2)+temp;
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

% This exports the matlab slices as a VTK file which is technically what we
% would use for paraview - an open source software developed by kitview for
% Python.
nanInd = isnan(summed);
summed2 = single(summed); %32 bit vs 64 bit - less memory used
summed2(nanInd) = -1000;
oldDir = cd;
cd(direc);
addpath(oldDir);

% summed_small=imresize3(summed,0.5);

% Mat2VTK('Construct.vtk',summed,'binary');
% Mat2VTK('Construct_small.vtk',summed_small,'binary');

% Mat2VTK('Construct_mask.vtk',maskSummed,'binary');

cd(oldDir);
convert=max(max(max(summed)))/255;
%% Creating X,Y and Z folders and putting the sliced images into those

disp('Outputting New Slices...')
mkdir(strcat(direc,'NewSlicesZ'));
for j=1:W
filename2 = strcat(direc,'NewSlicesZ\','Image_',num2str(j,'%03d'),'.bmp');
imwrite(permute(summed(:,j,:),[1 3 2]),filename2);
end

mkdir(strcat(direc,'NewSlicesX'));
for j=1:H
filename2 = strcat(direc,'NewSlicesX\','Image_',num2str(j,'%03d'),'.bmp');
imwrite(flipud(permute(summed(:,:,j),[2 1 3])),filename2);
end

mkdir(strcat(direc,'NewSlicesY'));
for j=1:L
filename2 = strcat(direc,'NewSlicesY\','Image_',num2str(j,'%03d'),'.bmp');
imwrite(flipud(permute(summed(j,:,:),[2 3 1])),filename2);
end

disp('Done!')
% profile viewer



