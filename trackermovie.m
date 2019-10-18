global FileNameImage ImageNumber FolderName
%trackerMovie.m by A. Smelser 12/12/12 for checking VST tracking accuracy
%MUST RUN GSE_Analysis.m on data prior to running this script; it requires
%M_sort variable in workspace.

%********************USER INPUTS*******************************************
loop            = 1; %0 for individual vesicle selection; 1 for loop over all vesicles

speed           = 10; %1 for slow (shows every frame); 10 for fast(shows every 10th 
%                      frame)
follow          = 0; %0 for stationary viewing field (ideal for random vesicle motion restricted to small area); 1 for 
%                      following vesicles (better for directed vesicle motion)
ROIsize         = 60;% size of viewing field. Use even number. Recommended = 50.

%**************************************************************************
clc;
close all;
% clear all;
pause on;


%Bring in images
fprintf('Select first image in image sequence for tracking\n');
pause(2);
[filename,FolderName] = uigetfile('*.tif','Select first image in image sequence');
% filename = '2peroxisomes1_000001.tif';
imagenum = char(regexp(filename,'\d{4}','match'));
ImageNumber = str2double(imagenum);
ImageNumber;
FileNameImage = char(regexp(filename,'\w+(?=\d{4}\.tif)','match'));
FileNameImage;

% find number of vesicles
numParticles = max(M_sort(:,3))+1; %number of tracked particles
if numParticles == 0
    fprintf('Check M_sort array');
end;

%*******************TRACK ALL PARTICLES SEQUENTIALLY:***********************
if loop == 1;
       
fprintf('**********    START LOOP over Particles         **************\n');
fprintf('*************************************************************\n');
fprintf('*************************************************************\n');
    

for i=1:numParticles % number of tracked particles
    close all
    clc
fprintf('**********    Inside i loop     *****************************\n');
fprintf('*************************************************************\n');
fprintf('*************************************************************\n');
i
% Pull in data for 1 bead
N = M_sort((M_sort(:,3) == i-1),:); %

A = N(:,4); % framenumber
B = N(:,3); % spotID
C = N(:,5); % x
D = N(:,6); % y
E = N(:,7); % r set in VST. We'll change col7 later.

C_offset =  0;  % adjusts center of ROI to center of bead; was -3
D_offset =  0;  %             "  was 4

C = C+C_offset; %             "
D = D+D_offset; %             "
frames = max(A) - min(A) +1;  % added 2012_10_11 gh

if speed == 1;
% loop over all frames to create tracker ROI movie
for k = ImageNumber:ImageNumber + frames -1  % *************************** k loop
    PICTNAME = sprintf('%s%s%04d.tif',FolderName,FileNameImage,k);
    img = imread(PICTNAME); % current image for loop
    cc = round(C(k-ImageNumber+1)); % center bead in box, col
    dd = round(D(k-ImageNumber+1)); % center bead in box, rows
   
  img(round(M_sort(i*frames - frames + k,6)),round(M_sort(i*frames - frames + k,5)+ 5):round(M_sort(i*frames - frames + k,5)+ 10))=[0 0 0 0 0 0];
  img(round(M_sort(i*frames - frames + k,6)),round(M_sort(i*frames - frames + k,5)- 10):round(M_sort(i*frames - frames + k,5)- 5))=[0 0 0 0 0 0];
  img(round(M_sort(i*frames - frames + k,6)+ 5):round(M_sort(i*frames - frames + k,6)+ 10),round(M_sort(i*frames - frames + k,5)))=[0 0 0 0 0 0];
  img(round(M_sort(i*frames - frames + k,6)- 10):round(M_sort(i*frames - frames + k,6)- 5),round(M_sort(i*frames - frames + k,5)))=[0 0 0 0 0 0];
    
    figure (2) % **********************************************  Fig 2
        if follow == 0;
             I_particle0 = img(round(D(1)-ROIsize/2):round(D(1)+ROIsize/2),...
                  round(C(1)-ROIsize/2):round(C(1)+ROIsize/2));
        elseif follow == 1;
             I_particle0 = img(round(D(k)-ROIsize/2):round(D(k)+ROIsize/2),...
                  round(C(k)-ROIsize/2):round(C(k)+ROIsize/2));
        end;
  
   imshow(I_particle0,[0 round(1.1*max(max(I_particle0)))],'InitialMagnification','fit')
       titletext2 = sprintf('Fig 2 ROI for ParticleID = %s, Frame %s', num2str(B(i)), num2str(k));
       title(titletext2);
       xlabel('x');
       ylabel('y');
       
    pause(0.05);
       
    count = 1;  % reset at start of each loop over j 
            
  count = count + 1;
  clearvars img I_particle0 cc dd PICTNAME
 end     % end of k loop
elseif speed == 10;
% loop over every 10th frame to create tracker ROI movie
for k = ImageNumber:10:ImageNumber + frames -1  % *************************** k loop
    PICTNAME = sprintf('%s%s%04d.tif',FolderName,FileNameImage,k);
    img = imread(PICTNAME); % current image for loop
        cc = round(C(k-ImageNumber+1)); % center bead in box, col
        dd = round(D(k-ImageNumber+1)); % center bead in box, rows
     
  img(round(M_sort(i*frames - frames + k,6)),round(M_sort(i*frames - frames + k,5)+ 5):round(M_sort(i*frames - frames + k,5)+ 10))=[0 0 0 0 0 0];
  img(round(M_sort(i*frames - frames + k,6)),round(M_sort(i*frames - frames + k,5)- 10):round(M_sort(i*frames - frames + k,5)- 5))=[0 0 0 0 0 0];
  img(round(M_sort(i*frames - frames + k,6)+ 5):round(M_sort(i*frames - frames + k,6)+ 10),round(M_sort(i*frames - frames + k,5)))=[0 0 0 0 0 0];
  img(round(M_sort(i*frames - frames + k,6)- 10):round(M_sort(i*frames - frames + k,6)- 5),round(M_sort(i*frames - frames + k,5)))=[0 0 0 0 0 0];
    
    figure (2) % **********************************************  Fig 2
    
        if follow == 0;
             I_particle0 = img(round(D(1)-ROIsize/2):round(D(1)+ROIsize/2),...
                  round(C(1)-ROIsize/2):round(C(1)+ROIsize/2));
        elseif follow == 1;
             I_particle0 = img(round(D(k)-ROIsize/2):round(D(k)+ROIsize/2),...
                  round(C(k)-ROIsize/2):round(C(k)+ROIsize/2));
        end;
   imshow(I_particle0,[0 round(1.1*max(max(I_particle0)))],'InitialMagnification','fit')
       titletext2 = sprintf('Fig 2 ROI for ParticleID = %s, Frame %s', num2str(B(i)), num2str(k));
       title(titletext2);
       xlabel('x');
       ylabel('y');
       
    pause(0.05);
       
    count = 1;  % reset at start of each loop over j 
            
  count = count + 1;
  clearvars img I_particle0 cc dd PICTNAME
end;     % end of k loop     
pause();
                    
end; % end of i loop
end;

%***************** TRACK SELECTED VESICLE ONLY:**************************
elseif loop == 0;
    
strResponse = input('Enter two-digit particle ID\n','s');
particleID    = str2num(strResponse)


% Pull in data for 1 particle
N = M_sort((M_sort(:,3) == particleID),:); %

A = N(:,4); % framenumber
B = N(:,3); % particleID
C = N(:,5); % x
D = N(:,6); % y
E = N(:,7); % r set in VST. We'll change col7 later.

C_offset =  0;  % adjusts center of ROI to center of bead; was -3
D_offset =  0;  %             "  was 4

C = C+C_offset; %             "
D = D+D_offset; %             "
frames = max(A) - min(A) +1;  % added 2012_10_11 gh
if speed == 1;
% loop over all frames to create tracker ROI movie
for k = ImageNumber:ImageNumber + frames -1  % *************************** k loop
    PICTNAME = sprintf('%s%s%04d.tif',FolderName,FileNameImage,k);
    img = imread(PICTNAME); % current image for loop
        cc = round(C(k-ImageNumber+1)); % center bead in box, col
        dd = round(D(k-ImageNumber+1)); % center bead in box, rows
     
  img(round(M_sort((particleID+1)*frames - frames + k,6)),round(M_sort((particleID+1)*frames - frames + k,5)+ 5):round(M_sort((particleID+1)*frames - frames + k,5)+ 10))=[0 0 0 0 0 0];
  img(round(M_sort((particleID+1)*frames - frames + k,6)),round(M_sort((particleID+1)*frames - frames + k,5)- 10):round(M_sort((particleID+1)*frames - frames + k,5)- 5))=[0 0 0 0 0 0];
  img(round(M_sort((particleID+1)*frames - frames + k,6)+ 5):round(M_sort((particleID+1)*frames - frames + k,6)+ 10),round(M_sort((particleID+1)*frames - frames + k,5)))=[0 0 0 0 0 0];
  img(round(M_sort((particleID+1)*frames - frames + k,6)- 10):round(M_sort((particleID+1)*frames - frames + k,6)- 5),round(M_sort((particleID+1)*frames - frames + k,5)))=[0 0 0 0 0 0];
    
    figure (2) % **********************************************  Fig 2
    
        if follow == 0;
             I_particle0 = img(round(D(1)-ROIsize/2):round(D(1)+ROIsize/2),...
                  round(C(1)-ROIsize/2):round(C(1)+ROIsize/2));
        elseif follow == 1;
             I_particle0 = img(round(D(k)-ROIsize/2):round(D(k)+ROIsize/2),...
                  round(C(k)-ROIsize/2):round(C(k)+ROIsize/2));
        end;
   imshow(I_particle0,[0 round(1.1*max(max(I_particle0)))],'InitialMagnification','fit')
       titletext2 = sprintf('Fig 2 ROI for ParticleID = %s, Frame %s', num2str(B(1)), num2str(k));
       title(titletext2);
       xlabel('x');
       ylabel('y');
       
    pause(0.05);
       
    count = 1;  % reset at start of each loop over j 
            
  count = count + 1;
  clearvars img I_particle0 cc dd PICTNAME
end;     % end of k loop

elseif speed == 10;
% loop over every 10th frame to create tracker ROI movie
for k = ImageNumber:10:ImageNumber + frames -1  % *************************** k loop
    PICTNAME = sprintf('%s%s%04d.tif',FolderName,FileNameImage,k);
    img = imread(PICTNAME); % current image for loop
        cc = round(C(k-ImageNumber+1)); % center bead in box, col
        dd = round(D(k-ImageNumber+1)); % center bead in box, rows
     
  img(round(M_sort((particleID+1)*frames - frames + k,6)),round(M_sort((particleID+1)*frames - frames + k,5)+ 5):round(M_sort((particleID+1)*frames - frames + k,5)+ 10))=[0 0 0 0 0 0];
  img(round(M_sort((particleID+1)*frames - frames + k,6)),round(M_sort((particleID+1)*frames - frames + k,5)- 10):round(M_sort((particleID+1)*frames - frames + k,5)- 5))=[0 0 0 0 0 0];
  img(round(M_sort((particleID+1)*frames - frames + k,6)+ 5):round(M_sort((particleID+1)*frames - frames + k,6)+ 10),round(M_sort((particleID+1)*frames - frames + k,5)))=[0 0 0 0 0 0];
  img(round(M_sort((particleID+1)*frames - frames + k,6)- 10):round(M_sort((particleID+1)*frames - frames + k,6)- 5),round(M_sort((particleID+1)*frames - frames + k,5)))=[0 0 0 0 0 0];
    
    figure (2) % **********************************************  Fig 2
    
        if follow == 0;
             I_particle0 = img(round(D(1)-ROIsize/2):round(D(1)+ROIsize/2),...
                  round(C(1)-ROIsize/2):round(C(1)+ROIsize/2));
        elseif follow == 1;
             I_particle0 = img(round(D(k)-ROIsize/2):round(D(k)+ROIsize/2),...
                  round(C(k)-ROIsize/2):round(C(k)+ROIsize/2));
        end;
   imshow(I_particle0,[0 round(1.1*max(max(I_particle0)))],'InitialMagnification','fit')
       titletext2 = sprintf('Fig 2 ROI for ParticleID = %s, Frame %s', num2str(B(1)), num2str(k));
       title(titletext2);
       xlabel('x');
       ylabel('y');
       
    pause(0.05);
       
    count = 1;  % reset at start of each loop over j 
            
  count = count + 1;
  clearvars img I_particle0 cc dd PICTNAME
end;     % end of k loop
end;
end
