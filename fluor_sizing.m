function [M_sort_r,numParticles] = fluor_sizing(M_sort_cut,pixeldistance,~,ROIsize,check)
global FolderName FileNameImage ImageNumber

%Measures the radius of fluorescent particles

%Nested function: myfun5

clc;
close all;


% Allocate, initialize M_sort_r
M_sort_r = M_sort_cut; 

% Generate 2D arrays to represent regions of interest around particles
xROI = linspace(round(-ROIsize/2),round(ROIsize/2),ROIsize+1);
yROI = xROI;
[XXROI,YYROI] = meshgrid(xROI,yROI); % xy space for ROI

%Bring in images, if haven't previously
if check == 0;
fprintf('Select first image in series');
pause(2);
[filename,FolderName] = uigetfile('*.tif','Select first image in series');
imagenum = char(regexp(filename,'\d{4}','match'));
ImageNumber = str2double(imagenum);
ImageNumber;
FileNameImage = char(regexp(filename,'\w+(?=\d{4}\.tif)','match'));
FileNameImage;
elseif check == 1;
fprintf('Sizing fluorescent particles in images\n')
pause(1);
end;



for m = ImageNumber:ImageNumber+max(M_sort_cut(:,4))      
    
%Add image numbers to 2nd column of M_sort_r
     M_sort_r((M_sort_r(:,4) == m - ImageNumber),2)= m;
end;
    
% Find number of particles
numParticles = max(M_sort_cut(:,3))+1 %number of tracked particles
if numParticles == 0
    fprintf('Check M_sort array');
end
      
for i=1:numParticles % number of tracked particles
    close all
    clc
fprintf('**********    Inside loop over Particles     ****************************\n');

i


%Designate columns within tracking data
N = M_sort_cut((M_sort_cut(:,3) == i-1),:);


A = N(:,4); %framenumber
B = N(:,3); %spotID
C = N(:,5); %x
D = N(:,6); %y
E = N(:,7); %r set in VST
C_offset =  2;  % move center of ROI to center of particle; was 3.5??? check this
D_offset = 0;  %             ";was -2.0
     
C = C+C_offset; %             "
D = D+D_offset; %             "
     


    
%Average the intensities of the first 20 frames

% For each particle, loop over all frames to average z at each pixel
% within the ROI that is centered on moving particle

    %Set z_sum at 0, set count at 1

        z_sum(ROIsize+1,ROIsize+1)=0;
        count = 1;  % reset at start of each loop over j 
                
 for j = ImageNumber:ImageNumber+19; 
        count
   
        PICTNAME = sprintf('%s%s%04d.tif',FolderName,FileNameImage,j);
        img = imread(PICTNAME); % first image in series
        
     if j == ImageNumber;    
        img1 = img; % set first image aside for later display
          
                      
                         
        figure (1) % **********************************************  Fig 1
            I_particle0 = img1(round(D(1)-ROIsize/2):round(D(1)+ROIsize/2),...
                           round(C(1)-ROIsize/2):round(C(1)+ROIsize/2));
            imshow(I_particle0,[0 round(1.1*max(max(I_particle0)))])
               titletext2 = sprintf('Fig 1 ROI for particleID = %s', num2str(B(1)));
               title(titletext2);
               xlabel('x');
               ylabel('y');
               pause(1);
   
     end % end of if j == ImageNumber
    cc = round(C(j-ImageNumber+1)); % x at center (columns) of box
    dd = round(D(j-ImageNumber+1)); % y at center (rows)    of box
  z_sum = z_sum + double(img(dd-ROIsize/2:dd+ROIsize/2,cc-ROIsize/2:cc+ROIsize/2));
  count = count + 1;
  
  
 end     % end of for loop, "j = ImageNumber:ImageNumber + 19" over first 20 frames.
    
% divide sum of z's by number of frames to get averaged image of vesicle
    z_ave = z_sum/20;
        

% scale z values into appropriate range for fitting; we will scale back later
scalingfactor = max(max(z_ave))- min(min(z_ave));

z_scaled = z_ave/scalingfactor;
 
  figure(2)% ******************************************************** Fig 2
   
    % mesh plot of ROI for averaged raw vesicle intensities
    % (0,0) at upper left of image;
                mesh(xROI,yROI,z_ave)
                set(gca,'YDir','reverse'); % flip y axis so 00 is at upper left.
                title('Fig 2. Mesh of averaged z');
                xlabel('x (0,0)topleft');
                ylabel('y ');
                pause(1);
                  
% *  *  *  *  *  *  *  finished with raw data   *  *  *  *  *  *  *  
  
% % *************   Use fmincon  to fit single Gaussian to raw data  ********* 
% % %initial guesses for the 10 variables which define 1d Gaussian
x0 = [ max(max(z_scaled))-min(min(z_scaled)); 0.22; 0; 0; -min(min(z_scaled))];
% x0 = [ -.25; 0.22;  0; 0; 0.1]

% constraints 
A = []; 
b = [];
Aeq = []; 
beq = [];
 
% set upper and lower bounds for variables in x
     lb = [-3; 0;-5;-5; -2]; %[-3; 0; -5;-5;-1]
     ub = [ 3; 6; 5; 5; 0]; %[3; 6 5; 5; 1]
             
myfunflag = 0; % set flag to print x on first use of myfun
 
% set options for fmincon
options = optimset('largescale','off','hessian','off',...
    'Display','iter','tolFun',1e-8,'maxfunEvals',15000); 
 
% % ***************  *****************  *************  ************
% % PRIMARY USE OF fmincon   ******************************************* L289
     
     [x,fvalfunction] = fmincon(@myfun5,x0,A,b,Aeq,beq,lb,ub)
     
% fprintf('*************finished fmincon*********************\n')



figure(3)% ***********************************************************  Fig 3
contour(xROI,yROI,z_fit,10);
        % mesh(xROI,yROI,z_fit);
    set(gca,'YDir','reverse'); % flip y axis so 00 is at upper left.
    set(gca,'XGrid','on','YGrid','on');
    title('Fig 3. zfit');
    xlabel('xROI (col)');
    ylabel('yROI (row)');
    axis([-ROIsize/2 ROIsize/2 -ROIsize/2 ROIsize/2 0 2000]);%was -0.3 0.3
    pause(1);


figure (4)  % *******************************************************   Fig 4    
    % mesh plot of ROI for fitted, unscaled particle 
    % (0,0) at upper left of image;
                z_fit_unscaled = (x(1)*scalingfactor)*exp(-((x(2)*(XXROI-x(3)).^2)+(x(2) *(YYROI-x(4)).^2)))-((x(5)*scalingfactor)); 
                mesh(z_fit_unscaled)
                set(gca,'YDir','reverse'); % flip y axis so 00 is at upper left.
                title('Fig 4. Mesh of fitted, unscaled z');
                xlabel('x (0,0)top_left');
                ylabel('y ');
                
%find fitted max

fprintf('fitted max= \n');
max_int = max(max(z_fit_unscaled))
distance_sigma = ((1/(2 * x(2)))^0.5)/(pixeldistance);
radius = distance_sigma;
radius % distance across half a single gaussian,in microns
pause(2);
       
 
z_sum = 0;

fprintf('**********    end of this particle    ***************\n');
fprintf('*************************************************\n');

% Pick lines in M_sort for particle i
M_sort_r((M_sort_cut(:,3) == i-1),7) = radius/(10^6);  % converts from microns to meters for GSE calculations
M_sort_r((M_sort_cut(:,3) == i-1),8) = max_int;




end % end of outermost loop over i=1:numberMatFiles


% write data to Excel spreadsheet
% OUTPUTPATH2 = strcat(FolderName,'sizedata.xls')
% xlswrite(OUTPUTPATH2,data)
 
% % *******************************************************************
% % **************************      Nested Functions       ************
function f = myfun5(x)% % x is a 1-d array of 4 parameters (for 1 2d Gaussian)
    % % x is optimized by fmincon
    % % f is chi-squared, the measure of the goodness of fit to the data
    % % parameters in x are
    % % Gaussian 1 definition
    %         % x(1)=A1
    %         % x(2)=1/(2sigmax1^2)
    %         % x(3)=x_center
    %         % x(4)=y_center
    %         % x(5)=baseline
f = 0;

if (myfunflag == 0)
     fprintf('first use of myfun5. x0, x ,x');
     x0
     x
     myfunflag = 1
end
 
%     integral = 0;

% construct z_fit, the 2D array of the FITTED SINGLE Gaussian         
  z_fit =x(1)*exp(-((x(2)*(XXROI-x(3)).^2)+(x(2) *(YYROI-x(4)).^2))) - x(5); 
  
% compute f, the sumsquare of the deviations between fitted and observed, scaled z
   f = sum(sum((z_fit-z_scaled).^2)); 
  
end % end of nested function myfun5
 
 
end % end of fluor_sizing 
