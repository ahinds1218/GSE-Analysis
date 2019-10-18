function M_sort_r_finalcut = onesize(M_sort_cut,check) %,pixeldistance,pixelvalues,ROIsize)

%dic_onesize adds specified radius (same for all particles) to an added
%column in M_sort_r.

global FolderName FileNameImage ImageNumber

clc;
close all;
% clear all;
fprintf('*************ADDING COLUMN OF KNOWN RADII TO M_SORT_CUT***********\n')
M_sort_r_finalcut = M_sort_cut; % allocate, initialize M_sort_r_finalcut


%Bring in images, if haven't previously
if check == 0;
fprintf('Select first image in series');
pause(2);
[filename,FolderName] = uigetfile('*.tif','Select first image in series');
imagenum = char(regexp(filename,'\d{4}','match'));
ImageNumber = str2double(imagenum);
FileNameImage = char(regexp(filename,'\w+(?=\d{4}\.tif)','match'));

elseif check == 1;
fprintf('Adding user-specified particle radius to M_sort array\n')
pause(1);
end;

% % find number of vesicles
numParticles = max(M_sort_cut(:,3))+1 %number of tracked vesicles
if numParticles == 0
    fprintf('Check M_sort array\n');
end
      
% fprintf('**********    START LOOP over Particles         **************\n');

for i=1:numParticles % number of tracked vesicles
% fprintf('**********    Inside i loop over particles *********************\n');

i

N = M_sort_cut((M_sort_cut(:,3) == i-1),:);
% 
% 
A = N(:,4) %framenumber
B = N(:,3) %spotID
C = N(:,5) %x
D = N(:,6) %y
E = N(:,7) %r set in VST
C_offset =  2;  % move center of ROI to center of bead; was 3.5
D_offset = 0;  %             ";was -2.0
     
C = C+C_offset; %             "
D = D+D_offset; %             "

fprintf('i=\n');
i

% --- put r data into col 7 of M_sort array to create M_sort_r_finalcut array

% Pick lines in M_sort for bead i
M_sort_r_finalcut((M_sort_cut(:,3) == i-1),7) = 1.025/(2*10^6); %distance2/(2 * 10^6);  %2 converts from diameter to radius, and from microns to meters for GSE calculations


end % end of outermost loop over i=1:numberMatFiles
