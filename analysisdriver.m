function [d,v]= analysisdriver(M_sort_r_finalcut,frame_rate,pixeldistance,window,T,imgscheme)
global FolderName FileNameImage 
pause on;
% 
% This program is a driver for msd and viscoelasticity calculations created by George Holzwarth 
% to set up input params for msd.m. Also writes output data to hard drive
% in *.txt files.


    
% calls:    
%            video_tracking_constants
%            reset_x0_and_y0
%            remove_drift
%            msd
%            viscoelasticity


    % video_tracking_constants.m sets
    %    TIME = 1    IMAGE = 2   ID = 3  
    %    FRAME = 4   X    = 5    Y  = 6    
    %    R = 7  SIGMA = 8    PEAK = 9   
    %    RADIAL=10 PULSE=11  SEQ   = 12


% Step 1. Load data from vrpn.mat file
 
[v] = M_sort_r_finalcut;

numParticles = max(v(:,3))+1;

% construct TIME column

if imgscheme == 0;
    v(:,1)= v(:,4)/frame_rate;  
elseif imgscheme == 1;
    for i = 1:numParticles
    v(((1:100)+(i-1)*999),1)    = v(1:100,4)/frame_rate;
    v(((101:200)+(i-1)*999),1)  = (1000:1099)/frame_rate;
    v(((201:300)+(i-1)*999),1)  = (2000:2099)/frame_rate;
    v(((301:400)+(i-1)*999),1)  = (3000:3099)/frame_rate;
    v(((401:500)+(i-1)*999),1)  = (4000:4099)/frame_rate;
    v(((501:600)+(i-1)*999),1)  = (5000:5099)/frame_rate;
    v(((601:700)+(i-1)*999),1)  = (6000:6099)/frame_rate;
    v(((701:800)+(i-1)*999),1)  = (7000:7099)/frame_rate;
    v(((801:900)+(i-1)*999),1)  = (8000:8099)/frame_rate;
    v(((901:999)+(i-1)*999),1) = (9000:9098)/frame_rate;
    end
    
end;
    
% global variables
global TIME IMAGE ID FRAME X Y R SIGMA PEAK;
video_tracking_constants;

%scale data to physical units
v(:,X:Y) = v(:,X:Y) / pixeldistance * 1e-6; %convert video coords from pixels to meters 

fprintf('reset_x0_and_y0 is next\n');
% **************************************************    
% Step 2. Offset each bead so x_mean=0, y_mean=0
 [v]=reset_x0_and_y0(v);
    fprintf('finished reset_x0_and_y0\n');
    fprintf('remove_drift is next\n');
    
% **************************************************    
% Step 3. Correct for COM drift
             
  [v]=remove_drift(v);
  fprintf('Finished remove_drift\n');
  fprintf('msd for individual particles is next\n');
      
% **************************************************    
% Step 4. Compute msd for each bead, plot msd vs tau
[d] = msd(v,window,imgscheme);
        fprintf('Finished msd\n');
% Output d is a struct   
    % d.tau  is matrix length(tau) x Nbeads
    % d.msd  is matrix length(tau) x Nbeads          "   x   "
    % d.n              length(tau)     1
    % d.mean_alpha is length(tau)-1
    % d.std_alpha  is length(tau)-1
%******************************************
% Determine mean msd, stdev_msd, dG', dG", deta', deta"
fprintf('printing tau, mean_msd, stdev(mean_msd)');
d.mean.tau       = mean(d.tau,2);
d.mean.msd       = mean(d.msd,2);
d.mean.stdev_msd = nanstd(d.msd,1,2);
d.mean.n         = d.n;

% Struct d now has fields for individual particles and for mean:
    % d.tau  is matrix length(tau) x Nbeads
    % d.msd  is matrix length(tau) x Nbeads          "   x   "
    % d.n              length(tau)     1
    % d.mean_alpha
    % d.std_alpha
    
    % d.mean.tau = tau
    % d.mean.msd = msd_mean. Arraylength = length(tau)
    % d.mean.stdev_msd = standard deviation of mean_msd at each tau
    % d.mean.n = sample_count % because particleID's are indexed by 0.

% **************************************************    
for p = 1:numParticles;
% Step 5 Construct msd_out,  one msd for each particle in matrix form.
% Units are um^2 instead of m^2    
msd_out=zeros(length(window),2);
msd_out(:,1) = d.mean.tau;
msd_out(:,2) = d.msd(:,p)*1e12;         % msd [um]


% show log-log plot of msds

figure (11);
% convert to log-log scale.
    xxx = log10(msd_out(:,1));
    yyy = log10(msd_out(:,2));
% fit straight line to data   
    plot(xxx,yyy, '-square');
        hold on
        xlabel('logbase10(\tau) [s]');
        ylabel('logbase10 msd [um^2]');
        title('Fig 11 Particle MSD vs Tau, log-log scale');
        grid on;
        hold on;
pause (2);      
% Save msd_out to HD. 

% construct filename based on raw data filename without overwriting previous analyses
% of same data

constrfilename = sprintf('%s%s%s%03d','msd_',FileNameImage,'ID',p-1);

listfiles = ls(fullfile(FolderName,strcat(constrfilename,'*.txt')));

 
if isempty(listfiles) == 1;

     filename_msd = strcat(constrfilename,'_01','.txt')
    
    
else
    [a,~] = size(listfiles);
    filecell = cellstr(listfiles);
    lastfilename_cell = filecell(a);
    lastfilename = char(lastfilename_cell);
    lastfileext = char(regexp(lastfilename,'\d{2}.txt','match'));
    lastfile = char(regexp(lastfileext,'\d{2}','match'));
    filesequence = str2double(lastfile)+1;
    FileNumber = num2str(filesequence,'%02.0f');
    filename_msd = strcat(constrfilename,'_',FileNumber,'.txt')
end

dlmwrite(fullfile(FolderName,filename_msd), msd_out,'delimiter','\t','newline','pc');
clear msd_out;
end; % end of p=1:numParticles ********************************************************** 
% --- END of loop over beads

fprintf('msd completed.\n');
fprintf('MSD data exported to excel files in "Current Folder"\n')
fprintf('viscoelasticity is next\n');
pause(2);

% **************************************************    
% Step 6 Compute G',G", eta', eta"
viscoelasticity(d,window,v,'f',T,numParticles);
return;
