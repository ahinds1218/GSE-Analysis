% GSE_Analysis.m by Amanda Smelser, created 6/28/12

% For GSE analysis of particles tracked by Video Spot Tracker (VST).

% VST output (*.vrpn file) must be converted to *.vrpn.mat file using
% vrpnLogToMat, provided by UNC, prior to running this program.



% This wrapper program calls the following functions:
%       [M_sort] = loadVSTdata
%       [check] = checktrackers(M_sort)
%       [M_sort_cut] = deleteparticles(M_sort)
%       [M_sort_r] = fluor_sizing(M_sort)
%       analysisdriver


clc;
clear all;
close all;

%************************USER INPUTS**************************************
%Shuttered or continuous imaging scheme?
imgscheme = 1; %0 for continuous, 1 for shuttered

% frame rate 
frame_rate =  100;

% pixels to microns conversion factor
pixeldistance = 9.23; %pixels/micron for PCO Edge camera at 60x (for magnification knob at 1x, use 9.23; for 1.5x, use 13.845)


%pixelvalues for 16-bit camera (PCO Edge)
pixelvalues = 65536;

%region of interest size around particle. Should be even number. Depends on particle size.
ROIsize = 30;

%temperature of particle environment (in Kelvins)
T = 310; % e.g. 298 (25 degrees C, room temperature) or 310 (37 degrees C, body temperature)

%automatic sizing of each particle? (0 for no, 1 for yes) If no, set fixed diameter in "onesize.m" script.
%(In onesize.m script, default diamter = 1.025 microns).
sizing = 0;

%**************************************************************************
% Define frame intervals for tau (tau = intervals over which MSDs are computed): 

    if imgscheme == 0;
        window     = [1 2 5 10 20 50 100 200 500 1000 2000 5000 9990]; %use for continuous data
   
    elseif imgscheme == 1;
        window     = [1 2 5 10 20 50 99 100 200 500 990]; %use for 10 burst shutter scheme
        
    end;

%Load video spot tracker data, and create sorted array in workspace named "M_sort":

    [M_sort] = loadVSTdata;
    
%Decide whether to check trackers
    [check] = checktrackers(M_sort);
 
%Provide the option of manually deleting certain particles/tracks prior to
%analysis
    [M_sort_cut] = deleteparticles(M_sort);

    
%Load images and determine radii of each tracked particle. Add column of
%radius data to M_sort, and create resulting "M_sort_r" array in workspace.

if sizing == 1;
    [M_sort_r_finalcut,numParticles] = fluor_sizing(M_sort_cut,pixeldistance,pixelvalues,ROIsize,check);
    
elseif sizing == 0;
    [M_sort_r_finalcut] = onesize(M_sort_cut,check); %for user-specified radius that is the same for all particles; no sizing
        
end;

%Analyze tracker data to determine MSDs. Incorporate radius data to
%perform GSE analysis.

    fprintf('BEGINNING COMPUTATION OF G* VALUES FROM MSDs OF INDIVIDUAL PARTICLES\n')
    fprintf('using analysisdriver, msd, and viscoelasticity');
[d,v] = analysisdriver(M_sort_r_finalcut,frame_rate,pixeldistance,window,T,imgscheme);

% 
% fprintf('MSD data has been saved to the folder containing the raw image data.\n');
% fprintf('Data is saved in columns as follows: 1)tau, 2)msd (um^2), 3)msd std dev (um^2), 4)number of particles \n')

