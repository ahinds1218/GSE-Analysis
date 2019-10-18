function [M_sort] = loadVSTdata()

%Loads data from Video Spot Tracker and sorts into analysis-friendly array.


clc;
close all;
pause on;
% 
fprintf('Select *.vrpn.mat file that contains tracker data')
pause(2);
S = uiimport('-file');
fieldnames(S)
fieldnames(S.tracking);
M2=S.tracking.spot3DSecUsecIndexFramenumXYZRPY;
deleteMsort_rows;
M_sort = sortrows(M,3);
end

