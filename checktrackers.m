function [check] = checktrackers(M_sort)
%Provides the option of visually inspecting movies of trackers with particles,
%so that misguided trackers can be deleted.

reply = input('Do you want to check Trackers with trackerMovie? (Don''t include single quotes)(y/n): ', 's');
if strcmp(reply, 'y')
   trackermovie
   check = 1;
elseif strcmp(reply, 'n')
    disp('trackerMovie skipped')
   check = 0;
end
end