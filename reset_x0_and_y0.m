function [v] = reset_x0_and_y0(v)
% subtract mean_x and mean_y from x,y for each bead
% g. holzwarth 01/21/08
% intended purpose: apply to raw tracking data 
% so (x,y)_bar(averaged over t)is at (0,0) for each bead.

% INPUT variable: v
%   data columns: TIME IMAGE ID FRAME X Y R ROLL PITCH YAW;
    data = v;
% OUTPUT: [v]
%   v = data corrected so that
%        x = x-x_mean and
%        y = y-y_mean
fprintf('Entered reset_x0_and_y0\n');
% global variables
    global TIME IMAGE ID FRAME X Y R ROLL PITCH YAW;     
    video_tracking_constants;  
% handle the argument list
    if ~exist('data');
        error('No data found. '); 
    end;
  
    if nargin < 1 || isempty(data); 
        error('No data found. Operator incompetence.'); 
    end;
   
% Select a bead, compute x_mean, subtract from each x, etc.    
for Nparticle = 0 : get_particlemax(data) % Nparticle is particle_number
   clc;
   thisparticle_idx = find(data(:,ID) == Nparticle);
   data_particle=(data(thisparticle_idx,:));
   % plot original track and shifted track
          figure(5);
          subplot(1,2,1);
          plot(data(thisparticle_idx,X),data(thisparticle_idx,Y));
                title(['Fig. 5 Original Track ', num2str(Nparticle)]);
                xlabel('X');
                ylabel('Y');
                  
   x_mean=mean(data(thisparticle_idx,X));
   y_mean=mean(data(thisparticle_idx,Y));

   data(thisparticle_idx,X)=data(thisparticle_idx,X)-x_mean;
   data(thisparticle_idx,Y)=data(thisparticle_idx,Y)-y_mean;
 
   % plot shifted track
            subplot(1,2,2);
            plot(data(thisparticle_idx,X),data(thisparticle_idx,Y));
                title(['Shifted by Mean ', num2str(Nparticle)]);
                xlabel('X');
                ylabel('Y');
                pause(2);
end    % end of for loop over Nparticle
close(5);
% Construct output matrix [v]
v = data;   % columns:TIME ID FRAME X Y R ROLL PITCH YAW;

return;     % end of function [v] = reset_x0_and_y0.


     
     
     

