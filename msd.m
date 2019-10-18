function [d] = msd(v, window, imgscheme)

% msd.m original version written 01/21/06 by Jeremy Cribb (UNC-CH CS), modified 
% by George Holzwarth (WFU Physics) for multiple particles, drift correction, 
% and alpha smoothed with spline, and by Amanda Smelser (WFU Biochemistry)for shuttered imaging scheme

% Calls:video_tracking_constants
%       slope_evaluator

%  Arguments 
            %%  "files"  is the filename containing video tracking data (wildcards ok) 
    %  "v"  is the matrix containing video tracking data (wildcards ok) 
    %  "window" is a vector containing window sizes of tau when computing MSD. 
    %  "calib_um" is the microns per pixel conversion unit

%  Fields of output structure "d" are for individual particles only
    % d.tau  double
    % d.msd  double
    % d.n = sample_count double because particleID's are indexed by 0.
% Output matrix for Excel: d_matrix
    % 3 columns per particle. Col 1=tau; Col 2 particleID; Col 3=msd.
    % rows: different values of tau
clc

% print out the input data table v, with x,y in micrometers
    % (need to convert X and Y to microns, print, and convert back to m
   
    video_tracking_constants;
    
  
for particleID = 0 : get_particlemax(v);% **** Primary loop one particle at a time  **
    
    b = get_particle(v, particleID);    % Select particular particle from data file
      
    % for every window size (or tau)
    for w = 1:length(window)
        clear r2;
        clear frame;
        clear A1;clear A2; clear B1; clear B2;
        
        if imgscheme == 0;
      
            A1 = b(1:(end-window(w)),X);
            A2 = b(1:(end-window(w)),Y);
       
            B1 = b((window(w)+1):end,X);
            B2 = b((window(w)+1):end,Y);

        
        elseif imgscheme == 1;
            
            
            if window(w) < 99;
                for l = 1:10;
                A1((1+(99-window(w))*(l-1)):(l*(99-window(w)))) = b((1+100*(l-1)):(99-window(w)+100*(l-1)),X);
                A2((1+(99-window(w))*(l-1)):(l*(99-window(w)))) = b((1+100*(l-1)):(99-window(w)+100*(l-1)),Y);
        
                
				B1((1+(99-window(w))*(l-1)):(l*(99-window(w)))) = b((1+window(w)+100*(l-1)):(99+100*(l-1)),X);
                B2((1+(99-window(w))*(l-1)):(l*(99-window(w)))) = b((1+window(w)+100*(l-1)):(99+100*(l-1)),Y);
                end;  
                
            elseif window(w) == 99;
                
                A1 = b(1:100:end-window(w),X);
                A2 = b(1:100:end-window(w),Y);
        
                B1 = b(window(w)+1:100:end,X);
                B2 = b(window(w)+1:100:end,Y);
                
                                              
            elseif window(w) > 99;
                A1 = b(1:end-window(w),X);
                A2 = b(1:end-window(w),Y);

                B1 = b(window(w)+1:end,X);
                B2 = b(window(w)+1:end,Y);
                
            end;
        end;   % ************************************** end of if shutter=1

% Compute MSDs

        r2 = ( B1 - A1 ).^2 + ( B2 - A2 ).^2;
        
        msd(w, particleID+1) = mean(r2);
        
            if imgscheme == 0;
             tau(w, particleID+1) = window(w) * mean(diff(b(:,TIME)));
        
       
            elseif imgscheme == 1;
             tau(w, particleID+1) = b(window(w)+1,TIME);
        
            end;
     end     % ends loop over window (tau values)
%    
      
     % compute alpha as point-by-point for 1 particle
     Aalpha(:,particleID+1) = log10(tau(1:end-1,particleID+1));
     Balpha(:,particleID+1) = log10(tau(2:end  ,particleID+1));
     Calpha(:,particleID+1) = log10(msd(1:end-1,particleID+1));
     Dalpha(:,particleID+1) = log10(msd(2:end  ,particleID+1));
               
     alpha(:,particleID+1) = (Dalpha(:,particleID+1) - Calpha(:,particleID+1))./(Balpha(:,particleID+1)-Aalpha(:,particleID+1));
         

end         % ends loop over a single particle *******************************

% ************************************************************************

% plot data for all particles, all tau on linear scale
figure(8);
    plot(tau,msd,'-or',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerSize',8);
        title('Fig 8 MSD vs Tau, LINEAR AXES');
        xlabel('Tau');
        ylabel('MSD');
        pause(2);
         
% Prepare data for log-log plots.
logtau = log10(tau);
logmsd = log10(msd);

sample_count = sum(~isnan(logmsd),2);

% remove window sizes that returned no data
idx = find(sample_count > 0);
logtau = logtau(idx,:);
logmsd = logmsd(idx,:);
sample_count = sample_count(idx);

       

% ************************************* Start insertion

% Compute MEAN values for msd, alpha, and alpha_from_msd_mean
msd_mean     = nanmean((msd)');  
msd_std      =  nanstd((msd)');   
msd_mean_log = log10(msd_mean);
    
% compute alpha_from_msd_mean    
[msd_mean_log_smooth,alpha_from_msd_mean] = slope_evaluator( logtau(:,1),msd_mean_log );
   
mean_alpha   = mean(alpha,2);       % alpha is array for each particle
                                    
std_alpha    = std(alpha,1,2);      % std particle-to-particle

% Plots to show benefits of smoothed spline 

figure(9);
  plot(logtau(:,1),msd_mean_log,'-ob',logtau,msd_mean_log_smooth,'--squarer');
       xlabel('log(tau)');
       ylabel('log(msd mean');
       title('Fig 9 Smoothed vs. Unaltered MSD');
       legend('msd-mean-log','msd-mean-log-smooth','Location','NorthWest');
  pause(1)
  
figure(10);  % *****     NEW     *****************************************
% plot alpha vs log tau and alpha_spline vs logtau
  plot(log10(tau(1:end-1,1)),mean_alpha,'-or',...
       'LineWidth',2,...
       'MarkerEdgeColor','k',...
       'MarkerFaceColor',[.49 1 .63],...
       'MarkerSize',8);
     title('Fig 10 mean-alpha and alpha-from-msd-mean');
     xlabel('log(tau)');
     ylabel('mean-alpha');
     hold on
  plot(log10(tau(:,1)),alpha_from_msd_mean,'--squareb');
     % errorbar(log10(tau(1:end-1,1)),mean_alpha,std_alpha)
     legend('point-by-point derivative','derivative using spline');
  pause(1);
00;
   

    
% output msd for each particle to an Excel-readable file:
% Excel Output matrix = d_matrix
    % 3 columns per particle. Col 1=particleID; Col 2=tau; Col 3=msd.
    % rows: different values of tau
    last_particle = get_particlemax(v);
    numParticles  = last_particle +1;
    numCols   = 4*numParticles;
    numRows   = length(window);
    d_matrix  = zeros(numRows,numCols);
       
    for particleID = 0:last_particle
           
        d_matrix(:,1+particleID*3)=particleID;
            
        d_matrix(:,2+particleID*3)=tau(:,particleID + 1);
           
        d_matrix(:,3+particleID*3)=msd(:,particleID+1);
           
        d_matrix(:,4+particleID*3)=alpha_from_msd_mean(:);
    end
    
    
% define struct d to be passed to ve_single
d.tau      = tau; % <11 x 5> for 5 particles . 5 cols are the same
d.msd      = msd;  % individual particles <11 x 5 for 5 particles>. 5 cols differ
d.alpha    = alpha;%   "
d.alpha_spline = repmat(alpha_from_msd_mean,1,numParticles); 
                        %  % <11 x 5> for 5 particles . 5 cols are the same
d.n            = sample_count;% <11 x 1 for 5 particles. 

d.mean.msd     = msd_mean; % msd averaged over all particles. <11 x 1>
d.mean.msd_std = msd_std;  % standard deviation of mean msd. <11 x 1>
d.mean.alpha   = alpha_from_msd_mean;  % alpha from msd-mean and slope-eval
clear alpha  % needed because alpha will be redefined in viscoelasticity
