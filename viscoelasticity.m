function viscoelasticity(d,window,v,freq_type,T,numParticles)    global FolderName FileNameImage; %based on function ve% Reference: Mason, Rheol Acta 39, 371-378(2000)% Original ve program written 10/06/05 by Jeremy Cribb, extensively% modified by George Holzwarth in 2008, and modified further by Amanda% Smelser in 2013.% Extensive modification of ve January-February, 2008 gholzwarth% ve_single generated from ve 2013 (Amanda)% computes the viscoelastic moduli G',G", eta',eta" from % mean-square displacement data of individual particles.% Inputs    % "d" is the input struct of msd's and alpha    % Struct d has fields for individual particles and for mean:    % FOR INDIVIDUAL particleS       % d.tau   is <tau x Nparticles>       % d.msd   is <tau x Nparticles>       % d.alpha is <tau-1 x Nparticles>       % d.alpha_spline is <tau x Nparticles>        % d.n is sample_count.    % FOR MEAN            % d.mean.msd     = msd_mean; % msd averaged over all particles. <11 x 1>        % d.mean.msd_std = msd_std;  % standard deviation of mean msd. <11 x 1>        % d.mean.alpha   = alpha_from_msd_mean using smoothed spline. 2013_08_22 gh                        % d.mean.tau;           % in [s]    % d.mean.msd;           % in [m^2]    % d.mean.error_msd;     % in [m^2]    % d.mean.n;    % particle_radius           in [m].    % "freq_type" is 'f' for [Hz] or 'w' for [rad/s]        % default is [Hz]% output struct v contains:    % particle  = G' and G" from msd of individual trackers    % mean  = G' and G" from mean msd    % error = stdev/sqrt(N)fprintf('Entering viscoelasticity function');% Allocate arrays to store data for all particlesgstar_spline_allVes = zeros(length(window),numParticles);gp_spline_allVes    = zeros(length(window),numParticles);gpp_spline_allVes   = zeros(length(window),numParticles);         for p = 1:numParticles;  % loop over particles. ****     % Call video_tracking_constants to select radii out of vvideo_tracking_constants;%Create array of radii (one per particle)particle_radius = v((v(:,4) == 0),R);     if (nargin < 3) || isempty(freq_type)         freq_type = 'f';   endif (nargin < 2) || isempty(particle_radius)       particle_radius = 0.5e-6; endif (nargin < 1) || isempty(d)        error('no data struct found'); end%Constants for calculationsk = 1.38e-23;% Multiply msd values for each particle by corresponding radiiparticle_radius2 = repmat(particle_radius,1,size(window,2));msd_x_radii = d.msd.*particle_radius2';% Determine the mean (across particles) of msd x radii for each window msd_radii_mean = mean(msd_x_radii,2);%Rename data for simplificationmsd     = d.msd;   % msd <tau x numParticles)tau     = d.mean.tau;  % <tau x 1>N       = d.mean.n(1:end); % number of each tracker at each tauclear alpha; % redefinition below    alpha = d.alpha;  %<tau-1 x numParticles>    alpha_spline = d.alpha_spline; % slope of splined log-log plot of MEAN msd.                                 % <tau x numParticles>fprintf('Print input msd*radii, msd, tau, N, alpha,alpha_spline\n');    msd_x_radii_um2 = msd_x_radii*1e12; %convert to [um]^2 so it's readable.    msd_um2=msd*1e12;                         %convert to[um]^2 so it's readable.%     std_msd_um2 = std_msd*1e12;               %      "    MYgamma_spline = gamma(1 + alpha_spline);  % alpha from smoothed mean msd, all particles    msd_x_radii = msd_x_radii(:,p); % size(msd_x_radii);msd_x_radii_um2 = msd_x_radii_um2(:,p);    % msd_x_radii_um2 = msd_x_radii_um2(1:end-1,p);  if p == 1;figure (12)    logtau   = log10(tau);    plot(logtau(1:length(window)-1,:),alpha(:,:),'-sk')         title('Fig 12. alpha vs log(tau)');         xlabel('tau (s)');         ylabel('alpha(radians)');         pause (1)end;        % compute f from tau       ****************************************f = 1 ./ tau;w = 2*pi*f;% ***********************                          ************************% ***********************                          ************************% ***********************      G*, G', G"  inside particle loop ***************% --- compute G*_spline, G'_spline, G"_spline for current particle.clcgstar_spline = (4*k*T)./((6)*(pi)*(msd_x_radii).*(MYgamma_spline(:,p))); % Note that we're     % inside the p loop but gstar does not have p as an explicit variable. gp_spline = gstar_spline .* cos(pi/2 .* alpha_spline(:,p));gpp_spline= gstar_spline .* sin(pi/2 .* alpha_spline(:,p));% store G*_spline, G'_spline, G"_spline for each particle for future averaging over all particlesgstar_spline_allVes(:,p) = gstar_spline(:);gp_spline_allVes(:,p)    = gp_spline(:);gpp_spline_allVes(:,p)   = gpp_spline(:);             loggstar_spline    = log10(gstar_spline);loggp_spline       = log10(gp_spline);loggpp_spline      = log10(gpp_spline);        figure (13)  % INSIDE loop over particles. Use alpha_spline   ******  Fig 13    logw = log10(w);plot(logw,loggstar_spline,'-or','LineWidth', 2); hold on;plot(logw,loggp_spline,'--sg'); hold on;plot(logw,loggpp_spline,'-.db');        title('Fig 13 Log G*-spline,Gp-spl,Gpp-spl ');        legend('G*','Gp(storage)','Gpp(loss)','Location','NorthWest');        xlabel('logw');         ylabel('logG');        hold on;        %  compute eta*, eta', eta" and their logs   ****************nstar_spline = gstar_spline .* tau;   % eta*-splinenp_spline    = gpp_spline   .* tau;   % eta'-spline  storagenpp_spline   = gp_spline    .* tau;   % eta"-spline  losslognstar_spline = log10(nstar_spline);lognp_spline    = log10(np_spline);lognpp_spline   = log10(npp_spline);logf     = log10(f);                figure (14)  % INSIDE LOOP over particles, n vs f using alpha  ***    Fig 41plot(logf,lognstar_spline,'-or','LineWidth', 2);hold on;plot(logf,lognp_spline,'--sk');hold on;plot(logf,lognpp_spline,'-.db');        title('Fig 14 Log eta*-spline,eta-p-spline,eta-pp-spline');        legend('eta*-spline','eta-p-spline (loss)','eta-pp-spline(storage)','Location','SouthEast');        xlabel('log f');         ylabel('log eta-spline');% Export logw,loggstar,logdgstar,loggp,logdgp,loggpp,logdgpp to excel file% Set up array G for export to Excelclc    G=zeros(length(window),4);   % G=zeros(length(window)-1,7);    size(G(:,1));    size(w);    size(gstar_spline);    size(gp_spline);    size(gpp_spline);     G(:,1) = w(1:length(window));       G(:,2) = gstar_spline;             G(:,3) = gp_spline;    G(:,4) = gpp_spline;        % Create filename        listfiles = ls(fullfile(FolderName,strcat('g_',FileNameImage,'ID',num2str(p-1),'*.xls')));    if isempty(listfiles) == 1        filename_g = strcat('g_',FileNameImage,'ID',num2str(p-1),'_01','.xls');        else        [a,~] = size(listfiles);        filecell = cellstr(listfiles);        lastfilename_cell = filecell(a);        lastfilename = char(lastfilename_cell);        lastfileext = char(regexp(lastfilename,'\d{2}.txt','match'));        lastfile = char(regexp(lastfileext,'\d{2}','match'));        filesequence = str2double(lastfile)+1;        FileNumber = num2str(filesequence,'%02.0f');        filename_g = strcat('g_',FileNameImage,'ID',num2str(p-1),'_',FileNumber,'.xls');            end% Export    xlswrite(filename_g,G);    fprintf('G data exported to excel file in "Current Folder"\n')% Export logf,lognstar,logdnstar,lognp,lognpp to excel file% Set up array Eta for export    Eta=zeros(length(window),4);% Eta=zeros(length(window)-1,7);    Eta(:,1) = f(1:length(window));    Eta(:,2) = nstar_spline;             Eta(:,3) = np_spline;    Eta(:,4) = npp_spline;    % Create filename for Eta        listfiles = ls(fullfile(FolderName,strcat('eta_',FileNameImage,'ID',num2str(p-1),'*.xls')));    if isempty(listfiles) == 1        filename_eta = strcat('eta_',FileNameImage,'ID',num2str(p-1),'_01','.xls');    else        [a,~] = size(listfiles);        filecell = cellstr(listfiles);        lastfilename_cell = filecell(a);        lastfilename = char(lastfilename_cell);        lastfile.ext = char(regexp(lastfilename,'\d{2}.txt','match'));        lastfile = char(regexp(lastfile.ext,'\d{2}','match'));        filesequence = str2double(lastfile)+1;        FileNumber = num2str(filesequence,'%02.0f');        filename_eta = strcat('eta_',FileNameImage,'ID',num2str(p-1),'_',FileNumber,'.xls');           end% Export    xlswrite(filename_eta,Eta);    fprintf('Eta data exported to excel file in "Current Folder"\n')     %  clear msd* g* n* log*;      clear gstar_spline gp_spline gpp_spline nstar_spline np_spline npp_spline loggstar_spline loggp_spline loggpp_spline lognstar_spline lognp_spline lognpp_spline;    clear msd_x_radii msd_x_radii_um2 logmsd_x_radii; end % end of loop over p = numVes. Loop started at L46  ***************** % *********************************************************************** % *********************************************************************** % **********************************************************************%--- compute AVERAGE G*-spline,G'-spline,G"-spline and their logsgstar_spline_ave    = nanmean(gstar_spline_allVes,2);gp_spline_ave        = nanmean(gp_spline_allVes,2);gpp_spline_ave       = nanmean(gpp_spline_allVes,2);gstar_spline_ave_log = log10(gstar_spline_ave);gp_spline_ave_log    = log10(gp_spline_ave);gpp_spline_ave_log   = log10(gpp_spline_ave); % Plot average G*-spline, average eta*-spline.figure (15) % G*-spline,G'-spline,G"-spline in log-log plot OUTSIDE LOOP ..Fig 51plot(logw,gstar_spline_ave_log,'-or','LineWidth', 2);hold onplot(logw,gp_spline_ave_log,'--sk');hold onplot(logw,gpp_spline_ave_log,'-.db');        title('Fig 15 Log G*,Gp,Gpp spline-ave vs log w');        legend('G*-spline-ave','Gp-spline(storage)','Gpp-spline(loss)','Location','NorthWest');        xlabel('logw'); ylabel('logG');pause (2)% --- compute AVERAGE eta*, eta',eta" and their logsnstar_spline_ave = gstar_spline_ave.*tau;   % eta*np_spline_ave    = gpp_spline_ave.*tau;     % eta' storagenpp_spline_ave   = gp_spline_ave.*tau;      % eta" losslognstar_spline_ave = log10(nstar_spline_ave);lognp_spline_ave    = log10(np_spline_ave);lognpp_spline_ave   = log10(npp_spline_ave);% Plot average eta-spline-ave*, eta-p-spline-ave, eta-pp-spline-ave as log-og figure (16) % eta*-spline-ave, eta', eta" in log-log plot  outside loop ******  51 **    plot(logf,lognstar_spline_ave,'-or','LineWidth', 2);hold on    plot(logf,lognp_spline_ave,'--sk');hold on    plot(logf,lognpp_spline_ave,'-.db');    title('Fig 16 Log eta*,etap,etapp spline-ave vs log f ');    legend('eta*-spline-ave','eta-p-spline-ave (loss)','eta-pp-spline-ave(storage)');    xlabel('log f'); ylabel('log eta'); pause (2) fprintf('Analysis complete') end