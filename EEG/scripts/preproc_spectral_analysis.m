%% ECSU 25.09.2024
clear all
cd '/home/nicolas/Documents/Okinawa/PCE/analyses_nico/';
restoredefaultpath; matlabrc


%% Data
pth_bh  = '/home/nicolas/Documents/Okinawa/PCE/';
pth_fg  =  '/home/nicolas/Documents/Okinawa/PCE/analyses_nico/figures/';
pth_raw = '/home/nicolas/Documents/Okinawa/PCE/analyses_nico/fornicholas/';
addpath(genpath('lbmap')); %
addpath(genpath('fooof_mat-main'));
addpath(genpath('/home/nicolas/Documents/Okinawa/PCE/analyses_nico/chronux_2_12.v03'));
addpath(genpath('/home/nicolas/Documents/Okinawa/PCE/analyses_nico/tensorlab_2016-03-28'));

dataset_folders = {'pce01230807', 'pce02230809', 'pce03230814', 'pce04230816', 'pce05230818', ...
    'pce06230821', 'pce07230823', 'pce08230825', 'pce09230828', 'pce10230830', ...
    'pce11230908', 'pce12230915', 'pce13230925', 'pce14230927', 'pce15230929', ...
    'pce16231006', 'pce17231010', 'pce18231016', 'pce19231020', 'pce20231023', ...
    'pce21231027', 'pce22231030', 'pce23231101', 'pce24231107', 'pce25231109', ...
    'pce26231121', 'pce27231201', 'pce28231204', 'pce29231208', 'pce30231211', ...
    'pce0x230804'};

%% Fixed parameters
ecg                             = 0;
chanlocs                    = importdata('actiCAP-10-20-Cap64.txt');
channels                    = chanlocs.textdata(:,2);
PCE_behaviour           = readtable([pth_bh 'PCE_Ultimate.csv']);
f_sampling                 = 1000;  % sampling frequency
dsample_factor          = 100;
dsampling_factor      = dsample_factor; % downsamplig factor
f_sampling_ds           = f_sampling/dsampling_factor;
n_trials                       = 2; % number of trials
n_subj                        = 2; % number of participants
subj_idx                     = vertcat([1:64],[65:128]);
n_experiments           = 31; % number of experiments
trial_duration             = 60; % duration of trial in seconds
normalization_type    = 1;
rest_trials                   = [2 3];
task_trials                   = [7 13];
foi                              = 1/20; % Frequency of interest
fpass                         =  [0.001 foi*4]; % Frqequency band
p_val                          = 0.05; % P-value for Jacknifing
ds_ecg                       = 5;
ds_skin                     = 50;
ds_resp                     = 25;

%% Variable parameters
poly_order                 = 8;           % Polynomial trend regression
poly_length               = 1055;     %  Length prior for polynomial trend
tapers_rs                    = [1 3];      % Multitaper configuration
tapers_task                = [1 3];      % Multitaper configuration
swin                           = [60 20];  % Sliding window for spectrum calculation
artifact_duration       = 10;          % Artifact length
tolerance                   = 1;            % Minimum number of artifacts to reject a channel
ica_components        = 8;           % Number of ICA components
do_ica                        = 1;           % Deploy ICA (optional but needed to have same number of signals for all subjects)
print_plots                 = 0;           % Only use it for a few subjects, else the computer wil choke making many plots :-)
print_spect                = 0;           % Only use it for a few subjects, else the computer wil choke making many plots :-)
factors                       = 3;           % Number of factors for matrix/tensor factorization of spectral results


%% Get spectral !!!
idx = 0; % Global index
for i_e =1:n_experiments  % there are only 31 experiments!
    for i_trial = 1:n_trials
        idx = idx+1;
        
        % For debugging:
        %i_e =1; i_t=1;
        
        %% Load EEG
        %         fname = strcat(pth_raw,dataset_folders{i_e},'_rest',...
        %             num2str(rest_trials(i_trial)) ,'trial', num2str(task_trials(i_trial)),  '.mat');
        %         load(fname);
        %         y_eeg = eval(['rest',num2str(rest_trials(i_trial))  'trial', num2str(task_trials(i_trial))]);
        %         clear  eval(strcat('rest',num2str(rest_trials(i_t)) ,'trial', num2str(task_trials(i_trial))'));
        
        
        fname = strcat([pth_raw  'combined/'],dataset_folders{i_e},'_rest',...
            num2str(rest_trials(i_trial)) ,'trial', num2str(task_trials(i_trial)),  '.mat');
        load(fname);
        y_eeg = eval(['rest',num2str(rest_trials(i_trial))  'trial', num2str(task_trials(i_trial))]);
        clear  eval(strcat('rest',num2str(rest_trials(i_t)) ,'trial', num2str(task_trials(i_trial))'));
        

        for i_p = 1:n_subj
            
            indices{idx} = [i_e; i_trial; i_p];
            
            rng(42); % to get ICA to work, must need to set the random generator to a fixed seed
            %close all
            
            %% Physiological signals
            if i_p ==1
                y_ecg   =  y_eeg(129,60000:(end-60000));
                y_skin  =  y_eeg(130,60000:(end-60000));
                y_resp  =  y_eeg(131,60000:(end-60000));
            end
            if i_p == 2
                y_ecg   =  y_eeg(132,60000:(end-60000));
                y_skin  =  y_eeg(133,60000:(end-60000));
                y_resp  =  y_eeg(134,60000:(end-60000));
            end
            
            y_ecg_ds   = downsample(y_ecg(:,(end- 239999):end)',ds_ecg); size(y_ecg_ds);
            y_skin_ds  = downsample(y_skin(:,(end- 239999):end)',ds_skin); size(y_skin_ds);
            y_resp_ds   = downsample(y_resp(:,(end- 239999):end)',ds_resp); size(y_resp_ds);
            
            A = zscore(y_ecg_ds);% (100:200);
            B = zscore(y_skin_ds); %(100:200);
            C = zscore(y_resp_ds); %(100:200);
            x_A = linspace(-180,60,240000 /ds_ecg);
            x_B = linspace(-180,60,240000 /ds_skin);
            x_C = linspace(-180,60,240000 /ds_resp);
            
            %print_plots =1
            if print_plots ==1
                figure,
                pos = get(gcf, 'Position');
                set(gcf, 'Position', [pos(1) pos(2)+200 1600, 800]); % Set plot size
                subplot 311
                plot(x_A,A(:,1),'b','LineWidth', 1.5);
                hold on
                plot(x_A,A(:,2),'r','LineWidth', 1.5);
                %line([0 0],[-0.5 0.5],'Color',[1 0 0],'LineStyle',':', 'LineWidth',2)
                grid on
                %ylim([-pi/6 +pi/6]);
                xlim([x_C(1)  x_C(end)]);
                xlabel('Time (seconds)');
                ylabel('ECG amplitude (a.u.)');
                set(gca, 'FontSize',12,'LineWidth', 1.5);
                set(gca, 'box', 'off');
                set(gcf, 'color', 'w');
                subplot 312
                plot(x_B,B(:,1),'b','LineWidth', 1.5);
                hold on
                plot(x_B,B(:,2),'r','LineWidth', 1.5);
                %line([0 0],[-3 3],'Color',[1 0 0],'LineStyle',':', 'LineWidth',2)
                grid on
                %ylim([-4 +4]);
                xlim([x_B(1)  x_B(end)]);
                ylabel('EDA amplitude (a.u.)');
                xlabel('Time (seconds)');
                set(gca, 'FontSize',12,'LineWidth', 1.5);
                set(gca, 'box', 'off');
                set(gcf, 'color', 'w');
                subplot 313
                plot(x_C,C(:,1),'b','LineWidth', 1.5);
                hold on
                plot(x_C,C(:,2),'r','LineWidth', 1.5);
                %line([0 0],[-3 3],'Color',[1 0 0],'LineStyle',':', 'LineWidth',2)
                grid on
                %ylim([-4 +4]);
                xlim([x_A(1)  x_A(end)]);
                ylabel('RP amplitude (a.u.)');
                xlabel('Time (seconds)');
                set(gca, 'FontSize',12,'LineWidth', 1.5);
                set(gca, 'box', 'off');
                set(gcf, 'color', 'w');
                fname = [ pth_fg  'physiological_signals.png']
                print(gcf, fname, '-dpng', '-r150', '-painters')
            end
            
          
            %% EEG  downsampling
            y        =  y_eeg(subj_idx(i_p,:),60000:(end-60000)); % original data
            dim_orig = size(y);
            y        = y(:,(end- 239999):end)'; size(y); % Extract from the end of the trial, 1 minute of task + 3 minutes of pre-task
            dim_epoch = size(y);
            rest_length = dim_orig(2)-60000*3;
            task_length = 60000;
            blank_length=  rest_length - task_length;
            task_onset = floor(dim_epoch(1)-task_length)/dsampling_factor;
            task_offset = floor(dim_epoch(1)-task_length)/dsampling_factor + 600;
            blank_onset = task_onset - floor(blank_length/dsampling_factor);
            
            y_ds   = downsample(y,dsampling_factor); size(y_ds);%
            
            
            %% Raw data:   full of sudden change artifacts due to movement
            if print_plots ==1
                figure,
                colormap(nawhimar)
                pos = get(gcf, 'Position');
                set(gcf, 'Position', [pos(1) pos(2) 1000, 250]); % Set plot size
                imagesc(y_ds'); colorbar
                hold on
                line([1800 1800],[64 0],'Color',[0 0 0],'LineWidth',1.5,'LineStyle','-.')
                x = linspace(1,240,2400);
                y_min = min(y_ds(:));
                y_max = max(y_ds(:));
                xlabel('Time point (a.u.)','FontSize',12);
                ylabel('Channel #','FontSize',12);
                set(gca, 'FontSize',12,'LineWidth', 1.5);
                set(gca, 'box', 'off');
                set(gcf, 'color', 'w');
                title('Raw EEG data')
                fname = [ pth_fg  'raw_timeSeries_exp-' num2str(i_e) '_trial-'  num2str(i_trial) '_P' num2str(i_p)  '.png']
                print(gcf, fname, '-dpng', '-r150', '-painters')
            end
            
            
            %% Slow local polynomial regression to remove sudden changes (but not spikes)
            %y_test= [y_ds(:,33), y_ds(:,44)]; size(y_test)
            %figure, plot(y_test);
            %sgf = sgolayfilt(y_test,8,155) -sgolayfilt(y_test,8,455);
            %figure, plot(sgf);
            %sgf = sgolayfilt(y_test,8,355); size(sgf)
            %y_dt =y_test - sgf;
            %figure, plot(y_dt);
            sgf = sgolayfilt(y_ds,poly_order,poly_length); size(sgf);
            %y_dt = sgolayfilt(y_ds,8,155) - sgolayfilt(y_ds,8,455);  % This works excellent but it is difficult to "prove" to the stubborn reviewer unfamiliar with signal processing
            %figure, imagesc(sgf'); colorbar
            y_dt =y_ds - sgf;
            %
            if print_plots ==1
                figure,
                colormap('gray')
                colormap(nawhimar)
                pos = get(gcf, 'Position');
                set(gcf, 'Position', [pos(1) pos(2) 1000, 250]); % Set plot size
                imagesc(y_dt'); colorbar
                hold on
                line([1800 1800],[64 0],'Color',[0 0 0],'LineWidth',1.5,'LineStyle','-.')
                x = linspace(1,240,2400);
                %y_min = min(y_dt(:))
                %y_max = max(y_dt(:))
                xlabel('Time point (a.u.)','FontSize',12);
                ylabel('Channel #','FontSize',12);
                set(gca, 'FontSize',12,'LineWidth', 1.5);
                set(gca, 'box', 'off');
                set(gcf, 'color', 'w');
                title('Detrended EEG data')
                fname = [ pth_fg  'detrended_timeSeries_exp-' num2str(i_e) '_trial-'  num2str(i_trial) '_P' num2str(i_p)  '.png']
                print(gcf, fname, '-dpng', '-r150', '-painters')
            end
            
            
            %% Channel rejection
            %y_filt = sosfilt(sos,y_ds);
            %figure, plot(y_filt);
            %outlier_idx = isoutlier(y_filt);  size(outlier_idx);
            outlier_idx = isoutlier(y_dt);  size(outlier_idx);%
            %figure, imagesc(outlier_idx'); colorbar
            %y_cleaned = filloutliers(y_ds ,'spline');
            %figure, plot(y_cleaned);
            if print_plots ==1
                figure,
                colormap('gray')
                pos = get(gcf, 'Position');
                set(gcf, 'Position', [pos(1) pos(2) 1000, 250]); % Set plot size
                imagesc(outlier_idx'); colorbar
                hold on
                line([1800 1800],[64 0],'Color',[0 0 0],'LineWidth',1.5,'LineStyle','-.')
                x = linspace(1,240,2400);
                %y_min = min(y_dt(:))
                %y_max = max(y_dt(:))
                xlabel('Time point (a.u.)','FontSize',12);
                ylabel('Channel #','FontSize',12);
                set(gca, 'FontSize',12,'LineWidth', 1.5);
                set(gca, 'box', 'off');
                set(gcf, 'color', 'w');
                title('Artifacts in detrended EEG data')
                fname = [ pth_fg  'flawed_timeSeries_exp-' num2str(i_e) '_trial-'  num2str(i_trial) '_P' num2str(i_p)  '.png']
                print(gcf, fname, '-dpng', '-r150', '-painters')
            end
            
            clear Pv nPv
            for n = 1:size(outlier_idx,2)
                f = find(diff([0,outlier_idx(:,n)',0]==1));
                p = f(1:2:end-1);  % Start indices
                y = f(2:2:end)-p;  % Consecutive ones’ counts
                Pv(n) = numel(find(y >= artifact_duration)); % Number of segments
                nPv(n) = sum(y(find(y >= artifact_duration))); % Length of segments
                %outliers_ep0chs()
            end
            
            
            %figure, imagesc(outlier_idx(:, Pv <  tolerance)'); colorbar
            
            %y_cleaned = filloutliers(y_dt ,'spline');
            %figure, plot(y_cleaned);
            
            if print_plots ==1
                figure,
                colormap('gray')
                pos = get(gcf, 'Position');
                set(gcf, 'Position', [pos(1) pos(2) 1000, 250]); % Set plot size
                subplot 121
                h=bar(Pv)
                h.EdgeColor = 'none';
                xlabel('Channel');
                ylabel('Artifact #');
                set(gca, 'FontSize',12,'LineWidth', 1.5);
                set(gca, 'box', 'off');
                set(gcf, 'color', 'w');
                xlim([0 65]);
                subplot 122
                h = bar(nPv)
                h.EdgeColor = 'none';
                xlabel('Channel');
                ylabel('time points');
                set(gca, 'FontSize',12,'LineWidth', 1.5);
                set(gca, 'box', 'off');
                set(gcf, 'color', 'w');
                xlim([0 65]);
                fname = [ pth_fg  'eeg_spikeArtifacts_exp-' num2str(i_e) '_trial-'  num2str(i_trial) '_P' num2str(i_p)  '.png']
                print(gcf, fname, '-dpng', '-r150', '-painters')
            end
            
            y_good = y_dt(:, Pv <=  tolerance);
            disp(strcat("Eligible channels :", num2str(size(y_good,2))))
            
            
            %% ICA (optional)
            if do_ica == 1
                
                if  size(y_good,2) >=  ica_components
                    
                    
                    %% Detrend
                    y_norm = zscore(zscore(y_good,0,2),0,1);
                    
                    %% ICA
                    
                    % PERFORM PCA
                    [coeff,Data_PCA,latent,tsquared,explained,mu] = pca(double(y_norm), 'NumComponents', ica_components); %this pca function need for run this code.
                    % compute and display explained variation
                    disp(strcat("Top ", string(ica_components), " principal components explain ", ...
                        string(sum(explained(1:ica_components))), " of variation"))
                    % compute independent components from principal components
                    Mdl = rica(Data_PCA,  ica_components); %this rica function need to run this code.
                    y_ICA = transform(Mdl, Data_PCA); size(y_ICA);
                    %             figure, imagesc(y_ICA'); colorbar
                    %             figure, plot(y_ICA);
                    %
                    %
                    %
                    %             %tSeries_cleaned = sosfilt(sos,y_ICA);
                    %             %figure, plot( tSeries_cleaned);
                    %
                    %
                    %             %% ICA spectra
                    %             %% Spectral analysis parameters
                    %             foi                 = 1/4; % Frequency of interest
                    %             fpass             =  [0.001 foi*4]; % Frqequency band
                    %             p_val             = 0.05; % P-value for Jacknifing
                    %             Slepians     = [3 5];
                    %             params = struct(...
                    %                 'Fs',f_sampling_ds ,...
                    %                 'tapers', Slepians ,...
                    %                 'pad', 1,...
                    %                 'fpass',fpass,...
                    %                 'err',[1 p_val],...
                    %                 'trialave',1);
                    %
                    %             %% Plot spectra
                    %             figure
                    %             c=0;
                    %             for i_a = 1:2
                    %                 for i_2 = 1:8
                    %                     c=c+1;
                    %                     subplot(8,2,c)
                    %                     title(num2str(c));
                    %                     %[pxx,f] = pwelch(y_ICA(:,c),f_sampling_ds );
                    %                     [S, f, Serr]     = mtspectrumc(y_ICA(:,c),params);
                    %                     plot(f,S ,'k-.','LineWidth', 1);
                    %                     hold on
                    %                     plot(f,Serr,'r','LineWidth', 1.5,'HandleVisibility','off')
                    %                     xlim([f(1) f(end)]);
                    %                     xlabel('Frequency (Hz)','FontSize',12);
                    %                     ylabel('Power (a.u.)','FontSize',12);
                    %                 end
                    %             end
                    y_clean = y_ICA; %y_norm ; %y_ica; %tSeries_cleaned;
                else
                    %% Bad experiment, do not save data
                    bad_experiment(idx,i_p,2)  = idx;
                end
                
            else
                
                %y_norm = zscore(zscore(y_good,0,2),0,1);
                y_norm = zscore(zscore(y_good,0,1),0,2);
                y_clean = y_norm;
                
                
            end
            
            
            %% Spectral analysis
            %tapers_rs                     = [1 3]; % Multitaper configuration
            %tapers_task                 = [1 3]; % Multitaper configuration
            %%
            params = struct(...
                'Fs',f_sampling_ds ,...
                'tapers', tapers_rs ,...
                'pad', 1,...
                'fpass',fpass,...
                'err',[1 p_val],...
                'trialave',1);
            
            total_minutes = size(y_dt,1)/f_sampling_ds/60;
            time_points    = size(y_dt,1)/total_minutes;
            
            y_pre  = y_clean(1:time_points*3,:); size(y_pre);
            y_task = y_clean((end- time_points+1):end,:); size(y_task);
            
            [S_pre, f_pre, Serr_pre]     = mtspectrumc(y_pre,params);
            dim_pre = size(S_pre);
            
            params = struct(...
                'Fs',f_sampling_ds ,...
                'tapers', tapers_task ,...
                'pad', 1,...
                'fpass',fpass,...
                'err',[1 p_val],...
                'trialave',1);
            
            [S_task_, f_task_, Serr_task_]     = mtspectrumc(y_task,params);
            dim_task = size(S_task_);
            
            S_task = interp1(1:dim_task(1), S_task_, linspace(1,dim_task(1),dim_pre(1)))';
            Serr_task = [interp1(1:dim_task(1), Serr_task_(1,:), linspace(1,dim_task(1),dim_pre(1)));...
                interp1(1:dim_task(1), Serr_task_(2,:), linspace(1,dim_task(1),dim_pre(1)))]; size( Serr_task );
            f = f_pre;
            
            
            % Min-max normalization
            dim = size(S_pre);
            % Pre-task (baseline): plot 1
            mindata = min(S_pre);
            maxdata = max(S_pre);
            %mindata = min(Serr_pre(:,1));
            %maxdata = max(Serr_pre(:,2));
            S_pre = bsxfun(@rdivide, bsxfun(@minus, S_pre, mindata), maxdata - mindata);
            Serr_pre = bsxfun(@rdivide, bsxfun(@minus, Serr_pre, mindata), maxdata - mindata);
            % Normalize task: plot 2
            mindata = min(S_task);
            maxdata = max(S_task);
            S_task_ = bsxfun(@rdivide, bsxfun(@minus, S_task_, mindata), maxdata - mindata);
            Serr_task_ = bsxfun(@rdivide, bsxfun(@minus, Serr_task_, mindata), maxdata - mindata);
            % interpolation of spectrum in task to baseline length for subtraction: plot 3
            mindata = min(S_task);
            maxdata = max(S_task);
            %mindata = min(Serr_task(:,1));
            %maxdata = max(Serr_task(:,2));
            S_task = bsxfun(@rdivide, bsxfun(@minus, S_task, mindata), maxdata - mindata);
            Serr_task = bsxfun(@rdivide, bsxfun(@minus, Serr_task, mindata), maxdata - mindata);
            
            
            params = struct(...
                'Fs',f_sampling_ds ,...
                'tapers', tapers_task ,...
                'pad', 1,...
                'fpass',fpass,...
                'err',[1 p_val],...
                'trialave',1);
            
            
            logScale        = 0; % TFR scale
            [S,tr,fr,Serr] = mtspecgramc(y_clean,swin,params); % spectrogram
            dim = size(S);
            P = zeros(dim); size(P);
            
            if logScale == 1
                S = real(10*log10(S));
            end
            
            
             spectras = 0;

            if print_spect ==1
                figure,
                %colormap(nawhimar)
                pos = get(gcf, 'Position');
                set(gcf, 'Position', [pos(1) pos(2) 1000, 250]); % Set plot size
                %pos = get(gcf, 'Position');
                %set(gcf, 'Position', [pos(1)-10000+idx*100 pos(2)+idx*100+idx*20 900, 200]);
                subplot 131
                spectrum = S_pre;
                spec_err = Serr_pre;
                plot(f,spectrum ,'k-.','LineWidth', 1);
                hold on
                plot(f,spec_err,'k','LineWidth', 1.5,'HandleVisibility','off')
                xlim([f(1) f(end)]);
                ylim([0 max(spec_err(:))]);
                xlabel('Frequency (Hz)','FontSize',12);
                ylabel('Power (a.u.)','FontSize',12);
                grid on
                set(gca, 'FontSize',12,'LineWidth', 1.5);
                set(gca, 'box', 'off');
                
                subplot 132
                spectrum = S_task_;
                spec_err = Serr_task_;
                plot(f_task_,spectrum ,'k-.','LineWidth', 1);
                hold on
                plot(f_task_,spec_err,'k','LineWidth', 1.5,'HandleVisibility','off')
                xlim([f_task_(1) f_task_(end)]);
                ylim([0 max(spec_err(:))]);
                xlabel('Frequency (Hz)','FontSize',12);
                ylabel('Power (a.u.)','FontSize',12);
                grid on
                set(gca, 'FontSize',12,'LineWidth', 1.5);
                set(gca, 'box', 'off');
                
                subplot 133
                %spectrum = (S_task-S_pre)./S_pre;
                spectrum = S_pre - S_task;
                %spec_err = (Serr_task-Serr_pre)./Serr_pre;
                plot(f,spectrum ,'r','LineWidth', 2);
                %hold on
                %plot(f,spec_err,'b','LineWidth', 1.5,'HandleVisibility','off')
                xlim([f(1) f(end)]);
                ylim([-1 1]);
                xlabel('Frequency (Hz)','FontSize',12);
                ylabel('task - baseline','FontSize',12);
                grid on
                set(gca, 'FontSize',12,'LineWidth', 1.5);
                set(gca, 'box', 'off');
                set(gcf, 'color', 'w');
                if do_ica == 1
                    fname = [ pth_fg  'slow_spectra_ICA_exp-' num2str(i_e) '_trial-'  num2str(i_trial) '_P' num2str(i_p)  '.png']
                else
                    fname = [ pth_fg  'slow_spectra_exp-' num2str(i_e) '_trial-'  num2str(i_trial) '_P' num2str(i_p)  '.png']
                end
                print(gcf, fname, '-dpng', '-r150', '-painters')
                
                
                
                %% Time frequency chart
                
                %% Plot spectrogram
                figure,
                %colormap(nawhimar)
                pos = get(gcf, 'Positio');
                set(gcf, 'Position', [pos(1) pos(2) 1000, 350]); % Set plot size
                %colormap(nawhimar)
                %colormap(lbmap(dim(2),'RedBlue'))
                colormap(lbmap(100,'BlueGray'))
                %colormap(lbmap(dim(2),'BrownBlue' ))
                switch spectras
                    case 0 % uncorrected
       
                        x = S(:,:);
                        mindata = min(x(:));
                        maxdata = max(x(:));
                        P_norm(:,:) = bsxfun(@rdivide, bsxfun(@minus, S(:,:), mindata), maxdata - mindata);
                        pcolor(tr,fr,P_norm');
                         pcolor(tr,fr,P_norm');
                        specType = '_uncorrrected';
                    case 1 % frequency corrected using global mean of frequency band
                        for i_t = 1:dim(1)
                            for i_f = 1:dim(2)
                                P(i_t,i_f)  =(S(i_t,i_f)-mean(S(:,i_f),1))./std(S(:,i_f),1);
                            end
                        end
                        pcolor(tr,fr,P');
                        specType = '_frequency_corrected';
                    case 2 % frequency corrected using baseline mean of frequency band
                        for i_t = 1:dim(1)
                            for i_f = 1:dim(2)
                                bs = S(1:round(dim(1)/4),i_f);
                                P(i_t,i_f) = (S(i_t,i_f)-mean(bs(:)))./std(bs(:));
                            end
                            
                        end
                        pcolor(tr,fr,P');
                        specType = '_baseline_corrected';
                        
                    case 3
                        
                         for i_t = 1:dim(1)
                            for i_f = 1:dim(2)
                                P(i_t,i_f)  =(S(i_t,i_f)-mean(S(:,i_f),1))./std(S(:,i_f),1);
                            end
                         end
                        
                            for i_f = 1:dim(2)
                                %P(i_s,i_t,i_f)  =(S(i_s,i_t,i_f)-mean(S(i_s,:,i_f),2))./std(S(i_s,:,i_f),1);
                                %P(i_s,i_t,i_f)  =S(i_s,i_t,i_f); %(S(i_s,i_t,i_f)-mean(S(i_s,:,i_f),2))./std(S(i_s,:,i_f),1);
                                P(:,i_f)  = (S(:,i_f)-mean(S(:,i_f),1))./std(S(:,i_f),1);
                            end
                        
                        x = P(:,:);
                        mindata = min(x(:));
                        maxdata = max(x(:));
                        P_norm(:,:) = bsxfun(@rdivide, bsxfun(@minus, P(:,:), mindata), maxdata - mindata);
                        pcolor(tr,fr,P_norm');
                        specType = '_frequency_corrected';
                        
                end
                shading  flat %interp;
                hold on
                line([180 180],[f(1) f(end)],'Color',[0 0 0],'LineWidth',1.5,'LineStyle','-.')
                xlim([tr(1) tr(end)]);
                %ylim([1 100]);
                c = colorbar;
                title('Spectrogram');
                xlabel('Time (sec)');
                legend('boxoff');
                set(gca, 'FontSize',18,'LineWidth', 2);
                if logScale == 1
                    set(get(c,'title'),'String','log(dB/Hz)');
                else
                    set(get(c,'title'),'String','Power');
                end
                ylabel('Frequency (Hz)');
                if do_ica == 1
                    fname = [ pth_fg  'slow_spectrograms_ICA' specType '_exp-' num2str(i_e) '_trial-'  num2str(i_trial) '_P' num2str(i_p)  '.png']
                else
                    fname = [ pth_fg  'slow_spectrograms' specType '_exp-' num2str(i_e) '_trial-'  num2str(i_trial) '_P' num2str(i_p)  '.png']
                end
                print(gcf, fname, '-dpng', '-r150', '-painters')
                
                
                
            end
            
            
            %% Save
            
            % Signals
            TFR.eeg_clean{idx,i_p} = y_good; % mean(y_norm,2);
            
            if do_ica == 1
                TFR.eeg_ica{idx,i_p} = y_ICA;
            else
                TFR.eeg_ica{idx,i_p} = [];
            end
            
            % Spectrum
            TFR.spectra_pret(idx,i_p,:) = S_pre;
            TFR.spectra_post(idx,i_p,:) = S_task_;
            TFR.spectra_error_pret(idx,i_p,:,:) = Serr_pre;
            TFR.spectra_error_post(idx,i_p,:,:) = Serr_task_;
            TFR.spectra_frequencies_pret = f;
            TFR.spectra_frequencies_post = f_task_;
            
            % Time frequency representation/chart
            TFR.spectrum(idx,i_p,:,:) = S;
            %TFR.spectrum(idx,i_p,:,:) = P_norm;
            TFR.spectrum_freqs = fr;
            TFR.spectrum_time = tr;
            TFR.blank_onset(idx,i_p) = blank_onset ;
            

            if print_plots ==1
                close all
            end
            
        end
        
        TFR.repetition_index = indices;
        PHYS.ecg(idx,i_p,:,:) = y_ecg_ds;
        PHYS.eda(idx,i_p,:,:) = y_skin_ds;
        PHYS.resp(idx,i_p,:,:) = y_resp_ds;
        PHYS.dsample = [ds_ecg, ds_skin, ds_resp]
        
        %%
        close all
    end
end
