%% Cross-decoding of Inter-Participant slow rhythms
% This script 1) decodes joint succes in two groups and  2) all other click types
% (including single success) in two groups. It uses Linear Discriminant
% Analysis with 5 fold crossvalidation across response-aligned time points.
%  ECSU 13.10.2024
% author: Nicolas Gravel, 13.20.2024 nicolas.gravel@gmail.com
clear all
%close all

preproc           = 0;
periclick          = 1;
ecg                  = 0;


cd '/home/nicolas/Documents/Okinawa/PCE/analyses_nico/';

addpath(genpath('reissuewiththematlabscript'));
%% Data
pth  = '/home/nicolas/Documents/Okinawa/PCE/';
pth_fg  =  '/home/nicolas/Documents/Okinawa/PCE/analyses_nico/figures/';

addpath(genpath('lbmap')); %
dataset_folders = {'pce01230807', 'pce02230809', 'pce03230814', 'pce04230816', 'pce05230818', ...
    'pce06230821', 'pce07230823', 'pce08230825', 'pce09230828', 'pce10230830', ...
    'pce11230908', 'pce12230915', 'pce13230925', 'pce14230927', 'pce15230929', ...
    'pce16231006', 'pce17231010', 'pce18231016', 'pce19231020', 'pce20231023', ...
    'pce21231027', 'pce22231030', 'pce23231101', 'pce24231107', 'pce25231109', ...
    'pce26231121', 'pce27231201', 'pce28231204', 'pce29231208', 'pce30231211', ...
    'pce32230804'};

%% Variables
chanlocs                    = importdata('actiCAP-10-20-Cap64.txt');
channels                    = chanlocs.textdata(:,2);
PCE_behaviour           = readtable([pth 'PCE_Ultimate.csv']);
f_sampling                 = 1000;  % sampling frequency
dsampling_factor       = 10; % downsamplig factor
f_sampling_ds            = f_sampling/dsampling_factor;
%freqs_low                   = [0.02,0.01,3,8, 12, 35];  % high pass
%freqs_high                  = [0.06,0.5,12,12, 30, 65]; % low pass
freqs_low                   = [0.03,0.08,8,30];  % high pass
freqs_high                 = [0.08,0.12,14,80]; % low pass
n_trials                       = 18; % number of trials
n_subj                        = 2; % number of participants
n_experiments            = 31; % number of experiments
trial_duration             = 60; % duration of trial in seconds
peri_click_time           = 15; % epoch (in seconds) around button click
click_time_stamp       = table2array(PCE_behaviour(:,7)); % button clicks
idx_experiment          = repelem(1:n_experiments,n_trials*2); % index to experiment
idx_P1                        = 1:n_subj:(length(click_time_stamp)); % P1
idx_P2                        = 2:n_subj:(length(click_time_stamp)); % P2
idx_subjects               = [idx_P1;idx_P2]; % index participants
idx_null                     = isnan(click_time_stamp);
idx_P1_early              = find(click_time_stamp(idx_P1)<peri_click_time);
idx_P2_early              = find(click_time_stamp(idx_P2)<peri_click_time);
idx_P1_edgy              = find(click_time_stamp(idx_P1)>(trial_duration-peri_click_time));
idx_P2_edgy              = find(click_time_stamp(idx_P2)>(trial_duration-peri_click_time));
null_clicks                 = union(find(idx_null(idx_P1)),find(idx_null(idx_P2)));
early_clicks               = union(idx_P1_early,idx_P2_early);
edgy_clicks               = union(idx_P1_edgy,idx_P2_edgy); intersect(early_clicks,edgy_clicks); % must be 0
bad_trials                  = sort([null_clicks;early_clicks;edgy_clicks]);
artifact_duration       = 5;          % Artifact length
tolerance                   = 3;          % Minimum number of artifacts to reject a channel
ica_components       = 10;        % Number of ICA components
do_ica                       = 0;           % Deploy ICA (optional but needed to have same number of signals for all subjects)
poly_order                 = 8;           % Polynomial trend regression
poly_length               = 455;     %  Length prior for polynomial trend


%% Get epochs around button click

idx_freq                       =1; % 0.05 Hz
f_low                           = freqs_low(idx_freq);
f_high                         = freqs_high(idx_freq);
f_nrm_low                  = f_low /(f_sampling_ds/2);
f_nrm_high                = f_high /(f_sampling_ds/2);
[z,p,k]                         = butter(4,[f_nrm_low f_nrm_high],'bandpass'); % determine filter coefficients:
filt_005                   = zp2sos(z,p,k); % Filter

idx_freq                       =2; % 0.1 Hz
f_low                           = freqs_low(idx_freq);
f_high                         = freqs_high(idx_freq);
f_nrm_low                  = f_low /(f_sampling_ds/2);
f_nrm_high                = f_high /(f_sampling_ds/2);
[z,p,k]                         = butter(4,[f_nrm_low f_nrm_high],'bandpass'); % determine filter coefficients:
filt_01                   = zp2sos(z,p,k); % Filter

idx_freq                       =3; % 0.1 Hz
f_low                           = freqs_low(idx_freq);
f_high                         = freqs_high(idx_freq);
f_nrm_low                  = f_low /(f_sampling_ds/2);
f_nrm_high                = f_high /(f_sampling_ds/2);
[z,p,k]                         = butter(4,[f_nrm_low f_nrm_high],'bandpass'); % determine filter coefficients:
filt_10                   = zp2sos(z,p,k); % Filter



if preproc == 1
    
    % Global index
    tic
    idx_exp = 0;
    idx = 0;
    
    %%
    for i_e =1:n_experiments  % there are only 31 experiments!
        
        idx_exp = idx_exp + 1;
        
        
        
        for i_t = 1:n_trials
            
            idx = idx+1;
            
            for i_p = 1:n_subj
                
                %%  Only keep eligible trials
                if ismember(idx,bad_trials) == 1
                    fprintf('%s\r\n',strcat('experiment: ',num2str(i_e),...
                        ' participant: ',num2str(i_p),' trial: ',num2str(i_t),'  skip trial'));
                else
                    
                    %% Click onset
                    click_onset = click_time_stamp(idx_subjects(i_p,idx));
                    preOnset  = round((click_onset-peri_click_time)*f_sampling);
                    postOnset = round((click_onset+peri_click_time)*f_sampling);
                    
                    
                    %% Respiration
                    fname = strcat(pth,'RESP/',dataset_folders{i_e},'/pce',dataset_folders{i_e}(4:5),'_P',...
                        num2str(i_p) ,'_Trial', num2str(i_t),  '.mat');
                    if exist(fname) == 0
                        fprintf('%s\r\n',strcat('experiment: ',num2str(i_e),' participant: ',num2str(i_p),' trial: ',num2str(i_t), '     no RESP file!'));
                        y_resp = [];
                    else
                        load(fname);
                        y_resp = eval(strcat('pce_P', num2str(i_p) ,'_Trial'));
                        if periclick == 1
                            epoch_resp{idx,i_p} =  downsample(zscore(y_resp(preOnset:postOnset))',dsampling_factor)';
                        else
                            epoch_resp{idx,i_p} =  downsample(zscore(y_resp)',dsampling_factor)';
                            
                        end
                        
                    end
                    
                    %% EDA
                    fname = strcat(pth,'EDA/',dataset_folders{i_e},'/pce',dataset_folders{i_e}(4:5),'_P',...
                        num2str(i_p) ,'_Trial', num2str(i_t),  '.mat');
                    if exist(fname) == 0
                        fprintf('%s\r\n',strcat('experiment: ',num2str(i_e),' participant: ',num2str(i_p),' trial: ',num2str(i_t), '     no EDA file!'));
                        y_eda = [];
                    else
                        load(fname);
                        y_eda = eval(strcat('pce_P', num2str(i_p) ,'_Trial'));
                        if periclick == 1
                            epoch_eda{idx,i_p} =  downsample(zscore(y_eda(preOnset:postOnset))',dsampling_factor)';
                        else
                            epoch_eda{idx,i_p} =  downsample(zscore(y_eda)',dsampling_factor)';
                        end
                    end
                    
                    %% Load ECG
                    fname = strcat(pth,'ECG/',dataset_folders{i_e},'/pce',dataset_folders{i_e}(4:5),'_P',...
                        num2str(i_p) ,'_Trial', num2str(i_t),  '.mat');
                    if exist(fname) == 0
                        fprintf('%s\r\n',strcat('experiment: ',num2str(i_e),' participant: ',num2str(i_p),' trial: ',num2str(i_t), '     no ECG file!'));
                        y_ecg = [];
                    else
                        load(fname);
                        y_ecg = eval(strcat('pce_P', num2str(i_p) ,'_Trial'));
                        if periclick == 1
                            epoch_ecg{idx,i_p} =  downsample(zscore(y_ecg(preOnset:postOnset))',dsampling_factor)';
                            %epoch_ecg{idx,i_p} =  zscore(y_ecg(preOnset:postOnset))';
                        else
                            epoch_ecg{idx,i_p} =  downsample(zscore(y_ecg)',dsampling_factor)';
                        end
                    end
                    
                    %% Compute HRV
                    if ecg == 1
                        fname = strcat(pth,'ECG/',dataset_folders{i_e},'/pce',dataset_folders{i_e}(4:5),'_P',...
                            num2str(i_p) ,'_Trial', num2str(i_t),  '.mat');
                        if exist(fname) == 0
                            fprintf('%s\r\n',strcat('experiment: ',num2str(i_e),' participant: ',num2str(i_p),' trial: ',num2str(i_t), '    no ECG file!'));
                            y_ecg = [];
                        else
                            load(fname);
                            y_ecg = eval(strcat('pce_P', num2str(i_p) ,'_Trial'));
                            
                            %epoch_ecg{idx,i_p} = zscore(y_ecg(preOnset:postOnset))';
                            y = zscore(y_ecg(preOnset:postOnset));
                            
                            %% Cardio
                            buffer = 3000;
                            overlap  = 250;
                            shift = buffer-overlap;
                            for ci=1:fix((size(y,2)-buffer)/shift +1)
                                start_index = 1+(ci-1)*shift;
                                stop_index = min(start_index+ buffer-1);
                                time_index(ci) = round((start_index+stop_index)/2);
                                chunk = y(start_index:stop_index);
                                rr =  HRV.rrx(chunk);
                                out_rrx(ci,:) = mean(rr(2:end));
                                out_entro(ci,:) = HRV.ApEn(chunk);
                                out_HRV(ci,:) = HRV.rrHRV(rr(2:end));
                                clear HRV
                            end
                            out_rrx(isnan(out_rrx))=0;
                            out_entro(isnan(out_entro))=0;
                            out_HRV(isnan(out_HRV))=0;
                            %% a little sanity check
                            %figure,
                            %pos = get(gcf, 'Position');
                            %set(gcf, 'Position', [pos(1) pos(2) 1400, 300]); % Set plot size
                            %subplot 411;plot(y);subplot 412;stairs(out_rrx);subplot 413;stairs(out_entro);subplot 414;stairs(out_HRV);
                            epoch_ecg_HRV{idx,i_p} = out_HRV;
                            epoch_ecg_entro{idx,i_p} = out_entro;
                            epoch_ecg_rrx{idx,i_p} = out_rrx;
                            clear out
                        end
                    end
                    
                    %% Load EEG
                    fname = strcat(pth,'EEG_raw/',dataset_folders{i_e},'/pce',dataset_folders{i_e}(4:5),'_P',...
                        num2str(i_p) ,'_Trial', num2str(i_t),  '.mat');
                    load(fname);
                    
                    y_eeg = eval(strcat('pce_P', num2str(i_p) ,'_Trial'))';
                    
                    y_ds = downsample(y_eeg,dsampling_factor);
                    
                    y_poly = sgolayfilt(double(y_ds),poly_order,poly_length); size(y_poly);
                    
                    %SSres = norm(y_ds-y_poly,2)^2
                    
%                     SSres = sum((y_ds-y_poly).^2);
%                     SStot = sum((y_ds - mean(y_ds,1)).^2);
%                     varex_poly = (1-(SSres/SStot))*100
%                     
%                     
%                     y_tot    = y_ds-mean(y_ds,1);
%                     y_res   = y_ds-y_poly;
%                     
%                     for i = 1:64
%                         varex(i) = corr(y_tot(:,i),y_res(:,i))^2;
%                     end
%                     
%                     varex_poly = varex; 
%                     
%                     y_dt = y_ds - y_poly ;
%                     
%                                         figure,
%                                         subplot 311
%                                         plot(y_ds)
%                                         subplot 312
%                                         plot(y_poly)
%                                         subplot 313
%                                         plot(y_dt)
                    
                    
                    
                    %% Click onset
                    
                    
                    click_onset = click_time_stamp(idx_subjects(i_p,idx));
                    
                    if periclick == 1
                        preOnset  = round((click_onset-peri_click_time)*f_sampling_ds);
                        postOnset = round((click_onset+peri_click_time)*f_sampling_ds);
                        
                        y_dt_ = y_dt(preOnset:postOnset,:);  % If there are artifcts within the desired range, then reject
                    else
                        y_dt_ = y_dt;
                    end
                    
                    %figure, imagesc(y_ds')
                    outlier_idx = isoutlier(y_dt_);  size(outlier_idx);%
                    clear Pv nPv
                    for n = 1:size(outlier_idx,2)
                        f = find(diff([0,outlier_idx(:,n)',0]==1));
                        p = f(1:2:end-1);  % Start indices
                        y = f(2:2:end)-p;  % Consecutive ones’ counts
                        Pv(n) = numel(find(y >= artifact_duration)); % Number of segments
                        nPv(n) = sum(y(find(y >= artifact_duration))); % Length of segments
                        %outliers_ep0chs()
                    end
                    y_good = y_dt(:,Pv <=  tolerance);
                    disp(strcat("Eligible channels :", num2str(size(y_good,2))))
                    
                    varex_poly = mean(varex_poly(Pv <=  tolerance));
                    disp(strcat("mean arhytmic variance (mean across eligible channels) :", num2str(varex_poly)));

                    
                    %y_norm = zscore(y_good,0,1);
                    %figure, imagesc(y_norm')
                    %figure,
                    %subplot 311
                    %plot(y_good(1:100,10:20))
                    %subplot 312
                    %y_norm = zscore(y_good,0,1);
                    %plot(y_norm(1:100,10:20))
                    %subplot 313
                    %y_norm = zscore(y_good,0,2);
                    %y_norm = (y_good - mean(y_good(:)))/std(y_good(:));
                    %plot(y_norm(1:100,10:20))
                    
                    y_norm = y_good;
                    
                    
                    
                    y_filt_005 = sosfilt(filt_005, y_norm);
                    y_mean_005 = nanmean(y_filt_005,2);
                    
                    z = hilbert(y_mean_005);
                    rho_005 = real(abs(z));
                    phi_005 = angle(z);
                    
                    
                    y_filt_01 = sosfilt(filt_01, y_norm);
                    y_mean_01 = nanmean(y_filt_01,2);
                    z = hilbert(y_mean_01);
                    rho_01 = real(abs(z));
                    phi_01 = angle(z);
                    
                    y_filt_10 = sosfilt(filt_10, y_norm);
                    y_mean_10 = nanmean(y_filt_10,2);
                    z = hilbert(y_mean_10);
                    rho_10 = real(abs(z));
                    phi_10 = angle(z);
                    
                    %figure, plot(y_mean_005)
                    %figure, plot(y_mean_01)
                    %figure, plot(rho_005)
                    %figure, plot(rho_01)
                    %figure, plot(phi_005)
                    %figure, plot(phi_01)
                    
                    
                    %figure, imagesc(y_filt')
                    %figure, plot(y_filt)
                    if do_ica ==1
                        % PERFORM PCA
                        [coeff,Data_PCA,latent,tsquared,explained,mu] = pca(y_good, 'NumComponents', ica_components); %this pca function need for run this code.
                        % compute and display explained variation
                        disp(strcat("Top ", string(ica_components), " principal components explain ", ...
                            string(sum(explained(1:ica_components))), " of variation"))
                        % compute independent components from principal components
                        Mdl = rica(Data_PCA,  ica_components); %this rica function need to run this code.
                        y_ICA = transform(Mdl, Data_PCA); size(y_ICA);
                        %figure, imagesc(y_ICA')
                        epoch_ica{idx,i_p} = y_ICA';
                    end
                    
                    if periclick == 1
                        preOnset  = round((click_onset-peri_click_time)*f_sampling_ds);
                        postOnset = round((click_onset+peri_click_time)*f_sampling_ds);
                        
                        y_mean_005 =  y_mean_005(preOnset:postOnset,:);
                        phi_005 = phi_005(preOnset:postOnset,:);
                        rho_005 = rho_005(preOnset:postOnset,:);
                        
                        y_mean_01 =  y_mean_01(preOnset:postOnset,:);
                        phi_01 = phi_01(preOnset:postOnset,:);
                        rho_01 = rho_01(preOnset:postOnset,:);
                        
                        y_mean_10 =  y_mean_10(preOnset:postOnset,:);
                        phi_10 = phi_10(preOnset:postOnset,:);
                        rho_10 = rho_10(preOnset:postOnset,:);
                    end
                    
                    
                    
                    
                    %hola
                    
                    %figure, plot(y_mean_005)
                    %figure, plot(y_mean_01)
                    %figure, plot(rho_005)
                    %figure, plot(rho_01)
                    %figure, plot(phi_005)
                    %figure, plot(phi_01)
                    
                    
                    epoch_phi_005{idx,i_p} = phi_005';
                    epoch_rho_005{idx,i_p} = rho_005';
                    epoch_filt_005{idx,i_p} = y_mean_005';
                    
                    
                    epoch_phi_01{idx,i_p} = phi_01';
                    epoch_rho_01{idx,i_p} = rho_01';
                    epoch_filt_01{idx,i_p} = y_mean_01';
                    
                    epoch_phi_10{idx,i_p} = phi_10';
                    epoch_rho_10{idx,i_p} = rho_10';
                    epoch_filt_10{idx,i_p} = y_mean_10';
                    
                    epoch_varex_poly{idx,i_p} = varex_poly;
                    epoch_eeg{idx,i_p} = nanmean(y_good,2)';
                    
                    indices{idx_exp,i_t,i_p,idx}= [idx_exp,i_t,i_p,idx];
                    
                end
            end
            
        end
    end
    
    
    
    %% Save ECSU physiological markers around click onset
    ECSU.varex_poly     = epoch_varex_poly;
    ECSU.eeg                = epoch_eeg;
    ECSU.filt_005         = epoch_filt_005;
    ECSU.theta_005     = epoch_phi_005;
    ECSU.rho_005        = epoch_rho_005;
    ECSU.filt_01         = epoch_filt_01;
    ECSU.theta_01     = epoch_phi_01;
    ECSU.rho_01        = epoch_rho_01;
    ECSU.filt_10         = epoch_filt_10;
    ECSU.theta_10     = epoch_phi_10;
    ECSU.rho_10        = epoch_rho_10;
    ECSU.resp      = epoch_resp;
    ECSU.eda       = epoch_eda;
    ECSU.ecg       = epoch_ecg;
    ECSU.indices = indices;
    
    if do_ica ==1
        ECSU.ica    = epoch_ica;
    end
    
    if ecg ==1
        ECSU.hrv     = epoch_ecg_HRV;
        ECSU.entro  = epoch_ecg_entro;
        ECSU.rrx      = epoch_ecg_rrx;
    end
    
    if periclick == 1
        save ECSU_aroundClick ECSU
    else
        save ECSU_wholeTask ECSU
    end
    
else
    
    if periclick == 1
        load ECSU_aroundClick ECSU
    else
        load ECSU_wholeTask ECSU
    end
    
    
    epoch_varex_poly = ECSU.varex_poly;
    epoch_eeg           = ECSU.eeg ;
    epoch_filt_005   = ECSU.filt_005 ;
    epoch_phi_005  = ECSU.theta_005;
    epoch_rho_005  = ECSU.rho_005;
    epoch_filt_01   = ECSU.filt_01 ;
    epoch_phi_01  = ECSU.theta_01;
    epoch_rho_01  = ECSU.rho_01;
    epoch_filt_10   = ECSU.filt_10 ;
    epoch_phi_10  = ECSU.theta_10;
    epoch_rho_10  = ECSU.rho_10;
    epoch_resp = ECSU.resp;
    epoch_eda  = ECSU.eda;
    epoch_ecg  = ECSU.ecg;
    indices        = ECSU.indices;
    
    if ecg ==1
        epoch_ecg_HRV   = ECSU.hrv;
        epoch_ecg_entro = ECSU.entro;
        epoch_ecg_rrx     = ECSU.rrx;
    end
    
    if do_ica ==1
        epoch_ica = ECSU.ica;
    end
    
end


%LDA_v31