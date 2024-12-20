%% Cross-validated regression of infra-slow EEG vs click types
% November 27, 2024
% nicolas.gravel@gmail.com

%close all
clear all

addpath(genpath('MVPA-Light-master')); %

clickResponse_classification

norm                = 1; % global normalization
step                 = 4; % time points to include in regression
steps               =  [1:step:3000]; % time span around task response (click onset)

% Regresion model 
cfg                   = [];
cfg.model        = 'ridge';
cfg.cv               = 'kfold'; % conservative and robust, yet a little noisy
cfg.metric        = {'mae' 'mse' 'r_squared'};
cfg.repeat        = 100;

cfg.hyperparameter = [];
cfg.hyperparameter.lambda = 5;
cfg.feedback   = 0;

cfg_stat                   = [];
cfg_stat.metric        = 'mse';
cfg_stat.test            = 'permutation';
cfg_stat.n_permutations  = 2;
cfg_stat.statistic       = 'wilcoxon';
cfg_stat.null            = [];
cfg_stat.alpha         = 0.05;
cfg_stat.tail            = -1;
cfg_stat.design          = 'between';
cfg_stat.feedback   = 0;


model_labels   = { 'js_xP' , 'njs_xP'}


clear rndperf perf acc

i_sig = 0;
for is_signal = 4:5
    
    i_sig = i_sig +1;
    
    i_label = 0;
    for label_type = 1:2%length(model_labels)
        
        i_label =i_label+1;
                
        %% Extract data
        
        % EEG
        idx_eeg_trials = ~cellfun(@isempty,epoch_eeg); % bad trials are removed from both participants
        epoch_eeg_trials_P1 =find(idx_eeg_trials(:,1));
        epoch_eeg_trials_P2 =find(idx_eeg_trials(:,2));
        epoch_eeg_trials = intersect(epoch_eeg_trials_P1 ,epoch_eeg_trials_P2);
        fprintf('%s\r\n',[num2str(length(epoch_eeg_trials)) '  eligible EEG trials']);
        
        % RESP
        idx_resp_trials = ~cellfun(@isempty,epoch_resp); % bad trials are removed from both participants
        epoch_resp_trials_P1 =find(idx_resp_trials(:,1));
        epoch_resp_trials_P2 =find(idx_resp_trials(:,2));
        epoch_resp_trials = intersect(epoch_resp_trials_P1 ,epoch_resp_trials_P2);
        fprintf('%s\r\n',[num2str(length(epoch_resp_trials)) '  actual RESP trials']);
        
        % EDA
        idx_eda_trials = ~cellfun(@isempty,epoch_eda); % bad trials are removed from both participants
        epoch_eda_trials_P1 =find(idx_eda_trials(:,1));
        epoch_eda_trials_P2 =find(idx_eda_trials(:,2));
        epoch_eda_trials = intersect(epoch_eda_trials_P1 ,epoch_eda_trials_P2);
        fprintf('%s\r\n',[num2str(length(epoch_eda_trials)) '  actual EDA trials']);
        
        % ECG
        epoch_ecg = epoch_ecg;
        idx_ecg_trials = ~cellfun(@isempty,epoch_ecg); % bad trials are removed from both participants
        epoch_ecg_trials_P1 =find(idx_ecg_trials(:,1));
        epoch_ecg_trials_P2 =find(idx_ecg_trials(:,2));
        epoch_ecg_trials = intersect(epoch_ecg_trials_P1 ,epoch_ecg_trials_P2);
        fprintf('%s\r\n',[num2str(length(epoch_ecg_trials)) '  actual ECG trials']);
        
        % Index data according to available recording modalities
        switch is_signal
            
            case 2
                [idx,P] = find(idx_eeg_trials & idx_eda_trials & idx_resp_trials);
                signal  = epoch_rho_10;
                titl = 'EEG (10 {\itHz}) ';
                fname = [ pth_fg  'linregmodel _alpha'];
                
            case 3
                [idx,P] = find(idx_eeg_trials & idx_eda_trials & idx_resp_trials);
                signal = epoch_resp;
                titl = 'Respiratory pressure';
                fname = [ pth_fg  'linregmodel _Respiration'];
                
            case 4
                [idx,P] = find(idx_eeg_trials & idx_eda_trials & idx_resp_trials);
                signal  = epoch_rho_005;
                titl = 'EEG (0.05 {\itHz}) ';
                fname = [ pth_fg  'linregmodel _005'];
                
            case 5
                [idx,P] = find(idx_eeg_trials & idx_eda_trials & idx_resp_trials);
                signal  = epoch_rho_01;
                titl = 'EEG (0.1 {\itHz}) ';
                fname = [ pth_fg  'linregmodel _01'];
                
            case 6
                [idx,P] = find(idx_eeg_trials & idx_eda_trials & idx_resp_trials);
                signal  = epoch_eda;
                titl = 'EDA ';
                fname = [ pth_fg  'linregmodel _EDA'];
                
            case 7
                [idx,P] = find(idx_eeg_trials & idx_eda_trials & idx_resp_trials);
                signal  = epoch_eeg;
                titl = 'EEG';
                fname = [ pth_fg  'linregmodel _EEG'];
                
            case 8
                [idx,P] = find(idx_ecg_trials);
                signal  = epoch_ecg_HRV;
                %minmax  = 1;
                titl = 'HRV';
                fname = [ pth_fg  'linregmodel _HRV'];
                
            case 9
                [idx,P] = find(idx_ecg_trials);
                signal  = epoch_ecg_rrx;
                titl = 'RR';
                fname = [ pth_fg  'linregmodel _RR'];
                
        end
        
        
        %% Index click types
        
        index = vertcat((indices{:}))
        
        eligibile_index = indices(:,:,:,P); size(eligibile_index)
        vertcat(squeeze(cell2mat(eligibile_index(:,1,1,1))));
        
        eligible_data = intersect(idx(P==1),idx(P==2));
        
        click_time_stamp   = table2array(PCE_behaviour(:,7)); % button clicks
        
        click_P1 = click_time_stamp(idx_P1);
        click_P2 = click_time_stamp(idx_P2);
        
        click_P1 = click_P1(eligible_data);
        click_P2 = click_P2(eligible_data);
        
        clickOnset    = [click_P1,click_P2]; size(clickOnset)
        
        click_type       = table2array(PCE_behaviour(:,10)); % click type
        
        click_type_P1 = click_type(idx_P1);
        click_type_P2 = click_type(idx_P2);
        
        click_type_P1 = click_type_P1(eligible_data);
        click_type_P2 = click_type_P2(eligible_data);
        
        clickType    = [click_type_P1,click_type_P2]; size(clickType)
        
        
        %% Select trials for each click type
        
        % Identify click types
        js_P1 = strmatch('JS',click_type_P1,'exact');
        ss_P1 = strmatch('SS',click_type_P1,'exact');
        w_P1 = strmatch('W',click_type_P1,'exact');
        
        js_P2 = strmatch('JS',click_type_P2,'exact');
        ss_P2 = strmatch('SS',click_type_P2,'exact');
        w_P2 = strmatch('W',click_type_P2,'exact');
        
        % Extract data
        P1 = cell2mat(signal(eligible_data,1)); size(P1)
        P2 = cell2mat(signal(eligible_data,2)); size(P2)
        
        % Variance explained by the polynomial regression of a-rhythmic
        P1_varex = epoch_varex_poly(eligible_data,1);
        P2_varex = epoch_varex_poly(eligible_data,2);
        varex = [P1_varex P1_varex];
        varex = vertcat(varex{:});
        
        %figure, hist(varex(:),10);
        
        
        %% Normalize each participant EEG data within experiment to minimize
        % arbitrary amplitude changes from experiment to experiment:
        
        % 1st) Find eligible EEG experiments:
        idx_experiment = repelem(1:n_experiments,n_trials); % index to experiment
        idx = 0
        for i_e =1:n_experiments  % there are only 31 experiments!
            for i_t = 1:n_trials
                idx = idx+1;
                for i_p = 1:n_subj
                    if ismember(idx,bad_trials) == 1  %%  Only keep eligible trials
                        fprintf('%s\r\n',strcat('experiment: ',num2str(i_e),...
                            ' participant: ',num2str(i_p),' trial: ',num2str(i_t),'  skip trial'));
                    else
                        experiments{idx} = idx_experiment(idx);
                    end
                end
            end
        end
        good_experiments  = idx_experiment(~cellfun(@isempty,experiments));
        within_good_experiments_trials = arrayfun(@(x)length(find(good_experiments == x)), unique(good_experiments), 'Uniform', false);
        
        % 2nd) Apply minmax normalization across good trials within experiment
        if is_signal == 4 | is_signal == 5
            [idx,P] = find(idx_eeg_trials);
            eeg_data = intersect(idx(P==1),idx(P==2));
            P1_ = cell2mat(signal(eeg_data,1)); size(P1);
            P2_ = cell2mat(signal(eeg_data,2)); size(P2);
            exp_index = [good_experiments',eeg_data]
            [~,ia,~] = unique(good_experiments)
            P1 =zeros(size(eeg_data,1),size(P1,2));
            P2 =zeros(size(eeg_data,1),size(P2,2));
            for i_exp = 1: 31-1
                % Normalize P1
                X_ = P1_(ia(i_exp):within_good_experiments_trials{i_exp},:);
                mindata = min(X_(:));
                maxdata = max(X_(:));
                X = bsxfun(@rdivide, bsxfun(@minus, X_, mindata), maxdata - mindata);
                P1(ia(i_exp):within_good_experiments_trials{i_exp},:) = X;
                % Normalize P2
                X_ = P2_(ia(i_exp):within_good_experiments_trials{i_exp},:);
                mindata = min(X_(:));
                maxdata = max(X_(:));
                X = bsxfun(@rdivide, bsxfun(@minus, X_, mindata), maxdata - mindata);
                P2(ia(i_exp):within_good_experiments_trials{i_exp},:) = X;
            end
        end
        
        % Select click type specific data
        P1_js = P1(js_P1,1:3000); size(P1_js)
        P1_ss = P1(ss_P1,1:3000); size(P1_ss)
        P1_w = P1(w_P1,1:3000); size(P1_w)
        
        P2_js = P2(js_P2,1:3000); size(P2_js)
        P2_ss = P2(ss_P2,1:3000); size(P2_ss)
        P2_w = P2(w_P2,1:3000); size(P2_w)
        
        
        %% Define independent and dependent variables for cross validated regression
        % Assign class labels accoding to classiffication target
        switch  model_labels{label_type}
            
            case 'js_xP'
                P1_class =  [zeros(size(P1_js,1),1)]; size(P1_class)
                P2_class =  [ones(size(P2_js,1),1)]; size(P2_class)
                PP = [P1_js; P2_js]; Ns1 = size(PP)
                
                cfg.k   = 36; % ~10% of sample size
                
            case 'njs_xP'
                P1_class =  [zeros(size(P1_ss,1),1); zeros(size(P1_w,1),1)]; size(P1_class)
                P2_class =  [ones(size(P2_ss,1),1); ones(size(P2_w,1),1)]; size(P2_class)
                PP = [P1_ss; P1_w; P2_ss; P2_w]; Ns2 = size(PP)      
                
                cfg.k   = 12; % ~10% of sample size    
        end
        
        % Concatenate class labels
        Y = [P1_class; P2_class]; size(Y);
        
        %% Cross-validated regression (optional: Multi-class decoding with stratified k-folds)
        if norm ==1
            PP = (PP - nanmean(PP,1))./nanstd(PP,[],1); size(X)
            find(isnan(PP))
        end        
        
        c = 0;
        tic
        for i_tp = 1:length(steps)
            c = c+1;
            X = PP(:,steps(i_tp):steps(i_tp)+step-1);
            y = Y;
            [perfm, result_reg] = mv_regress(cfg, X, y); % time series  to labels 
            %stat_permutation = mv_statistics(cfg_stat, result_reg, X, y);
            perf(1,i_sig,i_label,c) = result_reg.perf_std{2}/length(Y); % SEM
            perf(2,i_sig,i_label,c) = perfm{2}; %
            %perf(3,i_sig,i_label,c) = stat_permutation.mask;
        end
        toc   
    end
end


%% Plot results Mean Squared Error (MSE)
load crossValP1xP2RegH0

% JS
H0_005_js = rndperf(:,2,1,1,:);
H0_005_js = H0_005_js(:);
H0_01_js   = rndperf(:,2,2,1,:);
H0_01_js   = H0_01_js(:);
lcutoff_005_js = prctile(H0_005_js,25)
lcutoff_01_js   = prctile(H0_01_js,25)

% non-JS
H0_005_njs = rndperf(:,2,1,2,:);
H0_005_njs = H0_005_njs(:);
H0_01_njs   = rndperf(:,2,2,2,:);
H0_01_njs   = H0_01_njs(:);
lcutoff_005_njs = prctile(H0_005_njs,25)
lcutoff_01_njs   = prctile(H0_01_njs,25)

%% panels 1 and 2 are for participants as dependent variable

figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [(pos(1)+0-1)*1000 0 800,600]);
ca = [0 0.4470 0.7410]; 
cb = [0.6350 0.0780 0.1840]; 
tps   = linspace(-peri_click_time,peri_click_time,size(perf,4));

%% JS
subplot 211
yyaxis left
y = squeeze(perf(2,1,1,:));
ysem = squeeze(perf(1,1,1,:));
errorbar(tps,y,ysem*1.96,'Color',ca,'LineWidth', 2);
hold on

yh = 0.2523; 
noChance = y<=lcutoff_005_js;
f = find(diff([0,noChance',0]==1));
p = f(1:2:end-1); 
yy = f(2:2:end)-p;  

    for n = 1:length(p)
        line([tps(p(n)),tps((p(n)+yy(n))-1)],[yh,yh],'LineWidth', 3,'Color',ca);
    end

ylabel('MSE','FontSize', 12,'Color',ca);
yyaxis right
y = squeeze(perf(2,2,1,:));
ysem = squeeze(perf(1,2,1,:));
errorbar(tps,y,ysem*1.96,'Color',cb,'LineWidth', 2);
hold on

yh = 0.25213; 
noChance = y<=lcutoff_01_js;
f = find(diff([0,noChance',0]==1));
p = f(1:2:end-1); 
yy = f(2:2:end)-p;  

    for n = 1:length(p)
        line([tps(p(n)),tps((p(n)+yy(n))-1)],[yh,yh],'LineWidth', 3,'Color',cb);
    end

ylabel('MSE','FontSize', 12,'Color',cb);
set(findobj(gcf,'type','axes'),'FontName','Arial', 'FontSize', 12, 'LineWidth', 1.5);
set(gca, 'box', 'off');
set(gcf, 'color', 'w');
grid on
xlabel('Click onset ({\itsec})','FontSize', 12);
ax = gca;
ax.YAxis(1).Color ='k';
ax.YAxis(2).Visible = 'off';

%% non-JS
subplot 212
yyaxis left
y = squeeze(perf(2,1,2,:));
ysem = squeeze(perf(1,1,2,:));
errorbar(tps,y,ysem*1.96,'Color',ca,'LineWidth', 2);
hold on

yh = 0.248; 
noChance = y<=lcutoff_005_njs;
f = find(diff([0,noChance',0]==1));
f = find(diff([0,noChance',0]==1));

p = f(1:2:end-1);  
yy = f(2:2:end)-p; 
    for n = 1:length(p)
        line([tps(p(n)),tps((p(n)+yy(n))-1)],[yh,yh],'LineWidth', 3,'Color',ca);
    end

ylabel('MSE','FontSize', 12,'Color',ca);
yyaxis right
y = squeeze(perf(2,2,2,:));
ysem = squeeze(perf(1,2,2,:));
errorbar(tps,y,ysem*1.96,'Color',cb,'LineWidth', 2);
hold on

yh = 0.245; 
noChance = y<=lcutoff_01_njs;
f = find(diff([0,noChance',0]==1));
p = f(1:2:end-1);  
yy = f(2:2:end)-p; 

    for n = 1:length(p)
        line([tps(p(n)),tps((p(n)+yy(n))-1)],[yh,yh],'LineWidth', 3,'Color',cb);
    end
    
ylabel('MSE ','FontSize', 12,'Color',cb);
set(findobj(gcf,'type','axes'),'FontName','Arial', 'FontSize', 12, 'LineWidth', 1.5);
set(gca, 'box', 'off');
set(gcf, 'color', 'w');
grid on
xlabel('Click onset ({\itsec})','FontSize', 12);
ax = gca;
ax.YAxis(1).Color ='k';
ax.YAxis(2).Visible = 'off';
legend('boxoff')


%% Print figure
pth_fg  =  '/home/nicolas/Documents/Okinawa/PCE/analyses_nico/figures/';
print(gcf, [ fname '.svg'], '-dsvg', '-r150', '-painters')
print(gcf, [ fname '.png'], '-dpng', '-r150', '-painters')





