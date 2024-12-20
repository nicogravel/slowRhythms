%% Copuling of click types to the phase of infra-slow EEG rhythms
% December 2, 2024
% nicolas.gravel@gmail.com

%close all
clear all

addpath(genpath('circstat-matlab-master')); 
clickResponse_classification

fmax = 0.3;
nBins = 8;

i_sig = 0;
for is_signal = 4:5
    
    i_sig = i_sig +1;
    
    
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
            fname = [ pth_fg  'PCC_alpha'];
            
        case 3
            [idx,P] = find(idx_eeg_trials & idx_eda_trials & idx_resp_trials);
            signal = epoch_resp;
            titl = 'Respiratory pressure';
            fname = [ pth_fg  'PCC_Respiration'];
            
        case 4
            [idx,P] = find(idx_eeg_trials & idx_eda_trials );
            signal  = epoch_phi_005;
            titl = 'EEG (0.05 {\itHz}) ';
            fname = [ pth_fg  'PCC_005'];
            
        case 5
            [idx,P] = find(idx_eeg_trials & idx_eda_trials );
            signal  = epoch_phi_01;
            titl = 'EEG (0.1 {\itHz}) ';
            fname = [ pth_fg  'PCC_01'];
            
        case 6
            [idx,P] = find(idx_eeg_trials & idx_eda_trials & idx_resp_trials);
            signal  = epoch_eda;
            titl = 'EDA ';
            fname = [ pth_fg  'PCC_EDA'];
            
        case 7
            [idx,P] = find(idx_eeg_trials & idx_eda_trials & idx_resp_trials);
            signal  = epoch_eeg;
            titl = 'EEG';
            fname = [ pth_fg  'PCC_EEG'];
            
        case 8
            [idx,P] = find(idx_ecg_trials);
            signal  = epoch_ecg_HRV;
            %minmax  = 1;
            titl = 'HRV';
            fname = [ pth_fg  'PCC_HRV'];
            
        case 9
            [idx,P] = find(idx_ecg_trials);
            signal  = epoch_ecg_rrx;
            titl = 'RR';
            fname = [ pth_fg  'PCC_RR'];
            
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
    
    
    
    % Select click type specific data
    P1_js = P1(js_P1,1:3000); size(P1_js)
    P1_ss = P1(ss_P1,1:3000); size(P1_ss)
    P1_w = P1(w_P1,1:3000); size(P1_w)
    
    P2_js = P2(js_P2,1:3000); size(P2_js)
    P2_ss = P2(ss_P2,1:3000); size(P2_ss)
    P2_w = P2(w_P2,1:3000); size(P2_w)
    
    
    
    %%
    figure,
    
    step_length = 360/nBins;
    phase_bins = 0:step_length:360;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 1800, 350]); % Set plot size
    gray = [128 128 128]/255;
    
    %% JS
    subplot 151
    grid on; hold on;
    PP = [P1_js; P2_js]; Ns1 = size(PP)
    phase = PP(:,1500); %rad2deg(wrapToPi(PP(:,1500)))
    response = 180+phase*360/2/pi; % Phases in degrees
    [Ns]=histc(response,phase_bins);
    relativefreqNs = Ns ./ sum(Ns);
    b(1) = bar(phase_bins ,relativefreqNs,'histc');
    set(b(1),'FaceColor', gray,'EdgeColor', 'w','facealpha',1,'edgealpha',0.5);
    xlim([phase_bins(1)-5 phase_bins(end)+5]);
    ylim([0 fmax]);
    xlabel(['Phase ({\itdeg})'],'FontSize', 10);
    ylabel('Relative frequency','FontSize', 10);
    title('JS');
    set(gca, 'FontSize', 10);
    set(gca,'LineWidth',2)
    set(gca, 'box', 'off');
    axis square
    set(gcf, 'color', 'w');
    set(findobj(gcf,'type','axes'),'FontName','Arial', 'FontSize', 10, 'LineWidth', 1.5);
    
    %% SS
    subplot 152
    grid on; hold on;
    PP = [P1_ss; P2_ss;]; Ns2 = size(PP)
    phase = PP(:,1500); %rad2deg(wrapToPi(PP(:,1500)))
    response = 180+phase*360/2/pi; % Phases in degrees
    [Ns]=histc(response,phase_bins);
    relativefreqNs = Ns ./ sum(Ns);
    b(2) = bar(phase_bins ,relativefreqNs,'histc');
    set(b(2),'FaceColor', gray,'EdgeColor', 'w','facealpha',1,'edgealpha',0.5);
    xlim([phase_bins(1) phase_bins(end)]);
    ylim([0 fmax]);
    xlabel(['Phase ({\itdeg})'],'FontSize', 10);
    ylabel('Relative frequency','FontSize', 10);
    title('SS');
    set(gca, 'FontSize', 10);
    set(gca,'LineWidth',2)
    set(gca, 'box', 'off');
    axis square
    set(gcf, 'color', 'w');
    set(findobj(gcf,'type','axes'),'FontName','Arial', 'FontSize', 10, 'LineWidth', 1.5);
    
    %% W
    subplot 153
    grid on; hold on;
    PP = [P1_w; P2_w]; Ns3 = size(PP)
    phase = PP(:,1500); %rad2deg(wrapToPi(PP(:,1500)))
    response = 180+phase*360/2/pi; % Phases in degrees
    [Ns]=histc(response,phase_bins);
    relativefreqNs = Ns ./ sum(Ns);
    b(3) = bar(phase_bins ,relativefreqNs,'histc');
    set(b(3),'FaceColor', gray,'EdgeColor', 'w','facealpha',1,'edgealpha',0.5);
    xlim([phase_bins(1) phase_bins(end)]);
    ylim([0 fmax]);
    xlabel(['Phase ({\itdeg})'],'FontSize', 10);
    ylabel('Relative frequency','FontSize', 10);
    title('W');
    set(gca, 'FontSize', 10);
    set(gca,'LineWidth',2)
    set(gca, 'box', 'off');
    axis square
    set(gcf, 'color', 'w');
    set(findobj(gcf,'type','axes'),'FontName','Arial', 'FontSize', 10, 'LineWidth', 1.5);
    
    %% nJS (SS+W)
    subplot 154
    grid on; hold on;
    PP = [P1_ss; P1_w; P2_ss; P2_w]; Ns4 = size(PP)
    phase = PP(:,1500); %rad2deg(wrapToPi(PP(:,1500)))
    response = 180+phase*360/2/pi; % Phases in degrees
    [Ns]=histc(response,phase_bins);
    relativefreqNs = Ns ./ sum(Ns);
    b(4) = bar(phase_bins ,relativefreqNs,'histc');
    set(b(4),'FaceColor', gray,'EdgeColor', 'w','facealpha',1,'edgealpha',0.5);
    xlim([phase_bins(1) phase_bins(end)]);
    ylim([0 fmax]);
    xlabel(['Phase ({\itdeg})'],'FontSize', 10);
    ylabel('Relative frequency','FontSize', 10);
    title('non-JS');
    set(gca, 'FontSize', 10);
    set(gca,'LineWidth',2)
    set(gca, 'box', 'off');
    axis square
    set(gcf, 'color', 'w');
    set(findobj(gcf,'type','axes'),'FontName','Arial', 'FontSize', 10, 'LineWidth', 1.5);
    
    %% JS - nJS
    subplot 155
    grid on; hold on;
    PP_JS = [P1_js; P2_js]; Ns1 = size(PP_JS)
    phase = PP_JS(:,1500); %rad2deg(wrapToPi(PP(:,1500)))
    response = 180+phase*360/2/pi; % Phases in degrees
    JS_clicks = deg2rad(response);
    [Ns_JS]=histc(response,phase_bins);
    relativefreqNs_JS = Ns_JS ./ sum(Ns_JS); size(relativefreqNs_JS);
    
    PP_nJS = [P1_ss; P1_w; P2_ss; P2_w]; Ns4 = size(PP)
    phase = PP_nJS(:,1500); %rad2deg(wrapToPi(PP(:,1500)))
     nJS_clicks = deg2rad(response);
    response = 180+phase*360/2/pi; % Phases in degrees
    [Ns_nJS]=histc(response,phase_bins);
    relativefreqNs_nJS = Ns_nJS ./ sum(Ns_nJS); size(relativefreqNs_nJS);
    
    relativefreqNs = (Ns_JS - Ns_nJS)/sum(Ns_JS+Ns_nJS); % -  relativefreqNs_JS - relativefreqNs_nJS;
    

   
    b(5) = bar(phase_bins ,relativefreqNs,'histc');
    set(b(5),'FaceColor', gray,'EdgeColor', 'w','facealpha',1,'edgealpha',0.5);
    xlim([phase_bins(1) phase_bins(end)]);
    ylim([0 fmax]);
    xlabel(['Phase ({\itdeg})'],'FontSize', 10);
    ylabel('Relative frequency','FontSize', 10);
    title('JS {\itminus} non-JS');
    set(gca, 'FontSize', 10);
    set(gca,'LineWidth',2)
    set(gca, 'box', 'off');
    axis square
    set(gcf, 'color', 'w');
    set(findobj(gcf,'type','axes'),'FontName','Arial', 'FontSize', 10, 'LineWidth', 1.5);
    
    %% Print figure
    pth_fg  =  '/home/nicolas/Documents/Okinawa/PCE/analyses_nico/figures/';
    print(gcf, [ fname '.svg'], '-dsvg', '-r300', '-painters')
    print(gcf, [ fname '.png'], '-dpng', '-r300', '-painters')
    
     
    
    %%  Rayleigh test
    nBins = 8;
    step_length = 360/nBins;
    phase_bins = 0:step_length:360;
    PP_JS = [P1_js; P2_js]; Ns1 = size(PP_JS)
    phase = PP_JS(:,1500); %rad2deg(wrapToPi(PP(:,1500)))
    response = 180+phase*360/2/pi; % Phases in degrees
    JS_clicks = deg2rad(response);
    [Ns_JS]=histc(response,phase_bins);
     PP_nJS = [P1_ss; P1_w; P2_ss; P2_w]; Ns4 = size(PP)
    phase = PP_nJS(:,1500); %rad2deg(wrapToPi(PP(:,1500)))
     nJS_clicks = deg2rad(response);
    response = 180+phase*360/2/pi; % Phases in degrees
    [Ns_nJS]=histc(response,phase_bins);
    alpha = deg2rad((phase_bins(1:end-1) + phase_bins(2:end))/2);
    w       =  Ns_JS(1:end-1) - Ns_nJS(1:end-1);  
    [p, z] = circ_rtest(alpha,w);    
    % Display the results
    disp([' Rayleigh test : ', num2str(z)]);
    disp(['p-value : ', num2str(p)]);

    
end







