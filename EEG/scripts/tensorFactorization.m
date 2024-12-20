close all
rng(111) % to get NNTF to work, must need to set the random generator to a fixed seed

group = 0;
do_ica = 1;

if do_ica == 1
    load  ECSU_restTask_ica
else
    load  ECSU_restTask
end

%% Non-Negative Tensor Factorization (a low-rank approximation of tensor factors)
%% Define non-negative spectrogram tensor across subjects
%for group = 1:2
%%
switch group
    case 0
        S = cat(1,squeeze(TFR.spectrum(:,1,:,:)) , squeeze(TFR.spectrum(:,2,:,:)));
        label = 'all';
    case 1
        S = squeeze(TFR.spectrum(:,1,:,:));
        label = 'G1';
    case 2
        S = squeeze(TFR.spectrum(:,2,:,:));
        label = 'G2';
end

if do_ica == 1
   label = [label '_ica'];
end

tr = TFR.spectrum_time;
fr = TFR.spectrum_freqs;

dim = size(S)
clear P P_norm
for i_s = 1:dim(1)
    for i_t = 1:dim(2)
        for i_f = 1:dim(3)
            %P(i_s,i_t,i_f)  =(S(i_s,i_t,i_f)-mean(S(i_s,:,i_f),2))./std(S(i_s,:,i_f),1);
            %P(i_s,i_t,i_f)  =S(i_s,i_t,i_f); %(S(i_s,i_t,i_f)-mean(S(i_s,:,i_f),2))./std(S(i_s,:,i_f),1);
            P(i_s,i_t,i_f)  = (S(i_s,i_t,i_f)-mean(S(i_s,:,i_f),2))./std(S(i_s,:,i_f),1);          
        end
    end
    for i_f = 1:dim(1)
        x = P(i_s,:,:);
        mindata = min(x(:));
        maxdata = max(x(:));
        P_norm(i_s,:,:) = bsxfun(@rdivide, bsxfun(@minus, P(i_s,:,:), mindata), maxdata - mindata);
    end
end
T = fmt(P_norm);

%T = fmt(S);

%% Factor priors
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2)-800 800, 200]); % Set plot size
x = linspace(1,size(T, 3),size(T, 3));
mu = 8; %10.;
sig = 5;
y_gastric= gaussmf(x,[sig,mu]);
mu = 19; %20.;
sig = 5;
y_mayer= gaussmf(x,[sig,mu]);
plot(fr,y_gastric,'r-','LineWidth', 1.5);
hold on
plot(fr,y_mayer,'b-','LineWidth', 1.5);
plot(fr,[y_gastric + y_mayer],'k--','LineWidth', 2);
ylabel('Factor priors');
xlabel('Frequency (Hz)');
set(gca, 'FontSize',12,'LineWidth', 1.5);
set(gca, 'box', 'off');
set(gcf, 'color', 'w');
grid on
ylim([0 1.25]);
%title('Frequency tensor factorization')
legend( '0.05 {\itHz}','0.1 {\itHz}','both','Location','northeastoutside');
legend('boxoff');
fname = [ pth_fg  'NNTF_frequency_priors.png']
print(gcf, fname, '-dpng', '-r150', '-painters')
fname = [ pth_fg  'NNTF_frequency_priors.svg']
print(gcf, fname, '-dsvg', '-r150', '-painters')



%% NNTF model
factors = 3;
model = struct;
model.variables.A = rand(size(T, 1), factors);
model.variables.B = rand(size(T, 2), factors);
model.variables.C = [y_gastric; y_mayer; [y_gastric + y_mayer]]';% rand(size(T, 3), factors);
model.variables.C = [y_gastric; y_mayer;  0.5*ones(1,38)]';% rand(size(T, 3), factors);
%model.variables.C = [y_gastric; y_mayer;   rand(1,40)]';% rand(size(T, 3), factors);
%model.variables.C = rand(size(T, 3), factors);
model.transform = {1:factors, @struct_nonneg};
model.factorizations{1}.data = T;
model.factorizations{1}.cpd = 1:factors;
%[sol,output] = sdf_nls(model, 'Display', 10000, 'TolX', eps, 'TolFun', eps^2);
[sol,output] = sdf_nls(model, 'Algorithm',@nls_gndl,'JHasFullRank', true,'Display', 10, 'TolX', eps, 'TolFun', eps^2,...
    'CGMaxIter', 100); % Nonlinear least squares type NNTF algorithm
factor_subjects    = sol.factors{1};
factor_time           = sol.factors{2};
factor_frequency = sol.factors{3};

%% Subject/trial factors
figure,
colormap('Copper')
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2)-600 1200, 120]); % Set plot size
imagesc(factor_subjects');
c = colorbar
set(get(c,'title'),'String','Trial weight');
xlabel([ num2str(size(S,1)) ' repetitions ({\it2 rest-task trials across 31 experiments with 2 participants})'] );
ylabel('Factor');
set(gca, 'FontSize',12,'LineWidth', 1.5);
set(gca, 'box', 'off');
set(gcf, 'color', 'w');
%title('Subject tensor factorization')
fname = [ pth_fg  'NNTF_rep_' label '.svg']
print(gcf, fname, '-dsvg', '-r150', '-painters')
fname = [ pth_fg  'NNTF_rep_' label '.png']
print(gcf, fname, '-dpng', '-r150', '-painters')

%%
if group ==1
    factors_G1 = factor_subjects;
end
if group == 2
    factors_G2 = factor_subjects;
end

%% Frequency
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2)-400 300, 250]); % Set plot size
%set(gcf, 'Position', [pos(1) pos(2)-800 800, 200]); % Set plot size
plot(fr,factor_frequency(:,1),'r-.','LineWidth', 2);
hold on
plot(fr,factor_frequency(:,2),'b-.','LineWidth', 2);
plot(fr,factor_frequency(:,3),'k-.','LineWidth', 2);

ylabel('Weight');
xlabel('Frequency ({\it Hz})');
set(gca, 'FontSize',12,'LineWidth', 1.5);
set(gca, 'box', 'off');
set(gcf, 'color', 'w');
grid on
%ylim([0 0.5]);
%title('Frequency tensor factorization')
%legend('Peristaltic rhythm','Mayer waves','Noise','Location','northeastoutside');
%legend( '0.05 {\itHz}','0.1 {\itHz}','Broadband','Location','northeastoutside');
%legend('0.05 {\itHz}','0.1 {\itHz}','{\it Broadband}','Task onset','Orientation','horizontal', 'Location','SouthOutside');
%legend('0.05 {\itHz}','0.1 {\itHz}','{\it Broadband}','{\it Grand average}','Task onset','Orientation','horizontal', 'Location','SouthOutside');
legend('boxoff');
fname = [ pth_fg  'NNTF_freq_' label '.svg']
print(gcf, fname, '-dsvg', '-r150', '-painters')
fname = [ pth_fg  'NNTF_freq_' label '.png']
print(gcf, fname, '-dpng', '-r150', '-painters')

%% Time
figure,
colormap('Winter')
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2)+200 1200, 350]); % Set plot size
x = linspace(-180-tr(1), tr(end)-180, length(tr));
x = tr-180
plot(x,factor_time(:,1),'r-.','LineWidth', 1.5);
hold on
plot(x,factor_time(:,2),'b-.','LineWidth', 1.5);
plot(x,factor_time(:,3),'k-.','LineWidth', 2);
%line([180 180],[min(factor_time(:)) max(factor_time(:))],'Color',[0 0 0],'LineWidth',2,'LineStyle','-.','HandleVisibility','off')
switch group
    case 0
        onset = round(TFR.blank_onset(:)/10)
    case 1
        onset = round(TFR.blank_onset(:,1)/10)
    case 2
        onset = round(TFR.blank_onset(:,2)/10)
end
%histogram(onset ,100,'facealpha',.7,'edgecolor','none');
data = onset-180; % replace this with your own data
scalar = 0.1;%your_scalar; % replace this with your desired scalar
bin_edges = linspace(min(data), max(data), 50); % define bin edges
bin_counts = histcounts(data, bin_edges); % calculate bin counts
bin_counts_scaled = bin_counts * scalar; % scale the bin counts
b = bar(bin_edges(1:end-1), bin_counts_scaled, 'hist','facealpha',0.5,'edgecolor','none'); % plot the histogram
b.EdgeAlpha = 0;
b.FaceColor = [0.6290 0.3940 0.0250];
b.FaceAlpha = 0.5;
b.HandleVisibility ='off'
line([0 0],[-max(factor_time(:))*5 max(factor_time(:)*5)],'Color',[0.3 0.3 0.3],'LineStyle',':', 'LineWidth',2,'HandleVisibility','off')
ylim([0 5]);
ylabel('Weight');
xlabel('Time ({\it Seconds})');
yyaxis right
x = linspace(-180,60,2400);
clear y
do_ica = 0;
for i_rep =1:62
    switch group
        case 0
            %y(i_rep,:) =  sum([TFR.eeg_clean{i_rep,1} TFR.eeg_clean{i_rep,2}],2);
            if do_ica ==1
                y_1(i_rep,:) =   mean(TFR.eeg_ica{i_rep,1},2);
                y_2(i_rep,:) =   mean(TFR.eeg_ica{i_rep,2},2);
            else
                y_1(i_rep,:) =   mean(TFR.eeg_clean{i_rep,1},2);
                y_2(i_rep,:) =   mean(TFR.eeg_clean{i_rep,2},2);
            end
        case 1
            y(i_rep,:) =  sum(TFR.eeg_clean{i_rep,1},2);
        case 2
            y(i_rep,:) =  sum(TFR.eeg_clean{i_rep,2},2);
    end
end
if group ==0
    y = cat(1,y_1,y_2);
    size(y)
    %y_avg = -1*sum(y,1)/124; size(y_avg)
    y_avg = -1*mean(y,1); size(y_avg)
else
    size(y)
    y_avg = mean(y,1); size(y_avg)
end
y_avg = zscore(smoothdata(y_avg,'movmedian',10));
plot(x,y_avg,'k-','LineWidth', 1.5);
ylabel('Amplitude ({\it a.u.})')
set(gca,'YColor','k')
set(gca, 'FontSize',12,'LineWidth', 1.5);
set(gca, 'box', 'off');
set(gcf, 'color', 'w');
grid on
%ylim([0 3]);
if do_ica ==1
    ylim([-0.5 +0.5]);
    
    
end
xlim([-120 60]);
%title('Time tensor factorization')

legend('0.05 {\itHz}','0.1 {\itHz}','{\it Broadband}','{\it Grand average}','Task onset','Orientation','horizontal', 'Location','SouthOutside');
legend('boxoff')

fname = [ pth_fg  'NNTF_time_' label '.svg']
print(gcf, fname, '-dsvg', '-r150', '-painters')
fname = [ pth_fg  'NNTF_time_' label '.png']
print(gcf, fname, '-dpng', '-r150', '-painters')


%end

indices = zeros(62,3);
for i_rep =1:62
 indices(i_rep,:,:,:) = TFR.repetition_index{i_rep};
end