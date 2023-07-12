function [QA] = NWB_QA_graphs(nwbAll,units, is_sternberg)
% NWB_QA_graphs Generates behavioral and spike sorting metrics used in the SB
% data release. 
% nwbAll: A cell arrays containing all nwb files in the dataset 
% units: A units struct created using NWB_SB_extractUnits
% is_sternberg: 1 or 0 logical if sternberg is being processed. Will plot
% behavioral metrics if 1. 
% figHandle: The resulting figure that shows all relevant stat tests. 

QA = figure("Visible","on");
sessionCount = length(nwbAll);

if is_sternberg % If Sternberg, give behavioral plots
    %% Accuracy
    accTotals = zeros(sessionCount,1);
    for i = 1:length(accTotals)
        accTrials = nwbAll{i}.intervals_trials.vectordata.get('response_accuracy').data.load();
        accTotals(i) = sum(accTrials)/length(accTrials)*100;
    end
    accFig = subplot(2,5,5); %#ok<NASGU>
    plot(sort(accTotals),'k.','MarkerSize',10); hold on
    
    % SEM bars
    acc_sem = std(accTotals)/sqrt(length(accTotals));
    plot(1:length(accTotals), zeros(length(accTotals),1) + mean(accTotals),'b','LineWidth',1.5)
    plot(1:length(accTotals), zeros(length(accTotals),1) + mean(accTotals)-acc_sem,'b','LineWidth',1)
    plot(1:length(accTotals), zeros(length(accTotals),1) + mean(accTotals)+acc_sem,'b','LineWidth',1)
    
    title(sprintf('Accuracy: \\mu=%.2f SE=%.2f',mean(accTotals),acc_sem))
    xlabel('Session no.')
    ylabel('Accuracy (%)')
    xlim([0 length(accTotals)+1])
    ylim([min(accTotals)-5 100])
    if sessionCount>1; xticks([1 sessionCount]); end
    set(gca,'FontSize',13)
    hold off;
    
    %% RT vs Load
    filter_correct = 1;
    RTs_all = []; % Compiling raw RTs; {subject ID, load, RT}
    for i = 1:sessionCount
        session_Loads = double(nwbAll{i}.intervals_trials.vectordata.get('loads').data.load());
        tsResp = nwbAll{i}.intervals_trials.vectordata.get('timestamps_Response').data.load();
        tsProbe = nwbAll{i}.intervals_trials.vectordata.get('timestamps_Probe').data.load();
        if filter_correct % Filtering out incorrect trials
            accuracyBool = logical(nwbAll{i}.intervals_trials.vectordata.get('response_accuracy').data.load());
            session_Loads = session_Loads(accuracyBool);
            tsResp = tsResp(accuracyBool);
            tsProbe = tsProbe(accuracyBool);
        end
        session_RTs = tsResp-tsProbe; 
        session_Num = NaN(length(session_RTs),1); session_Num(:) = i;
        rtMat = [session_Num,session_Loads,session_RTs];
        RTs_all = vertcat(RTs_all,rtMat); %#ok<AGROW> % Grows, but doesn't scale drastically. 
    end
    
    % Median RTs
    anovanSet = double.empty(0,3);
    RTmedians_all = zeros(sessionCount,3);
    for i = 1:sessionCount % Accesses RT column in 3rd row
        RTs_session = RTs_all(RTs_all(:,1)==i,:);
        RT_1 = median(RTs_session(RTs_session(:,2)==1,3));
        RT_2 = median(RTs_session(RTs_session(:,2)==2,3));
        RT_3 = median(RTs_session(RTs_session(:,2)==3,3));
        RTmedians_all(i,:) = [RT_1, RT_2, RT_3];
        anovanSet = vertcat(anovanSet,horzcat([i i i]',[1 2 3]',RTmedians_all(i,:)')); %#ok<AGROW>
    end
    
    % ANOVAN 
    grp = { anovanSet(:,1), anovanSet(:,2)};
    [p_anova,tbl_anova_acc,~] = anovan( anovanSet(:,3), grp, 'random',1,'varnames',{'sessionID','Load'}, 'model','linear','display','off');
    
    % Plot RTs
    figRT = subplot(2,5,10); %#ok<NASGU>
    for i = 1:size(RTmedians_all,1)
        plot([1 2 3],RTmedians_all(i,:),'k.-.','MarkerSize',10); hold on
    end
    title('Median RT')
    xlabel(sprintf('Load | F_{%d,%d}: %.4f | p: %.4f',... 
        tbl_anova_acc{3,3},tbl_anova_acc{4,3},tbl_anova_acc{3,6},tbl_anova_acc{3, 7}))
    ylabel('Reaction Time (s)')
    xlim([.5 3.5])
    xticks([1 2 3])
    ylim([min(RTmedians_all(:))-0.25 max(RTmedians_all(:))+0.25])
    yticks([0.5 1 1.5 2])
    
    % SEM bars
    for i = 1:size(RTmedians_all,2)
        semRT = std(RTmedians_all(:,i))/sqrt(length(RTmedians_all(:,i)));
        offset = 0.25;
        plot([i-offset i+offset], [mean(RTmedians_all(:,i)) mean(RTmedians_all(:,i))],'b-','LineWidth',1.5)
        plot([i-offset i+offset], [mean(RTmedians_all(:,i)) mean(RTmedians_all(:,i))] + semRT,'b-','LineWidth',1)
        plot([i-offset i+offset], [mean(RTmedians_all(:,i)) mean(RTmedians_all(:,i))] - semRT,'b-','LineWidth',1)
    end
    set(gca,'FontSize',13)
    hold off;
end
%% Spike Sorting Metrics
% Aggregating timestamps across session types. 
total_ts = {units(:).spike_times}; 
total_chan = [units(:).subject_id;units(:).electrodes]';

ISI_sub3 = NaN(length(total_ts),1);
mean_rate = NaN(length(total_ts),1);
CV2s = NaN(length(total_ts),1);
% ts Metrics
for i = 1:length(total_ts) % Summary Over all sessions
    ts = total_ts{i}; % Import & Offset
    % ISI
    ts_isi = ts.*1e3; % Convert to ms
    ISIs = diff(ts_isi);
    ISI_sub3(i) = (sum(ISIs<3)/length(ISIs))*100;
    % Hz
    ts_hz = ts;
    mean_rate(i) = length(ts_hz)/(max(ts_hz)-min(ts_hz));  
    % CV2
    ts_cv2 = ts;
    isi_cv2 = diff(ts_cv2);
    ignoreMode = 1;
    [CV2, ~, ~, ~] = calcCV2(isi_cv2, ignoreMode);
    CV2s(i) = CV2;
end

% Channel Metrics
chanFreq = [];
for i = 1:length(unique(total_chan(:,1)))
    subj_cells = total_chan(total_chan(:,1)==i,:);
    [~, n] = RunLength(subj_cells(:,2));
    chanFreq = [chanFreq; n]; %#ok<AGROW> % For short term analysis
end

% Waveform Metrics
iso_dist = [];
mean_snr = [];
peak_snr = [];
proj_dist = [];
for i = 1:length(nwbAll)
    iso_dist = [iso_dist; nwbAll{i}.units.vectordata.get('waveforms_isolation_distance').data.load()]; %#ok<AGROW>
    mean_snr = [mean_snr; nwbAll{i}.units.vectordata.get('waveforms_mean_snr').data.load()]; %#ok<AGROW>
    peak_snr = [peak_snr; nwbAll{i}.units.vectordata.get('waveforms_peak_snr').data.load()]; %#ok<AGROW>
    proj_dist = [proj_dist; nwbAll{i}.units.vectordata.get('waveforms_mean_proj_dist').data.load()]; %#ok<AGROW>
end
iso_dist = rmmissing(iso_dist); % Removing invalid NaN values. 

%% Plotting

% colors: K2017 Blue, 
colors = {'#02009B'};
% 
% QA = figure(100);
%  % Doing some window wizardry to move the plot to the left screen. No idea
%  % how this works. Found  by trial and error. 
% QA.WindowState = 'maximized';
% movegui(QA,'west')
% QA.WindowState = 'maximized';


% Plot: Wires Per Num Units
wPuPlot = subplot(2,5,1);
uPw = histogram(chanFreq,'FaceColor',colors{1},'FaceAlpha',1);
% uPw = histogram(chanFreq,'BinMethod','integer','FaceColor',colors{1},'FaceAlpha',1);
uPw.BinWidth = 0.5;
uPw.BinEdges = [uPw.BinWidth uPw.BinEdges] + uPw.BinWidth/2;
title(sprintf('Units per wire: \\mu=%.2f \\sigma=%.2f',mean(chanFreq),std(chanFreq)))
xlabel('Nr of units recorded')
ylabel('Nr of wires')
set(gca,'FontSize',13)
% xlim([0 max(chanFreq)+1])
xticks(1:(max(chanFreq)))

% Plot: ISI
isiPlot = subplot(2,5,2); %#ok<NASGU>
thresh_isi = [4]; %#ok<NBRAK2>
if ~isempty(thresh_isi); ISI_sub3_filtered = ISI_sub3(ISI_sub3<thresh_isi); else; ISI_sub3_filtered = ISI_sub3; end
histogram(ISI_sub3_filtered,40,'FaceColor',colors{1},'FaceAlpha',1)
xlim([0 4])
title(sprintf('ISI<3%%: \\mu=%.2f \\sigma=%.2f',mean(ISI_sub3_filtered),std(ISI_sub3_filtered)))
xlabel('Percent of Interspike Intervals (ISI) < 3 ms')
ylabel('nr Units')
set(gca,'FontSize',13)

% Plot: Hz 
hzPlot = subplot(2,5,3); %#ok<NASGU>
thresh_Hz = [];
if ~isempty(thresh_Hz); mean_rate_filtered = ISI_sub3(mean_rate < thresh_Hz); else; mean_rate_filtered = mean_rate; end
histogram(mean_rate_filtered,40,'FaceColor',colors{1},'FaceAlpha',1)
xlim([0 40])
title(sprintf('Hz: \\mu=%.2f \\sigma=%.2f',mean(mean_rate_filtered),std(mean_rate_filtered)))
xlabel('Firing Rate (Hz)')
ylabel('nr Units')
set(gca,'FontSize',13)

% Plot: CV2
cv2Plot = subplot(2,5,4); %#ok<NASGU>
thresh_cv2 = [];
if ~isempty(thresh_cv2); CV2s_filtered = CV2s(CV2s < thresh_cv2); else; CV2s_filtered = CV2s; end
histogram(CV2s_filtered,40,'FaceColor',colors{1},'FaceAlpha',1)
% cdfplot(CV2s_filtered)
title(sprintf('CV2: \\mu=%.2f \\sigma=%.2f',mean(CV2s),std(CV2s)))
xlabel('CV2')
xlim([0 2])
ylabel('nr Units')
set(gca,'FontSize',13)
% set(gca,'XScale','log')

% Plot: Waveform Peak SNR
subplot(2,5,6)

histogram(peak_snr,40,'FaceColor',colors{1},'FaceAlpha',1)
xlim([0 max(peak_snr)])

title(sprintf('Peak SNR: \\mu:%.2f \\sigma:%.2f',mean(peak_snr),std(peak_snr)))
xlabel('Waveform peak SNR')
ylabel('nr Units')
set(gca,'FontSize',13)

% Plot: Waveform Mean SNR
subplot(2,5,7)

histogram(mean_snr,40,'FaceColor',colors{1},'FaceAlpha',1)
xlim([0 max(mean_snr)])

title(sprintf('Mean SNR: \\mu:%.2f \\sigma:%.2f',mean(mean_snr),std(mean_snr)))
xlabel('Waveform Mean SNR')
ylabel('nr Units')
set(gca,'FontSize',13)

% Plot: Pairwise Distance (Projection Test)
subplot(2,5,8)
proj_dist_plot = nonzeros(proj_dist);
histogram(proj_dist_plot,40,'FaceColor',colors{1},'FaceAlpha',1)
xlim([0 50])

title(sprintf('Mean Proj: \\mu:%.2f \\sigma:%.2f',mean(proj_dist_plot),std(proj_dist_plot)))
set(gca,'FontSize',13)
xlabel('Projection Distance (s.d)')
ylabel('nr Units')


% Plot: Isolation Distance (log10)
subplot(2,5,9)
iso_dist_log10 = log10(iso_dist);
iso_plot = iso_dist_log10;
histogram(iso_plot,30,'FaceColor',colors{1},'FaceAlpha',1)
xlim([0 4])
title(sprintf('Iso Dist: \\mu:%.2f \\sigma:%.2f',mean(iso_plot),std(iso_plot)))

xlabel('Isolation Distance (log 10)')
ylabel('nr Units')
set(gca,'FontSize',13)

% % Annotations
% annotation(gcf,'textbox',...
%     [0.100194150380022 0.93552636488386 0.0135833333333333 0.0209163346613544],...
%     'String',{'a)'},...
%     'FontSize',13,...
%     'FitBoxToText','off',...
%     'LineStyle','none');
% 
% annotation(gcf,'textbox',...
%     [0.312902483713355 0.941540933042034 0.0135833333333333 0.0209163346613544],...
%     'String',{'b)'},...
%     'FontSize',13,...
%     'FitBoxToText','off',...
%     'LineStyle','none');
% 
% annotation(gcf,'textbox',...
%     [0.518319150380021 0.938190256663265 0.0135833333333333 0.0209163346613544],...
%     'String',{'c)'},...
%     'FontSize',13,...
%     'FitBoxToText','off',...
%     'LineStyle','none');
% 
% annotation(gcf,'textbox',...
%     [0.720610817046688 0.932758414831843 0.0135833333333333 0.0209163346613544],...
%     'String',{'d)'},...
%     'FontSize',13,...
%     'FitBoxToText','off',...
%     'LineStyle','none');
% 
% annotation(gcf,'textbox',...
%     [0.103631650380021 0.456983180700732 0.0135833333333333 0.0209163346613544],...
%     'String',{'e)'},...
%     'FontSize',13,...
%     'FitBoxToText','off',...
%     'LineStyle','none');
% 
% annotation(gcf,'textbox',...
%     [0.309048317046688 0.456754252500942 0.0135833333333333 0.0209163346613544],...
%     'String',{'f)'},...
%     'FontSize',13,...
%     'FitBoxToText','off',...
%     'LineStyle','none');
% 
% annotation(gcf,'textbox',...
%     [0.520714983713354 0.459647072480132 0.0135833333333333 0.0209163346613544],...
%     'String',{'g)'},...
%     'FontSize',13,...
%     'FitBoxToText','off',...
%     'LineStyle','none');
% 
% annotation(gcf,'textbox',...
%     [0.720610817046688 0.459647072480132 0.0135833333333333 0.0209163346613544],...
%     'String',{'h)'},...
%     'FontSize',13,...
%     'FitBoxToText','off',...
%     'LineStyle','none');

%% Long-form stats

% test = struct([]);
% 
% for i=1:length(total_ts)
% [f,Pxxn,tvect,Cxx,edges1,n1,yGam1,edges2,n2,yGam2,mGlobal,m1,m2,percentageBelow,CV] = getStatsForCluster([],total_ts{i});
% 
% test(i).f = f;
% test(i).Pxxn = Pxxn;
% test(i).tvect = tvect;
% test(i).Cxx = Cxx;
% test(i).edges1 = edges1;
% test(i).n1 = n1;
% test(i).yGam1 = yGam1;
% test(i).edges2 = edges2;
% test(i).n2 = n2;
% test(i).yGam2 = yGam2;
% test(i).mGlobal = mGlobal;
% test(i).m1 = m1;
% test(i).m2 = m2;
% test(i).percentageBelow = percentageBelow;
% test(i).CV = CV;
% end
% 
% CV_test = [test(:).CV].^2;
% thresh_cv2_test = [2];
% if ~isempty(thresh_cv2_test); CV_test_filtered = CV_test(CV_test < thresh_cv2_test); else; CV_test_filtered = CV2s; end
% 
% figure(10)
% histogram(CV_test_filtered,100)
% title(sprintf('CV2 (getStatsForCluster): \\mu=%.2f \\sigma=%.2f',mean(CV_test),std(CV_test)))
% xlim([0 2])
end

























