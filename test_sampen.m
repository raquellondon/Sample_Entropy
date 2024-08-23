%% Calculate SampEn for each trial of Raquel's TOJ data (-1000 - 0 ms)
clear all
close all

readdir = 'C:\Users\rlondon\Dropbox\Projects_Ongoing\Jack\Londonetal_AlphaTOJ_revisionEJN\Revised_Analysis\Data\Prep\';            
writdir = 'C:\Users\rlondon\Dropbox\Projects_Ongoing\Transform_Entropy\Data\Sampen_Value_TOJ\';

cd(readdir);
filelist=dir('*.mat'); 
filelist={filelist.name};

m = 2;

for s = 1:length(filelist) 

    load([readdir filelist{1,s}])
    sampen_value = nan(length(FtEEG.trial),length(FtEEG.label));

    for t = 1:length(FtEEG.trial)
       
        for c = 1:length(FtEEG.label)
            
            data = FtEEG.trial{1,t}(c,1025:2048); %data = randn(1, 1000); 

            % SampEn
            r = 0.2 * std(data);
            sampen_value(t,c) = SampEn(data, m, r);
    
        end
        
    end
    
    save([writdir 'sampen_value_TOJ_' filelist{1,s}],"sampen_value", '-mat')

end

%%

% % SampEn over time averaged over participants and electrodes

clear all
close all

readdir = 'C:\Users\rlondon\Dropbox\Projects_Ongoing\Transform_Entropy\Data\Sampen_Value_TOJ\';
writdir = 'C:\Users\rlondon\Dropbox\Projects_Ongoing\Transform_Entropy\Data\Sampen_Value_TOJ_ALL\';

cd(readdir);
filelist=dir('*.mat'); 
filelist={filelist.name};

sampen_corr_ALL = nan(length(filelist),64);
sampen_mean_ALL = nan(length(filelist),64);

for s = 1:length(filelist) 

    load([readdir filelist{1,s}])

    for c = 1:size(sampen_value,2)
        
        sampen_corr_ALL(s,c) = corr(sampen_value(:,c),linspace(1,size(sampen_value,1),size(sampen_value,1))','type','spearman')';
        sampen_mean_ALL(s,c) = mean(sampen_value(:,c),1);

    end

end

save([writdir 'sampen_value_TOJ_Corr.mat'],"sampen_corr_ALL", '-mat')
save([writdir 'sampen_value_TOJ_Mean.mat'],"sampen_mean_ALL", '-mat')

% SampEn topography averaged over time

load channellocations.mat

% Topography of SampEn averaged over participants and time

figure
topovals = mean(sampen_mean_ALL,1);
cmax = max(max(topovals)) ;
cmin = 0;
clim = [cmin cmax];
topoplot(topovals(:),chanlocs,'maplimits', clim, 'emarker', {'.','k',8,2})
set(gca,'fontsize', 20)
mymap = flipud(brewermap(64,'RdBu'));
colormap(mymap)
colorbar
set(get(colorbar,'label'),'string','SampEn','FontSize',20); 


% Topography of the correlation between time elapsed and SampEn
% (How does SampEn evolve over time per electrode)
% Averaged over  participants

figure
topovals = mean(sampen_corr_ALL,1);
cmax = max(max(topovals)) ;
cmin = min(min(topovals));
lim = max(abs([cmax cmin]));
clim = [-lim lim];
topoplot(topovals(:),chanlocs,'maplimits', clim, 'emarker', {'.','k',8,2})
set(gca,'fontsize', 20)
mymap = flipud(brewermap(64,'RdBu'));
colormap(mymap)
colorbar

%% Run stats

clear all; close all; clc;

% load data
load C:\Users\rlondon\Dropbox\Projects_Ongoing\Transform_Entropy\Data\Sampen_Value_TOJ_ALL\sampen_value_TOJ_Corr.mat;
% The following three things can either be loaded in or can be specified below

load freq_elec_time_label.mat
freqlabel = 1;
timelabel = 1;%using only from -740 ms to -260 ms
% The following two things have to be loaded in 
% load in neighbourhood structure
load biosemi64_neighb.mat;
% load in neighbourhood matrix
load biosemi64_neighb_matrix.mat;
% data to retrieve for each comparison
entrydata = {'sampen_corr_ALL'};
% Names for saving output
savdat = {'sampen_corr_ALL_Stats'};
% Create vector of zeros (same size as entrydata matrices) to test real
% data against (one-sample t-test)
jackzeros = zeros(size(sampen_corr_ALL));

% Now do stats
for co = 1:length(entrydata)
    % create data entry and cfg structures
    powerent.label = eleclabel;
    powerent.dimord = 'subj_chan_freq_time';
    powerent.freq = freqlabel;
    powerent.time = timelabel;
    % !!!!! changes for each comparison !!!!!!!!!!
    powerent.powspctrm = [eval(entrydata{co}); jackzeros];
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cfg = [];
    % specify test statistic 
    %(see https://github.com/fieldtrip/fieldtrip/tree/master/statfun for
    %available statistical tests)
    cfg.statistic = 'ft_statfun_depsamplesT';
    % Specify design matrix
    cfg.design = [[ones(size(sampen_corr_ALL,1),1)', 2*ones(size(sampen_corr_ALL,1),1)'];[1:size(sampen_corr_ALL,1),1:size(sampen_corr_ALL,1)]];
    % Specify which rows of the design matrix contain independent variable
    % (condition) and unit variable (participants).
    cfg.ivar = 1;
    cfg.uvar = 2;
    % Which timepoints to include in analysis
    cfg.latency = 'all';
    % Keep single trials?
    cfg.keeptrials = 'no';
    % average over timepoints? 
    cfg.avgovertime = 'no';
    % Specify neighbourhood structure
    cfg.neighbours = neighbours;
    % method for test statistic calculation
    cfg.method = 'montecarlo';%'stats'
    % specify whether to correct for multiple comparisons and if so how?
    cfg.correctm = 'no';
    % This is the threshold p-value for inclusion of a single datapoint in a
    % cluster
    cfg.clusteralpha = 0.05;
    % statistic calculated from each cluster (maxsum simply sums all
    % t-values in the cluster)
    cfg.clusterstatistic = 'maxsum';
    % minimum number of neighbouring channels to a
    % single significant datapoint that must also be significant for
    % cluster inclusion
    cfg.minnbchan = 2;
    % Set alpha level for second level stats (testing significance of
    % clusters)
    cfg.alpha = 0.05;
    % specify one versus two-sided test (0 = two sided)
    cfg.tail = 0;
    % Corrects for the fact that positive and negative clusters are
    % tested separately. Basically diveds the alpha level by two at each
    % end of the distribution.
    cfg.correcttail = 'alpha';
    % Have critical value returned
    cfg.computecritval = 'yes';
    % Have p-values returned
    cfg.computeprob = 'yes';
    % Number of randomizations for permutation resampling
    cfg.numrandomization = 1000;
    % Specify neighbourhood matrix
    cfg.connectivity = neighb_mat;
    % There are lots of other parameters that can be specified in the cfg
    % structure:
    % see http://www.fieldtriptoolbox.org/reference/ft_freqstatistics
    % Perform statistical tests
    statsstrucjack = ft_freqstatistics(cfg,powerent);
    % Save output
    savename = [savdat{co},'_Corrected']; 
    save(savename, 'statsstrucjack');
    %clear statsstrucjack;
end

load 'channellocations.mat'; 

% Correct for multiple comparisons with FDR

tvals_sampen_corr = statsstrucjack.stat;
[p_fdr_sampen_corr, p_masked_sampen_corr] = fdr(statsstrucjack.prob, .05);

%Plot time-on-task effects for PLE

count = 0;
for chan = 1:size(p_masked_sampen_corr,1)
    if p_masked_sampen_corr(chan) == 1
        count = count + 1;
        markchans(count) = chan;
    end
end

%Plot group level topography of time on task effect on SampEn
figure
topoplot(tvals_sampen_corr,chanlocs,'emarker', {'.','k',8,2}, 'emarker2',{markchans,'o',[1 1 1],8,1})%'maplimits', clim, 
set(gca,'fontsize', 20)
ft_hastoolbox('brewermap', 1);% ensure this toolbox is on the path
mymap = flipud(brewermap(64,'RdBu'));
colormap(mymap)
colorbar
set(get(colorbar,'label'),'string','t-value','FontSize',20); 

%Plot scatter of time on task effect (correlation SampEn - trial order) for each participant averaged
%over electrodes in clusters
Ind_sampen_corr = nan(size(sampen_corr_ALL,1),1);
for s = 1:size(sampen_corr_ALL,1)
    Ind_sampen_corr(s) = mean(sampen_corr_ALL(s,markchans));
end

x = linspace(2,3,size(Ind_sampen_corr,1));
y = Ind_sampen_corr;
figure
h = scatterhist(x,y, 'Kernel', 'overlay', 'Direction','out','NBins',[10],'Color', [0,.5,.7],'Marker','.','MarkerSize',1);
axis([0 5 -0.5286    0.3140])
delete(h(2));
hold on
scatter(x,y,'filled', 'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',[1,1,1],'LineWidth', 1)
scatter(x,y,100,'filled', 'b','MarkerFaceColor',[0,.5,.7],'MarkerEdgeColor',[0,.5,.7],'MarkerFaceAlpha',.3,'LineWidth', 1)
scatter(2.5,mean(y),'k','filled')
xlabel([])
ylabel([])
set(gca,'XTick',[])
set(gca,'fontsize', 12)
set(gcf,'units','points','position',[0,0,230,330])
line([0,6],[0,0], 'color','k','LineWidth',1, 'LineStyle','--')

