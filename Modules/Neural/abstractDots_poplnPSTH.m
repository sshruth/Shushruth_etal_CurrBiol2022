function dataOut = abstractDots_poplnPSTH(eMat,respMat,plotpars)
% function dataOut = analyzePopuln(eMat,respMat,plotpars)
% plotpars: plotting parameters structure including the following fields
%     detrend - true or false : whether to detrend the target onset
%     grouping - all, correct: trial grouping
%     trlcutoff - (fraction) : fraction of trials that have dropped to cut off
%     filtsize - (positive integer): width of boxcar filter to smooth plots
%     monkID - AN or SM: monkey ID

%% Set up
E = eventCodes;
vGoodUnits = unique(eMat(:,E.unitID));
vGoodUnits = vGoodUnits(~isnan(vGoodUnits));
nUnits = length(vGoodUnits);
idx = ~isnan(eMat(:,E.coherence));
eMat = eMat(idx,:);
stiMatAll = respMat.spkMat.stimAligned(idx,:);
sacMatAll = respMat.spkMat.saccAligned(idx,:);
tarMatAll = respMat.spkMat.targAligned(idx,:);
vTStim = respMat.vTimes.tStim;
vTTarg = respMat.vTimes.tTarg;
vTSacc = respMat.vTimes.tSacc;

% Filter
lFilt = 80; % Length of boxcar filter
b = ones(lFilt,1)/lFilt; % boxcar filter
wFilt = plotpars.filtSize; % Smoothing boxcar filter size

stiMatF = 1e3*filter(b,1,stiMatAll')';
sacMatF = 1e3*filter(b,1,sacMatAll')';
tarMatF = 1e3*filter(b,1,tarMatAll')';

dataOut = struct();
dotDirs = [0 180];
vCohs = unique(eMat(:,E.coherence));
vCohs = vCohs(~isnan(vCohs));

% Find the distance of each target location to the RF
RFdists = nan(size(eMat,1),2);
RFdists(:,1) = sqrt((eMat(:,E.RF_x) - eMat(:,E.target1_x)).^2 + ...
    (eMat(:,E.RF_y) - eMat(:,E.target1_y)).^2);
RFdists(:,2) = sqrt((eMat(:,E.RF_x) - eMat(:,E.target2_x)).^2 + ...
    (eMat(:,E.RF_y) - eMat(:,E.target2_y)).^2);

% Remove unwanted time epochs
for f = 1:size(eMat,1)
    % Cut off 100 ms AFTER dots OFFset
    stiMatF(f,vTStim>eMat(f,E.dot_duration)+100) = NaN;
    
    % Cut off 100 ms before saccade onset for targ aligned and ...
    % ... 200/100 ms after target onset for saccade aligned
    if strcmp(plotpars.monkID,'SM') % For data from Monkey SM 
        cIdxT = vTTarg > (1e3*(eMat(f,E.time_saccade) - eMat(f,E.time_target_on))-100);
        cIdxS = vTSacc < (1e3*(eMat(f,E.time_target_on) - eMat(f,E.time_saccade))+200);
    else cIdxT = vTTarg > eMat(f,E.react_time)-100;
        cIdxS = vTSacc<(-200+eMat(f,E.react_time))*-1;
    end
    tarMatF(f,cIdxT) = nan;
    sacMatF(f,cIdxS) = nan;
end

% Define plot colors 
C_hi = [0.9569 0.3451 0.4000]; % Fiery Rose F45866
C_lo = [0.7804 0.8392 0.4275]; % June Bud C7D66D 
C_mid = [0.2549 0.8275 0.7412]; % Turquoise 41D3BD
C_0 = [0.65 0.65 0.65]; 

%% Normalize populn data and remove target onset response
dotMat = cell(nUnits,1);
tarMat = cell(nUnits,1);
sacMat = cell(nUnits,1);
cMat = cell(nUnits,1); % Event mats of each unit
dMat = cell(nUnits,1); % Mat of RF distances
tIdxD = vTStim>-100 & vTStim<800; % Dots
tIdxT = vTTarg>-300 & vTTarg<500; % Targ
tIdxT1 = vTTarg>100 & vTTarg<500; % Targ
tIdxS = vTSacc>-500 & vTSacc<100; % Sacc

for f = 1:nUnits
    % Normalize each unit to the max after stimulus onset
    idx = eMat(:,E.unitID) == vGoodUnits(f);
    cMax = nanmax(nanmean(tarMatF(idx,tIdxT1))); % Max of the MEANS
    cMin1 = nanmin(nanmin(tarMatF(idx,tIdxT1))); % Min during tar
    cMin2 = nanmin(nanmin(stiMatF(idx,tIdxD))); % Min during RDM
    cMin3 = nanmin(nanmin(sacMatF(idx,tIdxS))); % Min during sac
    cMin = min([cMin1 cMin2 cMin3]); % Absolute min
    tarMat{f} = (tarMatF(idx,tIdxT) - cMin)/(cMax-cMin);
    dotMat{f} = (stiMatF(idx,tIdxD) - cMin)/(cMax-cMin);
    sacMat{f} = (sacMatF(idx,tIdxS) - cMin)/(cMax-cMin);
    cMat{f} = eMat(idx,:);
    dMat{f} = RFdists(idx,:);
end

cMat = cell2mat(cMat);
dMat = cell2mat(dMat);
tarMat = cell2mat(tarMat);
tarMat_ = tarMat; % For non-detrended version

% Detrending of population means
if plotpars.detrend
    % Mean of all 0 & 4% coherence trials
    for f = 1:2
        if f ==1, xf = 2; else xf = 1; end
        idx1 = dMat(:,f)<4 & dMat(:,xf)>4;
        idx2 = idx1 & eMat(:,E.coherence)<=0.05;
        tarMat(idx1,:) = tarMat(idx1,:) - repmat(nanmean(tarMat(idx2,:)),sum(idx1),1);
    end
end

sacMat = cell2mat(sacMat);
dotMat = cell2mat(dotMat);

%% Population response during dots epoch
trCutoff = plotpars.trlcutoff; % fraction of trials dropping out
mnSpksStim = zeros(8,size(dotMat,2));
for f = 1:2 % dot direction
    % Low coherence (4 & 8%)
    switch plotpars.grouping
        case 'all'
            idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.01 & ...
                cMat(:,E.coherence)<0.1;
        case 'correct'
            idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.01 & ...
                cMat(:,E.coherence)<0.1 & cMat(:,E.isCorrect)==1;
    end
    mnSpksStim((f-1)*4+1,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(dotMat(idx,:)));
    nanIdx = sum(~isnan(dotMat(idx,:)))>sum(idx)*trCutoff;
    mnSpksStim((f-1)*4+1,~nanIdx) = nan;
    
    % Intermediate coherence (16%)
    switch plotpars.grouping
        case 'all'
            idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.1 & ...
                cMat(:,E.coherence)<0.2;
        case 'correct'
            idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.1 & ...
                cMat(:,E.coherence)<0.2 & cMat(:,E.isCorrect)==1;
    end
    mnSpksStim((f-1)*4+2,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(dotMat(idx,:)));
    nanIdx = sum(~isnan(dotMat(idx,:)))>sum(idx)*trCutoff;
    mnSpksStim((f-1)*4+2,~nanIdx) = nan;
    
    % High coherence (32% & 64%)
    switch plotpars.grouping
        case 'all'
            idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.2;
        case 'correct'
            idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.2 ...
                & cMat(:,E.isCorrect)==1;
    end
    mnSpksStim((f-1)*4+3,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(dotMat(idx,:)));
    nanIdx = sum(~isnan(dotMat(idx,:)))>sum(idx)*trCutoff;
    mnSpksStim((f-1)*4+3,~nanIdx) = nan;
  
    % 0% coherence trials
    switch plotpars.grouping
        case 'all'
            idx = cMat(:,E.coherence)==0;
        case {'correct','choice'}
            idx = cMat(:,E.coherence)==0 & cMat(:,E.target_choice)==f;
    end
    mnSpksStim((f-1)*4+4,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(dotMat(idx,:)));
    nanIdx = sum(~isnan(dotMat(idx,:)))>sum(idx)*trCutoff;
    mnSpksStim((f-1)*4+4,~nanIdx) = nan;
end

figure; hold all; subplot('Position',[0.1 0.1 0.47 0.8])
ax1 = plot(vTStim(tIdxD),mnSpksStim,'linewidth',2);
set(ax1,{'color'},{C_lo; C_mid; C_hi; C_0; C_lo;C_mid;C_hi;C_0}, ...
    {'LineStyle'},{'-';'-';'-';'-';'--';'--';'--';'--'}, ...
    'linewidth',3);
% legend({'C-lo','C-mid','C-hi','C-0','Y-lo','Y-mid','Y-hi','Y-0'}, ...
%     'location','northwest')
axis tight; 
if strcmp(plotpars.monkID,'AN')
    title ('Fig. 5A (AN)')
    set(gca,'xlim',[-100 599.99],'ylim',[0.25 1.1], ...
        'ytick',0.3:0.3:0.9)
else title ('Fig. 5D (SM)')
    set(gca,'xlim',[-100 599.99],'ylim',[0.1 1.4], ...
        'ytick',0.2:0.4:1.4)
end
line([0 0],ylim,'linewidth',2,'color','k')
xlabel 'Time after motion onset (ms)'
ylabel 'Mean of normalized responses'
axprefs(gca)

%% Population response during target epoch
trCutoff = plotpars.trlcutoff; % fraction of trials dropping out
mnSpksTar = cell(3,1);
mnSpksTarC = zeros(8,size(tarMat,2));
mnSpksTarY = zeros(8,size(tarMat,2));
mnSpksSacC = zeros(8,size(sacMat,2));
mnSpksSacY = zeros(8,size(sacMat,2));
mnSpksTarN = zeros(8,size(tarMat,2)); % No targets in RF
mnSpksSacN = zeros(8,size(sacMat,2));

for tid = 0:1:2 % target in RF ID
    mnSpksTar{tid+1} = zeros(8,size(tarMat,2));
    
    for f = 1:2 % dot direction
        % Low coherence (4 & 8%)
        switch plotpars.grouping
            case 'all'
                idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.01 & ...
                    cMat(:,E.coherence)<0.1 & cMat(:,E.object_in_RF)==tid;
            case 'correct'
                idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.01 & ...
                    cMat(:,E.coherence)<0.1 & cMat(:,E.isCorrect)==1 & cMat(:,E.object_in_RF)==tid;
        end
        % Always use correct trials for saccade mat
        idxS = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.01 & ...
            cMat(:,E.coherence)<0.1 & cMat(:,E.isCorrect)==1 & cMat(:,E.object_in_RF)==tid;
        
        mnSpksTar{tid+1}((f-1)*4+1,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat_(idx,:)));
        nanIdx = sum(~isnan(tarMat_(idx,:)))>sum(idx)*trCutoff;
        mnSpksTar{tid+1}((f-1)*4+1,~nanIdx) = nan;
        if tid == 1
            mnSpksTarC((f-1)*4+1,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarC((f-1)*4+1,~nanIdx) = nan;
            mnSpksSacC((f-1)*4+1,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacC((f-1)*4+1,~nanIdx) = nan;
        elseif tid == 2
            mnSpksTarY((f-1)*4+1,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarY((f-1)*4+1,~nanIdx) = nan;
            mnSpksSacY((f-1)*4+1,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacY((f-1)*4+1,~nanIdx) = nan;
        elseif tid == 0
            mnSpksTarN((f-1)*4+1,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarN((f-1)*4+1,~nanIdx) = nan;
            mnSpksSacN((f-1)*4+1,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacN((f-1)*4+1,~nanIdx) = nan;
            
        end
        
        % Intermediate coherence (16%)
        switch plotpars.grouping
            case 'all'
                idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.1 & ...
                    cMat(:,E.coherence)<0.2 & cMat(:,E.object_in_RF)==tid;
            case 'correct'
                idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.1 & ...
                    cMat(:,E.coherence)<0.2 & cMat(:,E.object_in_RF)==tid & cMat(:,E.isCorrect)==1;
        end
        % Always use correct trials for saccade mat
        idxS = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.1 & ...
            cMat(:,E.coherence)<0.2 & cMat(:,E.isCorrect)==1 & cMat(:,E.object_in_RF)==tid;
        
        mnSpksTar{tid+1}((f-1)*4+2,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat_(idx,:)));
        nanIdx = sum(~isnan(tarMat_(idx,:)))>sum(idx)*trCutoff;
        mnSpksTar{tid+1}((f-1)*4+2,~nanIdx) = nan;
        if tid==1
            mnSpksTarC((f-1)*4+2,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarC((f-1)*4+2,~nanIdx) = nan;
            mnSpksSacC((f-1)*4+2,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacC((f-1)*4+2,~nanIdx) = nan;
        elseif tid==2
            mnSpksTarY((f-1)*4+2,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarY((f-1)*4+2,~nanIdx) = nan;
            mnSpksSacY((f-1)*4+2,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacY((f-1)*4+2,~nanIdx) = nan;
        elseif tid==0
            mnSpksTarN((f-1)*4+2,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarN((f-1)*4+2,~nanIdx) = nan;
            mnSpksSacN((f-1)*4+2,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacN((f-1)*4+2,~nanIdx) = nan;
        end
        
        % High coherence (32 & 64%)
        switch plotpars.grouping
            case 'all'
                idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.2 & ...
                    cMat(:,E.object_in_RF)==tid;
            case 'correct'
                idx = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.2 & ...
                    cMat(:,E.object_in_RF)==tid & cMat(:,E.isCorrect)==1;
        end
        % Always use correct trials for saccade mat
        idxS = cMat(:,E.dot_dir)==dotDirs(f) & cMat(:,E.coherence)>0.02 & ...
            cMat(:,E.isCorrect)==1 & cMat(:,E.object_in_RF)==tid;
        
        mnSpksTar{tid+1}((f-1)*4+3,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat_(idx,:)));
        nanIdx = sum(~isnan(tarMat_(idx,:)))>sum(idx)*trCutoff;
        mnSpksTar{tid+1}((f-1)*4+3,~nanIdx) = nan;
        if tid==1
            mnSpksTarC((f-1)*4+3,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarC((f-1)*4+3,~nanIdx) = nan;
            mnSpksSacC((f-1)*4+3,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacC((f-1)*4+3,~nanIdx) = nan;
        elseif tid==2
            mnSpksTarY((f-1)*4+3,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarY((f-1)*4+3,~nanIdx) = nan;
            mnSpksSacY((f-1)*4+3,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacY((f-1)*4+3,~nanIdx) = nan;
        elseif tid==0
            mnSpksTarN((f-1)*4+3,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarN((f-1)*4+3,~nanIdx) = nan;
            mnSpksSacN((f-1)*4+3,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacN((f-1)*4+3,~nanIdx) = nan;
        end
        
        % 0% coherence trials
        idx = cMat(:,E.coherence)==0 & cMat(:,E.object_in_RF)==tid;
        switch plotpars.grouping
            case 'all'
                idx = cMat(:,E.coherence)==0 & cMat(:,E.object_in_RF)==tid;
            case {'correct','choice'}
                idx = cMat(:,E.coherence)==0 & cMat(:,E.object_in_RF)==tid & ...
                    cMat(:,E.target_choice)==f;
        end
        % Always use correct trials for saccade mat
        idxS = cMat(:,E.coherence)==0 & cMat(:,E.object_in_RF)==tid & cMat(:,E.target_choice)==f;
        
        mnSpksTar{tid+1}((f-1)*4+4,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat_(idx,:)));
        nanIdx = sum(~isnan(tarMat_(idx,:)))>sum(idx)*trCutoff;
        mnSpksTar{tid+1}((f-1)*4+4,~nanIdx) = nan;
        if tid==1
            mnSpksTarC((f-1)*4+4,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarC((f-1)*4+4,~nanIdx) = nan;
            mnSpksSacC((f-1)*4+4,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacC((f-1)*4+4,~nanIdx) = nan;
        elseif tid==2
            mnSpksTarY((f-1)*4+4,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarY((f-1)*4+4,~nanIdx) = nan;
            mnSpksSacY((f-1)*4+4,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacY((f-1)*4+4,~nanIdx) = nan;
        elseif tid==0
            mnSpksTarN((f-1)*4+4,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(tarMat(idx,:)));
            nanIdx = sum(~isnan(tarMat(idx,:)))>sum(idx)*trCutoff;
            mnSpksTarN((f-1)*4+4,~nanIdx) = nan;
            mnSpksSacN((f-1)*4+4,:) = filtfilt(ones(wFilt,1)/wFilt,1,nanmean(sacMat(idxS,:)));
            nanIdx = sum(~isnan(sacMat(idxS,:)))>sum(idxS)*trCutoff;
            mnSpksSacN((f-1)*4+4,~nanIdx) = nan;
        end
    end
end

% Blue in RF (main)
figure; hold all; 
subplot('Position',[0.1 0.1 0.47 0.8])
ax1 = plot(vTTarg(tIdxT),mnSpksTar{2},'linewidth',2);
set(ax1,{'color'},{C_lo; C_mid; C_hi; C_0; C_lo;C_mid;C_hi;C_0}, ...
    {'LineStyle'},{'-';'-';'-';'-';'--';'--';'--';'--'}, ...
    'linewidth',3);
% legend({'C-lo','C-mid','C-hi','C-0','Y-lo','Y-mid','Y-hi','Y-0'}, ...
%     'location','northwest')
axis tight;
if strcmp(plotpars.monkID,'AN')
    title ('Fig. 5B (AN)')
    set(gca,'xlim',[-100 599.99],'ylim',[0.25 1.1],'ytick',0.3:0.3:0.9)
else title ('Fig. 5E (SM)')
    set(gca,'xlim',[-100 599.99],'ylim',[0.1 1.4],'ytick',0.2:0.4:1.4)
end
line([0 0],ylim,'linewidth',2,'color','k')
xlabel 'Time after motion onset (ms)'
ylabel 'Mean of normalized responses'
axprefs(gca)

% Blue in RF (inset)
% figure; hold all; 
subplot('Position',[0.6 0.4 0.4 0.5])
ax1 = plot(vTTarg(tIdxT),mnSpksTarC,'linewidth',2);
set(ax1,{'color'},{C_lo; C_mid; C_hi; C_0; C_lo;C_mid;C_hi;C_0}, ...
    {'LineStyle'},{'-';'-';'-';'-';'--';'--';'--';'--'}, ...
    'linewidth',3);
% legend({'C-lo','C-mid','C-hi','C-0','Y-lo','Y-mid','Y-hi','Y-0'}, ...
%     'location','northwest')
axis tight; 
if strcmp(plotpars.monkID,'AN')
    set(gca,'xlim',[-100 599.99],'ylim',[-0.1 0.2], ...
        'ytick',-0.1:0.1:0.2,'yticklabel',{'',0,'',0.2})
else set(gca,'xlim',[-100 599.99],'ylim',[-0.2 0.25], ...
        'ytick',-0.2:0.1:0.2,'yticklabel',{'','',0,'',0.2})
end
line([0 0],ylim,'linewidth',2,'color','k')
axprefs(gca)


% Yellow in RF (main)
figure; hold all; 
subplot('Position',[0.1 0.1 0.47 0.8])
ax1 = plot(vTTarg(tIdxT),mnSpksTar{3},'linewidth',2);
set(ax1,{'color'},{C_lo; C_mid; C_hi; C_0; C_lo;C_mid;C_hi;C_0}, ...
    {'LineStyle'},{'-';'-';'-';'-';'--';'--';'--';'--'}, ...
    'linewidth',3);
% legend({'C-lo','C-mid','C-hi','C-0','Y-lo','Y-mid','Y-hi','Y-0'}, ...
%     'location','northwest')
axis tight; 
if strcmp(plotpars.monkID,'AN')
    title ('Fig. 5C (AN)')
    set(gca,'xlim',[-100 599.99],'ylim',[0.25 1.1],'ytick',0.3:0.3:0.9)
else title ('Fig. 5F (SM)')
    set(gca,'xlim',[-100 599.99],'ylim',[0.1 1.4],'ytick',0.2:0.4:1.4)
end
line([0 0],ylim,'linewidth',2,'color','k')
xlabel 'Time after motion onset (ms)'
ylabel 'Mean of normalized responses'
axprefs(gca)

% Yellow in RF (inset)
subplot('Position',[0.6 0.4 0.4 0.5])
ax1 = plot(vTTarg(tIdxT),mnSpksTarY,'linewidth',2);
set(ax1,{'color'},{C_lo; C_mid; C_hi; C_0; C_lo;C_mid;C_hi;C_0}, ...
    {'LineStyle'},{'-';'-';'-';'-';'--';'--';'--';'--'}, ...
    'linewidth',3);
% legend({'C-lo','C-mid','C-hi','C-0','Y-lo','Y-mid','Y-hi','Y-0'}, ...
%     'location','northwest')
axis tight; 
if strcmp(plotpars.monkID,'AN')
    set(gca,'xlim',[-100 599.99],'ylim',[-0.1 0.2], ...
        'ytick',-0.1:0.1:0.2,'yticklabel',{'',0,'',0.2})
else set(gca,'xlim',[-100 599.99],'ylim',[-0.2 0.25], ...
        'ytick',-0.2:0.1:0.2,'yticklabel',{'','',0,'',0.2})
end
line([0 0],ylim,'linewidth',2,'color','k')
axprefs(gca)
