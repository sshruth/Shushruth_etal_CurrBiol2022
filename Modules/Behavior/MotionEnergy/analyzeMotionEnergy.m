function dataOut = analyzeMotionEnergy(eMat,loCoh,ME,sTitle)
% loCoh: Coherence to use for ME calculation
% ME: Precomputed ME for use in plotting
% sMonk: Monkey string to add to figure title

%% Set up and process data

E = eventCodes;
scr_refresh = 1000/75; % Both monkeys have 75Hz screen refresh 
maxFrames = size(ME,2);

% Invert sign of 180° motion
idx = eMat(:,E.dot_dir)==180;
ME(idx,:) = ME(idx,:)*-1;

% Get relevant data
trlIdx = eMat(:,E.coherence)==loCoh & ~isnan(eMat(:,E.time_targ_acq));
cData = eMat(trlIdx,:);
ME = ME(trlIdx,:);

% If coh>0, remove mean coherence
if loCoh>0
    idx = cData(:,E.dot_dir)==180;
    ME(idx,:) = ME(idx,:)-repmat(nanmean(ME(idx,:)),sum(idx),1);
    ME(~idx,:) = ME(~idx,:)-repmat(nanmean(ME(~idx,:)),sum(~idx),1);
end

%% Plot
% Choice based MEs
dotTime = scr_refresh*(1:maxFrames); 
tIdx = dotTime<500;
fill_x = [dotTime(tIdx) fliplr(dotTime(tIdx))];
idx1 = cData(:,E.target_choice)==1;
idx2 = cData(:,E.target_choice)==2;

% Average both directions and plot
figure; hold all
cMEs = [ME(idx1,tIdx);ME(idx2,tIdx)*-1]; % Concatenate after sign correction
cMns = nanmean(cMEs);
cSD = nanstd(cMEs)./sqrt(sum([idx1;idx2]));
fill_y = [cMns+cSD fliplr(cMns)-fliplr(cSD)];
fill(fill_x,fill_y,[0.7 0.7 0.7],'facealpha',0.5,'edgecolor','none')
plot(dotTime(tIdx),cMns,'linewidth',3,'color','k')
line(xlim,[0 0],'color','k','linewidth',2)
set(gca,'ylim',[-0.12 0.22],'ytick',-0.1:0.05:0.2)
title(sTitle)
ylabel 'Motion energy (a.u.)'
xlabel 'Time from motion onset (ms)'
axprefs(gca)

dataOut.eMat = cData;
dataOut.MEs = ME;