function perfOut = abstractDots_psychometrics(eMat,nFig,figTitle)

% nFig - Figure no. to start
% figTitle - title for the figures (eg., Monkey name)

%% Common parameters
E = eventCodes;

% Remove uncompleted trials
idx = eMat(:,E.isCorrect)==0 | eMat(:,E.isCorrect)==1;
eMat = eMat(idx,:);
nTrials = size(eMat,1);

% Some useful fields
vCohs = unique(eMat(:,E.coherence));
dotDirs = unique(eMat(:,E.dot_dir));

% Signed coherences
cohList = eMat(:,E.coherence);
idx = eMat(:,E.dot_dir)==180;
cohList(idx) = cohList(idx)*-1;
signCohs = sort(unique(cohList));

%% Logistic fits
tChoice = eMat(:,E.isCorrect); % Use isCorrect and recode as choices
tChoice(eMat(:,E.dot_dir)==180 & eMat(:,E.isCorrect)==0) = 1;
tChoice(eMat(:,E.dot_dir)==180 & eMat(:,E.isCorrect)==1) = 0;
perfOut.fits = getLogistic(cohList,tChoice,figTitle,nFig);

end


%% SUBFUNCTIONS
function [dataOut] = getLogistic(vCoherences,vChoices,sTitle,figNo)

dataOut = struct();
vCohs = unique(vCoherences);
vCohs = vCohs(~isnan(vCohs)); % Remove NaNs
cohsForPlot = min(vCohs)-0.03:0.01:max(vCohs)+0.03;

perfOut = zeros(length(vCohs),2);
for f = 1:length(vCohs)
    idx = vCoherences==vCohs(f);
    perfOut(f,1) = sum(vChoices(idx)==1)/sum(idx);
    perfOut(f,2) = sqrt(perfOut(f,1)*(1-perfOut(f,1))/sum(idx));
end

dataForFits = [ones(length(vChoices),1) vCoherences vChoices];
dataForFits(dataForFits(:,end)==2,end) = 0;
[dataOut.logFits,~,~,dataOut.logSems] = logist_fit(dataForFits);%,0,0.1);

figure(figNo); hold all
errorbar(vCohs,perfOut(:,1),perfOut(:,2),'linestyle','none','marker','.',...
    'markerfacecolor','b','markersize',40)

vPs = dataOut.logFits(1) + (1-dataOut.logFits(1)-dataOut.logFits(2))./(1+exp(-(dataOut.logFits(3) ...
    + dataOut.logFits(4)*cohsForPlot)));

plot(cohsForPlot,vPs,'linewidth',2,'color','b')
text(-0.45,0.95,['\lambda = ' num2str(round(dataOut.logFits(1)*100)/100)])
text(-0.45,0.9,['\gamma = ' num2str(round(dataOut.logFits(2)*100)/100)])
text(-0.45,0.85,['\beta0 = ' num2str(round(dataOut.logFits(3)*100)/100)])
text(-0.45,0.8,['\beta1 ='  num2str(round(dataOut.logFits(4)*100)/100)])
text(-0.45,0.75,sprintf('nTrials: %g',length(vChoices)))

title (sTitle)

set(gca,'xlim',[min(cohsForPlot) max(cohsForPlot)],'ylim',[0 1])
line(xlim,[0.5 0.5],'color','k')
line([0 0],ylim,'color','k')

xlabel 'Motion strength (% coh)'
ylabel 'p (choose right)'
axprefs(gca)

end