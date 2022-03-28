function dataOut = logist_plot(vCoherences,vChoices,lnClr)
% function dataOut = logist_plot(vCoherences,vChoices,sTitle,lnClr)
% vCoherences - vector of signed coherences
% vChoices - preferred choice 1 other choice 0
% lnClr - Line property. Eg: 'r-' is red dots and red solid line

dataOut = struct();
vCohs = unique(vCoherences);
vCohs = vCohs(~isnan(vCohs)); % Remove NaNs

perfOut = zeros(length(vCohs),1);
for f = 1:length(vCohs)
    idx = vCoherences==vCohs(f);
    perfOut(f) = sum(vChoices(idx)==1)/sum(idx);
end

dataForFits = [ones(length(vChoices),1) vCoherences vChoices];
dataForFits(dataForFits(:,end)==2,end) = 0;
[dataOut.logFits,~,~,dataOut.logSems] = logist_fit(dataForFits,0);

scatter(vCohs,perfOut,'filled','MarkerFaceColor',lnClr(1))

cohsForPlot = min(vCohs)-0.03:0.01:max(vCohs)+0.03;
vPs = 1./(1+exp(-(dataOut.logFits(3) + dataOut.logFits(4)*cohsForPlot)));
plot(cohsForPlot,vPs,lnClr,'linewidth',2);

set(gca,'xlim',[min(cohsForPlot) max(cohsForPlot)],'ylim',[0 1])
xlabel 'Coherence'
ylabel 'p (right)'