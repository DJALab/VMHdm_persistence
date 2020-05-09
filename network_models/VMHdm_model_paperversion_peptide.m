
% makes the "peptide" model: all cells decay at same rate.

dt      = 0.005;
time = (1:10999)*dt - 5;
tauPop = 28; %approx. tau of the rat response
popMean = (time>0).*exp(-time/tauPop)*80;

figure(20); clf;
[activeMean,trMean,chanceMean,stimTriggeredAutocorr,manFix,manShuff,PCCStore,CDStore]=deal([]);
for i=1:10
    for stim = 1:2
        w = rand(1000,1).*(rand(1000,1)>.1); %input to a random X% of cells
        popRast = w*popMean;
        rstore_peptide{i,stim} = poissrnd(popRast*dt);
    end
    rIter1 = rstore_peptide{i,1}(:,time>0);
    rIter2 = rstore_peptide{i,2}(:,time>0);
    
    f = @(x) corr(x);
    f2 = @(x) zscore(x')';
    f3 = @(x,y) diag(corr(x,y))';
    [rMean,rFix,rShuff,PCCtemp,active,overlap,chance,timeDS] = simStatsAndPlots(rIter1,rIter2,time,f,f2,f3,1);
    
    stimTriggeredAutocorr     = cat(3,stimTriggeredAutocorr,rMean);
    manFix      = cat(3,manFix,rFix);
    manShuff    = cat(3,manShuff,rShuff);
    PCCStore    = [PCCStore; PCCtemp];
    activeMean  = [activeMean; active];
    trMean      = [trMean; overlap];
    chanceMean  = [chanceMean; chance];
end

figure(1);clf;
drawvar(timeDS(timeDS>-1),PCCStore,'b',1);
box off
ylim([-.25 1]);
xlim([-10 30])
set(gca,'xtick',-10:10:30)
set(gca,'ytick',-.25:.25:1)

figure(2);clf;
diagM = computeMetricDiagonal(stimTriggeredAutocorr,@(x)nanmean(x));
drawvar(timeDS(timeDS>-2)+2,diagM,'b',1);
box off;
set(gca,'xtick',0:15:45);
set(gca,'ytick',0:.2:1);
xlim([0 45])
ylim([0 1]);