function [rMean,rFix,rShuff,PCC,active,overlap,chance,timeDS] = simStatsAndPlots(rIter1,rIter2,time,f,f2,f3,simulated)

dt = time(2)-time(1);
keepCells = find(sum([rIter1(:,sum(time<5):end) rIter2(:,sum(time<5):end)],2)>10);

R1=rIter1;
if(simulated)
    r_iter  = rIter1(keepCells,:);
    temp    = [zeros(size(r_iter,1),1) r_iter(:,time<0) r_iter];
    R1      = smoothts(smoothts(temp,'e',0.5/dt),'e',1.5/dt);
    R1      = R1(:,(sum(time<0)):end);
end
R1 = smoothts(R1,'b',round(1/dt));
R1 = R1(:,1:round(1/dt):end);

R2=rIter2;
if(simulated)
    r_iter  = rIter2(keepCells,:);
    temp    = [zeros(size(r_iter,1),1) r_iter(:,time<0) r_iter];
    R2      = smoothts(smoothts(temp,'e',0.5/dt),'e',1.5/dt);
    R2      = R2(:,(sum(time<0)):end);
end
R2      = smoothts(R2,'b',round(1/dt));
R2      = R2(:,1:round(1/dt):end);

% compute within-trial distance measures
[rMean, rFix, rShuff] = applyMetricToRaster(R1,f,f2);

% compute between-stim similarity measure
PCC = f3(R1+rand(size(R1))*eps,R2+rand(size(R1))*eps);

% binarize data
timeDS = time(1:round(1/dt):end);
if(simulated)
    Rthr = 0.02;
else
    mu   = mean([R1(:,timeDS<-5) R2(:,timeDS<-5)],2);
    sig  = std([R1(:,timeDS<-5) R2(:,timeDS<-5)],[],2);
    Rthr = mu + 2*sig;
end
R1      = bsxfun(@gt,R1,Rthr);
R2      = bsxfun(@gt,R2,Rthr);

% and now compute overlap measures in the binarized representations
active = mean(R1|R2);
overlap = (sum(R1&R2)./sum(R1|R2));
overlap(isnan(overlap))=0;
chance = mean(R1).*mean(R2)./(mean(R1)+mean(R2)-mean(R1).*mean(R2));
chance(isnan(chance))=0;

if(simulated)
    subplot(3,1,1);hold on;
    plot(mean(R1),'b');
    plot(mean(R2),'r');
    xlim([0 90]);
    
    subplot(3,1,2);hold on;
    corrPlot = computeMetricDiagonal(rMean,@(x)nanmean(x));
    plot(corrPlot);
    plot([0 50],[0.5 0.5],'k');
    ylim([0 1]);
    xlim([0 45]);
    
    subplot(3,1,3);hold on;
    plot(PCC);
    plot([0 50],[0 0],'k');
    ylim([-.25 1]);
    xlim([0 30]);
    drawnow;
end