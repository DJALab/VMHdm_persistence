
N   = 1000;
[w1,w2] = deal(zeros(N,1));

% reset everything to be safe
doOnline     = 0;
plotTrial    = 0;
doOverlap    = 0;
doStructured = 0;
doBanded     = 0;

dt      = 0.005;
tau     = 0.1;
tauI    = 1;
tauS    = 2.5;
time    = dt:dt:120;
tOn     = 2.5;
tOff    = 4.5;

stimRaw = double(time>tOn & time<tOff);
stim1   = smoothts(stimRaw,'e',4/dt);stim1=stim1/max(eps+stim1);
stim1(time>tOn+10)=0;

time    = time - 60;
stimRaw = double(time>tOn & time<tOff);
stim2   = smoothts(stimRaw,'e',4/dt);stim2=stim2/max(eps+stim2);
stim2(time>tOn+10)=0;
time    = time + 60 - tOn;

% the version I use in the paper:
% gInh = 1.6; g=1; tauI=1; tauS=8; doDecorr=1; scPost = 0.65; beta = 2; tRun = 10000; scaleEvery = 11000;
% gInh = 1.5; g=1; tauI=1; tauS=6; doDecorr=1; scPost = 0.65; beta = 2; tRun = 10000; scaleEvery = 5000;

% parameters for the low-g and high-g examples of sRNN:
% gInh = 3.5; g=1; tauI=0.1; tauS=2; doDecorr=0;
% gInh = 3.5; g=4; tauI=0.1; tauS=2; doDecorr=0;

% gInh = 1.5; g=1; tauI=1; tauS=6; doDecorr=1; scPost = 0.65; beta = 2; tRun = 10000; scaleEvery = 11000;

% good banded performance:
% gInh = 3.75; g = 2; tauI = 0.05; tauE = 0.05; tauS = 6;
% doBanded = 1; doStructured = 1; Jsigma  = N*2/3; offset  = 0;

% good banded performance with fast + slow excitation:
% gInh = 2; g = 4; tauI = 0.05; tauE = 0.05; tauS = 6;
% doBanded = 1; doStructured = 1; Jsigma  = N*2/3; offset  = 0;


gInh = 2.4; g = 2; tauI = 0.05; tauE = 0.05; tauS = 6;
doBanded = 1; doStructured = 1; Jsigma  = N*.7; offset  = 0;

sp       = 0.1;
thr      = 0.1;
doOnline = 0;

% example of too much banding:
% gInh = 9; Jsigma = N/4; frO = 0;
%
% [stimTriggeredAutocorr, PCCStore, diagData] = deal({});
frOran = 0.5;%1:-0.2:0;
nreps = 10;
for OLiter = 1%:length(frOran)

    if(doOverlap)
        fr = round(N/2);
        frO = round(.7*fr);
        frI = fr - frO;
        overlap = randperm(N,frO);
        wshare = rand(frO,1);
        w1([randperm(N,frI) overlap]) = g*[rand(frI,1); wshare];
        w2([randperm(N,frI) overlap]) = g*[rand(frI,1); wshare];
    elseif(doStructured) % allow some overlap
        fr = round(N/2);
        frO = frOran(OLiter);
        w1(randperm(ceil(N*(0.5+frO/2)),fr)) = g * rand(fr,1);
        w2(randperm(ceil(N*(0.5+frO/2)),fr)+floor(N*(1-frO)/2)) = g * rand(fr,1);
    else
        fr = round(N/2);
        w1(randperm(N,fr)) = g * rand(fr,1);
        w2(randperm(N,fr)) = g * rand(fr,1);
    end

    figure(20);clf;
    [activeMean,trMean,chanceMean,stimTriggeredAutocorr,manFix,manShuff,PCCStore,CDStore]=deal([]);
    if(doOnline)
        figure(10);clf;
        h = plot(rand(1,10)');
        xlim(time([1 end]));
    end
    rstore = {};
    for rep = 1:nreps
        J = rand(N).*(rand(N)<sp)/sqrt(N*sp);
        if(doBanded)
            J = makeJBanded(N,sp,Jsigma,offset);
        elseif(doDecorr)
            J = decorrelateWeights(J,g/2,dt,tauS,tau,thr,beta,tRun,scaleEvery)*scPost;
        else
            J = rand(N).*(rand(N)<sp)/sqrt(N*sp);
        end
        [r,x,p,pS] = deal(zeros(N,length(stim2))); 
        I       = zeros(1,length(stim2));
        IS      = zeros(1,length(stim2));

        disp(rep);
        for t = 2:length(time)
            % synaptic conductances
            p(:,t)   = r(:,t-1) + (p(:,t-1) - r(:,t-1)) * exp(-dt/tauE);
            if(t>1/dt)
                smRates = r(:,t-1).*((sum(r(:,max(t-1/dt,1):t-1),2) - 20)>0);
                pS(:,t)  = smRates + (pS(:,t-1) - smRates) * exp(-dt/tauS);
            else
                pS(:,t) = pS(:,t-1);
            end

            % membrane potentials
            Iin     = max(  g/2*J*pS(:,t-1) ...
                          + g/2*J*p(:,t-1) ...
                          + w1*stim1(t) + w2*stim2(t)...
                          - (g*(g+gInh)/1000*I(t-1))... + g*(g+gInh)/1000*IS(t-1))...
                          + rand(N,1)*0.05,0);
            x(:,t)  = Iin + (x(:,t-1) - Iin) * exp(-dt/tau);

            % spikes
            r(:,t)  = (x(:,t)>=thr)*1/dt/100;
        %     x(any(r(:,t-4:t)~=0,2)',t) = 0;
            x(r(:,t)~=0,t) = 0;

            % inhibition
            I(t) = sum(r(:,t-1)) + (I(t-1) - sum(r(:,t-1))) * exp(-dt/tauI);
            
            if(mod(t,200)==0 && doOnline)
                disp([num2str(t) '/' num2str(length(time))]);
                set(h,'xdata',time(1:t),'ydata',smoothts(mean(r(:,1:t),1),'b',.1/dt));
                drawnow;
            end
            if(time(t)>56 && time(t)<60)
                r(:,t)=0;x(:,t)=0;I(t)=0;Iin=0;pS(:,t)=0;
            end
        end


    % compute R1 and R2 from model raster:
        rIter1 = r(:,time>0 & time<55)~=0;
        rIter2 = r(:,(time+tOn)>60 & (time+tOn)<(60+55))~=0;
        rstore{rep,1} = rIter1;
        rstore{rep,2} = rIter2;

        f = @(x) corr(x);
        f2 = @(x) zscore(x')';
        f3 = @(x,y) diag(corr(x,y))';
        figure(20);
        [rMean,rFix,rShuff,PCC,active,overlap,chance,timeDS] = simStatsAndPlots(rIter1,rIter2,time,f,f2,f3,1);
    %     f3 = @(x,y) diag(pdist2(x',y','cosine'))';
    %     [~,~,~,cosDist,~,~,~,~] = simStatsAndPlots(rIter1,rIter2,time,f,f2,f3,1);

        stimTriggeredAutocorr     = cat(3,stimTriggeredAutocorr,rMean);
        manFix      = cat(3,manFix,rFix);
        manShuff    = cat(3,manShuff,rShuff);
        PCCStore    = [PCCStore; PCC];
    %     CDStore    = [CDStore; cosDist];
        activeMean  = [activeMean; active];
        trMean      = [trMean; overlap];
        chanceMean  = [chanceMean; chance];
    end
    %
end

%%    
figure(331);clf
for OLiter = 1:length(frOran)
    timeDS=1:85;
    subplot(2,length(frOran),OLiter);
    hold on;
    drawvar(timeDS-tOn,PCCStore{OLiter},'k',1/sqrt(nreps));
    ylim([-.5 1]);
    xlim([0 30]);
    set(gca,'xtick',0:10:30,'ytick',-.25:.25:1);

    subplot(2,length(frOran),length(frOran)+OLiter);
    hold on;
    diagData{OLiter} = computeMetricDiagonal(stimTriggeredAutocorr{OLiter},@(x)nanmean(x));

    drawvar([0 timeDS]-tOn,diagData{OLiter},'r',1/sqrt(size(stimTriggeredAutocorr{OLiter},3)));
    xlim([0 45]);
    ylim([0 1]);
    set(gca,'xtick',0:15:45)
end


%%
rep = 1;
r1 = smoothts(smoothts(double(rstore{rep,1}),'e',.1/dt),'e',1.5/dt);
r2 = smoothts(smoothts(double(rstore{rep,2}),'e',.1/dt),'e',1.5/dt);
time=dt:dt:85-dt;

figure(50);clf;hold on;
n = 30;
plot(time-tOn,bsxfun(@plus,r1(1:n,:)',(1:n)*.1));
plot(time-tOn,bsxfun(@plus,r2(1:n,:)',(1:n)*.1));
axis tight;
xlim([-2 45]);

[pki,i] = max(r1,[],2);
[pkj,j] = max(r2,[],2);
thr = 0.04;
keepCells = (pki>thr)&(pkj>thr);
time = (dt:dt:85-dt) - tOn;

figure(60);clf;
plot(time(i(keepCells)),time(j(keepCells)),'.','markersize',20);
xlim([0 40]);
ylim([0 40]);
box off;
axis equal;














