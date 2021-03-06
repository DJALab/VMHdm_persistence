
N = 1000; sp = 0.1; thr = 0.1; [w1,w2] = deal(zeros(N,1));

% reset everything to be safe
rng(42);
[doOnline,plotTrial,doOverlap,doStructured,doDecorr,doBanded] = deal(false);

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

% good banded performance:
% gInh = 3.75; g = 2; tauI = 0.05; tauE = 0.05; tauS = 6;
% doBanded = 1; doStructured = 1; Jsigma  = N*2/3; offset  = 0;

% good banded performance with fast + slow excitation:
% gInh = 2; g = 4; tauI = 0.05; tauE = 0.05; tauS = 6;
% doBanded = 1; doStructured = 1; Jsigma  = N*2/3; offset  = 0;


gInh = 1.75; g = 6; tauI = 0.05; tauE = 0.05; tauS = 6;
doBanded = 0; doStructured = 0;

doOnline = 0;
nreps = 10;

fr = round(N/4);
w1(randperm(N,fr)) = g*2 * rand(fr,1);
w2(randperm(N,fr)) = g*2 * rand(fr,1);

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
        Iin     = max(  g/2*J*pS(:,t-1) + g/2*J*p(:,t-1) - (g*(g+gInh)/1000*I(t-1))...
                      + w1*stim1(t) + w2*stim2(t) + rand(N,1)*0.05,0);
        x(:,t)  = Iin + (x(:,t-1) - Iin) * exp(-dt/tau);

        % spikes
        r(:,t)  = (x(:,t)>=thr)*1/dt/100;
        x(r(:,t)~=0,t) = 0;

        % inhibition
        I(t) = sum(r(:,t-1)) + (I(t-1) - sum(r(:,t-1))) * exp(-dt/tauI);
        if(mod(t,200)==0 && doOnline)
            disp([num2str(t) '/' num2str(length(time))]);
            set(h,'xdata',time(1:t),'ydata',smoothts(mean(r(:,1:t),1),'b',.1/dt));
            drawnow;
        end
        if(time(t)>56 && time(t)<60)
            r(:,t)=0;x(:,t)=0;I(t)=0;Iin=0;pS(:,t)=0;p(:,t)=0; smRate=zeros(size(smRates));
        end
    end
    
    % compute R1 and R2 from model raster:
    rIter1 = r(:,time > 0  & time<55)~=0;
    rIter2 = r(:,time > 60 & time<(60+55))~=0;
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
    activeMean  = [activeMean; active];
    trMean      = [trMean; overlap];
    chanceMean  = [chanceMean; chance];
end

%%
thr         = 0.05;
sig         = smoothts(smoothts(double(rstore{rep,1}),'e',.1/dt),'e',2.5/dt);
response    = max(sig,[],2);
sig         = sig(response > thr,:);
plottime = (1:size(sig,2))*dt;

Ctime = plottime(100:100:end);
C = zeros(length(Ctime));
D = [];
for rep = 1:nreps
    sig     = smoothts(smoothts(double(rstore{rep,1}),'e',.1/dt),'e',2.5/dt);
    Ctr     = corr(sig(:,100:100:end));
    C = C + Ctr/nreps;
    D = [D; computeMetricDiagonal(Ctr,@(x)nanmean(x))];
end

figure(1);clf;
subplot(2,1,1);
image(C*32+32)
subplot(2,1,2);
drawvar([0 Ctime],D,'b',1);
ylim([-.25 1]);














