
N = 1000; sp = 0.1; thr = 0.1;

% reset everything to be safe
[doOnline,plotTrial,doOverlap,doStructured,doBanded,doDecorr] = deal(false);

dt      = 0.005;
tau     = 0.1;
time    = dt:dt:120;
tOn     = 2.5;
tOff    = 4.5;

stimRaw = double(time>tOn & time<tOff);
stim1   = smoothts(stimRaw,'e',4/dt);stim1=stim1/max(eps+stim1);
stim1(time>tOn+10)=0;
stim2 = circshift(stim1,sum(time<=60));
time    = time - tOn;


gInh = 2.8; g = 1; tauI = 0.05; tauE = 0.05; tauS = 0.2;
reroll = 1;
g*(g+gInh)
%%

doOnline = 0;
nreps = 10;

if(reroll)
    [w1,w2] = deal(zeros(N,1));
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
    if(reroll)
        J = rand(N).*(rand(N)<sp)/sqrt(N*sp);
        if(doBanded)
            J = makeJBanded(N,sp,Jsigma,offset);
        elseif(doDecorr)
            J = decorrelateWeights_fastSlow(J,g/2,dt,tauE,tauS,tau,thr,beta,tRun,doHebbian,w1,w2)*scPost;
        else
            J = rand(N).*(rand(N)<sp)/sqrt(N*sp);
        end
    end
    [r,x,p,pS] = deal(zeros(N,length(stim2))); 
    I       = zeros(1,length(stim2));
    IS      = zeros(1,length(stim2));

    disp(rep);
    for t = 2:length(time)
        % synaptic conductances
        p(:,t)   = r(:,t-1) + (p(:,t-1) - r(:,t-1)) * exp(-dt/tauE);
        pS(:,t)  = r(:,t-1) + (pS(:,t-1) - r(:,t-1)) * exp(-dt/tauS);

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

thr         = 0;
sig         = smoothts(smoothts(double(rstore{rep,1}),'e',.1/dt),'e',2.5/dt);
response    = max(sig,[],2);
sig         = sig(response > thr,:);

[pks,locs] = max(sig,[],2);
temp = sortrows([locs sig],1);
temp = temp(:,2:end);
temp = bsxfun(@minus,temp,median(temp(:,1:110),2));
temp = bsxfun(@times,temp,1./max(temp,[],2));
plottime = (1:size(temp,2))*dt;

figure(1);clf;
image(plottime,[],(temp')'*64*4/5+64/5);
xlim([0 45]);



%%
figure(1);clf;count=0;
for use = {rstore_sRNN_lowG, rstore_sRNN_highG, rstore_sRNN_banded}
    count=count+1;
    [Cearly,Cmid,Cearly_off,Cmid_off] = deal([]);
    for i=1:10
        sig     = smoothts(smoothts(double(use{:}{i,1}),'e',.1/dt),'e',1.5/dt);
        sig2    = smoothts(smoothts(double(use{:}{i,2}(:,1:length(sig))),'e',.1/dt),'e',1.5/dt);
        tEarly  = mean(sig(:,((dt/dt):(5/dt))),2);
        tMid    = mean(sig(:,((25/dt):(30/dt))),2);

        Cearly  = [Cearly; corr(sig,tEarly)'];
        Cmid    = [Cmid; corr(sig,tMid)'];

        Cearly_off  = [Cearly_off; corr(sig2,tEarly)'];
        Cmid_off    = [Cmid_off; corr(sig2,tMid)'];
    end
    time = (1:length(Cearly))*dt;

    data = traces.USS;
    tEarly      = mean(data(:,300+(((1/11)*11):(5*11))),2);
    tMid        = mean(data(:,300+((25*11):(30*11))),2);
    Cearly_data = corr(data,tEarly);
    Cmid_data   = corr(data,tMid);
    data = traces.rat;
    Cearly_off_data = corr(data,tEarly);
    Cmid_off_data = corr(data,tMid);
    dTime = (1:length(Cearly_data))/11 - 300/11;

    subplot(3,2,(count-1)*2+1);hold on;
    plot(dTime,Cearly_off_data,'k');
    plot(dTime,Cmid_off_data,'m');

    plot(dTime,Cearly_data,'b');
    plot(dTime,Cmid_data,'r');

    plot([0 45],[0 0],'k');
    xlim([-15 45]);
    ylim([-.1 1]);

    subplot(3,2,(count-1)*2+2);hold on;
    drawvar(time,Cearly_off,'k');
    drawvar(time,Cmid_off,'m')

    drawvar(time,Cearly,'b');
    drawvar(time,Cmid,'r')

    plot([0 45],[0 0],'k');
    box off;
    xlim([0 45]);
    ylim([-.1 1]);
    drawnow;
end
subplot(3,2,5);hold on;
legend('other early','other mid','early','mid','location','southeast');








