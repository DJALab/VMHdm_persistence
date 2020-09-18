
sess='session1';
rast = 'rast';
stimlist = {'male','rat','toy','USS','tone'};
mouseList = [1 3 4 5 7];
cmap = 'brgyc';
count=0;
figure(1);clf;hold on;
for s = stimlist
    count=count+1;
    stim = s{:};
    mu = {};
    x.(s{:}) = [];
    for ex = mouseList
        hits = find(strcmpi({expt(ex).(sess).stim},stim),2);
        
        win = -200:550;
        inds = expt(ex).(sess)(hits(1)).annot.stim.stim_on(1) + win;
        inds(inds>length(expt(ex).(sess)(hits(1)).(rast))) = length(expt(ex).(sess)(hits(1)).(rast));
        sig1 = expt(ex).(sess)(hits(1)).(rast)(:,inds);
        inds = expt(ex).(sess)(hits(2)).annot.stim.stim_on(1) + win;
        inds(inds>length(expt(ex).(sess)(hits(2)).(rast))) = length(expt(ex).(sess)(hits(2)).(rast));
        sig2 = expt(ex).(sess)(hits(2)).(rast)(:,inds);
        mu{ex} = mean([sig1]);
        mu{ex} = mu{ex} - mean(mu{ex}(win<-11));
        
        Kraw    = @(t1,del,T) (T>=0).*1./(-del).*(exp(-T/t1) - exp(-T/(t1+del)));
        
        pkHtRaw = @(t1,del) (-del).^-1*(exp(-(t1+del)./del.*log((t1+del)./t1)) - exp(-t1./del.*log((t1+del)./t1)));
        K       = @(x,T) x(1)*Kraw(abs(x(2)), abs(x(3)), T+x(4)/10)./pkHtRaw(x(2),x(3));

        cost = @(x) sum((mu{ex} - K(x,win/11)).^2);
        x.(stim)(ex,:) = fminunc(cost,[max(mu{ex}),10,10,0]);
    end
    means.(stim) = cat(1,mu{:});
    hold on;
    h(find(strcmpi(stimlist,s))) = drawvar(win/11,means.(stim),cmap(count),1/sqrt(5));
end
legend(h,stimlist);
xlim([-15 45]);
set(gca,'xtick',-15:15:45);
ylim([-5 42]);
set(gca,'ytick',0:10:40);

%%
clear tau_rise tau_fall;
count=0;
for s = stimlist(1:end-1)
    count=count+1;
    tau_rise(:,count) = abs(x.(s{:})(mouseList,2));
    tau_fall(:,count) = abs(x.(s{:})(mouseList,2)) + abs(x.(s{:})(mouseList,3));
end
figure(11);clf;
bar(mean(tau_fall));
hold on;
for s = 1:4
    plot(s,tau_fall(:,s),'k.','markersize',15);
    errorbar(s,mean(tau_fall(:,s)),std(tau_fall(:,s))/sqrt(5),'k');
end
set(gca,'xtick',1:4,'xticklabel',stimlist(1:4));
box off
