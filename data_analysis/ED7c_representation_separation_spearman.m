
stims = {'rat','male','toy','USS','tone'};
figure(1);clf;
allRast = [];
phi = [];
exran = [1 3 4 5 7];

sess = 'session3';sess2='session2';
type = 'rast';
fun = @(x) mean(x,2);
meas = @(s1,s2) corr(s1,s2,'type','pearson');
% meas = @(s1,s2) pdist2(s1',s2','correlation');

for ex = exran
    temp=zeros(5);
    for i = 1:5
        for j = setdiff(1:5,i)
            for tr = 1:2
                s1      = fun(expt(ex).resps.(sess).(type).(stims{i}){tr});
                s2      = fun(expt(ex).resps.(sess2).(type).(stims{j}){tr});
                temp(i,j) = temp(i,j) + meas(s1,s2)/4;
                s1      = fun(expt(ex).resps.(sess2).(type).(stims{i}){tr});
                s2      = fun(expt(ex).resps.(sess).(type).(stims{j}){tr});
                temp(i,j) = temp(i,j) + meas(s1,s2)/4;
            end
        end
        for tr = 1:2
            s1      = fun(expt(ex).resps.(sess).(type).(stims{i}){tr});
            s2      = fun(expt(ex).resps.(sess2).(type).(stims{i}){tr});
            temp(i,i) = temp(i,i) + meas(s1,s2)/2;
        end
    end
    phi     = [phi temp(:)];
end

bar(mean(phi,2));
hold on;
errorbar(mean(phi,2),std(phi,[],2)/sqrt(5),'k+');
plot(phi,'k.','markersize',15);
for i=1:4
    plot([i i]*5+.5,[-.5 1],'r');
end
xlim([0 26]);
ylim([-.5 1]);
box off;

