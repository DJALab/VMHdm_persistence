
stims   = {'rat','male','toy','USS','tone'};
labels  = 1:length(stims);
nstim   = length(stims);
exset   = [1 3 4 5 7];
nTr     = 6;
decodeFromPeak = 0;

clear guess;
for ex = exset
    disp(ex)
    means = [];rMat=[];
    clear resps;
    % for each mouse, compute stim-evoked responses for each trial on each
    % day of imaging
    for s = 1:length(stims)
        count=0;
        for d = 1:3
            day = ['session' num2str(d)];
            for tr = 1:2
                count   = count+1;
                r       = expt(ex).resps.(day).rast.(stims{s}){tr};
                sig     = bsxfun(@minus,r,mean(r(:,1:200),2)); % subtract pre-stimulus baseline
                if decodeFromPeak %this works pretty well too
                    [pks,locs]   = max(smoothts(sig,'g',5,1),[],2);
                    resp = pks;
                else
                    resp = mean(sig(:,201:end),2);
                end

                resps.(stims{s})(:,count) = resp;
            end
        end
        rMat = cat(2,rMat,resps.(stims{s}));
    end

    trialLabels = kron(labels,ones(1,nTr)); % the ground-truth labels
    guess{ex}   = zeros(1,size(rMat,2));    % where we'll store the decoder's guesses
    
    % fit Naive bayes classifier
    Mdl = fitcnb(rMat',trialLabels,'CrossVal','on','Leaveout','on');
    
    % predict (using cross-validation)
    guess{ex}(1,:) = kfoldPredict(Mdl);
end


% make a bar plot of decoder accuracy
gt = kron(labels,ones(1,nTr));
figure(2);clf;
for s = 1:nstim % loop computes the decoder accuracy for each stimulus class separately
    inds = (1:nTr) + (s-1)*nTr; %indices corresponding to stim s
    acc = [];
    for ex = exset % get accuracy for each imaged mouse for this stimulus
        acc = [acc mean(guess{ex}(:,inds)==gt(:,inds))];
    end
    
    bar(s,mean(acc))
    hold on;
    plot(s+randn(1,5)/10,acc,'k.','markersize',15)
    errorbar(s,mean(acc),std(acc)/sqrt(length(acc)),'k');
end
plot(get(gca,'xlim'),[1 1]/length(unique(gt(:))),'k--'); %dashed line = chance
set(gca,'xtick',labels,'xticklabel',stims);
ylim([0 1]);


% look at the confusion matrix between classes
confusion=zeros(length(stims));
for ex = exset
    for i=1:length(stims)
        confusion(i,:) = confusion(i,:) + hist(guess{ex}(gt==i),1:length(stims));
    end
end
figure(4);image(64-confusion/(length(exset)*nTr)*64);
set(gca,'ytick',labels,'yticklabel',stims);
colormap gray
colorbar
disp(confusion/length(exset)/nTr*100);

