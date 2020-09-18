%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User editable values
startTime=0;
stopTime=10^10;
xAxisLabel='Time (s)';
yAxisLabel='Fluorescent Intensity (% dF/F)';
colorList={[0, 0.8255, 0.7401], [0.1965, 0, 0.5650]}'; %RGB color specifications corresponding to 490 and 405 nm, respectively
legendLabels={'490 nm', '405 nm'}';

baselineStart=17.56; %Time to start counting as baseline
baselineStop=47.56; %Time to stop counting as baseline
startChar='r';
stopChar='s';

%Start and stop time for the part of the curve to find the area underneath:
% integralStartTime=202.0;
% integralStopTime=232.0;

photoCharacter='p';
rectFaceColors={[1,.9,.9];[.9,1,.9];[.9,.9,1];[1,1,.9];[1,.9,1];[.9,1,1]};
root = 'C:\Users\Ann\Desktop\Research\Prabhat\';
% root = 'C:\Users\R2\Desktop\MatlabTools\PrabhatStuff\PlotCSV\a\';
timeColumn=1;
controlColumn=3; %In this case, the column with the 405 nm intensity

stimTimeOffset=0; %Number of seconds stimulus time is ahead of recording
%time, e.g. if the stimulus started 17.1 seconds after the recording, then
%stimTimeOffset=-17.1. If the stimulus started 4.8 seconds before the
%luminosity recording, then stimTimeOffset=4.8




%% /User editable values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('pathName', 'var') || isnumeric(pathName) || ~exist(pathName, 'dir')
    pathName=root;
end

if ~strcmpi(pathName(end),filesep)
    pathName=[pathName,filesep];
end

[fileName, pathName]=uigetfile([pathName,'*.csv']);
if ~isnumeric(fileName)
    if ~strcmpi(pathName(end),filesep)
        pathName=[pathName,filesep];
    end
    fullPath=[pathName, fileName];
else
    return
end
b=csvread(fullPath);

[matFileName, pathName]=uigetfile('*.*');
if ~isnumeric(matFileName)
    if ~strcmpi(pathName(end),filesep)
        pathName=[pathName,filesep];
    end
    matFullPath=[pathName, matFileName];
else
    return
end

d=load(matFullPath,'-mat');
bChars=d.uc.bchar;
eChars=d.uc.echar;
behaviorLabels=d.uc.desc;
timeList=vertcat(d.fulllog.time{:})*60;
charList=d.rawlog;
bNumList=vertcat(d.fulllog.charnum{:});
bActionList=d.fulllog.action;

nActions=length(bActionList);

keepCharInds=true(nActions,1);
removeNextStop=false;


if ismember(photoCharacter, bChars)
    photoCharNum=find(strcmp(photoCharacter,bChars));
    if ismember(photoCharNum,bNumList)
        startCharInd=find(bNumList==photoCharNum,1,'first');
        timeList=timeList-timeList(startCharInd);
    end
    stopPhotoCharacter=eChars(photoCharNum);

    keepCharInds=bNumList~=photoCharNum;
end

timeList=timeList(keepCharInds);
%charList=charList(keepCharInds);
bNumList=bNumList(keepCharInds);
bActionList=bActionList(keepCharInds);

nActions=length(bActionList);




startTimeInds=find(strcmpi('start',bActionList));
stopTimeInds=find(strcmpi('stop', bActionList));

startTimes=timeList(startTimeInds);
stopTimes=timeList(stopTimeInds);

startNums=bNumList(startTimeInds);

nRects=length(startTimeInds);
%%



keepIndsLogical=b(:,1)>=startTime & b(:,1)<=stopTime;
keepInds=find(keepIndsLogical);

baselineIndsLogical=b(:,1) >= baselineStart & b(:,1) <= baselineStop;
baselineInds=find(baselineIndsLogical);

f=figure;
ax=axes;
hold all;

nColumns=size(b,2);
nColorsListed=length(colorList);

bCols=1:size(b,2);
dffCols=setdiff(bCols,[timeColumn, controlColumn]);
[dffVals, fitVals]=deal(b(:,dffCols));

%dffVals=b(:,2:end);
nColsRemoved=double(logical(timeColumn))+double(logical(controlColumn));

bMeans=mean(b(baselineInds,2:end),1);

nColumnsShown=length(dffCols);
fitMat=zeros(nColumnsShown,2);
p=pinv([b(baselineInds,controlColumn), ones(sum(baselineIndsLogical),1)]);
for i=1:nColumnsShown
    j=dffCols(i);
    [fitMat(i,:)]=p*b(baselineInds,j);
    fitVals(:,i)=fitMat(i,1)*b(:,controlColumn)+fitMat(i,2);
    dffVals(:,i)=(dffVals(:,i)-fitVals(:,i))./fitVals(:,i);
end
dffVals=dffVals*100;


if length(legendLabels) < nColumns-nColsRemoved
    legendLabels=vertcat(legendLabels,repmat({' '},[nColumns-nColsRemoved-length(legendLabels),1]));
end

nRectColors=length(rectFaceColors);

numsUsedSet=unique(bNumList);
nRectTypesUsed=length(numsUsedSet);

while nRectTypesUsed > nRectColors
    rectFaceColors=vertcat(rectFaceColors,{rand([1,3])});
    nRectColors=nRectColors+1;
end


rectMinY=min([min(min(dffVals))*2,0]);
rectMaxY=2*max(max(dffVals));
rectHeight=rectMaxY-rectMinY;
rectWidths=stopTimes-startTimes;
rectSet=zeros(nRects,1);

startTimes=startTimes-stimTimeOffset;

for i=1:nRects
    j=startNums(i);
    rectSet(i)=rectangle('Position',[startTimes(i),rectMinY,rectWidths(i),rectHeight],...
        'FaceColor',rectFaceColors{j},...
        'Parent', ax);
end



nLines=nColumns-1;
hLine=zeros(nLines,1);

for i=1:nColumnsShown
    if i <= nColorsListed
        hLine(i)=plot(b(keepInds,timeColumn),dffVals(keepInds,i),'Color',colorList{i});
    else
        %If not enough colors listed, let it generate them automatically
        hLine(i)=plot(b(keepInds,timeColumn),dffVals(keepInds,i));
    end
end

%Plot a horizontal dashed line at y=0
line('Parent', ax, 'XData',[min(b(:,timeColumn)),max(b(:,timeColumn))],'YData',[0,0],'LineStyle','--' ,'Color', [0,0,0]);

hLegend=legend(legendLabels{:});
xlabel(xAxisLabel);
ylabel(yAxisLabel);

yLim=get(ax,'YLim');
yLim(1)=max([yLim(1),rectMinY]);
yLim(2)=min([yLim(2),rectMaxY]);
set(ax,'YLim',yLim);

dff_integrator(f,ax,b(:,1),dffVals);


fullPathOut=[fullPath(1:end-4),'_dff.csv'];
csvwrite(fullPathOut,[b,fitVals,dffVals]);