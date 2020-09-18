%% User editable values
startTime=0;
stopTime=10^10;
xAxisLabel='Time (s)';
yAxisLabel='Fluorescent Intensity (dF/F)';
colorList={[0, 0.8255, 0.7401], [0.1965, 0, 0.5650]}'; %RGB color specifications corresponding to 490 and 405 nm, respectively
legendLabels={'490 nm', '405 nm'}';

baselineStart=19.00; %Time to start counting as baseline
baselineStop=49.00; %Time to stop counting as baseline
startChar='r';
stopChar='s';

photoCharacter='p';
rectFaceColors={[1,.9,.9];[.9,1,.9];[.9,.9,1];[1,1,.9];[1,.9,1];[.9,1,1]};

%% /User editable values



%b=csvread('C:\Users\R2\Desktop\Example_CSV_File.csv');
[fileName, pathName]=uigetfile('*.csv');
if ~isnumeric(fileName)
    if ~strcmpi(pathName(end),filesep)
        pathName=[pathName,filesep];
    end
    fullPath=[pathName, fileName];
else
    return
end
b=csvread(fullPath);






% [textFileName, pathName]=uigetfile('*.txt');
% if ~isnumeric(textFileName)
%     if strcmpi(pathName(end),filesep)
%         pathName=[pathName,filesep];
%     end
%     textFullPath=[pathName, textFileName];
% else
%     return
% end
% 
% 
% %c=textread('t.txt','%s');
% c=textread(textFullPath,'%s');
% c=reshape(c,[3,length(c)/3])';
% cTimes=c(:,2);
% cFlags=c(:,3);
% cMins=regexp(cTimes,'\d*(?=:)','match');
% cMins=cellfun(@str2num,vertcat(cMins{:}),'UniformOutput',true);
% cMinSecs=cMins*60;
% 
% 
% cSecs=regexp(cTimes,'(?<=:)\d*.\d*|(?<=:)\d*','match');
% cSecs=cellfun(@str2num,vertcat(cSecs{:}),'UniformOutput',true);
% cTimes=cSecs+cMinSecs;
% 
% ratStarts=~strcmpi(cFlags,stopChar);
% ratStops=strcmpi(cFlags,stopChar);
% startTimes=cTimes(ratStarts);
% stopTimes=cTimes(ratStops);
% 
% nRatTimes=min([sum(ratStarts),sum(ratStops)]);


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
    %startCharInd=find(strcmp(photoCharacter,charList),1,'first');
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


dffVals=b(:,2:end);


if length(legendLabels) < nColumns-1
    legendLabels=vertcat(legendLabels,repmat({' '},[nColumns-1-length(legendLabels),1]));
end

bMeans=mean(b(baselineInds,2:end),1);

for i=2:nColumns
    dffVals(:,i-1)=100*(b(:,i)-bMeans(i-1))/bMeans(i-1);
end



nRectColors=length(rectFaceColors);
%nRectTypes=length(unique(setdiff(bChars,photoCharacter)));

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

for i=1:nRects
    j=startNums(i);
    rectSet(i)=rectangle('Position',[startTimes(i),rectMinY,rectWidths(i),rectHeight],...
        'FaceColor',rectFaceColors{j},...
        'Parent', ax);
end


% rectMinY=min([min(min(dffVals))*2,0]);
% rectMaxY=2*max(max(b(:,2:end)));
% rectHeight=rectMaxY-rectMinY;
% rectWidths=stopTimes-startTimes;
% rectSet=zeros(nRatTimes,1);
% for i=1:nRatTimes
%     rectSet(i)=rectangle('Position',[startTimes(i),rectMinY,rectWidths(i),rectHeight],...
%         'FaceColor',[.8,.8,.8],...
%         'Parent', ax);
% end

nLines=nColumns-1;
hLine=zeros(nLines,1);

for i=2:nColumns
    j=i-1;
    if j <= nColorsListed
        %plot(b(keepInds,1),b(keepInds,i),'Color',colorList{i-1});
        hLine(j)=plot(b(keepInds,1),dffVals(keepInds,i-1),'Color',colorList{i-1});
    else
        %plot(b(keepInds,1),b(keepInds,i)); %If not enough colors listed, let it generate them automatically
        hLine(j)=plot(b(keepInds,1),dffVals(keepInds,i-1));
    end
end

%hLegend=legend(legendLabels{:});
xlabel(xAxisLabel);
ylabel(yAxisLabel);

yLim=get(ax,'YLim');
yLim(1)=max([yLim(1),rectMinY]);
yLim(2)=min([yLim(2),rectMaxY]);
set(ax,'YLim',yLim);





firstFigPos=get(f,'Position');
newFigPos=[firstFigPos(1)+firstFigPos(3),firstFigPos(2),firstFigPos(3)/3,firstFigPos(4)];

fLeg=figure('Position',newFigPos);
axLeg=axes('Parent',fLeg, 'Units', 'normalized','Position',[0,0,1,1]);
hold all;

nLineLabels=length(legendLabels);

nFullLegendLabels=nLines+nRectTypesUsed;
set(axLeg,'XLim',[0,1],'YLim',[0,nFullLegendLabels]);

[hText,legendHandles]=deal(zeros(nFullLegendLabels,1));

sxText=0.025;
sy=nFullLegendLabels+0.5;
for i=1:nFullLegendLabels
    sy=sy-1;
    if i <= nLines
        if i <= nLineLabels
            hText(i)=text(sxText, sy, legendLabels{i});
        else
            hText(i)=text(sxText, sy, ' ');
        end
    else
        hText(i)=text(sxText, sy, behaviorLabels(numsUsedSet(i-nLines)));
    end
end

textExtents=get(hText,'Extent');
textExtents=vertcat(textExtents{:});

sxLine=max(textExtents(:,1)+textExtents(:,3))+0.025;

rectHeight=0.48;
sy=nFullLegendLabels+0.5;
% sxLine=0.525;
lineMax=1;
rectWidth=lineMax-sxLine;

for i=1:nLines
    sy=sy-1;
    legendHandles(i)=line([sxLine,lineMax],[sy,sy],'Color',get(hLine(i),'Color'));
end

sy=sy-rectHeight/2;
for i=1:nRectTypesUsed
    j=i+nLines;
    sy=sy-1;
    if i <= nFullLegendLabels
        legendHandles(i)=rectangle('Position',[sxLine, sy, rectWidth, rectHeight],...
            'FaceColor',rectFaceColors{i});
    end
end

[pathIn, fileIn, extIn]=fileparts(fullPath);

if ~strcmpi(pathIn(end),filesep)
    pathIn=[pathIn,filesep];
end

fullPathOut=[pathIn,fileIn,'_DFF',extIn];

%     [ofn, outFilePath]=uiputfile({...
%         '*.csv', 'Comma Separated Value (*.csv)';
%         '*.*', 'All Files (*.*)'},...
%         'DF/F CSV file to export', pathName);
%    fullPathOut=[outFilePath,ofn];

csvwrite(fullPathOut,[b(:,1),dffVals]);