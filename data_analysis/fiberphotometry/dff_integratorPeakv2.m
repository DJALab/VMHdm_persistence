function dff_integrator(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User editable values
% global h;

h.defaultFrameRate=24; %Default frame rate (frames per second) to assume for annotation files from behavior annotator
h.frameRate=h.defaultFrameRate;

h.legendLabels={'490 nm', '405 nm'}';
h.colorList={[0, 0.8255, 0.7401], [0.1965, 0, 0.5650]}'; %RGB color specifications corresponding to 490 and 405 nm, respectively
h.rectFaceColors={[1,.9,.9];[.9,1,.9];[.9,.9,1];[1,1,.9];[1,.9,1];[.9,1,1]};
h.rectFaceColors2={[.7,.9,.9];[.9,.7,.9];[.9,.9,.7];[.7,.7,1];[.7,1,.7];[1,.7,.7]};
h.stimStartCharacter='p'; %Character in .mat indicating start of stimulus
h.stimBehaviorLabel='light';
h.timeColumn=1; %Column containing time data
h.rawColumn=2; %Main raw data column, in this case 490 nm
h.controlColumn=3; %In this case, the column with the 405 nm intensity
h.dffColumn=2*max([h.controlColumn,h.rawColumn])-1;
h.stimTimeOffset=0; %Number of seconds stimulus time is ahead of recording
%time, e.g. if the stimulus started 17.1 seconds after the recording, then
%h.stimTimeOffset=-17.1. If the stimulus started 4.8 seconds before the
%luminosity recording, then h.stimTimeOffset=4.8
h.defaultYLimInd=1; %Set this value to determine where the minimum of the
                  %y-axis should go: 1 - min at 0, 2 - min at minimum value
                  %of t-series, 3 - min at middle quarter of t-series
h.offsetVal=0;    %Default OffsetValue

h.rawColor=[0, 0.8255, 0.7401];
h.controlColor=[0.1965, 0, 0.5650];
h.dffColor=[0,0,1];

h.defaultBaselineRange=[40,75];
h.defaultIntegrationRange=[120,160];

h.defaultXLim=[0,250];
h.defaultDffYLim=[-50,50];
h.defaultVelocityYLim=[0,250];

h.velColor=[0,0.5,0]; %Color of velocity y-axis and curve
h.dffColor=[0,0,0.85]; %Color of DFF y-axis and curve

h.smoothingWindow=3;
h.velocityThreshold=50;
h.velocityMinSpan=3;
h.doUseVelocityThresholds=false;
h.doSmoothVelocity=false;

h.fitColor=[0.8,0.4,0];

%% /User editable values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if nargin >1 
    h.f=varargin{1};
    h.ax=varargin{2};
    h.times=varargin{3};
    h.dff=varargin{4};
else
    if ~exist('h.pathName', 'var')
        h.pathName=[pwd, filesep];
    end
    [fileName, h.pathName]=uigetfile([h.pathName,'*.csv']);
    if ~isnumeric(fileName)
        if ~strcmpi(h.pathName(end),filesep)
            h.pathName=[h.pathName,filesep];
        end
        fullPath=[h.pathName, fileName];
    else
        return
    end
    h.b=csvread(fullPath);
    
    h.baselineFig=figure();
    h.baselineAx=axes('Parent',h.baselineFig);
    set(h.baselineFig,'Visible','off');
    drawnow;    
    h.f=figure();
    h.ax=axes();
    
    
    h.times=h.b(:,h.timeColumn);
    if size(h.b,2)>=h.dffColumn
        h.dff=h.b(:,h.dffColumn);
    elseif size(h.b,2)==h.dffColumn
        %Allow loading of file without dff column, and generate dff column
        %with baseline set by user
        [h,useOldBaseline]=getDff(h);   
        
    end
    h.zeroLine=line('Parent', h.ax, 'XData',[min(h.times),max(h.times)],'YData',[0,0],'LineStyle','--' ,'Color', [0,0,0]);
    h.zeroLineXData=get(h.zeroLine,'XData');
    h.trueZeroLineXData=h.zeroLineXData;
    h.raw=h.b(:,h.rawColumn);
    h.control=h.b(:,h.controlColumn);
    %%%h.fRawPlot=plot(h.times,h.raw,'Parent',h.ax,'Color',h.rawColor);
    %%%h.fControlPlot=plot(h.times,h.raw,'Parent',h.ax,'Color',h.controlColor);
    drawnow;
    hold all;
    h.dffPlot=plot(h.times, h.dff, 'Parent', h.ax,'Color',h.dffColor);
    
    set(h.ax,'XLimMode','manual');
end
h.rectPosList=[];
h.stimRects=[];
h.integrationRects=[]; 
h.fitLines=[];

h.integralRects=[];
h.fitLine=[];


h.shiftPosVec=@(posVec,xShift) posVec+[xShift,0,0,0];

fPos=get(h.f,'Position');
nfWidth=350;
nfHeight=650;
nfX=fPos(1)+fPos(3)+10;
nfY=fPos(2)+fPos(4)-nfHeight; 

% set(h.integrationFig,'CloseRequestFcn', @closeIntegrationFig);


minTime=min(h.times);
maxTime=max(h.times);
minDff=min(h.dff);
maxDff=max(h.dff);


defaultBaselineMin=max([minTime,h.defaultBaselineRange(1)]);
defaultBaselineMax=min([maxTime,h.defaultBaselineRange(2)]);
defaultIntegrationMin=max([minTime,h.defaultIntegrationRange(1)]);
defaultIntegrationMax=min([maxTime,h.defaultIntegrationRange(2)]);

h.defaultBaselineRange=[defaultBaselineMin,defaultBaselineMax];
h.defaultIntegrationRange=[defaultIntegrationMin, defaultIntegrationMax];


h.peakFound=false;
h.peakVisible=true;
h.peakPoint=scatter(0,0,'d', 'filled');
set(h.peakPoint,'Visible','off');
h.xPeak=[];
h.yPeak=[];



%Create user interface control objects
h.integrationFig=figure(...
    'MenuBar', 'none',...
    'Name', 'Integrator',...
    'NumberTitle', 'off',...    
    'Resize', 'on',...
    'Toolbar', 'none',...
    'Units', 'pixels',...
    'Visible', 'on',...
    'Position', [nfX, nfY, nfWidth, nfHeight]);



sp=5;
ssp=2;
bw=70;
bh=20;
ew=80;

sx=sp;
sy=nfHeight-sp-bh;



sx=sx+sp+bw;
h.startTx=uicontrol('Parent', h.integrationFig, 'Style', 'text',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh], 'String', 'Start (s):');
sx=sx+sp+bw;
h.stopTx=uicontrol('Parent', h.integrationFig, 'Style', 'text',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh], 'String', 'Stop (s):');
sx=sx+sp+bw;
h.showLinesCb=uicontrol('Parent', h.integrationFig, 'Style', 'checkbox',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String', 'Show:', 'Value', true);

sy=sy-sp-bh;
sx=sp;

h.baselineTx=uicontrol('Parent', h.integrationFig, 'Style', 'text',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh], 'String', 'Baseline:');
sx=sx+sp+bw;
h.baselineStartEd=uicontrol('Parent', h.integrationFig, 'Style', 'edit',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String', num2str(h.defaultBaselineRange(1)),...
    'UserData',h.defaultBaselineRange(1));
sx=sx+sp+bw;
h.baselineStopEd=uicontrol('Parent', h.integrationFig, 'Style', 'edit',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String', num2str(h.defaultBaselineRange(2)),...
    'UserData',h.defaultBaselineRange(2));

sy=sy-ssp-bh;
sx=sp;
h.integrationTx=uicontrol('Parent', h.integrationFig, 'Style', 'text',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh], 'String', 'Integration:');
sx=sx+sp+bw;
h.integrationStartEd=uicontrol('Parent', h.integrationFig, 'Style', 'edit',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String', num2str(h.defaultIntegrationRange(1)),...
    'UserData',h.defaultIntegrationRange(1));
sx=sx+sp+bw;
h.integrationStopEd=uicontrol('Parent', h.integrationFig, 'Style', 'edit',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String', num2str(h.defaultIntegrationRange(2)),...
    'UserData',h.defaultIntegrationRange(2));


sx=sx+bw+sp;
h.baselineOrIntegrationBg=uibuttongroup('Parent', h.integrationFig,...
    'Units', 'pixels', 'Position', [sx, sy, bw, 2*bw+ssp],...
    'Title', '','BorderType', 'none');
h.baselineRb=uicontrol('Parent',h.baselineOrIntegrationBg,...
    'Style', 'radiobutton', 'Value', false, 'Units', 'pixels',...
    'String', 'Baseline', 'Position', [0,bh+ssp,bw,bh]);
h.integrationRb=uicontrol('Parent',h.baselineOrIntegrationBg,...
    'Style', 'radiobutton', 'Value', true, 'Units', 'pixels',...
    'String', 'Integration', 'Position', [0,0,bw,bh]);

if get(h.baselineRb,'Value')
    defaultStartLineXData=[h.defaultBaselineRange(1),h.defaultBaselineRange(1)];
    defaultStopLineXData=[h.defaultBaselineRange(2),h.defaultBaselineRange(2)];
else
    defaultStartLineXData=[h.defaultIntegrationRange(1),h.defaultIntegrationRange(1)];
    defaultStopLineXData=[h.defaultIntegrationRange(2),h.defaultIntegrationRange(2)];
end

h.startIndicatorLine=line('Parent',h.ax,'LineStyle', '--', 'Color', [0,0.7,0],...
    'XData',defaultStartLineXData, 'YData',[minDff, maxDff]);
h.stopIndicatorLine=line('Parent',h.ax,'LineStyle', '--', 'Color', [0.7,0,0] ,...
    'XData',defaultStopLineXData, 'YData',[minDff, maxDff]);


sx=sp;
sy=sy-2*sp-bh;

h.allXTx=uicontrol('Parent', h.integrationFig, 'Style', 'text',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh], 'String', 'X Limits:');
sx=sx+sp+bw;
h.allXLimStartEd=uicontrol('Parent', h.integrationFig, 'Style', 'edit',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String',num2str(h.defaultXLim(1)),...
    'UserData', h.defaultXLim(1));
sx=sx+sp+bw;
h.allXLimStopEd=uicontrol('Parent', h.integrationFig, 'Style', 'edit',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String',num2str(h.defaultXLim(2)),...
    'UserData', h.defaultXLim(2));


sy=sy-ssp-bh;
sx=sp;

h.dffYTx=uicontrol('Parent', h.integrationFig, 'Style', 'text',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh], 'String', 'Y (DFF):');
sx=sx+sp+bw;
h.dffYLimStartEd=uicontrol('Parent', h.integrationFig, 'Style', 'edit',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String',num2str(h.defaultDffYLim(1)),...
    'UserData', h.defaultDffYLim(1));
sx=sx+sp+bw;
h.dffYLimStopEd=uicontrol('Parent', h.integrationFig, 'Style', 'edit',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String',num2str(h.defaultDffYLim(2)),...
    'UserData', h.defaultDffYLim(2));



sy=sy-ssp-bh;
sx=sp;

h.velYTx=uicontrol('Parent', h.integrationFig, 'Style', 'text',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh], 'String', 'Y (Velocity):');
sx=sx+sp+bw;
h.velYLimStartEd=uicontrol('Parent', h.integrationFig, 'Style', 'edit',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String',num2str(h.defaultVelocityYLim(1)),...
    'UserData', h.defaultVelocityYLim(1));
sx=sx+sp+bw;
h.velYLimStopEd=uicontrol('Parent', h.integrationFig, 'Style', 'edit',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String',num2str(h.defaultVelocityYLim(2)),...
    'UserData', h.defaultVelocityYLim(2));


sy=sy-2*sp-bh;
sx=sp;

h.addVelocityPb=uicontrol('Parent', h.integrationFig,'Style', 'pushbutton',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String', 'Add Velocity');

sx=sx+sp+bw;

h.hideVelocityTb=uicontrol('Parent', h.integrationFig,...
    'Style', 'togglebutton', 'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
    'String', 'Hide Velocity');

sy=sy-sp-bh;
sx=sp;


h.useVelocityThresholdsCh=uicontrol('Parent', h.integrationFig,...
    'Style', 'checkbox', 'Units',' pixels', 'Position', [sx, sy, 1.6*bw, bh],...
    'String', 'Exclude Velocities',...
    'Value', h.doUseVelocityThresholds);

sx=sx+2*sp+1.6*bw;
h.velocityThresholdsTx=uicontrol('Parent', h.integrationFig,...
    'Style', 'text', 'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
    'String', 'greater than');

sx=sx+bw;
h.velocityThresholdEd=uicontrol('Parent', h.integrationFig,...
    'Style', 'edit', 'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
    'String', h.velocityThreshold);

% % % % % sy=sy-sp-bh;
% % % % % sx=sp+2*sp+1.6*bw;
% % % % % h.frameLengthThresholdTx=uicontrol('Parent', h.integrationFig,...
% % % % %     'Style', 'text', 'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
% % % % %     'String', 'For less than ');
% % % % % 
% % % % % sx=sx+bw;
% % % % % h.frameLengthThresholdEd=uicontrol('Parent', h.integrationFig,...
% % % % %     'Style', 'togglebutton', 'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
% % % % %     'String', h.velocityMinSpan);
% % % % % 
% % % % % sx=sx+bw;
% % % % % h.frameLengthThresholdTx2=uicontrol('Parent', h.integrationFig,...
% % % % %     'Style', 'text', 'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
% % % % %     'String', ' frames.');

sy=sy-sp-bh;
sx=sp;
h.smoothVelocityCh=uicontrol('Parent', h.integrationFig,...
    'Style', 'checkbox', 'Units', 'pixels', 'Position', [sx, sy, 1.6*bw, bh],...
    'String', 'Smooth Velocity',...
    'Value', h.doSmoothVelocity);

sx=sx+2*sp+1.6*bw;
h.smoothingWindowTx=uicontrol('Parent', h.integrationFig,...
    'Style', 'text', 'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
    'String', 'Window:');

sx=sx+bw;
h.smoothingWindowEd=uicontrol('Parent', h.integrationFig,...
    'Style', 'edit', 'Units', 'pixels', 'Position', [sx, sy, bw/2, bh],...
    'String', h.smoothingWindow);


sy=sy-2*sp-bh;
sy=sy-2*sp-bh;
sx=sp;

bw2=bw/2;

h.addStimuliPb=uicontrol('Parent', h.integrationFig, 'Style', 'pushbutton',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String', 'Add Stimuli');
% sx=sx+bw+sp;
sx=sx+bw+sp;
h.firstAnnotationYPosTx=uicontrol('Parent', h.integrationFig,...
    'Style', 'text', 'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
    'String', 'Y-position:');
sx=sx+bw+sp;
h.firstAnnotationYPosEd=uicontrol('Parent', h.integrationFig,...
    'Style', 'edit', 'Units', 'pixels', 'Position', [sx, sy, bw2, bh],...
    'String', '');

sx=sx+bw2+sp;
h.firstAnnotationHeightTx=uicontrol('Parent', h.integrationFig,...
    'Style', 'text', 'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
    'String', 'Height:');
sx=sx+bw+sp;
h.firstAnnotationHeightEd=uicontrol('Parent', h.integrationFig,...
    'Style', 'edit', 'Units', 'pixels', 'Position', [sx, sy, bw2, bh],...
    'String', '');

sx=sp;
sy=sy-sp-bh;


h.alignToStimuliCb=uicontrol('Parent', h.integrationFig,...
    'Style', 'checkbox', 'Units', 'pixels', 'Position', [sx,sy,ew,bh],...
    'String', 'Align to stim.', 'Value', true);


sx=sx+ew+sp;
h.offsetTx=uicontrol('Parent', h.integrationFig, 'Style', 'text', 'Units',...
    'pixels', 'Position', [sx, sy, bw, bh], 'String', 'Offset:');
sx=sx+bw+sp;
h.offsetEd=uicontrol('Parent', h.integrationFig, 'Style', 'edit', 'Units',...
    'pixels', 'Position', [sx, sy, bw, bh], 'String', '0');

sx=sp;
sy=sy-sp-bh;

h.annotationFrameRateTx=uicontrol('Parent', h.integrationFig, 'Style',...
    'text', 'Units', 'pixels', 'Position', [sx, sy, 2*bw, bh], 'String',...
    'Annotation Frame Rate:');
sx=sx+2*bw+sp;
h.annotationFrameRateEd=uicontrol('Parent', h.integrationFig, 'Style',...
    'edit', 'Units', 'pixels', 'Position', [sx, sy, bw, bh], 'String',...
    num2str(h.defaultFrameRate));

sx=sx+bw;
h.annotationFrameRateUnitsTx=uicontrol('Parent', h.integrationFig,...
    'Style', 'text', 'Units', 'pixels', 'Position', [sx, sy, 1.4*bw, bh],...
    'String', 'frames/second');

sy=sy-sp-bh;
sx=sp;

h.addSecondAnnotationPb=uicontrol('Parent', h.integrationFig, 'Style', 'pushbutton',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh],...
    'String', 'Add 2nd Ann.');


sx=sx+bw+sp;    
% sy=sy-sp-bh;

h.secondAnnotationYPosTx=uicontrol('Parent', h.integrationFig,...
    'Style', 'text', 'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
    'String', 'Y-position:');
sx=sx+bw+sp;
h.secondAnnotationYPosEd=uicontrol('Parent', h.integrationFig,...
    'Style', 'edit', 'Units', 'pixels', 'Position', [sx, sy, bw2, bh],...
    'String', '');

sx=sx+bw2+sp;
h.secondAnnotationHeightTx=uicontrol('Parent', h.integrationFig,...
    'Style', 'text', 'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
    'String', 'Height:');
sx=sx+bw+sp;
h.secondAnnotationHeightEd=uicontrol('Parent', h.integrationFig,...
    'Style', 'edit', 'Units', 'pixels', 'Position', [sx, sy, bw2, bh],...
    'String', '');


sx=sp;
sy=sy-2*sp-bh;

h.integralPb=uicontrol('Parent', h.integrationFig, 'Style', 'pushbutton',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh], 'String', 'Integrate');
sx=sx+bw+sp;
h.integralTx=uicontrol('Parent', h.integrationFig, 'Style', 'text', 'Units',...
    'pixels', 'Position', [sx,sy,ew,bh], 'String', '-');

sx=sx+ew+sp;
h.meanLabelTx=uicontrol('Parent', h.integrationFig, 'Style', 'text', 'Units',...
    'pixels', 'Position', [sx, sy, bw, bh], 'String', 'Mean:');
sx=sx+bw+sp;
h.meanTx=uicontrol('Parent', h.integrationFig, 'Style', 'text', 'Units',...
    'pixels', 'Position', [sx, sy, ew, bh], 'String', '-');
sy=sy-bh-sp;
sx=sx-bw-sp;
h.semLabelTx=uicontrol('Parent', h.integrationFig, 'Style', 'text', 'Units',...
    'pixels', 'Position', [sx, sy, bw, bh], 'String', 'SEM:');
sx=sx+bw+sp;
h.semTx=uicontrol('Parent', h.integrationFig, 'Style', 'text', 'Units',...
    'pixels', 'Position', [sx, sy, ew, bh], 'String', '-');

sx=sp;
sy=sy-bh-2*sp;

h.findPeakPb=uicontrol('Parent', h.integrationFig, 'Style', 'pushbutton',...
    'Units', 'pixels', 'Position', [sx,sy,bw,bh], 'String', 'Find peak');
sx=sx+bw+sp;
h.findPeakTx=uicontrol('Parent', h.integrationFig, 'Style', 'text', 'Units',...
    'pixels', 'Position', [sx,sy,2*bw,bh], 'String', '[?,?]');
sx=sx+2*bw+sp;
h.showPeakCb=uicontrol('Parent', h.integrationFig, 'Style', 'checkbox',...
    'Units', 'pixels', 'Position', [sx, sy, 1.8*bw, bh], 'String', 'Show peak mark',...
    'Value', true);

sx=sp;
sy=sy-2*sp-1.5*bh;
h.fitPb=uicontrol('Parent', h.integrationFig, 'Style', 'pushbutton',...
    'Units', 'pixels', 'Position', [sx, sy, bw, bh], 'String', 'Fit');
sx=sx+bw+sp;
h.fitTx=uicontrol('Parent', h.integrationFig, 'Style', 'pushbutton',...
    'Enable', 'off',...
    'Units', 'pixels', 'Position', [sx, sy, 3*ew, 1.5*bh],...
    'FontSize', 10,...
    'String', '<html>r(t)=(&tau;<sub>d</sub> - &tau;<sub>s</sub>)<sup>-1</sup>[exp(-t/&tau;<sub>d</sub>) - exp(-t/&tau;<sub>s</sub>)]</html>');
% h.fitTx=text(sx,sy,'$\frac{1}{\tau_d-\tau_s}$');

%,...
%    'BackgroundColor',get(h.integrationFig, 'Color'));
sy=sy-sp-bh;
sx=sp;
h.fitTxA=uicontrol('Parent', h.integrationFig, 'Style', 'text', 'Units',...
    'pixels', 'Position', [sx, sy, 2*ew, bh], 'String', 'T_r=');
sx=sx+2*ew+sp;
h.fitTxB=uicontrol('Parent', h.integrationFig, 'Style', 'text', 'Units',...
    'pixels', 'Position', [sx, sy, 2*ew, bh], 'String', 'T_f=');
% sx=sx+ew+sp;
% h.fitTxC=uicontrol('Parent', h.integrationFig, 'Style', 'text', 'Units',...
%     'pixels', 'Position', [sx, sy, ew, bh], 'String', 'c=');


sy=sy-2*sp-bh;
sx=sp;

h.clearPb=uicontrol('Parent', h.integrationFig, 'Style', 'pushbutton',...
    'Units', 'pixels', 'Position', [sx, sy, bw, bh],...
    'String', 'Clear');
sx=sx+bw+sp;
h.preserveStimuliCb=uicontrol('Parent', h.integrationFig, 'Style', 'checkbox',...
    'Units', 'pixels', 'Position', [sx, sy, 1.5*ew, bh],...
    'String', 'Keep Stimuli', 'Value', true);
sy=sy-sp-bh;
h.preserveVelocityCb=uicontrol('Parent', h.integrationFig, 'Style', 'checkbox',...
    'Units', 'pixels', 'Position', [sx, sy, 1.5*ew, bh],...
    'String', 'Keep Velocity', 'Value', true);




h.trueTimes=h.times; 
[h.trueRectXPos, h.offsetRectXPos, h.centeredRectXPos]=deal([]);
h.indicatorLines=[h.startIndicatorLine, h.stopIndicatorLine];

%Set callbacks for objects


% % set(h.showLinesCb, 'Callback', {@showLinesCb_callback,h});
% % set(h.integrationStopEd, 'Callback', {@callback_integrationEd,h.startIndicatorLine});
% % set(h.integrationStartEd, 'Callback', {@callback_integrationEd, h.stopIndicatorLine});
% % set(h.offsetEd, 'Callback', {@offsetEd_callback, h});
% % set(h.integralPb, 'Callback',{@callback_integralPb, h});
% % set(h.fitPb, 'Callback',{@callback_fitExponentialPb, h});
% % set(h.clearPb, 'Callback', {@callback_clearPb, h});
% % set(h.addStimuliPb, 'Callback', {@callback_addStimuliPb, h});
% % set(h.alignToStimuliCb, 'Callback', {@callback_alignToStimuliCb, h});
% % 
% % set(h.integrationFig, 'CloseRequestFcn', {@closeWindowFcn, h.f});
% % set(h.f, 'CloseRequestFcn', {@closeWindowFcn, h.integrationFig});

% if get(h.showLinesCb,'Value')
%     set(h.f, 'WindowButtonDownFcn',{},...
%         'WindowButtonMotionFcn',{@unclickedMotion, h.ax, h.indicatorLines},...
%         'WindowButtonUpFcn',{});
% else
%         set(h.f, 'WindowButtonDownFcn',{},...
%         'WindowButtonMotionFcn',{},...
%         'WindowButtonUpFcn',{});
% end

guidata(h.integrationFig,h);
h=setCallbacks(h);
guidata(h.integrationFig,h);



function h=setCallbacks(h)
h=guidata(h.integrationFig);
set(h.showLinesCb, 'Callback', {@showLinesCb_callback,h});
set(h.baselineOrIntegrationBg, 'SelectionChangeFcn', {@baselineOrIntegration_callback, h});
% set(h.baselineRb, 'Callback', {baselineOrIntegration_callback, h});
% set(h.integrationRb, 'Callback', {baselineOrIntegration_callback, h});

set(h.integrationStopEd, 'Callback', {@callback_integrationEd, h, h.stopIndicatorLine});
set(h.integrationStartEd, 'Callback', {@callback_integrationEd, h, h.startIndicatorLine});
set(h.baselineStopEd, 'Callback', {@callback_baselineEd, h, h.stopIndicatorLine});
set(h.baselineStartEd, 'Callback', {@callback_baselineEd, h, h.startIndicatorLine});

set(h.offsetEd, 'Callback', {@offsetEd_callback, h});
set(h.integralPb, 'Callback',{@callback_integralPb, h});
set(h.findPeakPb, 'Callback',{@callback_findPeakPb, h});
set(h.showPeakCb, 'Callback', {@callback_showPeakCb, h});
set(h.fitPb, 'Callback',{@callback_fitPb, h});
set(h.clearPb, 'Callback', {@callback_clearPb, h});
set(h.addStimuliPb, 'Callback', {@callback_addStimuliPb, h});
set(h.alignToStimuliCb, 'Callback', {@callback_alignToStimuliCb, h});


set(h.addSecondAnnotationPb, 'Callback', {@callback_addSecondAnnotationPb,h});
set(h.firstAnnotationHeightEd, 'Callback', {@callback_firstAnnotationHeightEd,h});
set(h.firstAnnotationYPosEd, 'Callback', {@callback_firstAnnotationYPosEd,h});
set(h.secondAnnotationHeightEd, 'Callback', {@callback_secondAnnotationHeightEd,h});
set(h.secondAnnotationYPosEd, 'Callback', {@callback_secondAnnotationYPosEd,h});
set(h.annotationFrameRateEd, 'Callback', {@callback_annotationFrameRateEd,h});

set(h.allXLimStartEd, 'Callback', {@allXLimEd_Callback, h});
set(h.allXLimStopEd, 'Callback', {@allXLimEd_Callback, h});
set(h.dffYLimStartEd, 'Callback', {@dffYLimEd_Callback, h});
set(h.dffYLimStopEd, 'Callback', {@dffYLimEd_Callback, h});
set(h.velYLimStartEd, 'Callback', {@velYLimEd_Callback, h});
set(h.velYLimStopEd, 'Callback', {@velYLimEd_Callback, h});

set(h.addVelocityPb, 'Callback', {@callback_addVelocityPb, h});
set(h.hideVelocityTb, 'Callback', {@callback_hideVelocityTb,h});


set(h.velocityThresholdEd, 'Callback', {@callback_fixVelocity,h});
set(h.useVelocityThresholdsCh, 'Callback', {@callback_fixVelocity,h});
set(h.smoothingWindowEd,'Callback', {@callback_fixVelocity,h});
set(h.smoothVelocityCh, 'Callback', {@callback_fixVelocity,h});



set(h.integrationFig, 'CloseRequestFcn', {@closeWindowFcn, h.f});
set(h.f, 'CloseRequestFcn', {@closeWindowFcn, h.integrationFig});

guidata(h.integrationFig,h);


%%

function callback_addSecondAnnotationPb(~, ~, h)
h=guidata(h.integrationFig);
if isempty(h.stimRects);
    %If no stimulus rectangles have been loaded, then load the second as if
    %they were stimulus rectangles
    callback_addStimuliPb(h.addStimuliPb, [], h);
    return
end

[matFileName, h.pathName]=uigetfile([h.pathName, '*.mat']);
if ~isnumeric(matFileName)
    if ~strcmpi(h.pathName(end),filesep)
        h.pathName=[h.pathName,filesep];
    end
    matFullPath=[h.pathName, matFileName];
else
    return
end

h.stimTimeOffset=str2double(get(h.offsetEd,'String'));
if isempty(h.stimTimeOffset) || isnan(h.stimTimeOffset)
    h.stimTimeOffset=0;
end

d=load(matFullPath, '-mat');

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


%%Remove old rectangles for second annotation
try
    delete(h.secondAnnotationRects);
    h.secondAnnotationRects=[];
catch ME
    %disp(ME.message);
end

%%/Remove old rectangles for second annotation

if ismember(h.stimStartCharacter, bChars)
    photoCharNum=find(strcmp(h.stimStartCharacter,bChars));
    if ismember(photoCharNum,bNumList)
        startCharInd=find(bNumList==photoCharNum,1,'first');
        timeList=timeList-timeList(startCharInd);
    end
    %h.stimEndChar=eChars(photoCharNum);
    keepCharInds=bNumList~=photoCharNum;
end

timeList=timeList(keepCharInds);


bNumList=bNumList(keepCharInds);
bActionList=bActionList(keepCharInds);

startTimeInds=strcmpi('start',bActionList);
stopTimeInds=strcmpi('stop',bActionList);

startTimes=timeList(startTimeInds);
stopTimes=timeList(stopTimeInds);
startNums=bNumList(startTimeInds);
stopNums=bNumList(stopTimeInds);
nRects=sum(startTimeInds);
    
nColumns=size(h.b,2);
% nColorsListed=length(h.colorList);

nColsRemoved=double(logical(h.timeColumn))+double(logical(h.controlColumn));
if length(h.legendLabels) < nColumns-nColsRemoved
    h.legendLabels=vertcat(h.legendLabels,repmat({' '}, [nColumns-nColsRemoved-length(h.legendLabels), 1]));
end

nRectColors=length(h.rectFaceColors2);
numsUsedSet=unique(bNumList);
nRectTypesUsed=length(numsUsedSet);

while nRectTypesUsed > nRectColors
    h.rectFaceColors2=vertcat(h.rectFaceColors2,{rand([1,3])});
    nRectColors=nRectColors+1;
end

h.rectWidths2=stopTimes-startTimes;
rectMaxY=max(max(h.dff));


% rectMinY=min([min(min(h.dff)),0]);


h.rectMinY2=str2double(get(h.secondAnnotationYPosEd,'String'));
h.rectHeight2=str2double(get(h.secondAnnotationHeightEd, 'String'));




if isBadNumericValue(h.rectMinY2)%isempty(h.rectMinY2) || ~isnumeric(h.rectMinY2) || isnan(h.rectMinY2)
    h.rectMinY2=h.rectMinY+h.rectHeight;
    %h.rectMinY2=min([min(min(h.dff)),0]);
    set(h.secondAnnotationYPosEd,'String',num2str(h.rectMinY2));
end

if isempty(h.rectHeight2) || ~isnumeric(h.rectHeight2) || isnan(h.rectHeight2)
    h.rectHeight2=(rectMaxY-h.rectMinY2);
    set(h.secondAnnotationHeightEd, 'String', num2str(h.rectHeight2));
end


% h.secondAnnotationHeightEd
% try
%     delete(get(hObject,'UserData'));
%     set(hObject,'UserData',[]);
% catch ME
% end


h.rectPosList2=zeros(nRects,4);
for i=1:nRects
    j=startNums(i);
%     %%%
%     disp([startTimes(i),h.rectMinY2,h.rectWidths2(i),h.rectHeight2])
%     %%%
    h.rectPosList2(i,:)=[startTimes(i),h.rectMinY2,h.rectWidths2(i),h.rectHeight2];
    h.r2(i)=rectangle('Position', h.rectPosList2(i,:),...
        'FaceColor', h.rectFaceColors2{j},...
        'Parent', h.ax);
end



doAlignToStims=get(h.alignToStimuliCb, 'Value');

if doAlignToStims
    xShift=-startTimes(1);    
    h.rectPosList2(:,1)=h.rectPosList2(:,1)+xShift;
    for i=1:nRects
        set(h.r2(i),'Position',h.rectPosList2(i,:));
    end
%     h=shiftAxes(h, xShift);
%     h=shiftPlottedLines(h, xShift);
%     h=shiftStimRectangles(h, xShift);
    %h.times=h.times-xShift;
% else
%     h=shiftStimRectangles(h, xShift);
%     h.times=h.trueTimes;
end



h.secondAnnotationLoaded=true;  
guidata(h.integrationFig,h);
 %%

function isBad=isBadNumericValue(n)
isBad=isempty(n) || ~isnumeric(n) || isnan(n);


function callback_firstAnnotationYPosEd(~,~,h)
h=guidata(h.integrationFig);

if isempty(h.stimRects) || ~any(ishandle(h.stimRects))
    return    
end

h.oldRectYPos=mode(h.rectPosList(:,2));

newYPos=str2double(get(h.firstAnnotationYPosEd, 'String'));

if isnumeric(newYPos) && ~isempty(newYPos) && ~isnan(newYPos)
    h.oldRectYPos=newYPos;
else
    newYPos=h.oldRectYPos;
end

set(h.firstAnnotationYPosEd, 'String', num2str(newYPos));
[h.rectPosList(:,2)]=deal(newYPos);

for i=1:length(h.stimRects)
    set(h.stimRects(i),'Position',h.rectPosList(i,:));
end


h.rectMinY=newYPos;


guidata(h.integrationFig,h);



function callback_showPeakCb(~, ~, h)
h=guidata(h.integrationFig);

h.peakVisible=get(h.showPeakCb,'Value');

if h.peakVisible
    if h.peakFound
        set(h.peakPoint,'Visible','on');
    end
    h.peakVisible=true;
else
    set(h.peakPoint,'Visible','off');
    h.peakVisible=false;
end


guidata(h.integrationFig, h);



function callback_findPeakPb(~, ~, h)
h=guidata(h.integrationFig);

startFrame=max([1,max(find(h.times<get(h.integrationStartEd,'UserData')))]);
stopFrame=min([length(h.times),min(find(h.times>get(h.integrationStopEd,'UserData')))]);
if stopFrame<startFrame
    return
end
    
    frameList=startFrame:stopFrame;
    dffList=h.dff(frameList);
    timeList=h.times(frameList);

[xPeak, yPeak]=findPeak(timeList, dffList);
set(h.peakPoint,'XData',xPeak);
set(h.peakPoint,'YData',yPeak);
h.peakFound=true;

set(h.findPeakTx, 'String', ['[',num2str(xPeak),',',num2str(yPeak),']']);

if h.peakVisible
    set(h.peakPoint,'Visible','on');
end

guidata(h.integrationFig, h);



function [xPeak, yPeak]=findPeak(xx, yy)


nPoints=length(xx);
if nPoints ~= length(yy)
    error('X and Y values must be of the same length.');
end



grossPeakVal=max(yy);
grossPeakInd=find(yy==grossPeakVal);

nGrossPeaks=length(grossPeakInd);
if nGrossPeaks>1
    localMeans=grossPeakVal*ones(nGrossPeaks,1);
    localMeanMaxInds=true(nGrossPeaks,1);
    curRadius=1;
    
    while sum(localMeanMaxInds)>1
        for i=1:nGrossPeaks
            startInd=max([1,grossPeakInd(i)-curRadius]);
            stopInd=min([nGrossPeaks,grossPeakInd(i)+curRadius]);
            localMeans(i)=mean(yy(startInd:stopInd));
        end
        if startInd==1 && stopInd==nGrossPeaks
            firstIndVal=find(localMeanMaxInds,1,'first');
            localMeanMaxInds=false(nGrossPeaks,1);
            localMeanMaxInds(firstIndVal)=true;
            break
        end
        curRadius=curRadius+1;
        localMeanMaxInds=localMeans==max(localMeans);
    end
    xPeak=xx(localMeanMaxInds);
    yPeak=yy(localMeanMaxInds);
else
    xPeak=xx(grossPeakInd);
    yPeak=grossPeakVal;
end


    








function callback_firstAnnotationHeightEd(~,~,h)
h=guidata(h.integrationFig);

if isempty(h.stimRects) || ~any(ishandle(h.stimRects))
    return
else
    h.stimRects=h.stimRects(ishandle(h.stimRects));
end

h.oldRectHeight=mode(h.rectPosList(:,4));
newRectHeight=str2double(get(h.firstAnnotationHeightEd, 'String'));
if isnumeric(newRectHeight) && ~isnan(newRectHeight)...
        && ~isempty(newRectHeight) && newRectHeight > 0
    h.oldRectHeight=newRectHeight;
else
    newRectHeight=h.oldRectHeight;
end

set(h.firstAnnotationHeightEd, 'String', num2str(newRectHeight));

[h.rectPosList(:,4)]=deal(newRectHeight);

for i=1:length(h.stimRects)
    set(h.stimRects(i),'Position',h.rectPosList(i,:));
end

h.rectHeight=newRectHeight;

guidata(h.integrationFig,h);


function callback_secondAnnotationYPosEd(~,~,h)
h=guidata(h.integrationFig);

if ~isfield(h,'r2') || isempty(h.r2) || ~any(ishandle(h.r2(1)))
    return
else
    h.stimRects=h.stimRects(ishandle(h.stimRects));
end

h.oldRectYPos2=mode(h.rectPosList2(:,2));

newYPos=str2double(get(h.secondAnnotationYPosEd,'String'));

if isnumeric(newYPos) && ~isempty(newYPos) && ~isnan(newYPos)
    h.oldRectYPos2=newYPos;
else
    newYPos=h.oldRectYPos2;
end

set(h.secondAnnotationYPosEd, 'String', num2str(newYPos));
[h.rectPosList2(:,2)]=deal(newYPos);

for i=1:length(h.r2)
    set(h.r2(i), 'Position', h.rectPosList2(i,:));
end
h.rectMinY2=newYPos;

guidata(h.integrationFig,h);


function callback_annotationFrameRateEd(~,~,h)
h=guidata(h.integrationFig);
    h.frameRate=str2double(get(h.annotationFrameRateEd,'String'));
    if ~isnumeric(h.frameRate) || isempty(h.frameRate) ||...
            ~all(size(h.frameRate)==[1,1]) || h.frameRate<=0
        h.frameRate=h.defaultFrameRate;
    end
    set(h.annotationFrameRateEd,'String',num2str(h.frameRate));
guidata(h.integrationFig,h);


function callback_secondAnnotationHeightEd(~,~,h)
h=guidata(h.integrationFig);

if ~isfield(h,'r2') || isempty(h.r2) || ~any(ishandle(h.r2(:)))
    return
else
    h.r2=h.r2(ishandle(h.r2));
end

h.oldRectHeight2=mode(h.rectPosList2(:,4));
newRectHeight=str2double(get(h.secondAnnotationHeightEd,'String'));

if isnumeric(newRectHeight) && ~isempty(newRectHeight)...
    && ~isnan(newRectHeight) && newRectHeight > 0
    h.oldRectHeight2=newRectHeight;
else
    newRectHeight=h.oldRectHeight2;
end

set(h.secondAnnotationHeightEd,'String',num2str(newRectHeight));

[h.rectPosList2(:,4)]=deal(newRectHeight);

for i=1:length(h.r2)
    set(h.r2(i),'Position',h.rectPosList2(i,:));
end

h.rectHeight2=newRectHeight;

guidata(h.integrationFig,h);




function baselineOrIntegration_callback(~, ~, h)
h=guidata(h.integrationFig);

if get(h.baselineRb,'Value')
        
    curBaselineStartVal=str2double(get(h.baselineStartEd, 'String'));
    curBaselineStopVal=str2double(get(h.baselineStopEd, 'String'));
    set(h.startIndicatorLine,'XData',[curBaselineStartVal,curBaselineStartVal]);
    set(h.stopIndicatorLine,'XData',[curBaselineStopVal,curBaselineStopVal]);
elseif get(h.integrationRb, 'Value')
    curIntegrationStartVal=str2double(get(h.integrationStartEd, 'String'));
    curIntegrationStopVal=str2double(get(h.integrationStopEd, 'String'));
    set(h.startIndicatorLine,'XData',[curIntegrationStartVal,curIntegrationStartVal]);
    set(h.stopIndicatorLine,'XData',[curIntegrationStopVal,curIntegrationStopVal]);
end

guidata(h.integrationFig,h);


function callback_fixVelocity(hObject, eventdata, h)
h=guidata(h.integrationFig);


h.doUseVelocityThreshold=get(h.useVelocityThresholdsCh,'Value');
h.doSmoothVelocity=get(h.smoothVelocityCh,'Value');


%Check velocity threshold for sensible entry:

newVelocityThreshold=abs(str2double(get(h.velocityThresholdEd,'String')));

if ~isnumeric(newVelocityThreshold) || isnan(newVelocityThreshold)
    newVelocityThreshold=h.velocityThreshold;
else
    h.velocityThreshold=newVelocityThreshold;
end


set(h.velocityThresholdEd, 'String', num2str(h.velocityThreshold));


%Check smoothing window for sensible entry:
newSmoothingWindow=abs(round(str2double(get(h.smoothingWindowEd,'String'))));
windowIsOdd=isOdd(newSmoothingWindow);
if ~exist('windowIsOdd','var') || isempty(windowIsOdd) || isnan(windowIsOdd)
    newSmoothingWindow=h.smoothingWindow;
elseif ~windowIsOdd
    newSmoothingWindow=newSmoothingWindow+1;
end
h.smoothingWindow=newSmoothingWindow;
set(h.smoothingWindowEd, 'String', num2str(h.smoothingWindow));

h.velocityVals=h.velocityRawVals;

if h.doUseVelocityThreshold
    %Cut out the velocities above the threshold
    h.velocityVals(h.velocityVals > h.velocityThreshold)=NaN;
end

if h.doSmoothVelocity
    %Smooth the velocity using the appropriate window
    h.velocityVals=smooth(h.velocityVals,h.smoothingWindow);
end

if isfield(h,'velCurve') && ishandle(h.velCurve)
    set(h.velCurve,'YData',h.velocityVals);
end

guidata(h.integrationFig,h);


function oddOrNot=isOdd(n)

if ~isnumeric(n) || n ~= round(n)
    oddOrNot=NaN;
else
    oddOrNot = floor(n/2)~=ceil(n/2);
end




function callback_hideVelocityTb(~, ~, h)
h=guidata(h.integrationFig);

if ishandle(h.velCurve)
    if get(h.hideVelocityTb, 'Value')
        set(h.velCurve,'Visible','off');
    else
        set(h.velCurve,'Visible','on');
    end
end

guidata(h.integrationFig,h);





function callback_addVelocityPb(~, ~, h)
h=guidata(h.integrationFig);
% filterSpec=[strcat(h.pathName,{'*.xls*';'*.csv'}),{'Microsoft Excel';'Comma Separated Value'}];
curDir=pwd;
try
    cd(h.pathName);
catch ME
    disp(ME.message);
end
[h.velocityFileName, h.pathName]=uigetfile({'*.xls*','Microsoft Excel';'*.csv','Comma Separated Value'});
cd(curDir);


% [h.velocityFileName,h.pathName]=uigetfile([h.pathName,{'*.xls','*.xls','*.csv'}]);
if ~isnumeric(h.velocityFileName)
    if ~strcmpi(h.pathName(end),filesep)
        h.pathName=[h.pathName,filesep];
    end
    velocityFullPath=[h.pathName, h.velocityFileName];
else
    return
end

h.xlsData=xlsread(velocityFullPath);


h.velTimes=h.xlsData(:,1);
h.velocityRawVals=h.xlsData(:,2);
h.velocityVals=h.velocityRawVals;

h.velAx=axes('Parent', get(h.ax,'Parent'),...
    'XAxisLocation', 'bottom',...
    'YAxisLocation', 'right',...
    'Color', 'none',...
    'XLim', get(h.ax,'XLim'),...
    'YColor', h.velColor);%For velocity graph's axes

h.velCurve=plot(h.velTimes,h.velocityVals,'Parent',h.velAx,'Color',h.velColor); %For velocity curve object

if get(h.hideVelocityTb,'Value')
    set(h.velCurve,'Visible','off');
end

h.velZeroLine=line('Parent', h.velAx, 'XData',[min(h.velTimes),max(h.velTimes)],'YData',[0,0],'LineStyle','--' ,'Color', [0,0,0]);


set(h.velAx,...
    'XAxisLocation', 'bottom',...
    'YAxisLocation', 'right',...
    'Color', 'none',...
    'XLim', get(h.ax,'XLim'),...
    'XLimMode', 'manual',...
    'YLimMode', 'manual',...
    'YColor', h.velColor);%For velocity graph's axes
set(h.ax,'YColor', h.dffColor);

h=setLineLengths(h);

guidata(h.integrationFig,h);


function h=setLineLengths(h)

xLim=get(h.ax,'XLim');
yLim=get(h.ax,'YLim');

horizontalLineFields={'zeroLine', 'velZeroLine'};
verticalLineFields={'startIndicatorLine','stopIndicatorLine'};

h=setLineToLimits(h, 'X', horizontalLineFields);
h=setLineToLimits(h, 'Y', verticalLineFields);
% h=setRectangleToLimits(h, 'Y', 'stimRects');

function h=setLineToLimits(h, axisLetter, lineFieldNames)
axisLetter=upper(axisLetter);
for i=1:length(lineFieldNames)
    curLineName=lineFieldNames{i};
    if isfield(h,curLineName) && ishandle(h.(curLineName))
        curLine=h.(curLineName);
        curAx=get(curLine,'Parent');
        curLim=get(curAx,[axisLetter,'Lim']);
        set(curLine, [axisLetter,'Data'], curLim);
    end
end

function h=setRectangleToLimits(h, axisLetter, rectFieldName)
axisLetter=upper(axisLetter);
if isfield(h,rectFieldName)
    stimRects=h.(rectFieldName);
    for i=1:length(stimRects)
        curRectangle=stimRects(i);
        if ishandle(curRectangle)
            curAx=get(curRectangle,'Parent');
            curPos=get(curRectangle,'Position');
            curLim=get(curAx,[axisLetter,'Lim']);
            if strcmpi(axisLetter,'X')
                set(curRectangle, 'Position',...
                    [curLim(1),curPos(2),curLim(2)-curLim(1),curPos(4)]);                
            elseif strcmpi(axisLetter,'Y')
                set(curRectangle, 'Position',...
                    [curPos(1),curLim(1),curPos(3),curLim(2)-curLim(1)]);
            end
        end
    end
end


function allXLimEd_Callback(hObject, eventdata, h)
h=guidata(h.integrationFig);

oldX=get(h.ax,'XLim');
oldXStart=oldX(1);
oldXStop=oldX(2);

newXStart=str2double(get(h.allXLimStartEd,'String'));
newXStop=str2double(get(h.allXLimStopEd,'String'));

if newXStart >= newXStop
    xDiff=oldXStop-oldXStart;
    if hObject==h.allXLimStartEd
        newXStop=newXStart+xDiff;
    elseif hObject==h.allXLimStopEd
        newXStart=newXStop-xDiff;
    end
end
set(h.allXLimStartEd,'String',num2str(newXStart));
set(h.allXLimStopEd,'String',num2str(newXStop));
set(h.ax,'XLim',[newXStart,newXStop], 'XLimMode', 'manual');
set(h.velAx,'XLim', [newXStart,newXStop], 'XLimMode', 'manual');

h=setLineLengths(h);
guidata(h.integrationFig,h);

function dffYLimEd_Callback(hObject, eventdata, h)
h=guidata(h.integrationFig);

oldY=get(h.ax,'YLim');
oldYStart=oldY(1);
oldYStop=oldY(2);
newYStart=str2double(get(h.dffYLimStartEd,'String'));
newYStop=str2double(get(h.dffYLimStopEd,'String'));

if newYStart >= newYStop
    yDiff=oldYStop-oldYStart;
    if hObject==h.dffYLimStartEd
        newYStop=newYStart+yDiff;
    elseif hObject==h.dffYLimStopEd
        newYStart=newYStop-yDiff;
    end
end
set(h.dffYLimStartEd,'String',num2str(newYStart));
set(h.dffYLimStopEd,'String',num2str(newYStop));
set(h.ax,'YLim', [newYStart,newYStop], 'YLimMode', 'manual');

h=setLineLengths(h);
guidata(h.integrationFig,h);

function velYLimEd_Callback(hObject, eventdata, h)
h=guidata(h.integrationFig);

oldY=get(h.velAx,'YLim');
oldYStart=oldY(1);
oldYStop=oldY(2);
newYStart=str2double(get(h.velYLimStartEd,'String'));
newYStop=str2double(get(h.velYLimStopEd,'String'));

if newYStart >= newYStop
    yDiff=oldYStop-oldYStart;
    if hObject==h.velYLimStartEd
        newYStop=newYStart+yDiff;
    elseif hObject==h.velYLimStopEd
        newYStart=newYStop-yDiff;
    end
end
set(h.velYLimStartEd,'String',num2str(newYStart));
set(h.velYLimStopEd,'String',num2str(newYStop));
set(h.velAx,'YLim', [newYStart,newYStop], 'YLimMode', 'manual');

h=setLineLengths(h);
guidata(h.integrationFig,h);



function callback_alignToStimuliCb(hObject, eventdata, h)
h=guidata(h.integrationFig);

newVal=get(hObject, 'Value');
% rects=get(h.addStimuliPb,'UserData');
nRects=length(h.stimRects);
if nRects==0
    guidata(h.integrationFig,h);
    return;
end
%rectPos=get(h.stimRects(1),'Position');

if newVal
    xShift=-h.offsetRectXPos(1);
else
    xShift=h.offsetRectXPos(1);
end


h=shiftStimRectangles(h, xShift);
h=shiftPlottedLines(h, xShift);
h=shiftAxes(h, xShift);




% h=shiftObjects(h,xShift);
guidata(h.integrationFig,h);




function h=shiftObjects(h,xShift)
curXLim=get(h.ax,'XLim');
set(h.dffPlot,'XData',h.times);
set(h.ax,'XLim',curXLim-xShift);

 
[h,curBaselineStartVal]=checkForNumericEntry(h,h.baselineStartEd,xShift);
[h,curBaselineStopVal]=checkForNumericEntry(h,h.baselineStopEd, xShift);
[h,curIntegrationStartVal]=checkForNumericEntry(h,h.integrationStartEd,xShift);
[h,curIntegrationStopVal]=checkForNumericEntry(h,h.integrationStopEd,xShift);

if get(h.baselineRb,'Value')
    set(h.startIndicatorLine,'XData',[curBaselineStartVal,curBaselineStartVal]);
    set(h.stopIndicatorLine,'XData',[curBaselineStopVal,curBaselineStopVal]);
elseif get(h.integrationRb, 'Value')
    set(h.startIndicatorLine,'XData',[curIntegrationStartVal,curIntegrationStartVal]);
    set(h.stopIndicatorLine,'XData',[curIntegrationStopVal,curIntegrationStopVal]);
end
h.zeroLineXData=h.zeroLineXData-xShift;


set(h.zeroLine,'XData',[min(h.times),max(h.times)]);

% set(h.startIndicatorLine, 'Position', h.startIndicatorLinePos);
% set(h.stopIndicatorLine, 'Position', h.stopIndicatorLinePos);


function [h,v]=checkForNumericEntry(h,hEdit, xShift)
v=str2double(get(hEdit,'String'));
if isempty(v) || isnan(v) || length(v)>1
    v=0;
else
    v=v-xShift;
end
set(hEdit,'String',num2str(v));



function closeWindowFcn(hObject, eventData, otherWindow)
if ishandle(otherWindow)
    try
        delete(h.baselineFig)
    catch ME
        %disp(ME.message);
    end
    
    try
        delete(otherWindow);
    catch ME
        disp(ME.message);
    end
    
    try
        delete(hObject);
    catch ME
        disp(ME.message);
    end
    

end

% yLim=get(h.ax, 'YLim');
% set(h.ax,'YLimMode', 'manual');
% if get(defaultYLimPu,'Value')==1
%     yLim(1)=0;
% elseif get(defaultYLimPu,'Value')==2
%     yLim(1)=min(h.dff);
% elseif get(defaultYLimPu,'Value')==3
%     yLim(1)=min(h.dff(floor(length(h.times)/4):ceil(3*length(h.times)/8)));
% end



function unclickedMotion(hObject, eventdata, h)
clickPoint=get(h.ax,'CurrentPoint');
x=clickPoint(1,1);
y=clickPoint(1,2);

xStart=get(h.indicatorLines(1),'XData'); xStart=xStart(1);
xStop=get(h.indicatorLines(2),'XData'); xStop=xStop(1);

if x <= xStart+1 && x >= xStart-1
    isOnStart=true;
else
    isOnStart=false;
end

if x <= xStop+1 && x >= xStop-1
    isOnStop=true;
else
    isOnStop=false;
end

if isOnStart || isOnStop
    set(hObject,'Pointer','crosshair');
else
    set(hObject,'Pointer','arrow');
end

guidata(h.integrationFig,h);





function clickedMotion(hObject, eventdata, h)
guidata(h.integrationFig,h);
h=guidata(h.integrationFig);

function buttonDown(hObject, eventdata, h)
guidata(h.integrationFig,h);
h=guidata(h.integrationFig);

function buttonUp(hObject, eventdata, h)
h=guidata(h.integrationFig);
set(hObject, 'WindowButtonMotionFcn', {@unclickedMotion, h.ax, h.indicatorLines});
guidata(h.integrationFig,h);


function callback_clearPb(hObject, eventdata, h)
h=guidata(h.integrationFig);

if ~get(h.preserveStimuliCb, 'Value')
    %h=clearUserDataObjects(h,h.addStimuliPb);
    if get(h.alignToStimuliCb,'Value')
        xShift=h.trueTimes(1)-h.times(1);
        h=shiftAxes(h, xShift);
        h=shiftPlottedLines(h, xShift);
        h=shiftStimRectangles(h, xShift);
    end
    
    try
        delete(h.stimRects);
        h.stimRects=[];
    catch ME
        disp(ME.message);
    end
    
    try
        delete(h.r2);
        h.r2=[];
    catch ME
        disp(ME.message);
    end
    %h.times=h.trueTimes;
    %set(h.dffPlot,'XData',h.times);
    
end

if ~get(h.preserveVelocityCb, 'Value')
    
    try
        delete(h.velCurve);
        delete(h.velZeroLine);
        delete(h.velAx);
    catch ME
        disp(ME.message);
    end
    
end

try
    delete(h.integrationRects);
    h.integrationRects=[];
catch ME
    disp(ME.message);
end
try
    delete(h.fitLine);
    h.fitLine=[];
catch ME
    disp(ME.message);
end


%h=clearUserDataObjects(h,h.integralPb);
%h=clearUserDataObjects(h,h.fitPb);

guidata(h.integrationFig,h);


function h=clearUserDataObjects(h,obj)
try
    udo=get(obj,'UserData');
    if ~isempty(udo)
        delete(udo);
        set(obj,'UserData',[]);
    end
catch ME
    disp(ME.message); 
end


function h=shiftAxes(h, xShift)
%Shifts axes and marker lines to match other shifts
curXLim=get(h.ax,'XLim')+xShift;
set(h.ax,'XLim',curXLim);

if isfield(h, 'velAx') && ishandle(h.velAx)
    set(h.velAx,'XLim',get(h.ax,'XLim'));
end

set(h.allXLimStartEd,'String',num2str(curXLim(1)));
set(h.allXLimStopEd, 'String', num2str(curXLim(2)));

h=setLineLengths(h);


 
function [h,v]=shiftNumericEntry(h,hEdit, xShift)
%Shifts the values on edit boxes to comport with new marker positions
v=str2double(get(hEdit,'String'));
if isempty(v) || isnan(v) || length(v)>1
    v=0;
else
    v=v+xShift;
end
set(hEdit,'String',num2str(v),'UserData',v);


function h=shiftPlottedLines(h, xShift)
%Shifts plotted lines for DF/F, fits, and integral
h.times=h.times+xShift;
set(h.dffPlot,'XData',h.times);

[h,curBaselineStartVal]=shiftNumericEntry(h,h.baselineStartEd,xShift);
[h,curBaselineStopVal]=shiftNumericEntry(h,h.baselineStopEd, xShift);
[h,curIntegrationStartVal]=shiftNumericEntry(h,h.integrationStartEd,xShift);
[h,curIntegrationStopVal]=shiftNumericEntry(h,h.integrationStopEd,xShift);
if get(h.baselineRb,'Value')
    set(h.startIndicatorLine,'XData',[curBaselineStartVal,curBaselineStartVal]);
    set(h.stopIndicatorLine,'XData',[curBaselineStopVal,curBaselineStopVal]);
elseif get(h.integrationRb, 'Value')
    set(h.startIndicatorLine,'XData',[curIntegrationStartVal,curIntegrationStartVal]);
    set(h.stopIndicatorLine,'XData',[curIntegrationStopVal,curIntegrationStopVal]);
end
zeroLineXData=get(h.zeroLine,'XData');
%zeroLineXData=zeroLineXData+xShift; %Commented March 2018

set(h.zeroLine,'XData',zeroLineXData);
%Shift integration rectangles and fit lines accordingly
if ~isempty(h.fitLine) && ishandle(h.fitLine)
    fitLineXData=get(h.fitLine,'XData');
    fitLineXData=fitLineXData+xShift;
    set(h.fitLine,'XData',fitLineXData);
end

if isfield(h,'velCurve') && ishandle(h.velCurve)
    h.velTimes=h.velTimes+xShift;
    set(h.velCurve,'XData',h.velTimes);
end


if ~isempty(h.integrationRects) && ishandle(h.integrationRects(1))
    
    nIntRects=length(h.integrationRects);
    for i=1:nIntRects
        integrationRectPos=get(h.integrationRects(i),'Position');
        integrationRectPos=integrationRectPos+[xShift,0,0,0];
        set(h.integrationRects(i),'Position',integrationRectPos);
    end
%     integrationRectPos=get(h.integrationRects,'Position');
%     if iscell(integrationRectPos)
%         integrationRectPos=cellfun(h.shiftPosVec,integrationRectPos,...
%             repmat({xShift},size(integrationRectPos)),...
%             'UniformOutput', false);
%     else
%         integrationRectPos=integrationRectPos+[xShift,0,0,0];
%     end
%     set(h.integrationRects,'Position',integrationRectPos);
end





function h=shiftStimRectangles(h, xShift, fieldList)
%Shifts stimulus rectangles
nRects=length(h.stimRects);
if nRects > 0
    
%     disp(h.rectPosList(:,:))
    
    
    h.rectPosList(:,1)=h.rectPosList(:,1)+xShift;
    for i=1:nRects
%         disp(i);
%         disp(h.stimRects(i));
%         disp(h.rectPosList(i,:));
        
        set(h.stimRects(i),'Position',h.rectPosList(i,:));
        
        %rectPos=get(h.stimRects(i),'Position');
        %rectPos(1)=h.offsetRectXPos(i)+xShift;
        %set(h.stimRects(i),'Position',rectPos);
    end
end

if isfield(h,'r2') && ~isempty(h.r2)
    nRects2=length(h.r2);
else
    nRects2=0;
end
if nRects2 > 0
    h.rectPosList2(:,1)=h.rectPosList2(:,1)+xShift;
    for i=1:nRects2
        set(h.r2(i),'Position',h.rectPosList2(i,:));
    end
    
end





function offsetEd_callback(hObject, eventdata, h)
h=guidata(h.integrationFig);

rectObjects=get(h.addStimuliPb,'UserData');
offsetStr=get(hObject,'String');
offsetVal=str2double(offsetStr);
basePositions=get(hObject, 'UserData');

if isempty(h.stimRects)
    return
end

nStimRects=length(h.stimRects);

if isempty(basePositions)
    for i=1:nStimRects
        nStimRects=length(h.stimRects);
    end
    set(hObject,'UserData',rectPosList);
end

if isempty(offsetVal) || isnan(offsetVal)
    offsetVal=h.offsetVal;
    set(hObject,'String',num2str(offsetVal));
end

offsetChange=offsetVal-h.offsetVal;
%If we are not centering the graph to the stimuli, we move them to the
%right for positive values. If we are centering the graph to the stimuli,
%then we move everything else to the left for positive values.



h.offsetRectXPos=h.trueRectXPos+offsetVal;
if ~get(h.alignToStimuliCb,'Value') %Shift rectangles by offset
    xShift=offsetChange;
    h=shiftStimRectangles(h, xShift);
    %h=shiftObjects(h,offsetChange);
else %Shift graph opposite direction by offset
    xShift=-offsetChange;
    h=shiftAxes(h, xShift);
    h=shiftPlottedLines(h, xShift);
end

% for i=1:nStimRects
%     try
%         curRectPos=get(rectObjects(i),'Position');
%         if get(h.alignToStimuliCb,'Value')
%             curRectPos(1)=h.centeredRectXPos(i,:);
%             
%         else
%             curRectPos(1)=h.offsetRectXPos(i,:);
%             
%         end
%         set(rectObjects(i),'Position', curRectPos);
%     catch ME
%         disp(ME.message);
%     end
%     
% end

h.offsetVal=offsetVal;
guidata(h.integrationFig,h);


function callback_addStimuliPb(hObject, ~, h)
h=guidata(h.integrationFig);
if ~exist('h.pathName', 'var')
    h.pathName='.\';
end

[matFileName, h.pathName]=uigetfile([h.pathName, '*.mat; *.txt']);
if ~isnumeric(matFileName)
    if ~strcmpi(h.pathName(end),filesep)
        h.pathName=[h.pathName,filesep];
    end
    matFullPath=[h.pathName, matFileName];
else
    return
end


h.stimTimeOffset=str2double(get(h.offsetEd,'String'));
if isempty(h.stimTimeOffset) || isnan(h.stimTimeOffset)
    h.stimTimeOffset=0;
end

[~, ~, fileExt]=fileparts(matFullPath);
if strcmpi(fileExt, '.txt')
    %s=textread(matFullPath, '%s');
    fid=fopen(matFullPath);
    s=textscan(fid,'%s');
    fclose(fid);
    s=s{1};
    
    
    startInd0=find(~cellfun(@isempty,regexp(s,'file:')),1,'first')+1;
    stopInd0=find(~cellfun(@isempty,regexp(s,'S1:')),1,'first')-1;
    s0=s(startInd0:stopInd0);
    s0=reshape(s0,[2,length(s0)/2])';
    
    bChars=s0(:,2);
    behaviorLabels=s0(:,1);
    %eChars=repmat({stopChar},size(bChars));
    
    startInd1=find(~cellfun(@isempty,regexp(s,'-------')),1,'first')+1;
    stopInd1=find(~cellfun(@isempty,regexp(s,'S2:')),1,'first')-1;
    s1=s(startInd1:stopInd1);
    s1=reshape(s1,[3,length(s1)/3])';
    startNums1=cellfun(@str2double, s1(:,1));
    stopNums1=cellfun(@str2double, s1(:,2))+1;
    s1Labels=s1(:,3);
    
    startInd2=find(~cellfun(@isempty,regexp(s,'-------')),1,'last')+1;
    stopInd2=length(s);
    s2=s(startInd2:stopInd2);
    s2=reshape(s2, [3,length(s2)/3])';
    startNums2=cellfun(@str2double, s2(:,1));
    stopNums2=cellfun(@str2double, s2(:,2))+1;
    s2Labels=s2(:,3);
    
    
    notOther=~strcmpi(s1Labels,'other');
    startNums1=startNums1(notOther);
    stopNums1=stopNums1(notOther);
    nRects=sum(notOther);
    rectLabels=s1Labels(notOther);
    
    framePeriod=1/h.frameRate;
    
    notOther2=~strcmpi(s2Labels,'other');
    startNums2=startNums2(notOther2);
    stopNums2=stopNums2(notOther2);
    
    %startNums=startNums1(notOther);
    startTimes=startNums1*framePeriod;
    stopTimes=stopNums1*framePeriod;
	
    labelsUsedSet=unique(s1Labels(notOther));
    nRectTypesUsed=length(labelsUsedSet);
    startNumList=1:nRectTypesUsed;
    
    startNums=ones(size(startTimes));
    
    startTimes2=startNums2*framePeriod;
    %stopTimes2=stopTimes2*framePeriod;
    
    if ~isempty(startNums2)
        zeroTime=startNums2(1);
        startTimes=startTimes-zeroTime;
        afterLightInds=startTimes>=0;
        startTimes=startTimes(afterLightInds);
        stopTimes=stopTimes(afterLightInds);
    else
        
    end
    
    
    for i=1:nRectTypesUsed
        curNumInds=strcmpi(rectLabels,labelsUsedSet{i});
        [startNums(curNumInds)]=deal(i);
    end
    
    
%     s1Labels
%     labelsUsedSet
elseif strcmpi(fileExt, '.mat')
    
    
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
    
    
    %%Remove old stimulus rectangles
    try
        if get(h.alignToStimuliCb,'Value')
            xShift=h.trueTimes(1)-h.times(1);
            h=shiftAxes(h, xShift);
            h=shiftPlottedLines(h, xShift);
            h=shiftStimRectangles(h, xShift);
        end
    catch ME
        disp(ME.message);
    end
    
    try
        delete(h.stimRects);
        h.stimRects=[];
    catch ME
        disp(ME.message);
    end
    %%/Remove old stimulus rectangles
    
    
    
    if ismember(h.stimStartCharacter, bChars)
        photoCharNum=find(strcmp(h.stimStartCharacter,bChars));
        if ismember(photoCharNum,bNumList)
            startCharInd=find(bNumList==photoCharNum,1,'first');
            timeList=timeList-timeList(startCharInd);
        end
        %h.stimEndChar=eChars(photoCharNum);
        keepCharInds=bNumList~=photoCharNum;
    end
    
    
    
    
    timeList=timeList(keepCharInds);
    %charList=charList(keepCharInds);
    bNumList=bNumList(keepCharInds);
    bActionList=bActionList(keepCharInds);
    
    
    nActions=length(bActionList);
    startTimeInds=strcmpi('start',bActionList);
    stopTimeInds=strcmpi('stop', bActionList);
    
    startTimes=timeList(startTimeInds);
    stopTimes=timeList(stopTimeInds);
    startNums=bNumList(startTimeInds);
    
    nRects=sum(startTimeInds);
    
    
    numsUsedSet=unique(bNumList);
    nRectTypesUsed=length(numsUsedSet);
    
    
end

nColumns=size(h.b,2);
nColorsListed=length(h.colorList);




%%%%


% [startTimes stopTimes]
% bNumList

nColsRemoved=double(logical(h.timeColumn))+double(logical(h.controlColumn));

if length(h.legendLabels) < nColumns-nColsRemoved
    h.legendLabels=vertcat(h.legendLabels,repmat({' '},[nColumns-nColsRemoved-length(h.legendLabels),1]));
end

nRectColors=length(h.rectFaceColors);

% numsUsedSet=unique(bNumList);
% nRectTypesUsed=length(numsUsedSet);


while nRectTypesUsed > nRectColors
    h.rectFaceColors=vertcat(h.rectFaceColors,{rand([1,3])});
    nRectColors=nRectColors+1;
end

rectMinY=min([min(min(h.dff)),0]);
rectMaxY=max(max(h.dff));
rectHeight=rectMaxY-rectMinY;
rectWidths=stopTimes-startTimes;


h.rectMinY=str2double(get(h.firstAnnotationYPosEd, 'String'));
h.rectHeight=str2double(get(h.firstAnnotationHeightEd, 'String'));

if isBadNumericValue(h.rectMinY)
    h.rectMinY=min([min(min(h.dff)),0]);
    set(h.firstAnnotationYPosEd, 'String', num2str(h.rectMinY));
end

if isBadNumericValue(h.rectHeight)
    h.rectHeight=(rectMaxY-h.rectMinY)/2;
    set(h.firstAnnotationHeightEd, 'String', num2str(h.rectHeight));
end

if isBadNumericValue(h.rectMinY)%isempty(h.rectMinY) || ~isnumeric(h.rectMinY) || isnan(h.rectMinY)
    h.rectMinY=min([min(min(h.dff)),0]);
    set(h.firstAnnotationYPosEd, 'String', num2str(h.rectMinY));
end
% set(h.rect(i), 'EdgeColor', h.r(i).col);

if isBadNumericValue(h.rectHeight)
    
end






h.rectWidths=rectWidths;

h.trueRectXPos=startTimes;
startTimes=startTimes-h.stimTimeOffset;
h.offsetRectXPos=startTimes;



try
    delete(get(hObject,'UserData'));
    set(hObject,'UserData',[]);
catch ME
    %disp(ME.message);
end

h.rectPosList=zeros(nRects,4);

doAlignToStims=get(h.alignToStimuliCb, 'Value');




% curXLim=get(h.ax,'XLim');
% set(h.dffPlot,'XData',h.times);
% set(h.ax,'XLim',curXLim-xShift);
% h=shiftObjects(h,xShift);


% try
% nRects
% size(h.stimRects)
% size(h.rectPosList)
% size(h.rectFaceColors)
% size(startNums)
% startNums
% h.rectFaceColors{:}

% size(startNums1)
% nRects

% [startNums1, stopNums1]

h.stimRects=zeros(nRects,1);
for i=1:nRects
    j=startNums(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h.rectPosList(i,:)=[startTimes(i),h.rectMinY,rectWidths(i),h.rectHeight];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         [startTimes(i), stopTimes(i)]
    %         h.rectPosList(i,:)
    
    h.stimRects(i)=rectangle('Position', h.rectPosList(i,:),...
        'FaceColor', h.rectFaceColors{j},...
        'Parent', h.ax);
    %h.offsetRectXPos(i)=rectPosList(i,1);
    %         if doAlignToStims
    %             rectPosList(i,1)=rectPosList(i,1)-xShift;
    %             set(rectSet(i),'Position',rectPosList(i,:));
    %         end
    
    uistack(h.stimRects(i),'bottom');
end
    %h.stimRects=rectSet;
    %set(hObject,'UserData',rectSet);
% catch ME
%     disp(ME.message);
% end



try
    set(h.offsetEd, 'UserData',h.rectPosList);
catch ME
    disp(ME.message);
end

if doAlignToStims
    xShift=-startTimes(1);
    h=shiftAxes(h, xShift);
    h=shiftPlottedLines(h, xShift);
    h=shiftStimRectangles(h, xShift);
    %h.times=h.times-xShift;
% else
%     h=shiftStimRectangles(h, xShift);
%     h.times=h.trueTimes;
end

guidata(h.integrationFig,h);

function ddff=differentiateData(dff,times)
ddff=dff-circshift(dff,1); ddff=ddff(2:end);
dt=times-circshift(times,1); dt=dt(2:end);
ddff=ddff./dt;


function [b,c]=fitDerivative(h)
ddff=differentiateData(h.dff,h.times);
t=h.times-h.times(1);
t=t(2:end); 
dddff=differentiateData(ddff,t);

length(ddff)

bc0=max(ddff(1:min([5,length(ddff)])));
bbc0=max(dddff(1:min([5,length(dddff)])));
b0=bbc0/bc0;
c0=-abs(bc0/b0);
beta0=[bc0, c0];


weightFactors=ones(size(ddff));
dexpc=@(t, bc, c) bc*exp(c*t);
residualSum=@(beta) dot(weightFactors,abs(dexpc(t, beta(1), beta(2))-ddff));
options = optimset('MaxFunEvals',10000);
beta=fminsearch(residualSum, beta0, options);
fitVals=dexpc(t,beta(1),beta(2));

bc=fitVals(1);
c=fitVals(2);
b=bc/c;
guidata(h.integrationFig,h);



function callback_fitPb(~, ~, h)
h=guidata(h.integrationFig);

baseOn=str2double(get(h.baselineStartEd,'String'));
baseOff=str2double(get(h.baselineStopEd,'String'));
baseLength=baseOff-baseOn;
time=h.times;



stimOn=str2double(get(h.integrationStartEd,'String'));
fitMax=str2double(get(h.integrationStopEd,'String'));

if isfield(h,'stimRects') && ~isempty(h.stimRects)
    stimRectPos=get(h.stimRects(1),'Position');
    stimOff=stimRectPos(3);
else
    stimOff=fitMax;
end


[p,h.fitVals] = getLinFit(time,baseOn,baseOff,stimOn,stimOff,fitMax,h.dff);


h.tau_r=p(1)^2;%abs(p(1));%(p(1)^2;
h.tau_f=(p(1)^2 + p(3)^2);%abs(p(1))+abs(p(3));%(p(1)^2 + p(3)^2);
h.fitLine=line('Parent', h.ax, 'XData', h.times, 'YData' , h.fitVals, 'Color',...
    h.fitColor, 'LineWidth', 2);

set(h.fitTxA,'String',['T_rise=',num2str(h.tau_r)]);
set(h.fitTxB,'String',['T_fall=',num2str(h.tau_f)]);

guidata(h.integrationFig,h);

function callback_fitExponentialPb(hObject, eventdata, h)
h=guidata(h.integrationFig);

startFrame=max([1,max(find(h.times<get(h.integrationStartEd,'UserData')))]);
stopFrame=min([length(h.times),min(find(h.times>get(h.integrationStopEd,'UserData')))]);
if stopFrame<startFrame
    return
end
fitInds=startFrame:stopFrame;

% goodIndsLogical=h.dff(fitInds)>0;
% fitInds=fitInds(goodIndsLogical);
beta0_a=min(h.dff(fitInds))-abs(0.1*min(h.dff(fitInds)));
beta0_b=h.dff(fitInds(1))-h.dff(fitInds(end));
beta0_c=-beta0_b/(h.times(fitInds(end))-h.times(fitInds(1)));


beta0=[beta0_a, beta0_b, beta0_c];
modelfun=@(bb,x)(bb(1)+bb(2)*exp(-bb(3)*x));

fitTimes=h.times(fitInds);
fitTimesAdjusted=fitTimes-fitTimes(1);
fitDff=h.dff(fitInds);
beta=nlinfit(fitTimesAdjusted,fitDff,modelfun,beta0);
fitVals=modelfun(beta,fitTimesAdjusted);


% p=polyfit(fitTimesAdjusted, log(fitDff), 1);
% beta(1)=p(1);
% beta(2)=exp(p(2));
% fitVals=beta(2)+exp(beta(1)*fitTimesAdjusted);

% global x y;
% x=fitTimesAdjusted;
% y=fitVals;
% [fitDff, log(fitDff), fitVals]

% weightFactors=fitDff-beta0_a;
% expc=@(t, a, h.b, c) a+h.b*exp(c*t);
% beta0=[beta0_a, beta0_b, beta0_c];
% residualSum=@(beta) dot(weightFactors,abs(expc(fitTimes, beta(1), beta(2), beta(3))-fitDff));
% options = optimset('MaxFunEvals',1000);
% beta=fminsearch(residualSum, beta0, options);
% fitVals=expc(fitTimesAdjusted,beta(1),beta(2),beta(3));
chiSquarePearson=sum(((fitDff-fitVals).^2)./fitVals);

[h.b,c]=fitDerivative(h);


set(h.fitTxA, 'String', ['a = ', num2str(beta(1))]);
set(h.fitTxB, 'String', ['h.b = ', num2str(beta(2))]);
set(h.fitTxC, 'String', ['c = ', num2str(beta(3))]);

hold(h.ax,'on');



try
    delete(h.fitLine);
    h.fitLine=[];

    %delete(get(hObject,'UserData'));
    %set(hObject,'UserData',[]);
    
catch ME
    %disp(ME.message);
end
xLim=get(h.ax,'XLim');    
yLim=get(h.ax,'YLim');

h.fitLine=line('Parent', h.ax, 'XData',fitTimes, 'YData' ,fitVals, 'Color',...
    [0.8,0.4,0], 'LineWidth', 2);
% fitObj=line('Parent', h.ax, 'XData',fitTimes, 'YData' ,fitVals, 'Color',...
%     [0.8,0.4,0], 'LineWidth', 2);
% set(hObject,'UserData', fitObj);

set(h.ax,'XLim',xLim,'YLim',yLim);
h=setLineLengths(h);
guidata(h.integrationFig,h);
% lineHandle=plot(fitTimes, fitCurve,'Color',[1,0.5,0]);
% set(hObject,'UserData',lineHandle);


function callback_integrationEd(hObject, eventdata, h, indicatorLine)
h=guidata(h.integrationFig);
newStr=get(hObject, 'String');
newNum=str2double(newStr);
if isnan(newNum)
    set(hObject,'String',num2str(get(hObject,'UserData')));
    return
else
    set(hObject,'UserData',newNum);
    if get(h.integrationRb,'Value')
        set(indicatorLine, 'XData', [newNum, newNum]);
    end
end
guidata(h.integrationFig,h);



function callback_baselineEd(hObject, eventdata, h, indicatorLine)
h=guidata(h.integrationFig);
newStr=get(hObject, 'String');
newNum=str2double(newStr);
if isnan(newNum)
    set(hObject,'String',num2str(get(hObject,'UserData')));
    return
else
    set(hObject,'UserData',newNum);
    if get(h.baselineRb,'Value')
        set(indicatorLine, 'XData', [newNum, newNum]);
    end
end
guidata(h.integrationFig,h);




function showLinesCb_callback(hObject, eventData, h)
if get(hObject,'Value')
    set(h.indicatorLines, 'Visible', 'on');
else
    set(h.indicatorLines, 'Visible', 'off');
end
guidata(h.integrationFig,h);

function mouseBoundSet()
rect = getrect();
guidata(h.integrationFig,h);


function moveBoundCallback()
guidata(h.integrationFig,h);


function callback_integralPb(hObject, eventdata, h)
h=guidata(h.integrationFig);
startFrame=max([1,max(find(h.times<get(h.integrationStartEd,'UserData')))]);
stopFrame=min([length(h.times),min(find(h.times>get(h.integrationStopEd,'UserData')))]);
if stopFrame<startFrame
    return
else
    
    frameList=startFrame:stopFrame;
    dffList=h.dff(frameList);
            
    frameSpan=stopFrame-startFrame;
    nFrames=frameSpan+1;
    totalVal=sum(dffList);
    meanVal=totalVal/nFrames;
    semVal=sqrt(sum((dffList-meanVal).^2)/(frameSpan*nFrames));
    
    set(h.integralTx,'String',num2str(totalVal));
    set(h.meanTx, 'String', num2str(meanVal));
    set(h.semTx, 'String', num2str(semVal));
    
    
    %get(h.integrationStartEd,'UserData')
    %get(h.integrationStopEd,'UserData')
    %startFrame
    %stopFrame
    %frameSpan
	
    try
        delete(h.integrationRects);
        h.integrationRects=[];
        %delete(get(hObject,'UserData'));
        %set(hObject,'UserData',[]);
    catch ME
        %disp(ME.message);
    end
    try
        rectHandles=zeros(nFrames,1);
        j=0;
        for i=startFrame:stopFrame
            j=j+1;
   
            if i<length(h.times)
                rectHandles(j)=rectangle('Parent', h.ax,...
                    'Position',[h.times(i),min([0,h.dff(i)]),h.times(i+1)-h.times(i),abs(h.dff(i))]);
            end
        end
        h.integrationRects=rectHandles;%(rectHandles~=0);
        %set(hObject,'UserData',rectHandles);
    catch ME
        disp(ME.message);
    end
end

guidata(h.integrationFig,h);


function [h,useOldBaseline]=getDff(h)
baselineStartString=get(h.baselineStartEd,'String');
baselineStopString=get(h.baselineStopEd,'String');
baselineStartTime=str2double(baselineStartString);
baselineStopTime=str2double(baselineStopString);

if ~isValidNum(baselineStartTime) || ~isValidNum(baselineStopTime)
    disp('Invalid baseline start and/or stop time.');
    useOldBaseline=true;
    return
end


postStartBaselineInds=h.times>=baselineStartTime;
preStopBaselineInds=h.times<=baselineStopTime;
baselineLogicalInds=postStartBaselineInds & preStopBaselineInds;

if sum(baselineInds)<=2
    disp('There must be AT LEAST two data points between the start and stop of the baseline.');
    useOldBaseline=true;
    return
end

h.baselineTimes=[baselineStartTime, baselineStopTime];
prevBaseline=h.baselineInds;
h.baselineInds=find(baselineLogicalInds);

h.rawFitVals=h.b(h.baselineInds, h.rawColumn);
h.controlFitVals=h.b(h.baselineInds, h.controlColumn); 

nColumnsShown=length(h.dffColumn);
fitMat=zeros(nColumnsShown,2);
p=pinv([h.b(h.baselineInds,h.controlColumn), ones(sum(h.baselineLogicalInds),1)]);
rawVals=h.b(:,h.rawColumn);
[fitVals, dffVals]=deal(zeros(size(h.b(:,h.rawColumn))));

for i=1:nColumnsShown
    j=h.rawColumn(i);
    [fitMat(i,:)]=p*h.b(h.baselineInds,j);
    fitVals(:,i)=fitMat(i,1)*h.b(:,controlColumn)+fitMat(i,2);
    dffVals(:,i)=(rawVals(:,i)-fitVals(:,i))./fitVals(:,i);
    h.b=[h.b,fitVals,dffVals];
end



if strcmpi(get(h.baselineFig,'Visible'),'on')
    
    
else
    
end

%If the figure with the control fit is visible, update it. If it doesn't
%exist, create it, and set its visibility accordingly.


%Try csvread raw, catch it if it fails. If the error is the one you get
%from trying to read non-numerical data with csvread, then try xlsread. 
%If successful, set program to require two input files, then make the user 
%give a second file containing control data, and average each one along
%with time data. Make sure time data matches.



% set(h.baselineFitLine,'XData',,'YData',);
% st(h.baselineFitScatter,'XData',,'YData',);




%Assuming legitimate baseline has been selected, now perform fit for
%display and show DFF values.










function isValidNum=isValidNumericEntry(numString)
n=str2double(numString);

isValidNum=~isempty(n) && isnumeric(n) && ~isnan(n);






function [p,resp] = getLinFit(time,baseOn,baseOff,stimOn,stimOff,fitmax,dfof)
%% Fitting function by Ann Kennedy:
% The model I'm fitting is as follows:
% first, when the stimulus is turned on, the effective stimulus strength
% decays exponentially (think of this as some sort of adaptation/
% habituation to the stimulus):
%           tau_s ds/dt = -s + delta(t-stimOn)
% (where delta == Dirac delta function, so delta(t-stimOn) = 1 at t=stimOn
% and 0 everywhere else, and s(t) == the stimulus strength at time t.)
% This equation has the solution s(t) = (tau_s^-1) exp(-t/tau_s)
% So this describes the strength of the stimulus over time: we next need to
% describe how it feeds into the VMHdm circuit. If we make the fairly
% general assumption that VMHdm works like a leaky integrator with some
% time constant tau_d, then the equation for the VMHdm dynamics is
%           tau_d dr/dt = -r + s(t)
% This has the solution r(T) = (tau_d^-1)\int_0^T exp((t-T)/tau_d) s(t) dt
% where the equation for s(t) is given above. When we plug in this equation
% and evaluate the integral (which is a convolution of two exponentials),
% we find the solution
%           r(t) = (tau_d - tau_s)^-1 (exp(-t/tau_d) - exp(-t/tau_s))
% I make one further simplifying solution here, which is that tau_d >
% tau_s; this allows me to reparameterize as: tau_s = p(1), tau_d = p(1) +
% p(3) in the equations below (which makes the fits a little better
% behaved.)

dt      = 1/mean(time(2:end)-time(1:end-1));
frOn    = round(stimOn*dt  - time(1)*dt); %frame on which the stimulus turns on
bOn     = round(baseOn*dt  - time(1)*dt);
bOff    = round(baseOff*dt - time(1)*dt);

bOn=max([bOn,1]);
bOff=max([bOff,1]); %Ensure that these indices are greater than 0


inds = time<fitmax;
t = time(inds)-stimOn; %perform fit over the range specified by fitmax
% the equation for fit (below) is derived in my big comment above. I don't
% bother fitting the baseline, but rather assume it to be the mean value of
% dfof up to the time the stimulus turns on.
fit     = @(p) mean(dfof(bOn:bOff)) + p(2)^2/p(3)^2 * (t>0).*(exp(-t/(p(1)^2+p(3)^2)) - exp(-t/p(1)^2));
f       = @(p) mean((dfof(inds) - fit(p)).^2); %function to minimize == MSE between dfof and fit

p0 = [1 sqrt(max(dfof) - mean(dfof(bOn:bOff))) 3]; %initial conditions guess
p = fminunc(f,p0); %fminunc == "minimize the value of this function"

t = time - stimOn; %after fitting, reevaluate the same function over the entire time range
fit     = @(p) mean(dfof(bOn:bOff)) + p(2)^2/p(3)^2 * (t>0).*(exp(-t/(p(1)^2+p(3)^2)) - exp(-t/p(1)^2));
resp = fit(p);
