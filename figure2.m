%%%%%%%%%%%%%%%%%%%%%%%%%%____Figure 2____%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program computes the  time frequency spectra of a LFP signal using
% the STFT, multi-taper (MTM) and Wavelet Transform (WT). It shows the
% tiling diagram and some basis functions also.

% STFT and MTM Plots (plot A-C) are analyzed and plotted first. For these,
% Chronux software package is needed and can be downloaded from
% http://chronux.org/.

% Wavelet analysis (plot D) requires Matlab Wavelet toolbox. The program
% saves the analysis results in an intermediate file. In case the toolbox
% is not available, results can still be obtained by loading the
% intermediate file (IntermediateFileWavelet.mat) which is provided. 
% Delete or rename the intermediate file if you want to run the 
% Wavelet analysis from scratch

% The tiling diagram uses the length of window (T) and its inverse
% (Rayleigh frequency: 1/T) as the sides of the tile. Therefore each tile
% has an area of 1, irrespective of the window shape.

% Subhash Chandran and Supratim Ray

%%%%%%%%%%%%%%%%%%%%%______Input Signal_____%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
electrodeNum=81;
load(['elec' num2str(electrodeNum) '.mat']);
load('lfpInfo.mat'); % Loads timing information about the signal

%%%%%%%%%%%%%%%%%%%%% STFT and MTM Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqLimsHz = [0 150];        % frequencies to display in Hz.
winLenS    = [0.05 0.2 0.2]; % window lengths for STFT short, STFT long and MTM
BWList     = [1 1 2];        % Bandwidths. BW=1 for STFT
winStepS   = 0.001;          % Window step size in seconds
baselineS  = [-0.3 -0.1];    % Baseline period for computing change
tStartTile = 0;              % Assume the tiles have an edge here.

Ts=timeVals(2)-timeVals(1); % Sampling period
Fs = round(1/Ts);           % Hz

params.pad = -1;
params.Fs = Fs;
params.fpass =[freqLimsHz(1) freqLimsHz(2)+50] ;
params.trialave = 0;

for i=1:3
    params.tapers = [BWList(i) 2*BWList(i)-1];
    movingWin  = [winLenS(i) winStepS];
    disp(['Running MTM analysis with tapers: ' num2str(params.tapers) ', WindowLength: ' num2str(movingWin(1))]);
    
    [S{i},tList{i},fList{i}] = mtspecgramc(data',movingWin,params); %#ok<*SAGROW>
    logMeanS{i} = log10(mean(S{i},3));  % Average across trials and take log
    tList{i} = tList{i}+timeVals(1)-1/Fs; % Center the times with respect to the stimulus onset time
    
    % Change in power
    blPos = intersect(find(tList{i}>=baselineS(1)),find(tList{i}<baselineS(2)));
    baseline = mean(logMeanS{i}(blPos,:),1);
    dPowerdB{i} = 10*(logMeanS{i}'-repmat(baseline,length(tList{i}),1)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%______Display ______%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeLimsS = [-0.1 0.5];  % time interval to be displayed in seconds. Stimulus onset is at 0.
cLims = [-1 3];          % colormap limits for the spectrogram
cLimsWavelet=[2 6];      % colormap limits for the scalogram 
cLimsDiff = [-5 12.5];    % colormap limits for change in power

trialNum = 72;           % Trial number to be used for analysis
bfFrequency = [50 100];  % Display these basis functions
bfColors       = 'krb';  % Display basis functions using these colors
bfWindowColors = {[0.7 0.7 0.7],'m','c'};  % Display windows using these colors
fontSizeLarge=30;
fontSizeSmall=13;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure Handles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:4
    hBasis1{i} = subplot('position',[0.05  0.975-0.225*i 0.15 0.09]);
    hBasis2{i} = subplot('position',[0.05  0.975-0.225*i+0.11 0.15 0.09]);
    hTile{i}   = subplot('position',[0.3  0.975-0.225*i 0.15 0.2]); %#ok<*SAGROW>
    hSingle{i} = subplot('position',[0.475 0.975-0.225*i 0.15 0.2]);
    hAll{i}    = subplot('position',[0.65  0.975-0.225*i 0.15 0.2]);
    hDiff{i}   = subplot('position',[0.825 0.975-0.225*i 0.15 0.2]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot A-C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:3
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Column 1: Basis Fn %%%%%%%%%%%%%%%%%%%%%%%%%
    t=-winLenS(i)/2:Ts:(winLenS(i)/2)-Ts;
    windowSignalAll=dpsschk([BWList(i) 2*BWList(i)-1],round(Fs*winLenS(i)),Fs);
    
    for j=1:min(size(windowSignalAll,2),2)
        windowSignal = windowSignalAll(:,j);
        basis1=windowSignal'.*cos(2*pi*bfFrequency(1)*t);
        basis2=windowSignal'.*cos(2*pi*bfFrequency(2)*t);
        
        plot(hBasis1{i},t,basis1,'color',bfColors(j));
        hold(hBasis1{i},'on');
        plot(hBasis1{i},t,windowSignal,'color',bfWindowColors{j},'linestyle','--');
        plot(hBasis1{i},t,-windowSignal,'color',bfWindowColors{j},'linestyle','--');
        xlim(hBasis1{i},[-max(winLenS)/2 max(winLenS)/2]);
        
        plot(hBasis2{i},t,basis2,'color',bfColors(j));
        hold(hBasis2{i},'on');
        plot(hBasis2{i},t,windowSignal,'color',bfWindowColors{j},'linestyle','--');
        plot(hBasis2{i},t,-windowSignal,'color',bfWindowColors{j},'linestyle','--');
        xlim(hBasis2{i},[-max(winLenS)/2 max(winLenS)/2]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Column 2: Tiling %%%%%%%%%%%%%%%%%%%%%%%%%%%
    tListTiling = unique([flipdim((tStartTile:-winLenS(i):timeVals(1)),2) tStartTile:winLenS(i):timeVals(end)]);
    for j=1:length(tListTiling)        % vertical line
        x=[tListTiling(j) tListTiling(j)];
        line(x,freqLimsHz,'color','k','Parent',hTile{i});
        hold(hTile{i},'on');
    end
    
    fRes = (BWList(i)/winLenS(i));
    for j=freqLimsHz(1)-fRes/2:fRes:freqLimsHz(2)+fRes/2 %horizontal line
        y=[j j];
        line(timeLimsS,y,'color','k','Parent',hTile{i});
        hold(hTile{i},'on');
    end
    axis(hTile{i},[timeLimsS freqLimsHz]);
    
    %%%%%%%%%%%%%%%%% Column 3: TF plot for single trial %%%%%%%%%%%%%%%%%%
    pcolor(tList{i},fList{i},log10(squeeze(S{i}(:,:,trialNum)')),'Parent',hSingle{i});
    shading(hSingle{i},'interp');
    axis(hSingle{i},[timeLimsS freqLimsHz]);
    caxis(hSingle{i},cLims);
    
    %%%%%%%%%%%%%%%%% Column 4: TF plot for all trials %%%%%%%%%%%%%%%%%%%%
    pcolor(tList{i},fList{i},logMeanS{i}','Parent',hAll{i});
    shading(hAll{i},'interp');
    axis(hAll{i},[timeLimsS freqLimsHz]);
    caxis(hAll{i},cLims);
    
    %%%%%%%%%%%%%%%%%%% Column 5: Change power plot  %%%%%%%%%%%%%%%%%%%%%%
    pcolor(tList{i},fList{i},dPowerdB{i},'Parent',hDiff{i});
    shading(hDiff{i},'interp');
    axis(hDiff{i},[timeLimsS freqLimsHz]);
    caxis(hDiff{i},cLimsDiff);
end

%%%%%%%%%%%%%%%___________ Wavelet Analysis_____________%%%%%%%%%%%%%%%%%%%
if exist('IntermediateFileWavelet.mat','file')
    disp('Loading saved Wavelet analysis results...');
    load('IntermediateFileWavelet.mat');  % If intermediate file is present, load that
else
    disp('Saved Wavelet analysis results not found. Running wavelet analysis (needs wavelet toolbox)...');

    % We make some approximations to represent the wavelet decomposition
    % (which is a time-scale decomposition) into its corresponding
    % time-frequency representation. The command cwt in matlab computes the
    % wavelet coefficients at any positive real scale factor by
    % stretching/compressing the mother wavelet by that factor. Therefore,
    % we need to choose an appropriate set of scales.  
    
    % For time-frequency tiling, we assume a dyadic decomposition where the
    % frequency space is divided into two segments and the process is
    % repeated for the lower half. This gives us the desired center
    % frequencies of the wavelet decomposition.
    
    maxPower=10;
    fListTilingWavelet = Fs./(2.^(1:maxPower));  % desired centre frequencies of the wavelet tiles
    fListWavelet = Fs./(2.^(1:1/4:maxPower));    % more dense sampling in the frequency domain (leads to smoothing in frequency domain) 
    
    % We choose the morlet wavelet
    wname = 'morl';
    [psi,xVal]=wavefun(wname);
    cf = centfrq(wname);                                % centre frequency of morlet wavelet
    
    % We compute the scales for which the wavelet has the desired centre
    % frequency.
    scaleList = cf./(Ts*fListWavelet);                  % list of scales
    tListWavelet = timeVals;

    scAll=zeros(length(fListWavelet),length(tListWavelet));
    numTrials=size(data,1);
    
    for i=1:numTrials                                    % Wavelet coefficients of all trials
       cw1 = cwt(data(i,:),scaleList,wname);             % Continuous Wavelet Transform
       scAll=scAll+cw1.^2;                               % Squaring to get the scalogram
       
       if i==trialNum
           scSingle = cw1.^2;
       end
    end
    scAll = scAll/numTrials;
    
    % Save scalograms
    save('IntermediateFileWavelet.mat','scSingle','scAll','fListWavelet','fListTilingWavelet','tListWavelet','scaleList','psi','xVal','cf','Ts','Fs');
end

logMeanSc = log10(scAll);  % Average across trials and take log

% Change in power
blPos = intersect(find(tListWavelet>=baselineS(1)),find(tListWavelet<baselineS(2)));
baselineSc = mean(logMeanSc(:,blPos),2);
dScdB = 10*(logMeanSc-repmat(baselineSc,1,length(tListWavelet)));

%%%%%%%%%%%%%%%____________ Plot D: Wavelet_____________%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Column 1: Tiling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(fListTilingWavelet)
    
    f=fListTilingWavelet(i);
    plot(hTile{4},timeVals,f+zeros(1,length(timeVals)),'k'); % For simplicity, take these as edges, not as centers
    hold(hTile{4},'on');
    
    if i>1
        fLine = fListTilingWavelet(i):0.1:fListTilingWavelet(i-1);
        tListTilingWavelet = unique([flipdim((tStartTile:-1/f:timeVals(1)),2) tStartTile:1/f:timeVals(end)]);
        for j=1:length(tListTilingWavelet)
            plot(hTile{4},tListTilingWavelet(j)+zeros(1,length(fLine)),fLine,'k');
        end
    end
end
axis(hTile{4},[timeLimsS freqLimsHz]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Column 2: Basis Fn %%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:2 % Basis Function
    % Find the scale at which center frequency becomes bfFrequency
    s = cf./(Ts*bfFrequency(i));
    xValBasis{i} = xVal*s/Fs;
end
plot(hBasis1{4},xValBasis{1},psi/sum(psi.^2),bfColors(1)); 
xlim(hBasis1{4},[-max(winLenS)/2 max(winLenS)/2]);
plot(hBasis2{4},xValBasis{2},psi/sum(psi.^2),bfColors(1));
xlim(hBasis2{4},[-max(winLenS)/2 max(winLenS)/2]);

%%%%%%%%%%%%%%%%%%%% Column 3: TF plot for single trial %%%%%%%%%%%%%%%%%%%
pcolor(tListWavelet,fListWavelet,log10(scSingle),'Parent',hSingle{4});
shading(hSingle{4},'interp');
axis(hSingle{4},[timeLimsS freqLimsHz]);
caxis(hSingle{4},cLimsWavelet);

%%%%%%%%%%%%%%%%%%%% Column 4: TF plot for all trials %%%%%%%%%%%%%%%%%%%%%
pcolor(tListWavelet,fListWavelet,logMeanSc,'Parent',hAll{4});
shading(hAll{4},'interp');
axis(hAll{4},[timeLimsS freqLimsHz]);
caxis(hAll{4},cLimsWavelet);

%%%%%%%%%%%%%%%%%%%%%% Column 5: Change power plot  %%%%%%%%%%%%%%%%%%%%%%%
pcolor(tListWavelet,fListWavelet,dScdB,'Parent',hDiff{4}); 
shading(hDiff{4},'interp');
axis(hDiff{4},[timeLimsS freqLimsHz]);
caxis(hDiff{4},cLimsDiff);

%%%%%%%%%%%%%%%%%%%%%%%%%% Add Figure Details %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Titles
title(hBasis2{1},'Basis Functions','FontSize',fontSizeSmall);
title(hTile{1},'Tiling Diagram','FontSize',fontSizeSmall);
title(hSingle{1},'Single-trial Spectrum','FontSize',fontSizeSmall);
title(hAll{1},'Average Spectrum','FontSize',fontSizeSmall);
title(hDiff{1},'Change in Power','FontSize',fontSizeSmall);

% Labels and Fonts
for i=1:3
    set(hBasis1{i},'XTicklabel',[],'FontSize',fontSizeSmall);
    set(hBasis2{i},'XTicklabel',[],'FontSize',fontSizeSmall);
    set(hTile{i},'XTick',[0 0.4],'XTicklabel',[],'FontSize',fontSizeSmall);
    set(hSingle{i},'XTick',[0 0.4],'XTicklabel',[],'YTicklabel',[],'FontSize',fontSizeSmall);
    set(hAll{i},'XTick',[0 0.4],'XTicklabel',[],'YTicklabel',[],'FontSize',fontSizeSmall);
    set(hDiff{i},'XTick',[0 0.4],'XTicklabel',[],'YTicklabel',[],'FontSize',fontSizeSmall);
    
    ylabel(hTile{i},'Frequency (Hz)','FontSize',fontSizeSmall);
end

set(hBasis1{4},'FontSize',fontSizeSmall);
set(hBasis2{4},'XTicklabel',[],'FontSize',fontSizeSmall);
set(hTile{4},'XTick',[0 0.4],'FontSize',fontSizeSmall);
set(hSingle{4},'XTick',[0 0.4],'YTicklabel',[],'FontSize',fontSizeSmall);
set(hAll{4},'XTick',[0 0.4],'YTicklabel',[],'FontSize',fontSizeSmall);
set(hDiff{4},'XTick',[0 0.4],'YTicklabel',[],'FontSize',fontSizeSmall);

ylabel(hTile{4},'Frequency (Hz)','FontSize',fontSizeSmall);

xlabel(hBasis1{4},'Time (s)','FontSize',fontSizeSmall);
xlabel(hTile{4},'Time (s)','FontSize',fontSizeSmall);
xlabel(hSingle{4},'Time (s)','FontSize',fontSizeSmall);
xlabel(hAll{4},'Time (s)','FontSize',fontSizeSmall);
xlabel(hDiff{4},'Time (s)','FontSize',fontSizeSmall);

% Boxes
for i=1:4
    box(hBasis1{i},'on');
    box(hBasis2{i},'on');
    box(hTile{i},'on');
    box(hSingle{i},'on');
    box(hAll{i},'on');
    box(hDiff{i},'on');
end

% Plot Numbers
plotNumberList = {'A','B','C','D'};
for i=1:4
    text(-0.3,1.05, plotNumberList{i}, 'unit','normalized','fontsize',fontSizeLarge,'Parent',hBasis2{i});
end