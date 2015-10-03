%%%%%%%%%%%%%%%%%%%%%%%%%%%%____Figure 4____%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program computes time-frequency spectrum using MP algorithm
% Run these codes from the MP/matlab folder or add this folder to your
% Matlab path.

% The code saves intermediate MP files in folder MPData. Remove this folder
% if you want to run MPdecomposition from scratch.


%%%%%%%%%%%%%%%%%%%%%%%%%%______Input Signal_____%%%%%%%%%%%%%%%%%%%%%%%%%%
electrodeNum=81;
trialNum = 72;              % Trial number to be used for analysis
load(['elec' num2str(electrodeNum) '.mat']);
load('lfpInfo.mat'); % Loads timing information about the signal

Fs=1/(timeVals(2)-timeVals(1));       % Sampling Frequency
[numTrials,L]=size(data);             % length of signal

%%%%%%%%%%%%%%%%%%%%%%%%%%______Perform MP_____%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folderName = 'MPData/';       % code places intermediate data in this folder.
tag = ['elec' num2str(electrodeNum) '/'];


if ~exist(platformSpecificNameMPP([folderName 'elec' num2str(electrodeNum) '/GaborMP/mp0.bok.000']),'file')
    disp('Running MP decomposition...');
    
    %Import data
    X=data';
    signalRange = [1 L];                  % full range
    importData(X,folderName,tag,signalRange,Fs);
    
    % Perform Gabor decomposition
    Numb_points = L;                      % length of the signal
    Max_iterations = 500;                 % number of iterations for MP algorithm
    runGabor(folderName,tag,Numb_points,Max_iterations);
else
    disp('MP decomposition data already exists. Skipping MP decomposition...');
end

% Retrieve information
gaborInfo = getGaborData(folderName,tag,1);
    
% Reconstruct energy
wrap=1;
atomList=[];                          % all atoms
frequency = 0:Fs/L:Fs/2;              % frequency axis
recEnAll=zeros(length(frequency),L);  % initialising energy matrix

disp('Reconstructing energy...');
for i=1:numTrials
    disp(['trial ' num2str(i) ' of ' num2str(numTrials)]);
    recEn=reconstructEnergyFromAtomsMPP(gaborInfo{i}.gaborData,L,wrap,atomList);
    if i==trialNum
        recEnSing=recEn;
        recEnAll =recEnAll + recEn;
    else
        recEnAll =recEnAll + recEn;
    end
end

logMeanEn=conv2LogMPP(recEnAll/numTrials);      % log mean TF energy matrix for all trials

% Change in power from baseline
baselineS = [-0.3 -0.1];    % define baseline
blL = find(timeVals>=baselineS(1),1);           % lower index of baseline
blU = find(timeVals<baselineS(2),1,'last');     % upper index of baseline
baselineEn=mean(logMeanEn(:,blL:blU),2);        % baseline TF Energy Matrix
dEndB = 10*(logMeanEn-repmat(baselineEn,1,size(logMeanEn,2))); % change in dB
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%____PLOTS____%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%______Display parameters______%%%%%%%%%%%%%%%%%%%%%%%
timeLimsS = [-0.1 0.5];     % time interval to be displayed in seconds. Stimulus onset is at 0.
freqLimsHz = [0 150];       % frequencies to display in Hz.   
cLims = [-1 3];             % colormap limits for the spectrogram
cLimsDiff = [-5 12.5];      % colormap limits for change in power
fontSizeLarge=30;
fontSizeSmall=13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfData{1} = conv2LogMPP(recEnSing); % single trial
tfData{2} = logMeanEn;           % average of all trials
tfData{3} = dEndB;               % change in power

hTF=cell(1,3);
for i=1:3
    hTF{i} = subplot(2,3,i);  % figure handle
    imagesc(timeVals,frequency,tfData{i},'Parent',hTF{i});
    set(hTF{i},'YDir','normal');
    axis(hTF{i},[timeLimsS freqLimsHz]);
    if i==3
        caxis(hTF{i},cLimsDiff);
    else
        caxis(hTF{i},cLims);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Add Figure Details %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
titleList =[{'Single-trial Spectrum'} {'Average Spectrum'} {'Change in Power'}];
plotNumberList = {'A','B','C'};

for i=1:3
    % Titles
    title(hTF{i},titleList{i},'FontSize',fontSizeSmall);
    % Labels and Fonts
    set(hTF{i},'XTick',[0 0.4],'YTick',(0:1/3:1)*freqLimsHz(2),'FontSize',fontSizeSmall);
    xlabel(hTF{i},'Time (s)','FontSize',fontSizeSmall);
    % Boxes
    box(hTF{i},'on');
    %colorbars
    colorbar('peer',hTF{i},'location','northoutside');
    % Plot Numbers
    text(-0.1,1.5, plotNumberList{i}, 'unit','normalized','fontsize',fontSizeLarge,'Parent',hTF{i});
end
ylabel(hTF{1},'Frequency (Hz)','FontSize',fontSizeSmall);