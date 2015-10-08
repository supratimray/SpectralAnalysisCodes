%%%%%%%%%%%%%%%%%%%%%%%%%%%%____Figure 5____%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program computes time-frequency spectrum using the HHT method
% Add to path 'Matlab runcode' available at -
% http://rcada.ncu.edu.tw/research1.htm
% Add to path the folders containing the input files:
% elec81.mat and lfpInfo.mat
%
% Vinay Shirhatti, Ashutosh Mishra
% Dr. Supratim Ray's Lab, IISc, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%______Input Signal_____%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

electrodeNum=81;
load(['elec' num2str(electrodeNum) '.mat']); % Loads the signal (data)
load('lfpInfo.mat'); % Loads timing information about the signal

bl=[-0.3 -0.1]; % baseline period
blL=find(timeVals>=bl(1),1);            % lower baseline index
blU=find(timeVals<bl(2),1,'last');      % upper baseline index

% Choose the method - EMD or EEMD
runEMDFlag = 0; % set this 1 for EMD, anything else to run EEMD

%%
if runEMDFlag
    hilbertFileName = 'hhtSpectraEMD.mat';
else
    hilbertFileName = 'hhtSpectraEEMD.mat'; %#ok<*UNRCH>
end

if exist(hilbertFileName,'file')
    disp('Loading saved Hilbert Data...');
    load(hilbertFileName);
else
    
    disp('Generating Hilbert Data...');
    %%%%%%%%%%%%%%%%%%%%_______HHT Parameters____________%%%%%%%%%%%%%%%%%%%%%%
    %---------------Parameters for running the nnspe command-------------------
    
    % Required parameters for implementation
    t0=timeVals(1);         % starting time value of the signal
    t1=timeVals(end);       % last time value of the signal
    fres=201;               % frequency resolution i.e. no. of frequencies between the maximum and minimum frequencies
    tres=length(timeVals);  % time resolution i.e. no. of time points between the maximum and minimum time values
    fw0=0;                  % minimum frequency
    fw1=200;                % maximum frequency
    tw0=t0;                 % minimum time for zooming
    tw1=t1;                 % maximum time for zooming
    lscale=[];              % => linear axis
    nfilter=0;              % number of points to be median filtered
    
    if runEMDFlag == 1
        ifmethod=[];        % instantaneous frequency method
        normmethod=[];      % normalization method
        Nstd = 0;
        NE = 1;
    else
        ifmethod = 'hilbert';   % instantaneous frequency method
        normmethod = 'hilbert'; % normalization method
        % to perform eemd, chosing Nstd = 0.1, NE = 100
        Nstd = 0.1;             % ratio of added noise
        NE = 100;               % ensemble number
    end
    % Refer nnspe.m for more information on these parameters
    
    
    %%%%%%%%%%%%%%%%________HHT Implementation__________%%%%%%%%%%%%%%%%%%%%%%%
    
    if runEMDFlag == 1
        disp('Implementing EMD');
    else
        disp(['Implementing EEMD, with Nstd = ' num2str(Nstd) ', NE = ' num2str(NE)]);
    end
    
    % Initialize the energy distribution matrices
    E = zeros(fres,tres);
    
    %------------------For a single trial----------------
    % Eg. Trial no. 72
    trialNum = 72;
    disp(['Computing for a single trial, trial number: ' num2str(trialNum)]);
    emdList = eemd(data(trialNum,:),Nstd,NE);
    % Compute the HHT energy spectrum
    [E1,t,f] = nnspe(emdList(:,2:10),t0,t1,fres,tres,fw0,fw1,tw0,tw1,ifmethod,normmethod,nfilter,lscale);
    % 1st IMF is the signal itself in this list, so ignoring this
    % And taking the first 9 IMFs
    
    %------------------For multiple trials----------------
    numTrials = size(data,1);
    
    for i=1:numTrials
        disp(['trial ' num2str(i) ' of ' num2str(numTrials)]);
        emdList = eemd(data(i,:),Nstd,NE);
        [P,t,f] = nnspe(emdList(:,2:10),t0,t1,fres,tres,fw0,fw1,tw0,tw1,ifmethod,normmethod,nfilter,lscale);
        E=E+P;
    end
    E = E/numTrials; % averaged across trials
    
    if runEMDFlag
        save(hilbertFileName,'E1','E','t','f'); % save the HHT spectrum matrices
    else
        save(hilbertFileName,'E1','E','t','f'); % save the HHT spectrum matrices
    end
end

%%
%------------------Change in power from Baseline-----------------

% limiting/trunctating the lowest values of for meaningful log calculations
% otherwise the log of extremely small values gives '-infinity' and
% discolours the time-frequency spectrum
tE = E;
thresh = 10^(-16);
tE(tE<thresh)=thresh;

% log of average spectrum across trials at every frequency and time
logtE = log10(tE);

% log of average spectrum across trials at every frequency - in the baseline period
blogtE = logtE(:,blL:blU);

% mean baseline power at every frequency
meanblogtE = mean(blogtE,2);

% calculate the change in power from the baseline
meanblogtE = repmat(meanblogtE,1,length(timeVals));
difftEn = 10*(log10(tE) - meanblogtE); % in dB

%%

%%%%%%%%%%%%%%%%%%%%%___________Plots_________________%%%%%%%%%%%%%%%%%%%%%
% Plots without smoothing of the HHT spectra
%%%%%%%%%%%%%%%%%%%%%%______Display parameters______%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
timeLimsS = [-0.1 0.5]; % time interval to be displayed in seconds. Stimulus onset is at 0.
freqLimsHz = [0 150];   % Range of frequencies in Hz to be displayed
cLims = [-1 3];         % color axis limits for spectra
cLimsDiff = [-5 12.5];  % color axis limits for difference spectrum
fontSizeSmall=13;       % for Plots
fontSizeLarge=30;       % for Plots

% compute a thresholded form of E1 for display
tE1 = E1;
thresh = 10^(-16);
tE1(tE1<thresh)=thresh;

tfData{1} = tE1;        % single trial
tfData{2} = logtE;      % average of all trials
tfData{3} = difftEn;    % change in power from baseline


hTF=cell(1,3);
for i=1:3
    if runEMDFlag
        hTF{i} = subplot(2,3,i);  % figure handle
    else
        hTF{i} = subplot(2,3,3+i);  % figure handle
    end
    imagesc(t,f,tfData{i},'Parent',hTF{i});
    set(hTF{i},'YDir','normal');
    axis(hTF{i},[timeLimsS freqLimsHz]);
    if i==3
        caxis(hTF{i},cLimsDiff);
    else
        caxis(hTF{i},cLims);
    end
end

%------------------------- Add Figure Details -----------------------------
titleList =[{'Single-trial Spectrum'} {'Average Spectrum'} {'Change in Power'}];
if runEMDFlag
    plotNumberList = {'A','B','C'};
else
    plotNumberList = {'D','E','F'};
end

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


%%
%%%%%%%%%%%%%%%%%____Smoothing of the HHT spectra____%%%%%%%%%%%%%%%%%%%%%%

gaussFilter = fspecial('gaussian', 8, 1); % a Gaussian low-pass filter

% filter the single-trial and the average energy distributions, E1 and E
E1f = filter2(gaussFilter, E1);
Ef = filter2(gaussFilter, E);

% limiting/trunctating the lowest values of for meaningful log calculations
% otherwise the log of extremely small values gives '-infinity' and
% discolours the time-frequency spectrum
tEf = Ef;
thresh = 10^(-16);
tEf(tEf<thresh)=thresh;

% log of average spectrum across trials at every frequency and time
logtEf = log10(tEf);

% log of average spectrum across trials at every frequency - in the baseline period
blogtEf = logtEf(:,blL:blU);

% mean baseline power at every frequency
meanblogtEf = mean(blogtEf,2);

% calculate the change in power from the baseline
meanblogtEf = repmat(meanblogtEf,1,length(timeVals));
difftEnf = 10*(log10(tEf) - meanblogtEf); % in dB

%%
%%%%%%%%%%%%%%%%%%%%%_______Plots after smoothing_________%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%______Display parameters______%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
timeLimsS = [-0.1 0.5]; % time interval to be displayed in seconds. Stimulus onset is at 0.
freqLimsHz = [0 150];   % Range of frequencies in Hz to be displayed
cLims = [-1 3];         % color axis limits for spectra
cLimsDiff = [-5 12.5];  % color axis limits for difference spectrum
fontSizeSmall=13;       % for Plots
fontSizeLarge=30;       % for Plots

% compute a thresholded form of E1f for display
tE1f = E1f;
thresh = 10^(-16);
tE1f(tE1f<thresh)=thresh;

tfDataf{1} = tE1f;        % single trial
tfDataf{2} = logtEf;      % average of all trials
tfDataf{3} = difftEnf;    % change in power from baseline

hTF2=cell(1,3);
for i=1:3
    if runEMDFlag
        hTF2{i} = subplot(2,3,i);  % figure handle
    else
        hTF2{i} = subplot(2,3,3+i);  % figure handle
    end
    imagesc(t,f,tfDataf{i},'Parent',hTF2{i});
    set(hTF2{i},'YDir','normal');
    axis(hTF2{i},[timeLimsS freqLimsHz]);
    if i==3
        caxis(hTF2{i},cLimsDiff);
    else
        caxis(hTF2{i},cLims);
    end
end

%------------------------- Add Figure Details -----------------------------
titleList =[{'Single-trial Spectrum'} {'Average Spectrum'} {'Change in Power'}];
if runEMDFlag
    plotNumberList = {'A','B','C'};
else
    plotNumberList = {'D','E','F'};
end

for i=1:3
    % Titles
    title(hTF2{i},titleList{i},'FontSize',fontSizeSmall);
    % Labels and Fonts
    set(hTF2{i},'XTick',[0 0.4],'YTick',(0:1/3:1)*freqLimsHz(2),'FontSize',fontSizeSmall);
    xlabel(hTF2{i},'Time (s)','FontSize',fontSizeSmall);
    % Boxes
    box(hTF2{i},'on');
    %colorbars
    colorbar('peer',hTF2{i},'location','northoutside');
    % Plot Numbers
    text(-0.1,1.5, plotNumberList{i}, 'unit','normalized','fontsize',fontSizeLarge,'Parent',hTF2{i});
end
ylabel(hTF2{1},'Frequency (Hz)','FontSize',fontSizeSmall);
