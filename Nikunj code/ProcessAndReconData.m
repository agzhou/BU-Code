% Process AZ001 Data
clearvars
close all

%% Initialize

% Extract Filepath relative to current file location
currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

% Define a few additional filepath variables. Uncomment Desired folder to
% process
% dataFolder = "Allen Data\01-29-2025 AZ001 ULM Data\run 0 baseline\";
dataFolder = "Allen Data\01-29-2025 AZ001 ULM Data\run 1 left eye\";
% dataFolder = "Allen Data\01-29-2025 AZ001 ULM Data\run 2 right eye\";


% Load parameters
tic; p = initMultiFileParams(fullfile(dataFilePath,dataFolder,"params.mat")); toc


%% Test Beamform
% Extract Filenames
% fileNames = extractFileList(fullfile(dataFilePath,dataFolder));
% tic; load(fullfile(dataFilePath,dataFolder,fileNames{1})); toc

pN = computeNewGrid3D(p,[1,p.szX],[1,p.szY],[1,p.szZ],p.szX*2,p.szY*2,p.szZ);
tic; cRF = computeCRF(RcvData(:,:,1),pN);
Recon = RCA_DAS_CPP(pN,cRF);
Recon = reshape(Recon,[pN.szY,pN.szX,pN.szZ]); toc

plotImgs(Recon,pN)
figure
volshow(abs(Recon))

%% Process baseline Dataset
dataFolder = "Allen Data\01-29-2025 AZ001 ULM Data\run 0 baseline\";

% Insert desired code here

%% Process run 1 left eye
dataFolder = "Allen Data\01-29-2025 AZ001 ULM Data\run 1 left eye\";
fileNames = extractFileList(fullfile(dataFilePath,dataFolder));
ReconFolder = ""; % NOTE INSERT FOLDER PATH HERE WHERE YOU WANT RECON DATA SAVED

for nFile = 1:length(fileNames)
    load(fullfile(dataFilePath,dataFolder,fileNames{1}));
    cRF = computeCRF(RcvData,p);
    clearvars RcvData
    Recon = zeros(p.szX*p.szY,p.szZ,size(cRF,3));
    for nFrame = 1:size(cRF,3)
        Recon(:,:,nFrame) = RCA_DAS_CPP(p,cRF(:,:,nFrame));
    end
    clearvars cRF
    filename = ReconFolder + "ReconDat_" + num2str(nFile) + ".mat";
    save(filename,'Recon','-v7.3');
    clearvars Recon
end

%% Process run 2 right eye
dataFolder = "Allen Data\01-29-2025 AZ001 ULM Data\run 2 right eye\";

% Insert desired code here

%% Helper Functions
function [p] = initMultiFileParams(dataFile)

    load(dataFile);
    struct2vars(P); % Stolen from the Matlab File Exchange
    
    % Define the default PData from the setup script since PData doesn't
    % appear to be saved
    PData.PDelta = [Trans.spacing, Trans.spacing, 0.5];
    PData.Size(1) = 80;
    PData.Size(2) = 80;
    PData.Size(3) = ceil((endDepth-startDepth)/PData.PDelta(3)); % startDepth, endDepth and pdelta set PData.Size(3).
    PData.Origin = [-Trans.spacing*39.5,Trans.spacing*39.5,startDepth]; % x,y,z of upper lft crnr of page 1.
    
    % Pre define a few more parameters
    lambda = Resource.Parameters.speedOfSound/(Trans.frequency*1e6);
    TXangle = reshape([TX(:).Steer],2,[]).';    % Extracts all TX angles
    xCoord = (PData(1).Origin(1) + (0:PData(1).Size(2)-1)*PData(1).PDelta(1))*lambda;
    yCoord = (PData(1).Origin(2) - (0:PData(1).Size(1)-1)*PData(1).PDelta(2))*lambda;
    zCoord = (PData(1).Origin(3) + (0:PData(1).Size(3)-1)*PData(1).PDelta(3))*lambda;

    % Define TXApod (column 1 indicates transmitting on X, column 2 indicates
    % transmitting on Y)
    TXApod = zeros(length(TX),2);
    for i = 1:length(TX)
        if (any(TX(i).Apod(1:end/2)))
            TXApod(i,1) = 1;
        elseif (any(TX(i).Apod(end/2+1:end)))
            TXApod(i,2) = 1;
        end
    end
    
    p = struct('fs',Receive(1).decimSampleRate*1e6,... % sampling frequency [Hz]
        'pitch', Trans.spacingMm*1e-3,... % Element Spacing [m]
        'fc', Trans.frequency*1e6,... % center frequency [Hz]
        'c', Resource.Parameters.speedOfSound,... % speed of sound [m/s]
        'fnumber', [0.1,0.1],...cot(2*asin(sin(range(TXangle)/2))),... % angular aperature ratio [ul]
        't0',(-Receive(1).startDepth + TW(1).peak)/(Trans.frequency*1e6),... % start of transmit [s]
        'TXangle',TXangle,... % vector of transmit angles where col 1 is tilt along x and col 2 is tilt along y [rad] 
        'ElemPos',Trans.ElementPos(:,1:2)*lambda,... % element position [m]
        'xCoord',xCoord,... % x-coordinates of grid [m]
        'yCoord',yCoord,...
        'zCoord',zCoord,... % z-coordinates of grid [m]
        'numEl',int32(Trans.numelements),... % Number of elements [ul]
        'szRF',int32(Resource.RcvBuffer(1).rowsPerFrame),... %
        'szAcq',int32(Receive(1).endSample),... % Number of time samples in dataset [ul]
        'szX',int32(length(xCoord)),... % length of x-coordinate [ul]
        'szY',int32(length(yCoord)),...
        'szZ',int32(length(zCoord)),... % length of z-coordinate [ul]
        'na',int32(length(TX)),... % Number of beams [ul]
        'nc',int32(Resource.Parameters.numRcvChannels),... % Number of channels [ul]
        'ConnMap',int32(Trans.Connector),... % Element connector mapping [el]
        'startSample',[Receive(1:length(TX)).startSample],...
        'endSample',[Receive(1:length(TX)).endSample],...
        'TXApod', int32(TXApod));

    p.nPoints = p.szX*p.szY*p.szZ;
    p.L = double(p.numEl)/2*p.pitch;
end


function [sortedFileNames] = extractFileList(folderPath)

    % Extract filenames (excluding '.' and '..')
    files = dir(folderPath); % Get directory contents
    fileNames = {files(~[files.isdir]).name};
    
    % Define regex pattern to match filenames
    pattern = "RF-5-11-500-200-1-(\d+).mat";
    
    % Filter filenames that match the pattern and extract numbers
    filteredFiles = [];
    fileNumbers = [];
    
    for i = 1:length(fileNames)
        token = regexp(fileNames{i}, pattern, 'tokens');
        if ~isempty(token)
            fileNumbers(end+1) = str2double(token{1}{1}); % Extract filenumber
            filteredFiles{end+1} = fileNames{i}; % Store valid filenames
        end
    end
    
    % Sort files by filenumber
    [~, sortIdx] = sort(fileNumbers);
    sortedFileNames = filteredFiles(sortIdx).';

end

function plotImgs(IQb,p)
    figure('Position',[151,111,1234,653]);
    tiledlayout("flow");

    nexttile
    plotGammaScaleImage(p.xCoord*1e3,p.zCoord*1e3,abs(squeeze(IQb(:,20,:))).',0.25)
    axis image
    title("XZ Plane DAS")
    % colormap jet

    nexttile
    plotGammaScaleImage(p.yCoord*1e3,p.zCoord*1e3,abs(squeeze(IQb(20,:,:))).',0.25)
    axis image
    title("YZ Plane DAS")
    % colormap jet

    nexttile
    plotGammaScaleImage(p.xCoord*1e3,p.yCoord*1e3,abs(squeeze(IQb(:,:,20))),0.25)
    axis image
    title("XY Plane DAS")
    % colormap jet

end

function plotImgs2(IQb,varargin)

    if (nargin == 2)
        p = varargin{1};
        [xi,yi,zi] = meshgrid(p.xCoord,p.yCoord,p.zCoord);
    end

    env = abs(IQb);
%     I = 20*log10(env/max(env(:)));
    I = env;

    % Display
    figure
    I(1:round(size(I,1)/2),1:round(size(I,2)/2),:) = NaN;
    for k = [-40:10:-10 -5 -1]
        isosurface(xi*1e2,yi*1e2,zi*1e2,I,k)
    end
    view(-60,40)
    colormap([1-hot;hot])
    c = colorbar;
    c.Label.String = 'dB';
    box on, grid on
    zlabel('[cm]')
    title('PSF at the focal point [dB]')
%     xlim([-1,1])
%     ylim([-1,1])
%     zlim([2.4,3.2])
end