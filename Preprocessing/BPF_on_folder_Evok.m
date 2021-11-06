function BPF_on_folder_Evok(eFolder,meFolder,varargin)

Rotate = 0;
theta = 0;
Glut = 0;
if nargin > 2
    Rotate = varargin{1};
end
if Rotate && nargin > 3
    theta = varargin{2};
end
if nargin > 4
    Glut = varargin{3};
end
Fs = 150; % sampling rate
if Glut
   Fs = 150; 
end
Fstop1 = 0.25; Fpass1 = 0.5; Fpass2 = 6; Fstop2 = 6.5; dur = 5;
if nargin > 5
    filtOptions = varargin{4};
    Fstop1 = filtOptions.Fstop1; Fpass1 = filtOptions.Fpass1; Fpass2 = filtOptions.Fpass2;
    Fstop2 = filtOptions.Fstop2; dur = filtOptions.dur;
end

try
    badTrials = load(makeName('badTrials.txt',meFolder));
catch
    badTrials = [];
end
osFolders = dir(sprintf('%s\\*.tif',eFolder));
goodTrials = true(1,length(osFolders));
goodTrials(badTrials) = 0;
sFolders = osFolders(goodTrials);
ImgSeq_all = [];
meanImg_all = [];
fstart = zeros(1,length(sFolders));
fend = zeros(1,length(sFolders));
for jj = 1:length(sFolders)
    trialFolder = makeName(sFolders(jj).name,eFolder);
    ImgSeq = imreadalltiff(trialFolder,'slow'); ImgSeq = double(ImgSeq);
    if ~isreal(ImgSeq)
        ImgSeq = abs(ImgSeq);
    end
    ImgSeq(isnan(ImgSeq)) = 0;
    ImgSeq(isinf(ImgSeq)) = 0;
    ImgSeq_all = cat(3,ImgSeq_all,ImgSeq);
    meanImg = mean(ImgSeq,3); meanImg = repmat(meanImg,1,1,size(ImgSeq,3));
    meanImg_all = cat(3,meanImg_all,meanImg);  
    
    if jj == 1
        fstart(jj) = 1;
    else
        fstart(jj) = fend(jj-1)+1;
    end
    fend(jj) = size(ImgSeq_all,3);
end

% BPF
if Glut
    Hd = BPF_LS(Fs,Fstop1,Fpass1,Fpass2,Fstop2,dur);
    gBPF = Hd.Numerator;
    DF_F0 = tempFilt_imgStack(ImgSeq_all,gBPF);
    DF_F0 = 100*DF_F0./meanImg_all;
end

% Gaussian spatial filtering
% 2D Gaussian
param.GaussWindow = 5;
param.Sigma       = 1;
Center  = fix([param.GaussWindow/2,param.GaussWindow/2])+1;
[R,C]   = ndgrid(1:param.GaussWindow, 1:param.GaussWindow);
Gauss2d = gauss2dC(R,C,param.Sigma,Center);
[~,~,nFrames] = size(ImgSeq_all);
for idx = 1:nFrames
    DF_F0(:,:,idx) = conv2(DF_F0(:,:,idx),Gauss2d, 'same');
    DF_F0(:,:,idx) = conv2(DF_F0(:,:,idx),Gauss2d, 'same');
end
if Rotate
    n90 = round(theta/90);
    for idx = 1:size(DF_F0,3)
        im = DF_F0(:,:,idx);
        im = rot90(im,n90);
        DF_F0(:,:,idx) = im;
    end
end

% saving
Fpass1 = round(Fpass1*100);
Fpass2 = round(Fpass2);
mosFolders = dir(meFolder);
msFolders = cleanFolderList(mosFolders,badTrials);
for jj = 1:length(msFolders)
    trialFolder = makeName(msFolders(jj).name,meFolder);
    f1 = fstart(jj);
    f2 = fend(jj);
    ImgSeq = DF_F0(:,:,f1:f2);
    ImgSeq(isnan(ImgSeq)) = 0;
    ImgSeq(isinf(ImgSeq)) = 0;  
    fileName = sprintf('%s\\ImgSeq_BPF_0%dHz_%dHz.mat',trialFolder, Fpass1, Fpass2);
    display('Saving file ... ');
    save(fileName,'ImgSeq','-v7.3');
%     fileName = sprintf('%s\\ImgSeq_BPF_0%dHz_%dHz.raw',trialFolder, Fpass1, Fpass2);
%     imwriteallraw(fileName,ImgSeq,[])
end

n = 0;

