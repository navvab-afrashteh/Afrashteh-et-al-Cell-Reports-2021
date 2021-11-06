function PreProcessEV_HiLo(file,filepath,varargin)
trial = file(3:end);

savePath = filepath;
if nargin > 2
    savePath = varargin{1};
end

Rotate = 0;
theta = 0;
Glut = 0;
AC = 0;
Tw = 0.3;
Ts = 0.1;
if nargin > 3
    Rotate = varargin{2};
end
if Rotate && nargin > 4
    theta = varargin{3};
end

if nargin > 5
    Glut = varargin{4};
end
if nargin > 6
    AC = varargin{5};
end

Fs = 150; % sampling rate
if AC
   % get Fs from file name
   IndexHz = max(strfind(lower(filepath),'hz'));
   Fs = filepath(1:IndexHz-1);
   IndexU = max(strfind(Fs,'_'));
   Fs = Fs(IndexU+1:end);
   Fs = str2double(Fs);
end
if Glut
    Tw = 2;
    Ts = 0.5;
end
% 2D Gaussian
param.GaussWindow = 5;
param.Sigma       = 1;
Center  = fix([param.GaussWindow/2,param.GaussWindow/2])+1;
[R,C]   = ndgrid(1:param.GaussWindow, 1:param.GaussWindow);
Gauss2d = gauss2dC(R,C,param.Sigma,Center);
% Load ImgSeq_raw
filetype = '.tif';
filenameImgSeq = [filepath, file, filetype];
ImgSeq_raw = imreadalltiff(filenameImgSeq); ImgSeq_raw = double(ImgSeq_raw);
if ~AC && ~Glut
    ImgSeq_raw(:,:,32) = (ImgSeq_raw(:,:,31)+ImgSeq_raw(:,:,33))/2;
end
[width,height,nFrames] = size(ImgSeq_raw);
% Load no trials and use it to calc F0
if ~AC && ~Glut
    filetype = '.tif';
    filenameNo = [filepath, 'no',trial, filetype];
    F0 = imreadalltiff(filenameNo); F0 = double(F0);
else
    F0 = ImgSeq_raw;
end
for r = 1:width
    parfor c = 1:height
        sig = F0(r,c,:);sig =sig(:);
        yhat = sig - locdetrend(sig,Fs,[Tw,Ts]);
        F0(r,c,:) = yhat;
    end
end
% Calc DF/F0
ImgSeq = 100*(ImgSeq_raw-F0)./F0; ImgSeq(isnan(ImgSeq))=0; ImgSeq(isinf(ImgSeq))=0;
% LPF
ImgSeq = reshape(ImgSeq,width*height,nFrames);
if Glut
    N = 11; sigma = 3; Gauss1d = gausswin(N,N/2/sigma);
    N1 = floor(N*Fs/100);
    t = 0:N-1; t = t/100;
    t1 = 0:N1-1; t1 = t1/Fs;
    Gauss1d = interp1(t,Gauss1d,t1);
    Gauss1d = Gauss1d/sum(Gauss1d);
    parfor p = 1:width*height
        sig = ImgSeq(p,:); sig = sig(:);
        ImgSeq(p,:) = filtfilt(Gauss1d,1,sig);
    end
elseif AC
    Gauss1d = LPF_Gauss_25Hz(Fs);
    parfor p = 1:width*height
        sig = ImgSeq(p,:); sig = sig(:);
        ImgSeq(p,:) = filtfilt(Gauss1d,1,sig);
    end
else
    Gauss1d = LPF_Gauss_25Hz(Fs); Lhalf = round((length(Gauss1d)-1)/2);
    parfor p = 1:width*height
        sig = ImgSeq(p,:); sig = sig(:);
        ImgSeq(p,:) = conv(sig,Gauss1d,'same');
    end
end
ImgSeq = reshape(ImgSeq,width,height,nFrames);

% Gaussian spatial filtering
for idx = 1:nFrames
    ImgSeq(:,:,idx) = conv2(ImgSeq(:,:,idx),Gauss2d, 'same');
    ImgSeq(:,:,idx) = conv2(ImgSeq(:,:,idx),Gauss2d, 'same');
end
% Trimm the filtering edge effect
if ~AC && ~Glut
    ImgSeq = ImgSeq(:,:,Lhalf+1:nFrames-Lhalf);
end
if Rotate
    n90 = round(theta/90);
    for idx = 1:size(ImgSeq,3)
        im = ImgSeq(:,:,idx);
        im = rot90(im,n90);
        ImgSeq(:,:,idx) = im;
    end
end

% save the results
% imwriteallraw([savePath, file, '_DF_F0_LPF_25Hz_G1_G1.raw'],ImgSeq,'*float32');
save([savePath, 'ImgSeq.mat'],'ImgSeq','-v7.3');
end
