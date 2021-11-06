function PreProcessSpon_BPF(filepath, file, varargin)

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
   Fs = 100; 
end
Fstop1 = 0.25; Fpass1 = 0.5; Fpass2 = 6; Fstop2 = 6.5; dur = 5;
if nargin > 5
    filtOptions = varargin{4};
    Fstop1 = filtOptions.Fstop1; Fpass1 = filtOptions.Fpass1; Fpass2 = filtOptions.Fpass2;
    Fstop2 = filtOptions.Fstop2; dur = filtOptions.dur;
end
% 2D Gaussian
param.GaussWindow = 5;
param.Sigma       = 1;
Center  = fix([param.GaussWindow/2,param.GaussWindow/2])+1;
[R,C]   = ndgrid(1:param.GaussWindow, 1:param.GaussWindow);
Gauss2d = gauss2dC(R,C,param.Sigma,Center);
% Load ImgSeq
filetype = '.tif';
filename = [filepath, file, filetype];
ImgSeq = imreadalltiff(filename); ImgSeq = double(ImgSeq); 
[~,~,nFrames] = size(ImgSeq);

% BPF
if Glut
    meanImg = mean(ImgSeq,3);
    Hd = BPF_LS(Fs,Fstop1,Fpass1,Fpass2,Fstop2,dur);
    gBPF = Hd.Numerator;
    DF_F0 = tempFilt_imgStack(ImgSeq,gBPF);
    DF_F0 = 100*DF_F0./repmat(meanImg,1,1,nFrames);
end

% Gaussian spatial filtering
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

% save the results
Fpass1 = round(Fpass1*100);
Fpass2 = round(Fpass2);
fileName = sprintf('_DF_F0_BPF_0%dHz_%dHz_G1_G1.mat', Fpass1, Fpass2);
save([filepath, file, fileName],'DF_F0','-v7.3');

% fileName = sprintf('_DF_F0_BPF_0%dHz_%dHz_G1_G1.raw', Fpass1, Fpass2);
% imwriteallraw([filepath, file, fileName],DF_F0,'*float32');

