function HPF_on_folder_Evok(meFolder,fileName, varName)

try
    badTrials = load(makeName('badTrials.txt',meFolder));
catch
    badTrials = [];
end
osFolders = dir(meFolder);
sFolders = cleanFolderList(osFolders,badTrials);
ImgSeq_all = [];
fstart = zeros(1,length(sFolders));
fend = zeros(1,length(sFolders));
for jj = 1:length(sFolders)
    trialFolder = makeName(sFolders(jj).name,meFolder);
    ImgSeq = load(makeName(fileName,trialFolder));
    cmdText = sprintf('ImgSeq = ImgSeq.%s;',varName);
    eval(cmdText);
    if ~isreal(ImgSeq)
        ImgSeq = abs(ImgSeq);
    end
    ImgSeq(isnan(ImgSeq)) = 0;
    ImgSeq(isinf(ImgSeq)) = 0;
    ImgSeq_all = cat(3,ImgSeq_all,ImgSeq);    
    if jj == 1
        fstart(jj) = 1;
    else
        fstart(jj) = fend(jj-1)+1;
    end
    fend(jj) = size(ImgSeq_all,3);
end

[dim1,dim2,nFrames] = size(ImgSeq_all);
ImgSeq_all = reshape(ImgSeq_all,dim1*dim2,nFrames);
nPix = dim1*dim2;

disp('Running parfor ... plz wait');
Fs = 150;
Hd = HPF_LS_05Hz(Fs);
Hd = Hd.Numerator;  
parfor p = 1:nPix
    sig = ImgSeq_all(p,:);
    sigHPF = filtfilt(Hd,1,sig);
    ImgSeq_all(p,:) = sigHPF;
end
disp('Done !!!');
ImgSeq_all = reshape(ImgSeq_all,dim1,dim2,nFrames);

for jj = 1:length(sFolders)
    trialFolder = makeName(sFolders(jj).name,meFolder);
    f1 = fstart(jj);
    f2 = fend(jj);
    ImgSeq = ImgSeq_all(:,:,f1:f2);
    ImgSeq(isnan(ImgSeq)) = 0;
    ImgSeq(isinf(ImgSeq)) = 0;    
    fileName = sprintf('%s\\ImgSeq_HPF.mat',trialFolder);
    display('Saving file ... ');
    save(fileName,'ImgSeq','-v7.3');
end

n = 0;

