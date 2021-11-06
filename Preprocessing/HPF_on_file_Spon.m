function ImgSeq = HPF_on_file_Spon(filePath,fileName, varName, varargin)

temp = load(makeName(fileName,filePath));
cmdText = sprintf('ImgSeq = temp.%s;',varName);
eval(cmdText);
ImgSeq(isnan(ImgSeq)) = 0;
ImgSeq(isinf(ImgSeq)) = 0;

[dim1,dim2,nFrames] = size(ImgSeq);
ImgSeq = reshape(ImgSeq,dim1*dim2,nFrames);
nPix = dim1*dim2;

disp('Running parfor ... plz wait');
Hd = HPF_LS_05Hz;
Hd = Hd.Numerator;
parfor p = 1:nPix
    sig = ImgSeq(p,:);
    sigHPF = filtfilt(Hd,1,sig);
    ImgSeq(p,:) = sigHPF;
end
disp('Done !!!');
ImgSeq = reshape(ImgSeq,dim1,dim2,nFrames);

