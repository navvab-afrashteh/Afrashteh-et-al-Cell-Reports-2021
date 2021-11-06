function CCcoeff = MotifsCorr(TemplateSeq,ImgSeq,Mask)
% gets template frames, spontaneous image sequence, and mask and calculates
% correlation between template and image sequence. Only inside mask pixels
% are used to calculate correlation.

TemplateSeq = getMaskedValues(TemplateSeq,Mask);
ImgSeq = getMaskedValues(ImgSeq,Mask);
nFrames = size(ImgSeq,1);
nFramesTemplate = size(TemplateSeq,1);
TemplateSeq = reshape(TemplateSeq,1,numel(TemplateSeq));
allThisFrame = cell(1,nFrames-nFramesTemplate+1);
for cnt = 1:(nFrames-nFramesTemplate+1)
      allThisFrame{cnt} = ImgSeq(cnt:(cnt+nFramesTemplate-1),:);
end
parfor cnt = 1:(nFrames-nFramesTemplate+1)
      thisFrame = allThisFrame{cnt};
      thisFrame = reshape(thisFrame,1,numel(thisFrame));
      CC = corrcoef(TemplateSeq,thisFrame);
      CCcoeff(cnt) = CC(1,2);
end