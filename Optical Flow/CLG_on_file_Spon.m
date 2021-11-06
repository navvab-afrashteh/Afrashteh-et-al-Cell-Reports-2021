function varargout = CLG_on_file_Spon(filePath,fileName, varName, masks, varargin)
[r, c] = find(masks.mask);
idxMask = sub2ind(size(masks.mask),r,c);
makeName(fileName,filePath)
temp = load(makeName(fileName,filePath));
cmdText = sprintf('ImgSeq = temp.%s;',varName);
eval(cmdText);
ImgSeq = applyMask(ImgSeq,double(masks.bigMask));
if ~isreal(ImgSeq)
    ImgSeq = abs(ImgSeq);
end
ImgSeq(isnan(ImgSeq)) = 0;
ImgSeq(isinf(ImgSeq)) = 0;

[dim1,dim2,totalFramesP] = size(ImgSeq);
mW = round(min(dim1,dim2)*0.5/2);
para = [0.04,0.5,mW,15,3,30];
% hWaitBar = waitbar(0,sprintf(' '));
% Performing CLG method
uv = zeros(size(ImgSeq,1), size(ImgSeq,2), totalFramesP-1);
display('Running parfor ... plz wait');
parfor idx = 1:(totalFramesP-1)
%         waitbar(idx/totalFramesP,hWaitBar,sprintf('Processing frames ... %d of %d with CLG',idx,totalFramesP));
    im1 = ImgSeq(:,:,idx);
    im2 = ImgSeq(:,:,idx+1);
%     [im1, im2] = normalize_two_consequtive_frames(im1,im2,idxMask);
    [u, v, ~] = Coarse2FineTwoFrames(im1,im2,para);
    uv(:,:,idx) = (u +1i*v);
end
display('Done !!!');
% close(hWaitBar);
if nargin == 4
    display('Saving file now ... plz wait');
    uvFile = makeName(sprintf('uv%s',fileName),filePath);
    save(uvFile,'uv','-v7.3');
    display('Done !!!');
else
    varargout{1} = uv;
end


