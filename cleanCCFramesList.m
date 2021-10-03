function frames = cleanCCFramesList (CCcoeff,th,varargin)
% frames = cleanCCFramesList(CCcoeff,threshold,CrossOrMax);
% Gets the output from 'MotifsCorr' function and a threshold value and find
% frames that correlation coefficient is greater than the threshold value.
% If sevelar nearby frames are greater than threshold then you can use 
% CrossOrMax input to specify if you like to choose the first frame that
% crosses the threshold or the frame with the maximum value.

frames = [];
CrossOrMax = 'Max'; % default
if nargin >=3
    CrossOrMax = varargin{1};
end
CCcoeff_th = CCcoeff>th;
[~,idx] = find(CCcoeff_th > 0);
% If you want to select the first point that crosses the threshold
if strcmp(CrossOrMax,'Cross')
    frameNumbers = idx;
    allNums = 1:length(frameNumbers);
    for ii = 1:30
        dFN = diff(frameNumbers);
        inds = find(dFN == ii);
        frameNumbers(inds+1) = [];
        allNums(inds+1) = [];
    end
    frames = frameNumbers;
    flz = frames <= 0;
    frames(flz) = [];
end
% If you want to select the local maximas
if strcmp(CrossOrMax,'Max')
    diff_idx = find(diff(idx)>1);
    if isempty(diff_idx)
        return
    end
    idx_new = zeros(1,length(diff_idx)+1);
    for p = 1:length(diff_idx)+1
        if p == 1
            temp = idx(1:diff_idx(p));
        elseif p == length(diff_idx)+1
            temp = idx(diff_idx(p-1)+1:end);
        else
            temp = idx(diff_idx(p-1)+1:diff_idx(p));
        end
        [~,maxcorr_temp] = max(CCcoeff(temp));
        idx_new(p) = temp(maxcorr_temp);
    end
    diff_idx_new = find(diff(idx_new)>30);
    if isempty(diff_idx_new)
        return
    end
    idx_new2 = zeros(1,length(diff_idx_new)+1);
    for p = 1:length(diff_idx_new)+1
        if p == 1
            temp = idx_new(1:diff_idx_new(p));
        elseif p == length(diff_idx_new)+1
            temp = idx_new(diff_idx_new(p-1)+1:end);
        else
            temp = idx_new(diff_idx_new(p-1)+1:diff_idx_new(p));
        end
        [~,maxcorr_temp] = max(CCcoeff(temp));
        idx_new2(p) = temp(maxcorr_temp);
    end
    frames = idx_new2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;plot(CCcoeff);
% hold on;stem(idx,CCcoeff(idx),'.')
% hold on;stem(idx_new,CCcoeff(idx_new),'ko')
% hold on;stem(idx_new2,CCcoeff(idx_new2),'g*')
% hold on;stem(frameNumbers,CCcoeff(frameNumbers),'ms')

% n=0;
