function TrajInfo = calcAllCenterOfMassTrajectory(ImgSeq,mask,x0,y0,r)

% ImgSeq contains frames from the start to the end of trajectory. For
% example for evoked trials I input the onset of response to peak frame.
% After peak frame most of the time there is no clear movement in traveling
% waves. Waves just fade away.

% x0 and y0 are x and y location of expected onset of trajectory. I input
% center of roi for stim. For example center of FLS1.

% r is to make sure that in current frame we capture main activity that is
% normally closer to previous point in the trajectory. I set r = 45.

[dim1,dim2,~] = size(ImgSeq);
xc = x0; yc = y0;
ImgSeq = applyMask(ImgSeq,mask,0);
idxS = 1; idxE = size(ImgSeq,3);
idxs = idxS:idxE;
for ii = idxs
    ROI = CircROI(xc,yc,r,dim1,dim2).*mask;
    ii1 = max(ii-1,1);
    im = ImgSeq(:,:,ii1);
    maskVals = getMaskedValues(im,ROI);
    threshold = prctile(maskVals,40);
    if isnan(threshold); threshold = 0; end
    if isinf(threshold); threshold = 1; end
    if threshold > 1; threshold = 1; end
    threshold = max(0,threshold);
    [BlargestP, xc, yc] = getlargestROI(ImgSeq(:,:,ii1).*ROI,threshold);
    ii2 = ii;
    im = ImgSeq(:,:,ii2);
    maskVals = getMaskedValues(im,ROI);
    threshold_s = prctile(maskVals,25);
    if isnan(threshold_s); threshold_s = 0; end
    if isinf(threshold_s); threshold_s = 1; end
    if threshold_s > 1; threshold_s = 1; end
    threshold_s = max(0,threshold_s);
    [Blargest_s, x1, y1] = getlargestROI(ImgSeq(:,:,ii).*ROI,threshold_s);
    nDir = 4*5;
    flowdir = getFlowDir(Blargest_s{1}, x1, y1,BlargestP{1}, xc, yc,nDir);
    ii3 = ii;
    im = ImgSeq(:,:,ii3);
    maskVals = getMaskedValues(im,ROI);
    threshold_s = prctile(maskVals,50);
    if isnan(threshold_s); threshold_s = 0; end
    if isinf(threshold_s); threshold_s = 1; end
    if threshold_s > 1; threshold_s = 1; end
    threshold_s = max(0,threshold_s);
    [Blargest, xc, yc] = getlargestROI(ImgSeq(:,:,ii).*ROI,threshold_s);
    [xe,ye] = FlowPoint(Blargest{1},x0,y0,flowdir);
    x0 = xc; y0 = yc;
    if ii == 1
        TrajInfo.x(ii) = xc;
        TrajInfo.y(ii) = yc;
    end
    if isempty(xe) || isempty(ye)
        xe = nan; ye = nan;
    end
    TrajInfo.x(ii+1) = xe;
    TrajInfo.y(ii+1) = ye;
end
idx = find(~isnan(TrajInfo.x));
if length(idx)>1
    TrajInfo.x = interp1(idx,TrajInfo.x(idx),1:max(idxs)+1);
end
idx = find(~isnan(TrajInfo.y));
if length(idx)>1
    TrajInfo.y = interp1(idx,TrajInfo.y(idx),1:max(idxs)+1);
end
end

function flowdir = getFlowDir(B1,x1,y1,B2,x2,y2,nDir)
alltheta = linspace(0,2*pi,nDir);
alldir = cos(alltheta)+1i*sin(alltheta);
xe1 = zeros(nDir-1,1); ye1 = zeros(nDir-1,1); idx1 = zeros(nDir-1,1);
xe2 = zeros(nDir-1,1); ye2 = zeros(nDir-1,1); idx2 = zeros(nDir-1,1);
for k = 1:nDir-1
    [xe1(k),ye1(k),idx1(k)] = FlowPoint(B1,x1,y1,alldir(k));
    [xe2(k),ye2(k),idx2(k)] = FlowPoint(B2,x2,y2,alldir(k));
end
[idx1, I1] = sort(idx1);
[idx2, I2] = sort(idx2);

x1dir = zeros(nDir-1,1); y1dir = zeros(nDir-1,1); x2dir = zeros(nDir-1,1); y2dir = zeros(nDir-1,1);
for k = 1:nDir-2
    if ~isnan(idx1(k)) && ~isnan(idx1(k+1))
        x1dir(I1(k)) = mean(B1(idx1(k):idx1(k+1),2));
        y1dir(I1(k)) = mean(B1(idx1(k):idx1(k+1),1));
    end
    if ~isnan(idx2(k)) && ~isnan(idx2(k+1))
        x2dir(I2(k)) = mean(B2(idx2(k):idx2(k+1),2));
        y2dir(I2(k)) = mean(B2(idx2(k):idx2(k+1),1));
    end
end
for k = nDir-1
    if ~isnan(idx1(k)) && ~isnan(idx1(1))
        x1dir(I1(k)) = mean([B1(1:idx1(1),2); B1(idx1(k):end,2)]);
        y1dir(I1(k)) = mean([B1(1:idx1(1),1); B1(idx1(k):end,1)]);
    end
    if ~isnan(idx2(k)) && ~isnan(idx2(1))
        x2dir(I2(k)) = mean([B2(1:idx2(1),2); B2(idx2(k):end,2)]);
        y2dir(I2(k)) = mean([B2(1:idx2(1),1); B2(idx2(k):end,1)]);
    end
end
orig1 = x1 - 1i*y1;
vec1 = x1dir + 1i*y1dir - orig1;
orig2 = x2 - 1i*y2;
vec2 = x2dir + 1i*y2dir - orig2;
vec = vec1 - vec2;
flowdir = sum(vec);
flowdir = flowdir/abs(flowdir);

% plot(B1(:,2),B1(:,1),'w','linewidth',2);
% plot(x1,y1,'*w');
% plot(B2(:,2),B2(:,1),'k','linewidth',2);
% plot(x2,y2,'*k');
% plot(x1dir,y1dir,'r.')
% plot(x2dir,y2dir,'g.')
% quiver(x2,y2,real(flowdir),imag(flowdir),50)
% quiver(x2dir,y2dir,real(vec),imag(vec),5)
% n=0;
end


function [Blargest, xc, yc] = getlargestROI(im,threshold)
BW = im2bw(im,threshold);
CC = bwconncomp(BW,4);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idxM] = max(numPixels);
if isempty(idxM)
    centerOfMass = [];
    Blargest{1} = [];
    Blargest{2} = [];
else
    largestROI = zeros(size(BW));
    largestROI(CC.PixelIdxList{idxM}) = 1;
    largestROI = imfill(largestROI,'holes');
    Blargest = bwboundaries(largestROI,4);
    grayImage = im.*largestROI;
    centerOfMass = regionprops(largestROI, grayImage, 'WeightedCentroid');
end
if isempty(centerOfMass)
    centerOfMass = [nan,nan];
else
    centerOfMass = centerOfMass.WeightedCentroid;
end
xc = centerOfMass(1); yc = centerOfMass(2);
end


function [xe,ye,idx] = FlowPoint(B,xe,ye,avguv)
if isempty(B)
    idx = nan;
    xe = nan;
    ye = nan;
    return;
end
m = imag(avguv)/real(avguv);
y = m*(B(:,2)-xe)+ye;
b = B(:,1)-y;
b = b(1:end-1).*b(2:end);
idx = find(b<=0);
b = [B(idx,2), B(idx,1)];
b = b - repmat([xe,ye],size(b,1),1);
b = b * [real(avguv); imag(avguv)];
b = b > 0;
idx = idx(b);
B = B(idx,:);
d = B - repmat([ye,xe],size(B,1),1);
d = sqrt(sum(d.^2,2));
[~,idx1] = max(d);
ye = B(idx1,1);
xe = B(idx1,2);
idx = idx(idx1);
if isempty(idx)
    idx = nan;
    xe = nan;
    ye = nan;
end
end

function ROI = CircROI(xc,yc,r,dim1,dim2)

hf = figure; h = viscircles([xc,yc],r);
ch = get(h,'Children');
xv = get(ch(1),'XData');
yv = get(ch(1),'YData');
close(hf);

[X,Y] = meshgrid(1:dim2, 1:dim1);
[in,on] = inpolygon(X,Y,xv,yv);
ROI = or(in,on);
end

function data = applyMask (data,mask,varargin)
mask = double(mask);

p = inputParser;
addRequired(p,'data',@isnumeric);
addRequired(p,'mask',@isnumeric);
addOptional(p,'shiftToZero',0,@isnumeric);
addOptional(p,'shiftToZeroAndNormalize',0,@isnumeric);
addOptional(p,'maskValue',0,@isnumeric);
parse(p,data,mask,varargin{:});

shiftToZero = p.Results.shiftToZero;
shiftToZeroAndNormalize = p.Results.shiftToZeroAndNormalize;
maskValue = p.Results.maskValue;

for ii = 1:size(data,3)
    data(:,:,ii) = data(:,:,ii) .* mask;
end
mdata = min(data(:));

if maskValue ~=0
    Imin = maskValue;
else
    Imin = mdata;
end

imask = ~mask * Imin;
for ii = 1:size(data,3)
    data(:,:,ii) = data(:,:,ii) + imask;
end

if shiftToZero
    data = data - mdata;
end

if shiftToZeroAndNormalize
    data = data - mdata;
    data = data./max(data(:));
end
end

function values = getMaskedValues (data,mask)
[r, c] = find(mask==1);
maskI = sub2ind(size(mask),r,c);
values = zeros(size(data,3),length(maskI));
for ii = 1:size(data,3)
    temp = data(:,:,ii);
    values(ii,:) = transpose(temp(maskI));
end
end