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




