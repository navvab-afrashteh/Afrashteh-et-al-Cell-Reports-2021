function ImgSeq = imreadallraw(filename,x,y,nFrames,precision,type)

% '*uint8' 8 bit imaging, raw data, behavioral camera
% '*uint16' 16 bit imaging, raw data, VSD camera
% '*float32' 32 bit, filtered data, VSD camera
if ~exist('nFrames','var')
    nFrames = [];
end
if ~exist('precision','var')
    precision = [];
end
if ~exist('type','var')
    type = 'slow';
end
if isempty(nFrames)
    nFrames = 1000000;
end
if isempty(precision)
    precision = '*float32';
end
% switch between 'slow' and 'fast' (using parallel loop) loading
switch type
    case 'slow'
        fid0 = fopen(filename, 'r', 'b');
        ImgSeq = fread(fid0,[x*y nFrames],precision);
        fclose(fid0);
        nFrames = size(ImgSeq,2);
        
    case 'fast'
        switch precision
            case '*uint8'
                unitSize = 1;
            case '*uint16'
                unitSize = 2;
            case '*float32'
                unitSize = 4;
            otherwise
                throw('precision not defined')
        end
        numThread = 10; % number of parallel threads
        ImgSeq = {};
        batchsize = nFrames / numThread;
        parfor i = 1:numThread
            fid0 = fopen(filename, 'r', 'b');
            fseek(fid0, unitSize*x*y*(i-1)*batchsize, -1);
            ImgSeq{i} = fread(fid0, [x*y batchsize], precision);
            ImgSeq{i} = single(ImgSeq{i});
            fclose(fid0);
        end
        ImgSeq = cat(2, ImgSeq{:});
end
ImgSeq = reshape(ImgSeq,x,y,nFrames);
ImgSeq = permute(ImgSeq, [2, 1, 3]);
