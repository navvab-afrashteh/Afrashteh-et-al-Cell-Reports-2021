function [ImgSeq,nFrames] = imreadalltiff(filename, varargin)
% IMREADALLTIFF Utility method to read all frames from a TIFF file.
% imreadalltiff(filename): all frames with fast loading
% imreadalltiff(filename, nFrames): loads first nFrames frames with fast loading
% imreadalltiff(filename, 'info'): gives the file info only. no file loading.
% imreadalltiff(filename, type). type can be 'slow' or 'fast' (default). 
% 'slow' uses feval and 'fast' uses fread. It loads all frames.
% You can also use imreadalltiff(filename, nFrames, type) or 
% imreadalltiff(filename, type, nFrames)

infoOnly = 0; % if 0 only returns the file info
if ~isempty(varargin)
    try
        if contains([varargin{:}],'info')
            infoOnly=1;
        end
    end
end
type = 'fast'; % type of loading the files. 
for n=1:nargin-1
    if isnumeric(varargin{n})
        nFrames = varargin{n};
    elseif strcmpi(varargin{n},'slow') || strcmpi(varargin{n},'fast')
        type = varargin{n};
    end
end
% file info
info = imfinfo(filename);
if infoOnly
    ImgSeq = info;
    return;
end
if ~exist('nFrames','var')
    nFrames = numel(info);
end
% load the data
switch type
    case 'fast'
        bitDepth = info(1).BitDepth;
        unitSize = bitDepth/8;
        if bitDepth==8
            precision = '*uint8';
            cl = 'uint16';
        elseif bitDepth==12
            precision='ubit12=>uint16'; % convert 12bit to 16bit on reading
            cl = 'uint16';
        elseif bitDepth==16
            precision = '*uint16';
            cl = 'uint16';
        elseif bitDepth==32
            precision = '*float32';
            cl = 'single';
        else
            throw('precision not defined')
        end
        machinefmt='b';
        % file info
        x = info(1).Width;
        y = info(1).Height;
        d = info(2).Offset-info(1).Offset; % these are header bytes
        d = d - unitSize*x*y;
        ImgSeq = zeros(y*x,nFrames,cl);
        fid0 = fopen(filename, 'r', machinefmt);
        for f=1:nFrames
            fseek(fid0, info(f).Offset+d, 'bof');
            ImgSeq(:,f)=fread(fid0,x*y,precision,machinefmt);
        end
        fclose(fid0);
        ImgSeq=permute(reshape(ImgSeq,x,y,nFrames),[2 1 3]);
        
    case 'slow'
        % read the first frame
        tf = imformats('tif');
        I = feval(tf.read, filename, 1);
        % preallocate output
        ImgSeq = zeros(size(I,1), size(I,2), nFrames, class(I));
        % copy in first image
        ImgSeq(:,:,1) = I;
        % try to read all frames
        frame = 1;
        try
            while frame < nFrames
                frame = frame + 1;
                ImgSeq(:,:,frame) = feval(tf.read, filename, frame);
            end
        catch
            frame = frame - 1;
            disp(strcat('Frames read: ', int2str(frame)));
            % optional line to trim off any extra frames read
            % comment out if not needed
            ImgSeq = ImgSeq(:,:,1:frame);
        end
end

