function motifs = refineMotifs_AC(motifs,dfMatch,uvMatch,varargin)
Fs = 150;
if nargin > 3
   Fs = varargin{1};
end
dur = round(Fs*0.1);

nFrames = length(dfMatch);
frames = motifs.frames;
frames = frames(frames>400);
frames = frames(frames<(nFrames-(3*dur)));
cc = motifs.allCCs;

% here is where we find the closest df peak to CC peak
motifs.framesdf = [];
if ~isempty(frames)
    for f = 1:length(frames)
        a = dfMatch(frames(f)-dur:min(frames(f)+dur,nFrames));
        [~,idx] = max(a);
%         c = cc(frames(f)-dur:min(frames(f)+dur,nFrames));
%         [~,idxc] = max(c);
%         figure(1000000);clf;plot(a);hold on;
%         plot(idx,a(idx),'*b');
%         plot(c,'r');
%         plot(idxc,c(idxc),'*r');
        
        frames(f) = frames(f) + idx - dur -1;
    end
    frames = frames(frames>400);
    motifs.framesdf = frames;
end

nFrames = length(uvMatch);
frames = motifs.framesdf;
% frames = frames(frames>400);
motifs.framesuv = [];
if ~isempty(frames)
    for f = 1:length(frames)
        a = uvMatch(frames(f)-dur:min(frames(f)+dur,nFrames));
        [~,idx] = max(a);
        frames(f) = frames(f) + idx - dur -1;
    end
    motifs.framesuv = frames;
end

frames = motifs.frames;
frames = frames(frames>400);
motifs.frames = frames;


%%%%%%%%%%%%%%%%%%% old%%%%%%%%%%%%%%%%%
% thePeaks = detectPeaks(dfMatch);
% frames = motifs.frames;
% frames = frames(frames>400);
% motifs.framesdf = [];
% if ~isempty(frames)
%     for f = 1:length(frames)
%         a = thePeaks(:,1) - frames(f)>=0;
%         [~,idx] = max(a);
%         [~,idx] = max(dfMatch(frames(f):thePeaks(idx,1)));
%         frames(f) = frames(f) + idx -1;
%     end
%     motifs.framesdf = frames;
% end
% 
% thePeaks = detectPeaks(uvMatch);
% frames = motifs.frames;
% frames = frames(frames>400);
% motifs.framesuv = [];
% if ~isempty(frames)
%     for f = 1:length(frames)
%         a = thePeaks(:,1) - frames(f)>=0;
%         [~,idx] = max(a);
%         [~,idx] = max(uvMatch(frames(f):thePeaks(idx,1)));
%         frames(f) = frames(f) + idx -1;
%     end
%     motifs.framesuv = frames;
% end

