function output = getDistributions_simpleAmps_artur_2_glut (d,anum,roi,ss,amplitude_threshold,varargin)
multLevels = 1;
if nargin > 5
   multLevels = varargin{1}; 
end
pfuvall = [];
onsets = [];
pfs = [];
for kk = 1:length(anum)
    an = anum(kk);
    mask = d.masks{an}{1};
    for ii = 1:size(d.d,2)
        if isempty(d.d{an,ii})
            continue;
        end
        tdf = [];        tuv = [];        tang = [];        tangdur = [];
        trT = [];        trspT = [];        tuvAreaUC = [];
        tdfd = [];
        tfd = [];
        pre_frames = 25;
        post_frames = 25;
        jj = 1;
        jjc = 0;
        jjcd = 0;
        while jj < length(d.pEvokedFolderTrials{an,ii})    
            thisdf = d.d{an,ii}{jj}{roi};
            whole_trace = thisdf;
            thisuv = d.d1{an,ii}{jj}{roi}; 
            whole_trace_uv = thisuv;
            idx1 = 27;
            idx = 27:27+10;
%             idx = 27:27+30;
            baseline = mean(thisdf(12:26));
            stdbaseline = std(thisdf(12:26));
            thisuv1ang = thisuv(27:end);
            respThreshold = baseline + 2 *stdbaseline;
            rinds = find(thisdf>respThreshold);
            responseFrame = rinds(find(rinds>27,1,'first'));

            thisdf = thisdf(idx); thisuv = thisuv(idx);
            pf = find(thisdf == max(thisdf));
            pfuv = find(abs(thisuv) == max(abs(thisuv)));
            pfuvall = [pfuvall pfuv];
            
            motif_frames_N = ((pf+26)-pre_frames):((pf+26)+post_frames);
            motif_onset = 27;
            if motif_onset < motif_frames_N(1)
                motif_onset = motif_frames_N(1);
            end
            if isempty(motif_onset)
                onsets(ii,jj) = NaN;
                pfs(ii,jj) = NaN;
                jj = jj + 1;
                continue;
            end
            onsets(ii,jj) = motif_onset;
            pfs(ii,jj) = pf+idx1-1;
            try
                baseline_mean = mean(whole_trace((motif_onset-15):motif_onset));
                baseline_std = std(whole_trace((motif_onset-15):motif_onset));
                dfval = thisdf(pf)-baseline_mean;
                uvval = thisuv(pfuv);
                pfang = pfuv; % max of Speed
                rT = pf*6.67;
                fper = 0;
                fpost = 2;
                uvdur = thisuv1ang(pfang-fper:pfang+fpost);
                uvdur = uvdur./abs(uvdur);
                uvdur = mean(uvdur);
                uvdur =  uvdur./abs(uvdur);
                
                % flow directionality as a measure of the similarity betwwn
                % velocity vector at max speed frame and 6 frames before and after
                thisuv = whole_trace_uv(27:27+30);
                m = min(pfang-1,6);
                fd = [thisuv(pfang); thisuv(pfang-m:pfang-1); thisuv(pfang+1:pfang+16)];
                NN = length(fd);
                fd1 = 0;
                u = fd(1); % velocity vector at max speed frame
                for nn = 1:NN-1
                    v = fd(nn+1);
                    Theta = angleUV(u,v); % find angle between the two vector
                    fd1 = fd1 + cosd(Theta)/(NN-1);
                end
                fd = fd1;
                thisuv = whole_trace_uv(idx);
                
%                 trialFolder = makeName(trials{jj},makeName(names{ss}{ii},makeName('pEvoked',makeName(d.animalNumber{an},getMainDataFolderGlut))));
                trialFolder = makeName(d.pEvokedFolderTrials{an,ii}{jj},getMainDataFolderGlut);
                fileName = makeName('ImgSeq.mat',trialFolder);
                maguv = abs(thisuv);
                uvAreaUC = trapz(maguv);                
            catch
                jj = jj + 1;
                continue;
            end
            mean_post_onset = mean(whole_trace(motif_onset:(pf+26+15)));
            if mean_post_onset < (baseline_mean + amplitude_threshold*baseline_std)
                jj = jj + 1;
                jjcd = jjcd+1;
                allThisdf_discarded{ii}(jjcd,:) = whole_trace((motif_onset-15):(motif_onset+25))-baseline_mean;
                tdfd = [tdfd dfval];
                continue;
            end
            tfd = [tfd fd];
            tdf = [tdf dfval];
            tuv = [tuv 10*abs(uvval)];
            tang = [tang uvval];
            trT = [trT rT];
            tuvAreaUC = [tuvAreaUC uvAreaUC];
            trspT = [trspT (responseFrame-27)*6.67];
            tangdur = [tangdur uvdur];
            jjc = jjc + 1;
            preNFrames = 0.1*150;
            postNFrames = 0.4*150;
            allThisdf{ii}(jjc,:) = whole_trace((motif_onset-preNFrames):(motif_onset+postNFrames))-baseline_mean;
            allThisuv{ii}(jjc,:) = 10*abs(whole_trace_uv((motif_onset-preNFrames):(motif_onset+postNFrames)));
            allWholeTrace{ii}(jjc,:) = whole_trace;
            allWholeTraceuv{ii}(jjc,:) = 10*abs(whole_trace_uv);
            Trials(ii,jjc) = jj; % trials that made it to the end
            riseTime(jjc) = (pf+26-motif_onset)*6.67;
            display_img_seq = 0;
            if ~display_img_seq
                jj = jj + 1;
            end
        end
        allfd{kk,ii} = tfd;
        mdf(kk,ii) = mean(tdf);
        edf(kk,ii) = std(tdf);
        muv(kk,ii) = mean(tuv);
        euv(kk,ii) = std(tuv);
        mrT(kk,ii) = mean(trT);
        erT(kk,ii) = std(trT);
        mrspT(kk,ii) = mean(trspT);
        erspT(kk,ii) = std(trspT);
        mrspTST(kk,ii) = (getResponseFrame(d.animalNumber(an),d.pEvokedFolders{an,ii})-27)*6.67;
        muvAreaUC(kk,ii) = mean(tuvAreaUC);
        if an == 4
            mang(kk,ii) = angle(mean(tang)) + pi;
            eang(kk,ii) = angle(std(tang)/sqrt(length(tang))) + pi;
        else
            mang(kk,ii) = angle(mean(tang));
            eang(kk,ii) = angle(std(tang)/sqrt(length(tang)));
        end
        allang{kk,ii} = tangdur;
        alldf{kk,ii} = tdf;
        allRiseTime{kk,ii} = riseTime;
        alldfd{kk,ii} = tdfd;
        alluv{kk,ii} = tang;
        allrT{kk,ii} = trT;
        allrspT{kk,ii} = trspT;
        avgAllThisdf{ii} = mean(allThisdf{ii});
        avgAllThisuv{ii} = mean(allThisuv{ii});
        avgAmp(ii) = max(avgAllThisdf{ii}) - min(avgAllThisdf{ii});
        avgSpeed(ii) = max(avgAllThisuv{ii}) - min(avgAllThisuv{ii});
    end
end
output.Trials = Trials;
output.allRiseTime = allRiseTime;
output.allThisuv = allThisuv;
output.allThisdf = allThisdf;
output.allWholeTrace = allWholeTrace;
output.allWholeTraceuv = allWholeTraceuv;
output.avgAllThisuv = avgAllThisuv;
output.avgAllThisdf = avgAllThisdf;
output.avgAmp = avgAmp;
output.avgSpeed = avgSpeed;
if exist('allThisdf_discarded','var')
    output.allThisdf_discarded = allThisdf_discarded;
else
    output.allThisdf_discarded = [];
end
output.onsets = onsets;
output.pfs = pfs;
output.mdf = mdf;
output.allfd = allfd;
output.edf = edf;

output.muv = muv;
output.euv = euv;

output.mang = mang;
output.eang = eang;
output.allang = allang;

output.mrT = mrT;
output.erT = erT;

output.mrspT = mrspT;
output.erspT = erspT;
output.mrspTST = mrspTST;

output.alldf = alldf;
output.alldfd = alldfd;
output.alluv = alluv;
for ii = 1:length(alluv)
    output.absalluv{ii} = 10*abs(alluv{ii});
end
output.allrT = allrT;
output.allrspT = allrspT;

output.muvAreaUC = muvAreaUC;


%% Average files
dataVis = 0;
for kk = 1:length(anum)
    an = anum(kk);
    for ii = 1:size(d.davg,2)
        if isempty(d.davg{an,ii})
            continue;
        end
        thisdf = d.davg{an,ii}{roi};
        uthisdf = interp(thisdf,10);
        thisuv = d.d1avg{an,ii}{roi};
        whole_trace = thisdf;
        whole_trace_uv = thisuv;
        if ss == 3
            idx = 31:31+25;
            baseline = mean(thisdf(12:26));
            thisuv1ang = thisuv(31:end);
        else
            idx = 27:27+10;
            baseline = mean(uthisdf(120:269));
            thisuv1ang = thisuv(27:end);
            stdbaseline = std(uthisdf(120:269));
            respThreshold = baseline + 2 *stdbaseline;
            rinds = find(uthisdf(270:end)>respThreshold);
            responseFrame = rinds(find(rinds>=1,1,'first'));
            n = 0;
        end
        thisdf = thisdf(idx); thisuv = thisuv(idx);
        pf = find(thisdf == max(thisdf));
        pfuv = find(abs(thisuv) == max(abs(thisuv)));

%         motif_frames_N = ((pf+26)-pre_frames):((pf+26)+post_frames);
        pf1 = find(whole_trace == max(whole_trace));
        motif_frames_N = 1:100;
        thisMotif = whole_trace(motif_frames_N);
        baseline_thisMotif = findBleachingTrend(thisMotif);
        intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
%         motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < (pf+26),1,'last')));
        motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < (pf1),1,'last')));
%         figure(100);clf;plot(thisMotif);hold on;plot(baseline_thisMotif)
        motif_onset = 27;
        AvgallThisdf{kk}{ii} = whole_trace;
        AvgallThisuv{kk}{ii} = whole_trace_uv;
        preNFrames = 0.1*150;
        postNFrames = 0.4*150;
        AvgallThisdf{kk}{ii} = whole_trace((motif_onset-preNFrames):(motif_onset+postNFrames))-mean(whole_trace((motif_onset-preNFrames):motif_onset));
        AvgallThisuv{kk}{ii} = 10*abs(whole_trace_uv((motif_onset-preNFrames):(motif_onset+postNFrames)));%-mean(10*abs(whole_trace_uv((motif_onset-preNFrames):(motif_onset))));
        if dataVis
            figure(200);clf;
            subplot 211;
            plot(AvgallThisdf{kk}{ii});
            subplot 212;
            plot(AvgallThisuv{kk}{ii});
            if isempty(d.suffix)
                ImgSeq = load(makeName(sprintf('ImgSeq.mat'),makeName(d.pEvokedFolders{an,ii},getMainDataFolderGlut)));
            else
                ImgSeq = load(makeName(sprintf('ImgSeq_%s.mat',d.suffix),makeName(d.pEvokedFolders{an,ii},getMainDataFolderGlut)));
            end
            playdf(ImgSeq.ImgSeq,d.masks{an}{1},d.masks{roi}{2},'figureNumber',300,'title',d.pEvokedFolders{an,ii});
%             winopen(makeName(d.pEvokedFolders{1,ii},getMainDataFolderGlut))
%             playdf(ImgSeq.ImgSeq,d.masks{1}{1},d.masks{roi}{2},'figureNumber',300)
            n = 0;
        end
        if isempty(motif_onset)
            continue;
        end
%         if motif_onset <= 15
%             continue;
%         end
        baseline_mean = mean(whole_trace(1:motif_onset));
        baseline_std = std(whole_trace(1:motif_onset));
        
        try
            dfval = thisdf(pf)-baseline_mean;
            uvval = thisuv(pfuv);
            pfang = pfuv; % max of Speed
            fper = 0;
            fpost = 2;
            uvdur = thisuv1ang(pfang-fper:pfang+fpost);
            uvdur = uvdur./abs(uvdur);
            uvdur = mean(uvdur);
            uvdur =  uvdur./abs(uvdur);
            
            ttp = pf * 6.67;
            rlat = (responseFrame) * 0.667;
        catch
            continue;
        end
        Avgallang{kk,ii} = uvdur;
        Avgalldf{kk,ii} = dfval;
        Avgalluv{kk,ii} = max(AvgallThisuv{kk}{ii}-min(AvgallThisuv{kk}{ii}));
%         Avgalluv{kk,ii} = uvval;
        Avgallttp(kk,ii) = ttp;
        avg_onsets(kk,ii) = motif_onset;
        try
            Avgallrlat(kk,ii) = rlat;
        catch
        end
%         AvgallThisdf{kk}{ii} = whole_trace;
    end
end
output.Avgallang = Avgallang;
output.Avgalldf = Avgalldf;
output.Avgalluv = Avgalluv;
output.Avgallttp = Avgallttp;
output.avg_onsets = avg_onsets;
output.AvgallThisdf = AvgallThisdf;
output.AvgallThisuv = AvgallThisuv;
try
    output.Avgallrlat = Avgallrlat;
catch
end

output.avg = getValueFromAverageFiles(d,anum,roi,ss,amplitude_threshold,150,10,varargin);

function avg = getValueFromAverageFiles(d,anum,roi,ss,amplitude_threshold,Fs,pf2mms,varargin)
mainFolder = getMainDataFolderGlut;
for ii = 1:3
    pEvokedFolder = d.pEvokedFolders{anum,ii};
    ISF = makeName(sprintf('ImgSeq_BPF_050Hz_6Hz.mat',anum),pEvokedFolder);
    ISF = makeName(ISF,mainFolder);
    uvISF = makeName(sprintf('uv_ImgSeq_BPF_050Hz_6Hz.mat'),pEvokedFolder);
    uvISF = makeName(uvISF,mainFolder);
    temp = load(ISF);
    imgseq = temp.ImgSeq;
    temp = load(uvISF);
    uv = temp.uv;
    mask = double(d.masks{anum}{roi});
    vals = getMaskedValues(imgseq,mask);
    avg.df = mean(vals,2);
    avg.df = d.davg{anum,ii}{roi}; %%%
    vals = getMaskedValues(uv,mask);
    avg.uv = mean(vals,2);
    avg.uvabs = mean(abs(vals),2);
    motif_onset = 27;%round(0.3*Fs);
    avg.Fs = Fs;
    avg.Ts = 1000/Fs;
    thisdf = avg.df;
    dur = round(Fs*0.1);
    baseline_mean = mean(thisdf((motif_onset-dur):motif_onset));
    preDur = 0.1; postDur = 0.4;
    postF = round(postDur*Fs);
    preF = round(preDur*Fs);

    if (motif_onset-preF) <= 0
        preF = motif_onset - 1;
    end
    [~,pf] = max(thisdf(motif_onset:(motif_onset+round(0.1*Fs))));
    pf = pf + motif_onset;
    avg.df_chopped{ii} = avg.df((motif_onset-preF):(motif_onset+postF))-baseline_mean;
    avg.uv_chopped{ii} = pf2mms*abs(avg.uv((motif_onset-preF):(motif_onset+postF)));
    avg.uv_chopped{ii} = avg.uv_chopped{ii} - min(avg.uv_chopped{ii}(1:round(0.1*Fs)));
    avg.uvabs_chopped{ii} = pf2mms*abs(avg.uvabs((motif_onset-preF):(motif_onset+postF)));
    avg.riseTime(ii) = (pf-motif_onset)*avg.Ts;
    avg.ts = linspace(0,(preDur+postDur),length(avg.df_chopped));

    avg.Amp(ii) = max(avg.df_chopped{ii});
%     avg.Speed(ii) = max(avg.uv_chopped{ii}) - mean(avg.uv_chopped{ii}(1:round(0.1*Fs)));
    avg.Speed(ii) = max(avg.uv_chopped{ii});
    avg.Speedabs(ii) = max(avg.uvabs_chopped{ii});
end





function arr = getArea(fileName,pf,mask)
temp = load(fileName);
temp.ImgSeq = applyMask(temp.ImgSeq,mask,'same');
minI = min(temp.ImgSeq(:)); 
maxI = max(temp.ImgSeq(:)); 
for ii = 25:40
    thisFrame = temp.ImgSeq(:,:,ii);
    vals = getMaskedValues(thisFrame,mask);
    figure(2);clf;imagesc(thisFrame,[minI maxI]);hold on;
    BW = im2bw(thisFrame,mean(vals)+1.5*std(vals));
    s = regionprops(BW,'centroid');
    centroids = cat(1, s.Centroid);
    contour(BW);
    plot(centroids(:,1),centroids(:,2), 'b*')
    pause(0.1);
end
n = 0;

function responseFrame = getResponseFrame(animalNumber,folderName)
temp = makeName(animalNumber{1},getMainDataFolderGlut);
temp = makeName('pEvoked',temp);
folderName = makeName(folderName,temp);
fileName = makeName('stimTemplate.mat',folderName);
try
    temp = load(fileName);
    responseFrame = temp.stimTemplate.frameNumbers(1);
catch
    responseFrame = 28;
end


