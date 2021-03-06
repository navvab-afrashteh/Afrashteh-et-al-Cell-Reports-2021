function output = getDistributions_simpleAmps_artur_2_AC (d,anum,roi,ss,amplitude_threshold,varargin)
multLevels = 1;
dataType = 'AC';
if nargin > 5
   multLevels = varargin{1}; 
end
pfuvall = [];
onsets = [];
pfs = []
% amplitude_threshold = 0.1;
for kk = 1:length(anum)
    an = anum(kk);
    mask = d.masks{an}{1};
    [lists, names]= getEvokedListsDiffStimAmps_glut(d.animalNumber(an),multLevels,dataType);
    for ii = 1:size(d.d,2)
        if isempty(d.d{an,ii})
            continue;
        end
        file = names{4}{ii};
        % get Fs from file name
        IndexHz = max(strfind(lower(file),'hz'));
        Fs = file(1:IndexHz-1);
        IndexU = max(strfind(Fs,'_'));
        Fs = Fs(IndexU+1:end);
        Fs = str2double(Fs);
        idx1 = round(Fs*0.3)+1;
        idxR = idx1:idx1+round(Fs*0.1)-1; % response indeces
        idxB = idx1-round(Fs*0.1):idx1-1; % baseline indeces
        pf2mms = Fs/15; % pixel/frame to mm/sec
        Ts = 1000/Fs; % frame-rate time duration
            
        tdf = [];        tuv = [];        tang = [];        tangdur = [];
        trT = [];        trspT = [];        tuvAreaUC = [];
        tdfd = [];
%         for jj = 1:length(trialNums)
        pre_frames = round(Fs*0.1667);
        post_frames = round(Fs*0.1667);
        jj = 1;
        jjc = 0;
        jjcd = 0;
        while jj < length(d.pEvokedFolderTrials{an,ii})    
            thisdf = d.d{an,ii}{jj}{roi};
            whole_trace = thisdf;
            thisuv = d.d1{an,ii}{jj}{roi}; 
            thisuvabs = d.d1abs{an,ii}{jj}{roi}; 
            whole_trace_uv = thisuv;
            whole_trace_uvabs = thisuvabs;
%             figure(12);clf;plot(whole_trace);

            baseline = mean(thisdf(idxB));
            stdbaseline = std(thisdf(idxB));
            thisuv1ang = thisuv(idx1:end);
            respThreshold = baseline + 2 *stdbaseline;
            rinds = find(thisdf>respThreshold);
            responseFrame = rinds(find(rinds>idx1,1,'first'));

            thisdf = thisdf(idxR); thisuv = thisuv(idxR);
            pf = find(thisdf == max(thisdf));
            pfuv = find(abs(thisuv) == max(abs(thisuv)));
            pfuvall = [pfuvall pfuv];
            
            motif_frames_N = ((pf+idx1-1)-pre_frames):((pf+idx1-1)+post_frames);
            thisMotif = whole_trace(motif_frames_N);
            baseline_thisMotif = findBleachingTrend(thisMotif);
            intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
            motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < (pf+idx1-1),1,'last')))-1;
            motif_onset = round(0.3*Fs);
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
            
            fig_for_paper = 0;
            if fig_for_paper
                ff = makeFigureWindow__one_axes_only(4,[1 1 3.5 2],[0.15 0.25 0.8 0.6]);
                axes(ff.ha);
                peak_vals = detectPeaksN(thisMotif);
                plot(Ts*(motif_frames_N-pfs(ii,jj)),thisMotif,'b');hold on;
                plot(Ts*(motif_frames_N-pfs(ii,jj)),baseline_thisMotif,'m');
                plot((motif_frames_N(peak_vals(:,1))-pfs(ii,jj))*Ts,peak_vals(:,2),'k*');
                plot(Ts*(motif_onset-pfs(ii,jj)),thisMotif(motif_frames_N == motif_onset),'*r');
                xlim([-170 170]);
                ylim([min([thisMotif;baseline_thisMotif])-0.1 max([thisMotif;baseline_thisMotif])+0.1]);
                set(gca,'TickDir','out','FontSize',8,'FontWeight','bold');
                xlabel('Time (ms)');
                hy = ylabel('VSD Signal (\DeltaF/F_0 %)');
                thisCols = {'r','m','k','b'};
                x1 = -150; x2 = x1+25; y1 = (0.2:0.07:7); y1 = y1(1:4); y2 = y1;
                legendFontSize = 6-1;
                legs = {'Estimated onset','Fitted curve','Local minimas','Signal'};
                for ii = 1:4
                    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii});
                    text(x2+10,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize);
                end
                
                save2pdf(makeName('motif_onset_calculation.pdf',getpdffolder),ff.hf,600);
                n = 0;
            end
            
            display_img_seq = 0;
            if display_img_seq
                ImgSeq_file_name = makeName(makeName('ImgSeq.mat',d.pEvokedFolderTrials{anum,ii}{jj}),getMainDataFolderAC);
                msf = matfile(ImgSeq_file_name);
                thisImgSeq = msf.ImgSeq;
                motif_frames = thisImgSeq(:,:,motif_frames_N);

                thisROI = d.masks{an}{2};
                bo = bwboundaries(thisROI);
                xs = bo{1}(:,1);        ys = bo{1}(:,2);
                figure(11);clf;
                for mfn = 1:length(motif_frames_N)
                    this_motif_frame = applyMask(motif_frames(:,:,mfn),mask);
                    subplot 211
                    cla
                    imagesc(this_motif_frame);
                    colorbar;
                    hold on;
                    axis equal; axis off;
                    plot(ys,xs,'m','linewidth',0.5);
                    set(gca,'Ydir','reverse');
                    titleT = sprintf('Frame %d - %d',mfn,motif_frames_N(mfn));
                    title(titleT,'Interpreter','none');
%                     if length(SoSi.source) > 0
%                         if mfn == SoSi.source(3)
%                             plot(SoSi.source(1),SoSi.source(2),'*','color','r');
%                             pause;
%                         end
%                     end
                    subplot 212
                    cla
                    plot(motif_frames_N,thisMotif);hold on;
                    plot(motif_frames_N,baseline_thisMotif,'m');
%                     baseline_thisMotif = findBleachingTrend(thisMotif,5);
%                     plot(motif_frames_N,baseline_thisMotif,'k');
%                     baseline_thisMotif = findBleachingTrend(thisMotif,7);
%                     plot(motif_frames_N,baseline_thisMotif,'c');
%                     plot(motif_frames_N(3:end),DthisMotif);
%                     plot(motif_frames_N,ones(size(thisMotif))*dfThreshold);
                    titleT = sprintf('%s - %s - %s    %d',d.animalNumber{anum},d.pEvokedFolders{anum,ii},d.pEvokedFolderTrials{anum,ii}{jj},motif_onset-27);
                    title(titleT,'Interpreter','none');
                    plot(motif_onset,thisMotif(motif_frames_N==motif_onset),'*g');
                    ylims = ylim;
                    plot([motif_frames_N(mfn) motif_frames_N(mfn)],[ylims(1) ylims(2)],'r');
                    if motif_frames_N(mfn) == motif_onset
                        subplot 211
                        text(10,10,'Estimated Onset','Color','w');
                    end
%                     if motif_frames_N(mfn) > (pf+26-7) && motif_frames_N(mfn) < (pf+26+7)                    
                    if motif_frames_N(mfn) > (motif_onset-3) && motif_frames_N(mfn) < (motif_onset+3)
                        pause(0.01);
                    else
                        pause(0.01);
                    end
                end
                disp('press a key');
                keyVal = getkey;%pause;
                if keyVal == 27
                   the_error_creating
                end
                if keyVal == 30 % Up arrow
                    jj = jj + 1;
                    continue;
                end
                if keyVal == 31 % down arrow
                    if jj > 1
                        jj = jj - 1;
                    end
                    continue;
                end
                if keyVal == 113
                    break;
                end
            end
            
            
            
            
            try
                if roi == 2
                    Tpre = round(Fs*0.1);
                else
                    Tpre = round(Fs/30);
                end
                baseline_mean = mean(whole_trace((motif_onset-Tpre):motif_onset));
                baseline_std = std(whole_trace((motif_onset-Tpre):motif_onset));
                dfval = thisdf(pf)-baseline_mean;
%                 dfval = thisdf(pf)-mean(baseline);
                uvval = thisuv(pfuv);
                pfang = pfuv; % max of Speed
%                 pfang = 4;
                rT = pf*Ts;
                fper = 0;
                fpost = 2;
                uvdur = thisuv1ang(pfang-fper:pfang+fpost);
                uvdur = uvdur./abs(uvdur);
                uvdur = mean(uvdur);
                uvdur =  uvdur./abs(uvdur);
%                 trialFolder = makeName(trials{jj},makeName(names{ss}{ii},makeName('pEvoked',makeName(d.animalNumber{an},getMainDataFolderGlut))));
                trialFolder = makeName(d.pEvokedFolderTrials{an,ii}{jj},getMainDataFolderAC);
                fileName = makeName('ImgSeq.mat',trialFolder);
                maguv = abs(thisuv);
                uvAreaUC = trapz(maguv);
%                 areaval = getArea(fileName,pf+26,mask);
                
            catch
                jj = jj + 1;
                continue;
            end
            Tpost = round(Fs*0.1);
            mean_post_onset = mean(whole_trace(motif_onset:(pf+idx1-1+Tpost)));
            if mean_post_onset < (baseline_mean + amplitude_threshold*baseline_std)
                jj = jj + 1;
                jjcd = jjcd+1;
%                 allThisdf_discarded{ii}(jjcd,:) = whole_trace((motif_onset-15):(motif_onset+25))-baseline_mean;
                tdfd = [tdfd dfval];
                continue;
            end
            tdf = [tdf dfval];
            tuv = [tuv pf2mms*abs(uvval)];
            tang = [tang uvval];
            trT = [trT rT];
            tuvAreaUC = [tuvAreaUC uvAreaUC];
            trspT = [trspT (responseFrame-idx1)*Ts];
            tangdur = [tangdur uvdur];
            jjc = jjc + 1;
            if roi == 2
                allThisdf{ii}(jjc,:) = whole_trace((motif_onset-round(0.1*Fs)):(motif_onset+round(0.4*Fs)))-baseline_mean;
                allThisuv{ii}(jjc,:) = pf2mms*abs(whole_trace_uv((motif_onset-round(0.1*Fs)):(motif_onset+round(0.4*Fs))));
                allThisuvabs{ii}(jjc,:) = pf2mms*abs(whole_trace_uvabs((motif_onset-round(0.1*Fs)):(motif_onset+round(0.4*Fs))));
            else
%                 allThisdf{ii}(jjc,:) = whole_trace((motif_onset-round(Fs/30)):(motif_onset+round(0.4*Fs)))-baseline_mean;
%                 allThisuv{ii}(jjc,:) = pf2mms*abs(whole_trace_uv((motif_onset-round(Fs/30)):(motif_onset+round(0.4*Fs))));
            end
            riseTime(jjc) = (pf+idx1-1-motif_onset)*Ts;
            if ~display_img_seq
                jj = jj + 1;
            end
        end
        mdf(kk,ii) = mean(tdf);
        edf(kk,ii) = std(tdf);
        muv(kk,ii) = mean(tuv);
        euv(kk,ii) = std(tuv);
        mrT(kk,ii) = mean(trT);
        erT(kk,ii) = std(trT);
        mrspT(kk,ii) = mean(trspT);
        erspT(kk,ii) = std(trspT);
        mrspTST(kk,ii) = (getResponseFrame(d.animalNumber(an),d.pEvokedFolders{an,ii})-idx1)*Ts;
        muvAreaUC(kk,ii) = mean(tuvAreaUC);
        if an == 6
            mang(kk,ii) = angle(-ctranspose(mean(tang)));
            eang(kk,ii) = angle(-ctranspose(std(tang)/sqrt(length(tang))));
        else
            mang(kk,ii) = angle(mean(tang));
            eang(kk,ii) = angle(std(tang)/sqrt(length(tang)));
        end
        %         allang{kk,ii} = tang;
        allang{kk,ii} = tangdur;
        alldf{kk,ii} = tdf;
        allRiseTime{kk,ii} = riseTime;
        alldfd{kk,ii} = tdfd;
        alluv{kk,ii} = tang;
        allrT{kk,ii} = trT;
        allrspT{kk,ii} = trspT;
        
        %         allangDur{kk,ii} = tangdur;
         avgAllThisdf{ii} = mean(allThisdf{ii});
        avgAllThisuv{ii} = mean(allThisuv{ii});
        avgAllThisuvabs{ii} = mean(allThisuvabs{ii});
        avgAmp(ii) = max(avgAllThisdf{ii}) - mean(avgAllThisdf{ii}(1:round(0.1*Fs)));
        avgSpeed(ii) = max(avgAllThisuv{ii}) - mean(avgAllThisuv{ii}(1:round(0.1*Fs)));
        avgSpeedabs(ii) = max(avgAllThisuvabs{ii}) - mean(avgAllThisuvabs{ii}(1:round(0.1*Fs)));
%         avgAmp(ii) = max(avgAllThisdf{ii}) - min(avgAllThisdf{ii});
%         avgSpeed(ii) = max(avgAllThisuv{ii}) - min(avgAllThisuv{ii});
    end
end
output.allRiseTime = allRiseTime;
output.allThisuv = allThisuv;
output.allThisdf = allThisdf;

output.avgAllThisuv = avgAllThisuv;
output.avgAllThisdf = avgAllThisdf;
output.avgAmp = avgAmp;
output.avgSpeed = avgSpeed;
output.avgSpeedabs = avgSpeedabs;

if exist('allThisdf_discarded','var')
    output.allThisdf_discarded = allThisdf_discarded;
else
    output.allThisdf_discarded = [];
end
output.onsets = onsets;
output.pfs = pfs;
output.mdf = mdf;
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
    output.absalluv{ii} = pf2mms*abs(alluv{ii});
end
output.allrT = allrT;
output.allrspT = allrspT;

output.muvAreaUC = muvAreaUC;
output.Fs = Fs;

%% Average files
for kk = 1:length(anum)
    an = anum(kk);
    for ii = 1:size(d.davg,2)
        if isempty(d.davg{an,ii})
            continue;
        end
        [lists, names]= getEvokedListsDiffStimAmps_glut(d.animalNumber(an),multLevels,dataType);
        inrpFactor = 10;
        file = names{4}{ii};
        % get Fs from file name
        IndexHz = max(strfind(lower(file),'hz'));
        Fs = file(1:IndexHz-1);
        IndexU = max(strfind(Fs,'_'));
        Fs = Fs(IndexU+1:end);
        Fs = str2double(Fs);
        idx1 = round(Fs*0.3)+1;
        idxR = idx1:idx1+round(Fs*0.1)-1; % response indeces
        idx1i = inrpFactor*round(Fs*0.3)+1;
        idxBi = idx1i-round(Fs*0.1)*inrpFactor:idx1i-1; % baseline indeces
        pf2mms = Fs/15; % pixel/frame to mm/sec
        
        thisdf = d.davg{an,ii}{roi};
        uthisdf = interp(thisdf,inrpFactor);
        thisuv = d.d1avg{an,ii}{roi};
        whole_trace = thisdf;
        whole_trace_uv = thisuv;
        if ss == 3
            idxR = 31:31+25;
            baseline = mean(thisdf(12:26));
            thisuv1ang = thisuv(31:end);
        else
            baseline = mean(uthisdf(idxBi));
            thisuv1ang = thisuv(idx1:end);
            stdbaseline = std(uthisdf(idxBi));
            respThreshold = baseline + 2 *stdbaseline;
            rinds = find(uthisdf(idx1i:end)>respThreshold);
            responseFrame = rinds(find(rinds>=1,1,'first'));
            n = 0;
        end
        thisdf = thisdf(idxR); thisuv = thisuv(idxR);
        pf = find(thisdf == max(thisdf));
        pfuv = find(abs(thisuv) == max(abs(thisuv)));
        
        motif_frames_N = ((pf+idx1-1)-pre_frames):((pf+idx1-1)+post_frames);
        thisMotif = whole_trace(motif_frames_N);
        baseline_thisMotif = findBleachingTrend(thisMotif);
        intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
        motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < (pf+idx1-1),1,'last')));
        
        AvgallThisdf{kk}{ii} = whole_trace;
        AvgallThisuv{kk}{ii} = whole_trace_uv;
        preNFrames = round(0.1*Fs);
        postNFrames = round(0.4*Fs);
        AvgallThisdf{kk}{ii} = whole_trace((motif_onset-preNFrames):(motif_onset+postNFrames))-mean(whole_trace((motif_onset-preNFrames):motif_onset));
        AvgallThisuv{kk}{ii} = pf2mms*abs(whole_trace_uv((motif_onset-preNFrames):(motif_onset+postNFrames)));%-mean(10*abs(whole_trace_uv((motif_onset-preNFrames):(motif_onset))));
        
        if isempty(motif_onset)
            continue;
        end
        Tmin = round(Fs*0.1);
        if motif_onset <= Tmin
            continue;
        end
        baseline_mean = mean(whole_trace((motif_onset-Tmin):motif_onset));
        baseline_std = std(whole_trace((motif_onset-Tmin):motif_onset));
        
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
            
            ttp = pf * Ts;
            rlat = (responseFrame) * Ts/inrpFactor;
        catch
            continue;
        end        
        Avgallang{kk,ii} = uvdur;
        Avgalldf{kk,ii} = dfval;
        
%         Avgalluv{kk,ii} = max(AvgallThisuv{kk}{ii}-min(AvgallThisuv{kk}{ii}));
        Avgalluv{kk,ii} = max(AvgallThisuv{kk}{ii}-mean(AvgallThisuv{kk}{ii}(1:preNFrames)));
%         Avgalluv{kk,ii} = uvval;
        Avgallttp(kk,ii) = ttp;
        avg_onsets(kk,ii) = motif_onset;
        try
            Avgallrlat(kk,ii) = rlat;
        catch
        end
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

output.avg = getValueFromAverageFiles(d,anum,roi,ss,amplitude_threshold,Fs,pf2mms,varargin);

function avg = getValueFromAverageFiles(d,anum,roi,ss,amplitude_threshold,Fs,pf2mms,varargin)
anum
mainFolder = getMainDataFolderAC;
for ii = 1:3
    pEvokedFolder = d.pEvokedFolders{anum,ii};
    ISF = makeName(sprintf('ImgSeq.mat',anum),pEvokedFolder);
    ISF = makeName(ISF,mainFolder);
    uvISF = makeName(sprintf('uv_ImgSeq.mat'),pEvokedFolder);
    uvISF = makeName(uvISF,mainFolder);
    temp = load(ISF);
    imgseq = temp.ImgSeq;
    temp = load(uvISF);
    uv = temp.uv;
    mask = double(d.masks{anum}{roi});
    vals = getMaskedValues(imgseq,mask);
    avg.df = mean(vals,2);
    vals = getMaskedValues(uv,mask);
    avg.uv = mean(vals,2);
    avg.uvabs = mean(abs(vals),2);
    
%     pf = round(length(avg.uv)/2)+1;
%     pre_frames = round(Fs*0.17);
%     post_frames = round(Fs*0.17);
%     motif_frames_N = (pf-pre_frames):(pf+post_frames);
%     thisMotif = avg.df(motif_frames_N);
%     % motif_frames_N = 1:length(thisMotif);
%     baseline_thisMotif = findBleachingTrend(thisMotif);
%     intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
%     motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < pf,1,'last')))-1;
    motif_onset = round(0.3*Fs);
%     avg.baseline_df = baseline_thisMotif;
%     figure(10000);clf;plot(motif_frames_N,thisMotif);hold on;plot(motif_frames_N,baseline_thisMotif);
    avg.Fs = Fs;
    avg.Ts = 1000/Fs;
    thisdf = avg.df;
    % thisuv = avg.uv;
    % thisuvabs = avg.uvabs;
    dur = round(Fs*0.1);
%     if (motif_onset-dur) <= 0
%         dur = motif_onset - 1;
%     end
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
    avg.uvabs_chopped{ii} = avg.uvabs_chopped{ii} - min(avg.uvabs_chopped{ii}(1:round(0.1*Fs)));
    avg.riseTime(ii) = (pf-motif_onset)*avg.Ts;
    avg.ts = linspace(0,(preDur+postDur),length(avg.df_chopped));

    avg.Amp(ii) = max(avg.df_chopped{ii});% - mean(avg.df_chopped{ii}(1:round(0.1*Fs)));
%     avg.Speed(ii) = max(avg.uv_chopped{ii}) - mean(avg.uv_chopped{ii}(1:round(0.1*Fs)));
    avg.Speed(ii) = max(avg.uv_chopped{ii});% - mean(avg.uv_chopped{ii}(1:round(0.1*Fs)));
    avg.Speedabs(ii) = max(avg.uvabs_chopped{ii});% - mean(avg.uvabs_chopped{ii}(1:round(0.1*Fs)));
    avg.auc(ii) = sum(avg.uv_chopped{ii}((round(0.1*Fs):(round(0.1*Fs)+round(0.2*Fs)))));
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
temp = makeName(folderName,getMainDataFolderAC);
fileName = makeName('stimTemplate.mat',temp);
try
    temp = load(fileName);
    responseFrame = temp.stimTemplate.frameNumbers(1);
catch
    responseFrame = 28;
end


