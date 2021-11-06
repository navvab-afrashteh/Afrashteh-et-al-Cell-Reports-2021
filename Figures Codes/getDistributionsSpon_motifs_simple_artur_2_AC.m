function dists = getDistributionsSpon_motifs_simple_artur_2_AC (fde,d,anum,roi,match,CrossOrMax,thn,LMH,amplitude_threshold,varargin)
% CrossOrMax is either Cross or Max
% matches is either FL or HL
multLevels = 1;
if nargin > 9
    multLevels = varargin{1};
end
dataType = 'AC';
if nargin > 10
    dataType = varargin{2};
end

ss = 1;
animalNumber = d.animalNumber;

% amplitude_threshold = 0.1;
alldf = [];
alldf1 = [];alldf2 = [];
alldfd = [];
alluv = [];
allang = [];
allFrames = [];
allFramesuv = [];
stimAmpNum = [];
for an = anum
    mask = double(d.masks{anum}{1});
    [lists, names] = getEvokedListsDiffStimAmps_glut(fde{ss}.animalNumber(an),multLevels,dataType);
    if multLevels
        ssLMH = cell2mat(fde{ss}.anGroup{an});
    else
        ssLMH = 1;
        LMH = 1;
    end
    LMHname = names{4}{ssLMH(LMH)};
    
    listR = listOfRecordingsSpon_glut(animalNumber(an),dataType);
    jjc = 0;
    jjcd = 0;
    for ii = 1:length(listR)
        % [an ii]
        motifs = getMotifsFromCCResultsSam_simple_AC(animalNumber(an),listR,ii,thn,match,CrossOrMax,LMHname,multLevels);
        
        if isempty(motifs)
            continue;
        end
        file = listR{ii};
        % get Fs from file name
        IndexHz = min(strfind(lower(file),'hz'));
        Fs = file(1:IndexHz-1);
        IndexU = max(strfind(Fs,'_'));
        Fs = Fs(IndexU+1:end);
        Fs = str2double(Fs);
        pf2mms = Fs/15; % pixel/frame to mm/sec
        Ts = 1000/Fs; % frame-rate time duration
        
        if ii == 1
            stimAmpNum = [stimAmpNum motifs.stimAmpNum];
        end
        %         dfMatch = d.d{an,ii}{matchROI};
        %         uvMatch = abs(d.d1{an,ii}{matchROI});
        dfMatch = d.d{an,ii}{roi};
        uvMatch = abs(d.d1{an,ii}{roi});
        
        motifs = refineMotifs_AC(motifs,dfMatch,uvMatch,Fs);
        framesdf = motifs.framesdf;
        framesuv = motifs.framesuv;
        
        thisdf = d.d{an,ii}{roi};
        thisuv = d.d1{an,ii}{roi};
        thisuvabs = d.d1abs{an,ii}{roi};
%         abs(thisuv)-thisuvabs
        % dfth = mean(dfMatch)+std(dfMatch);
        % baseline_mean = mean(thisdf);%(dfMatch<dfth));
        mthisuv1 = 0*mean(thisuv);%((dfMatch(1:end-1)+dfMatch(2:end))/2 < dfth));
        
        % [an ii length(frames)]
        % for jj = 1:length(framesdf)
        spon_file_name = makeName(sprintf('%s.mat',listR{ii}),makeName(d.animalNumber{anum},getMainDataFolderAC));
        msf = matfile(spon_file_name);
        spon_file_name_uv = makeName(sprintf('uv_%s.mat',listR{ii}),makeName(sprintf('%s',listR{ii}),makeName('pSpon',makeName(d.animalNumber{anum},getMainDataFolderAC))));
        msf_uv = matfile(spon_file_name_uv);
        pre_frames = round(Fs*0.17);
        post_frames = round(Fs*0.17);
        jj = 1;
        stdn = 0.5;
        while jj <= length(framesdf)
            pf = framesdf(jj);
            pfcc = motifs.frames(jj);
            pfuv = framesuv(jj);
            try
                dur = round(Fs);
                thisdfChopped = thisdf((pf-dur):(pf+dur));
            catch
                jj = jj + 1;
                continue;
            end
            mWholeTrace = mean(thisdf((pf-dur):(pf+dur)));
            stdWholeTrace = std(thisdf((pf-dur):(pf+dur)));
            thresh = mWholeTrace + (stdn * stdWholeTrace);
            
            indOnset = pf-dur-1+find(thisdfChopped(1:dur)>thresh,1,'last');
            dur2 = round(dur/2);
            mb = mean(thisdf((indOnset-dur2):(indOnset)));
            dfval2 = thisdf(pf)-mb;
            alldf2 = [alldf2 dfval2];
            
            dur1 = round(Fs*0.1);
            mthisuv = mean(thisdf((pf-2*dur1):(pf-dur1)));
            dfval1 = thisdf(pf)-mthisuv;
            alldf1 = [alldf1 dfval1];
            
            motif_frames_N = (pf-pre_frames):(pf+post_frames);
            thisMotif = thisdf(motif_frames_N);
            % DthisMotif = diff(diff(thisMotif));
            baseline_thisMotif = findBleachingTrend(thisMotif);
            intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
            motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < pf,1,'last')))-1;
            dur = round(Fs*0.1);
            baseline_mean = mean(thisdf((motif_onset-dur):motif_onset));
            baseline_std = std(thisdf((motif_onset-dur):motif_onset));
            if motif_onset < motif_frames_N(1)
                motif_onset = motif_frames_N(1);
            end
            if isempty(motif_onset)
                jj = jj + 1;
                continue;
            end
            %             baseline_mean = mean(thisdf((motif_onset-10):motif_onset));
            %             baseline_std = std(thisdf((motif_onset-10):motif_onset));
            if isempty(motif_onset)
                %                 onsets(ii,jj) = NaN;
                jj = jj + 1;
                continue;
            end
            %% for displaying image sequence ... comment out later
            display_img_seq = 0;
            if display_img_seq
                motif_frames = msf.DF_F0(:,:,motif_frames_N);
                motif_frames = applyMask(motif_frames,mask);
                Imin = min(motif_frames(:));
                Imax = max(motif_frames(:));
                thisROI = double(d.masks{anum}{2});
                
                bo = bwboundaries(thisROI);
                xs = bo{1}(:,1);        ys = bo{1}(:,2);
                figure(11);clf;
                dfval = thisdf(pf)-baseline_mean;
                for mfn = 1:length(motif_frames_N)
                    this_motif_frame = motif_frames(:,:,mfn);%applyMask(motif_frames(:,:,mfn),mask);
                    subplot 211
                    cla
                    imagesc(this_motif_frame,[Imin Imax]);
                    colorbar;
                    hold on;
                    axis equal; axis off;
                    plot(ys,xs,'m','linewidth',0.5);
                    set(gca,'Ydir','reverse');
                    titleT = sprintf('Frame %d - %d',mfn,motif_frames_N(mfn));
                    title(titleT,'Interpreter','none');

                    subplot 212
                    cla
                    plot(motif_frames_N,thisMotif);hold on;
                    plot(motif_frames_N,baseline_thisMotif,'m');
                    titleT = sprintf('Animal Number %s - Recording %s - Onset Frame %d - A_old %.3f - A %.3f',d.animalNumber{anum},listR{ii},motif_onset,dfval1,dfval);
                    title(titleT,'Interpreter','none');
                    ylims = ylim;
                    plot([motif_frames_N(mfn) motif_frames_N(mfn)],[ylims(1) ylims(2)],'r');
                    plot(pf,thisMotif(pre_frames+1),'*r');
                    plot(motif_onset,thisMotif(find(motif_frames_N==motif_onset)),'*g');
                    plot([motif_frames_N(1) motif_frames_N(end)],[mthisuv mthisuv],'k');
                    plot([motif_frames_N(1) motif_frames_N(end)],[baseline_mean baseline_mean],'c')
                    ylim([ylims(1) 2]);
                    if motif_frames_N(mfn)==pf
                        pause(0.1);
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
            n = 0;
            try
                dfval = thisdf(pf)-baseline_mean;
                uvval = thisuv(pfuv)-mthisuv1;
                
                pfang = pfuv;
                fper = 0;
                fpost = 2;
                uvdur = thisuv(pfang-fper:pfang+fpost);
                uvdur = uvdur./abs(uvdur);
                uvdur = mean(uvdur);
                uvdur =  uvdur./abs(uvdur);
            catch
                jj = jj + 1;
                continue;
            end
            % if dfval < amplitude_threshold
            durR = round(Fs*0.15);
            mean_post_onset = mean(thisdf(motif_onset:(pf+durR)));
            if mean_post_onset < (baseline_mean + amplitude_threshold*baseline_std)
                jj = jj + 1;
                jjcd = jjcd+1;
                preD = round(Fs*0.15);
                postD = round(Fs*0.25);
                allThisdf_discarded(jjcd,:) = thisdf((motif_onset-preD):(motif_onset+postD))-baseline_mean;
                alldfd = [alldfd dfval];
                continue;
            end
            allang = [allang uvdur];
            %             if ~(ss==3) && an == 4
            % %                 allang = [allang angle(conj(-uvval))];
            %                 allang = [allang angle(conj(-uvdur))];
            %             else
            % %                 allang = [allang angle(conj(uvval))];
            %                 allang = [allang angle(conj(uvdur))];
            %             end
            alldf = [alldf dfval]; alluv = [alluv (uvval)];
            allFrames = [allFrames;[an ii pf motif_onset]];
            allFramesuv = [allFramesuv;[an ii pfuv]];
            jjc = jjc + 1;
            allThisdf(jjc,:) = thisdf((motif_onset-round(0.1*Fs)):(motif_onset+round(0.4*Fs)))-baseline_mean;
            allThisuv(jjc,:) = pf2mms*abs(thisuv((motif_onset-round(0.1*Fs)):(motif_onset+round(0.4*Fs))));
            allThisuvabs(jjc,:) = pf2mms*abs(thisuvabs((motif_onset-round(0.1*Fs)):(motif_onset+round(0.4*Fs))));
            riseTime(jjc) = (pf-motif_onset)*Ts;
            %             allra = [allra dfval-baseline_mean];
            if ~display_img_seq
                jj = jj + 1;
            end
        end
    end
    avgAllThisdf = mean(allThisdf);
    avgAllThisuv = mean(allThisuv);
    avgAllThisuvabs = mean(allThisuvabs);
    avgAmp = max(avgAllThisdf) - mean(avgAllThisdf(1:round(0.1*Fs)));
    avgSpeed = max(avgAllThisuv) - mean(avgAllThisuv(1:round(0.1*Fs)));
    avgSpeedabs = max(avgAllThisuvabs) - mean(avgAllThisuvabs(1:round(0.1*Fs)));
%     avgAmp = max(avgAllThisdf) - min(avgAllThisdf);
%     avgSpeed = max(avgAllThisuv) - min(avgAllThisuv);
end
dists.riseTime = riseTime;
alldf1(alldf1<0)=[];
alldf2(alldf2<0) = [];
alldf2(isnan(alldf2)) = [];
dists.alldf = alldf;
dists.alldfd = alldfd;
dists.alluv = alluv;

dists.avgAllThisuv = avgAllThisuv;
dists.avgAllThisdf = avgAllThisdf;
dists.avgAmp = avgAmp;
dists.avgSpeed = avgSpeed;
dists.avgSpeedabs = avgSpeedabs;

dists.absalluv = pf2mms*abs(alluv);
dists.allang = allang;
dists.allFrames = allFrames;
dists.allFramesuv = allFramesuv;
dists.stimAmpNum = stimAmpNum;
if ~exist('allThisdf','var')
    allThisdf = [];
end
dists.allThisdf = allThisdf;
dists.allThisuv = allThisuv;
if exist('allThisdf_discarded','var')
    dists.allThisdf_discarded = allThisdf_discarded;
else
    dists.allThisdf_discarded = [];
end
dists.ss = ss;
dists.Fs = Fs;
dists.avg = getValueFromAverageFiles(fde,d,anum,roi,match,CrossOrMax,thn,LMH,amplitude_threshold,Fs,pf2mms,varargin);

function avg = getValueFromAverageFiles(fde,d,anum,roi,match,CrossOrMax,thn,LMH,amplitude_threshold,Fs,pf2mms,varargin)
allmotifonsets = [140,204,208,140,141,145,135];
anum
mainFolder = '\\mohajerani-nas.uleth.ca\storage2\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test';
ISF = makeName(sprintf('ImgSeq_avg_an%d_AC_anaesthetics_all_04corr.mat',anum),mainFolder);
uvISF = makeName(sprintf('uv_ImgSeq_avg_an%d_AC_anaesthetics_all_04corr.mat',anum),mainFolder);
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

pf = round(length(avg.uv)/2)+1;
pre_frames = round(Fs*0.17);
post_frames = round(Fs*0.17);
motif_frames_N = (pf-pre_frames):(pf+post_frames);
thisMotif = avg.df(motif_frames_N);
% motif_frames_N = 1:length(thisMotif);
baseline_thisMotif = findBleachingTrend(thisMotif);
intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < pf,1,'last')))-1;
motif_onset = allmotifonsets(anum);
avg.baseline_df = baseline_thisMotif;
figure(10000);clf;plot(motif_frames_N,thisMotif);hold on;plot(motif_frames_N,baseline_thisMotif);
avg.Fs = Fs;
avg.Ts = 1000/Fs;
thisdf = avg.df;
% thisuv = avg.uv;
% thisuvabs = avg.uvabs;
dur = round(Fs*0.1);
if (motif_onset-dur) <= 0
    dur = motif_onset - 1;
end
baseline_mean = mean(thisdf((motif_onset-dur):motif_onset));
preDur = 0.1; postDur = 0.4;
postF = round(postDur*Fs);
preF = round(preDur*Fs);

if (motif_onset-preF) <= 0
    preF = motif_onset - 1;
end

avg.df_chopped = avg.df((motif_onset-preF):(motif_onset+postF))-baseline_mean;
avg.uv_chopped = pf2mms*abs(avg.uv((motif_onset-preF):(motif_onset+postF)));
avg.uv_chopped = avg.uv_chopped - min(avg.uv_chopped(1:round(0.1*Fs)));
avg.uvabs_chopped = pf2mms*abs(avg.uvabs((motif_onset-preF):(motif_onset+postF)));
avg.uvabs_chopped = avg.uvabs_chopped - min(avg.uvabs_chopped);
riseTime = (pf-motif_onset)*avg.Ts;
avg.ts = linspace(0,(preDur+postDur),length(avg.df_chopped));

avg.Amp = max(avg.df_chopped);% - mean(avg.df_chopped(1:round(0.1*Fs)));
avg.Speed = max(avg.uv_chopped);% - mean(avg.uv_chopped(1:round(0.1*Fs)));
avg.Speedabs = max(avg.uvabs_chopped);% - mean(avg.uvabs_chopped(1:round(0.1*Fs)));

avg.auc = sum(avg.uv_chopped(round(0.1*Fs):(round(0.1*Fs)+round(0.2*Fs))));

