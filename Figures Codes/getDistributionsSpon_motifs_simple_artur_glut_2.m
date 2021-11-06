function dists = getDistributionsSpon_motifs_simple_artur_glut_2 (fde,d,anum,roi,match,CrossOrMax,thn,LMH,amplitude_threshold,varargin)
% CrossOrMax is either Cross or Max
% matches is either FL or HL
multLevels = 1;
if nargin > 9
   multLevels = varargin{1}; 
end
ss = 1;
animalNumber = d.animalNumber;

alldf = [];
alldf1 = [];alldf2 = [];
alldfd = [];
alluv = [];
alluvpp = [];
allang = [];
allfd = [];
allFrames = [];
allFramesuv = [];
stimAmpNum = [];
for an = anum
    mask = double(d.masks{anum}{1});
    [lists, names] = getEvokedListsDiffStimAmps_glut(fde{ss}.animalNumber(an),multLevels);
    if multLevels
        ssLMH = cell2mat(fde{ss}.anGroup{an});
    else
        ssLMH = 1;
        LMH = 1;
    end
    LMHname = names{4}{ssLMH(LMH)};
    listR = listOfRecordingsSpon_glut(animalNumber(an));
    jjc = 0;
    jjcd = 0;
    for ii = 1%:length(listR)
        motifs = getMotifsFromCCResultsSam_simple_glut(animalNumber(an),listR,ii,thn,match,CrossOrMax,LMHname,fde{1},multLevels);
        if isempty(motifs)
            continue;
        end
        if ii == 1
            stimAmpNum = [stimAmpNum motifs.stimAmpNum];
        end
        dfMatch = d.d{an,ii}{roi};
        uvMatch = abs(d.d1{an,ii}{roi});
        motifs = refineMotifs_glut_1(motifs,dfMatch,uvMatch);
        framesdf = motifs.framesdf;
        framesuv = motifs.framesuv;
        thisdf = d.d{an,ii}{roi};
        thisuv = d.d1{an,ii}{roi};
        mthisuv1 = 0*mean(thisuv);%((dfMatch(1:end-1)+dfMatch(2:end))/2 < dfth));
        spon_file_name = makeName(sprintf('%s.mat',listR{ii}),makeName(d.animalNumber{anum},getMainDataFolderGlut));
        msf = matfile(spon_file_name);
        spon_file_name_uv = makeName(sprintf('uv_%s.mat',listR{ii}),makeName(sprintf('%s',listR{ii}),makeName('pSpon',makeName(d.animalNumber{anum},getMainDataFolderGlut))));
        pre_frames = 22;
        post_frames = 22;
        jj = 1;
        stdn = 0.5;
        while jj <= length(framesdf)
            pf = framesdf(jj);
            pfcc = motifs.frames(jj);
            pfuv = framesuv(jj);
            try
                thisdfChopped = thisdf((pf-100):(pf+100));
            catch
                jj = jj + 1;
                continue;
            end
            mWholeTrace = mean(thisdf((pf-100):(pf+100)));
            stdWholeTrace = std(thisdf((pf-100):(pf+100)));
            thresh = mWholeTrace + (stdn * stdWholeTrace);
            
            indOnset = pf-101+find(thisdfChopped(1:100)>thresh,1,'last');
            mb = mean(thisdf((indOnset-50):(indOnset)));
            dfval2 = thisdf(pf)-mb;
            alldf2 = [alldf2 dfval2];
            
            mthisuv = mean(thisdf((pf-19):(pf-10)));
            dfval1 = thisdf(pf)-mthisuv;
            alldf1 = [alldf1 dfval1];
            motif_frames_N = (pf-pre_frames):(pf+post_frames);
            thisMotif = thisdf(motif_frames_N);
            baseline_thisMotif = findBleachingTrend(thisMotif);
            intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
            motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < pf,1,'last')))-1;
            baseline_mean = mean(thisdf((motif_onset-10):motif_onset));
            baseline_std = std(thisdf((motif_onset-10):motif_onset));
            if motif_onset < motif_frames_N(1)
                motif_onset = motif_frames_N(1);
            end
            if isempty(motif_onset)
                jj = jj + 1;
                continue;
            end
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
                % flow directionality as a measure of the similarity betwwn
                % velocity vector at max speed frame and 6 frames before and after
                fd = [thisuv(pfang); thisuv(pfang-10:pfang-1); thisuv(pfang+1:pfang+10)];
                NN = length(fd);
                fd1 = 0;
                u = fd(1); % velocity vector at max speed frame
                ccc = [];
                for nn = 1:NN-1
                    v = fd(nn+1);
                    Theta = angleUV(u,v); % find angle between the two vector
                    ccc(nn) = cosd(Theta);
                    fd1 = fd1 + cosd(Theta)/(NN-1);
                end
                fd = fd1;
                
            catch
                jj = jj + 1;
                continue;
            end
            mean_post_onset = mean(thisdf(motif_onset:(pf+15)));
            if mean_post_onset < (baseline_mean + amplitude_threshold*baseline_std)
                jj = jj + 1;
                jjcd = jjcd+1;
                allThisdf_discarded(jjcd,:) = thisdf((motif_onset-15):(motif_onset+25))-baseline_mean;
                alldfd = [alldfd dfval];
                continue;
            end
            allang = [allang uvdur];
            allfd = [allfd,fd];
            alldf = [alldf dfval]; alluv = [alluv (uvval)];
            allFrames = [allFrames;[an ii pf motif_onset motifs.frames(jj)]];
            allFramesuv = [allFramesuv;[an ii pfuv]];
            jjc = jjc + 1;
            preNFrames = 0.1*100;
            postNFrames = 0.4*100;
            allThisdf(jjc,:) = thisdf((motif_onset-preNFrames):(motif_onset+postNFrames))-baseline_mean;
            allThisuv(jjc,:) = 6.67*abs(thisuv((motif_onset-preNFrames):(motif_onset+postNFrames)));
            riseTime(jjc) = (pf-motif_onset)*10;
            jj = jj + 1;
        end
        avgAllThisdf{ii} = mean(allThisdf);
        avgAllThisuv{ii} = mean(allThisuv);
        avgAmp(ii) = max(avgAllThisdf{ii}) - min(avgAllThisdf{ii});
        avgSpeed(ii) = max(avgAllThisuv{ii}) - min(avgAllThisuv{ii});
    end
end
dists.riseTime = riseTime;
alldf1(alldf1<0)=[];
dists.alldf = alldf;
dists.alldfd = alldfd;
dists.alluv = alluv;
dists.avgAllThisuv = avgAllThisuv;
dists.avgAllThisdf = avgAllThisdf;
dists.avgAmp = avgAmp;
dists.avgSpeed = avgSpeed;

mUV = mean(allThisuv);
dists.absalluv = max(mUV - min(mUV));
dists.absalluv = abs(alluv)*6.67;
dists.allang = allang;
dists.allfd = allfd;
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
dists.avg = getValueFromAverageFiles(fde,d,anum,roi,match,CrossOrMax,thn,LMH,amplitude_threshold,100,6.667,varargin);

function avg = getValueFromAverageFiles(fde,d,anum,roi,match,CrossOrMax,thn,LMH,amplitude_threshold,Fs,pf2mms,varargin)
mainFolder = '\\mohajerani-nas.uleth.ca\storage2\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test';
ISF = makeName(sprintf('ImgSeq_BPF_05_6Hz_avg_an%d_AC_glut_all_04corr.mat',anum),mainFolder);
uvISF = makeName(sprintf('uv_ImgSeq_BPF_05_6Hz_avg_an%d_AC_glut_all_04corr.mat',anum),mainFolder);
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
baseline_thisMotif = findBleachingTrend(thisMotif);
intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < pf,1,'last')))-1;
avg.baseline_df = baseline_thisMotif;
avg.Fs = Fs;
avg.Ts = 1000/Fs;
thisdf = avg.df;
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
avg.riseTime = (pf-motif_onset)*avg.Ts;
avg.ts = linspace(0,(preDur+postDur),length(avg.df_chopped));
avg.Amp = max(avg.df_chopped);
avg.Speed = max(avg.uv_chopped);
avg.Speedabs = max(avg.uvabs_chopped);


