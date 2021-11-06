function dists = getDistributionsSpon_motifs_simple_artur (fde,d,anum,roi,match,CrossOrMax,thn,LMH,amplitude_threshold,varargin)
% CrossOrMax is either Cross or Max
% matches is either FL or HL
multLevels = 1;
if nargin > 9
   multLevels = varargin{1}; 
end

if strcmp(match,'FL')
    matchROI = 2;
    ss = 1;
elseif strcmp(match,'HL')
    matchROI = 3;
    ss = 2;
elseif strcmp(match,'VC')
    matchROI = 5;
    ss = 3;
end
animalNumber = {'M100511';'M101111';'M081712_Cnt_Rig4';'M082212_Cnt_Rig4';'M082412_Cnt_Rig4'};
if ss == 3
% animals with alternative FL stim, VC stim
animalNumber = {'M011211_motif_infraslow_Allen'; 'M012111_motif_infraslow_Allen';
        'M101911_motif_infraslow_alt_stim'; 'M102411_motif_infraslow_alt_stim'};
end
% amplitude_threshold = 0.1;
alldf = [];
alldfd = [];
alluv = [];
allang = [];
allfd = [];
allFrames = [];
allFramesuv = [];
stimAmpNum = [];
for an = anum
    mask = getMask(d.animalNumber(anum));
    [lists, names] = getEvokedListsDiffStimAmps(fde{ss}.animalNumber(an),multLevels);
    if multLevels
        ssLMH = cell2mat(fde{ss}.anGroup{an});
    else
        ssLMH = 1;
        LMH = 1;
    end
    LMHname = names{ss}{ssLMH(LMH)};
    
    listR = listOfRecordingsSpon(animalNumber(an));
    jjc = 0;
    jjcd = 0;
    for ii = 1:length(listR)
        %         [an ii]
        motifs = getMotifsFromCCResultsSam_simple(animalNumber(an),listR,ii,thn,match,CrossOrMax,LMHname,multLevels);
        
        if isempty(motifs)
            continue;
        end
        if ii == 1
            stimAmpNum = [stimAmpNum motifs.stimAmpNum];
        end
        %         dfMatch = d.d{an,ii}{matchROI};
        %         uvMatch = abs(d.d1{an,ii}{matchROI});
        dfMatch = d.d{an,ii}{roi};
        uvMatch = abs(d.d1{an,ii}{roi});
        
        motifs = refineMotifscc(motifs,dfMatch,uvMatch);
        framesdf = motifs.framesdf;
        framesuv = motifs.framesuv;
        
        thisdf = d.d{an,ii}{roi};
        thisuv = d.d1{an,ii}{roi};
        thisfd = d.fd{an,ii}{roi};
        %         dfth = mean(dfMatch)+std(dfMatch);
%         baseline_mean = mean(thisdf);%(dfMatch<dfth));
        mthisuv1 = 0*mean(thisuv);%((dfMatch(1:end-1)+dfMatch(2:end))/2 < dfth));
        
        %         [an ii length(frames)]
%         for jj = 1:length(framesdf)
        spon_file_name = makeName(sprintf('%s.mat',listR{ii}),makeName(d.animalNumber{anum},getMainDataFolder));
        msf = matfile(spon_file_name);
        spon_file_name_uv = makeName(sprintf('uv_%s.mat',listR{ii}),makeName(sprintf('%s',listR{ii}),makeName('pSpon',makeName(d.animalNumber{anum},getMainDataFolder))));
        msf_uv = matfile(spon_file_name_uv);
        pre_frames = 25;
        post_frames = 25;
        jj = 1;
        while jj < length(framesdf)
            pf = framesdf(jj);
            pfuv = framesuv(jj);
            motif_frames_N = (pf-pre_frames):(pf+post_frames);
            thisMotif = thisdf(motif_frames_N);
            DthisMotif = diff(diff(thisMotif));
            baseline_thisMotif = findBleachingTrend(thisMotif);
            intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
            motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < pf,1,'last')))-1;
            if isempty(motif_onset)
                jj = jj + 1;
                continue;
            end
            baseline_mean = mean(thisdf((motif_onset-15):motif_onset));
            baseline_std = std(thisdf((motif_onset-15):motif_onset));
            % for applying method to find spatiotemporal patterns 
%             m_motif_frames = applyMask(motif_frames,mask);
%             st_demo_vsdi (m_motif_frames);


%             uv_motif_frames = msf_uv.uv(:,:,motif_frames_N);
%             uv_motif_frames(:,:,end) = [];
%             SoSi = SourceSink_Analysis(motif_frames,uv_motif_frames,mask);
%             if length(SoSi.source) > 0
%                 n = 0;
%             end
            % for displaying image sequence ... comment out later
%             dfval = thisdf(pf)-baseline_mean;
%             if dfval < (baseline_mean + amplitude_threshold*baseline_std)
%                 jj = jj + 1;
%                 continue;
%             end
            display_img_seq = 0;
            if display_img_seq
                motif_frames = msf.DF_F0(:,:,motif_frames_N);
                motif_frames = applyMask(motif_frames,mask);
                Imin = min(motif_frames(:));
                Imax = max(motif_frames(:));
%                 motif_frames = motif_frames - Imin;
                thisROI = getROI(d.animalNumber(anum),'ROI_FL');

                bo = bwboundaries(thisROI);
                xs = bo{1}(:,1);        ys = bo{1}(:,2);
                figure(11);clf;
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
%                     if length(SoSi.source) > 0
%                         if mfn == SoSi.source(3)
%                             plot(SoSi.source(1),SoSi.source(2),'*','color','r');
%                             pause;
%                         end
%                     end
                    subplot 212
                    cla
                    plot(motif_frames_N,thisMotif);hold on;
                    baseline_thisMotif = findBleachingTrend(thisMotif);
                    plot(motif_frames_N,baseline_thisMotif,'m');
                    intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
%                     baseline_thisMotif = findBleachingTrend(thisMotif,5);
%                     plot(motif_frames_N,baseline_thisMotif,'k');
%                     baseline_thisMotif = findBleachingTrend(thisMotif,7);
%                     plot(motif_frames_N,baseline_thisMotif,'c');
%                     plot(motif_frames_N(3:end),DthisMotif);
%                     plot(motif_frames_N,ones(size(thisMotif))*dfThreshold);
                    motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < pf,1,'last')));
                    titleT = sprintf('%s - %s - %d',d.animalNumber{anum},listR{ii},motif_onset);
                    title(titleT,'Interpreter','none');
                    ylims = ylim;
                    plot([motif_frames_N(mfn) motif_frames_N(mfn)],[ylims(1) ylims(2)],'r');
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

                % flow directionality as a measure of the similarity betwwn
                % velocity vector at max speed frame and 6 frames before and after
                fd = [thisuv(pfang); thisuv(pfang-6:pfang-1); thisuv(pfang+1:pfang+6)];
                NN = length(fd);
                fd1 = 0;
                u = fd(1); % velocity vector at max speed frame
                for nn = 1:NN-1
                    v = fd(nn+1);
                    Theta = angleUV(u,v); % find angle between the two vector
                    fd1 = fd1 + cosd(Theta)/(NN-1);
                end
                fd = fd1;
            catch
                jj = jj + 1;
                continue;
            end
%             if dfval < amplitude_threshold
%             if dfval < (baseline_mean + amplitude_threshold*baseline_std)
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
            allThisdf(jjc,:) = thisdf((motif_onset-15):(motif_onset+65))-baseline_mean;
            allThisuv(jjc,:) = 10*abs(thisuv((motif_onset-15):(motif_onset+60)));
            riseTime(jjc) = (pf-motif_onset)*6.67;
%             allra = [allra dfval-baseline_mean];
            if ~display_img_seq
                jj = jj + 1;
            end
        end
    end
    
end
dists.riseTime = riseTime;
dists.alldf = alldf;
dists.alldfd = alldfd;
dists.alluv = alluv;
dists.absalluv = 10*abs(alluv);
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
dists.ss = ss;