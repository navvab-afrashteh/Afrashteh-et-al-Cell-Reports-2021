function dists = getDistributionsSpon_simple_artur_glut (d,anum,roi,ss,stdn,amplitude_threshold)
alldf = [];
alldf1 = [];
alldfd = [];
alluv = [];
allang = [];
allFrames = [];
allFramesuv = [];
riseTime = [];
allThisdf = [];
allThisuv = [];
% amplitude_threshold = 0.1;
for an = anum
    listR = listOfRecordingsSpon_glut(d.animalNumber(anum));
    jjc = 0;
    jjcd = 0;
    for ii = 1:size(d.d,2)

        if isempty(d.d{an,ii})
            continue;
        end
        if isfield(d,'DFthFrames')
            frames = d.DFthFrames{an,ii}{ss};
            fdf = d.d{an,ii}{roi};
            mfdf = mean(fdf);
            stdfdf = std(fdf);
            dfThreshold = mfdf+ 1 * stdfdf;
        else
            if ss == 1
                fdf = d.d{an,ii}{2};
            else
                fdf = d.d{an,ii}{3};
            end
            fdf = d.d{an,ii}{roi};
            % early and late part of signal is not cnosidered to remove 
            % possible outliers due to filtering
            mfdf = mean(fdf(800:end-800));
            stdfdf = std(fdf(800:end-800));
            dfThreshold = mfdf+ stdn * stdfdf;
            thePeaks = detectPeaks(fdf); % finds peaks ... returns nx2 matrix ... column 1 is peak positions
            thePeaks = thePeaks(thePeaks(:,1)>400,:);
            % column 2 is peak values
            ppos = find(thePeaks(:,2) > (dfThreshold)); % apply threshold
            frames = thePeaks(ppos,1); % find frame of interest
        end

        if isempty(frames)
            continue;
        end

        dfMatch = d.d{an,ii}{roi};
        uvMatch = abs(d.d1{an,ii}{roi});
        
        dur = 15;
        diff_idx = find(diff(frames)>dur);
        frames_15 = zeros(1,length(diff_idx)+1);
        if ~isempty(diff_idx)
            for p = 1:length(diff_idx)+1
                if p == 1
                    temp = frames(1:diff_idx(p));
                elseif p == length(diff_idx)+1
                    temp = frames(diff_idx(p-1)+1:end);
                else
                    temp = frames(diff_idx(p-1)+1:diff_idx(p));
                end
                [~,maxdf_temp] = max(dfMatch(temp));
                frames_15(p) = temp(maxdf_temp);
            end
        else
            [~,maxdf_temp] = max(dfMatch(frames));
                frames_15 = frames(maxdf_temp);
        end
    
        motifs.frames = frames_15;
        motifs = refineMotifs(motifs,dfMatch,uvMatch);
        framesdf = motifs.framesdf;
        framesuv = motifs.framesuv;
        
        thisdf = d.d{an,ii}{roi};
%         figure(12);clf;hist(thisdf(400:end));
        thisuv = d.d1{an,ii}{roi};
%         baseline_mean = mean(thisdf);
        mthisuv1 = 0*mean(thisuv);
%         for jj = 1:length(framesdf) % for each frame find df, uv, angle, and activity amplitude
        pre_frames = 17;
        post_frames = 17;
        jj = 1;
        while jj < length(framesdf)
            % activity amplitude is peak value - mean of the whole signal
            pf = framesdf(jj);
            pfuv = framesuv(jj);
            mthisuv = mean(thisdf((pf-19):(pf-10)));
            dfval1 = thisdf(pf)-mthisuv;
            alldf1 = [alldf1 dfval1];
%             baseline_mean = mean(thisdf((pf-29):(pf-15)));
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
            
            % for displaying image sequence ... comment out later
            display_img_seq = 0;
            if display_img_seq
                spon_file_name = makeName(sprintf('%s.mat',listR{ii}),makeName(d.animalNumber{anum},getMainDataFolder));
                msf = matfile(spon_file_name);
                motif_frames = msf.DF_F0(:,:,motif_frames_N);
                mask = getMask(d.animalNumber(anum),'ec',0.9);
                thisROI = getROI(d.animalNumber(anum),'ROI_FL');
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
                    titleT = sprintf('Frame %d',motif_frames_N(mfn));
                    title(titleT);
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
                    plot(motif_frames_N,ones(size(thisMotif))*dfThreshold);
                    motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < pf,1,'last')));
                    titleT = sprintf('%d',motif_onset);
                    title(titleT);
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
                uvval = thisuv(pfuv);
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
%             if dfval < amplitude_threshold
            if dfval < (baseline_mean + amplitude_threshold*baseline_std)
                jj = jj + 1;
                jjcd = jjcd+1;
                allThisdf_discarded(jjcd,:) = thisdf((motif_onset-15):(motif_onset+25))-baseline_mean;
                alldfd = [alldfd dfval];
                continue;
            end
            allang = [allang uvdur];
%             if an == 4
% %                 allang = [allang angle(conj(-uvval))];
%                 allang = [allang angle(conj(-uvdur))];
%             else
% %                 allang = [allang angle(conj(uvval))];
%                 allang = [allang angle(conj(uvdur))];
%             end
            alldf = [alldf dfval]; alluv = [alluv (uvval)];
            allFrames = [allFrames;[an ii pf]];
            allFramesuv = [allFramesuv;[an ii pfuv]];
            jjc = jjc + 1;
            allThisdf(jjc,:) = thisdf((motif_onset-10):(motif_onset+40))-baseline_mean;
            allThisuv(jjc,:) = 10*abs(thisuv((motif_onset-10):(motif_onset+40)));
            riseTime(jjc) = (pf-motif_onset)*10;
%             allra = [allra dfval-baseline_mean];
            if ~display_img_seq
                jj = jj + 1;
            end
        end
        n = 0;
    end
    n = 0;
end
% the return values are the same as for evoked file ... dfs, uvs, angles,
% ras from all animals
dists.riseTime = riseTime;
alldf1(alldf1<0)=[];
dists.alldf = alldf;
dists.alldfd = alldfd;
dists.alluv = alluv;
dists.absalluv = 10*abs(alluv);
dists.allang = allang;
dists.allFrames = allFrames;
dists.allFramesuv = allFramesuv;
dists.allThisdf = allThisdf;
dists.allThisuv = allThisuv;
if exist('allThisdf_discarded','var')
    dists.allThisdf_discarded = allThisdf_discarded;
else
    dists.allThisdf_discarded = [];
end
% dists.allra = allra;
