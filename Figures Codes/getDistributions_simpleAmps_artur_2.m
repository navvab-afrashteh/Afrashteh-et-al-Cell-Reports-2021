function output = getDistributions_simpleAmps_artur_2 (d,anum,roi,ss,amplitude_threshold,varargin)
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

            idx = 27:27+10;
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
            thisMotif = whole_trace(motif_frames_N);
            baseline_thisMotif = findBleachingTrend(thisMotif);
            intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
            motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < (pf+26),1,'last')))-1;
            if isempty(motif_onset)
                onsets(ii,jj) = NaN;
                pfs(ii,jj) = NaN;
                jj = jj + 1;
                continue;
            end
            onsets(ii,jj) = motif_onset;
            pfs(ii,jj) = pf+26;
            
            fig_for_paper = 0;
            if fig_for_paper
                ff = makeFigureWindow__one_axes_only(4,[1 1 3.5 2],[0.15 0.25 0.8 0.6]);
                axes(ff.ha);
                peak_vals = detectPeaksN(thisMotif);
                plot(6.67*(motif_frames_N-pfs(ii,jj)),thisMotif,'b');hold on;
                plot(6.67*(motif_frames_N-pfs(ii,jj)),baseline_thisMotif,'m');
                plot((motif_frames_N(peak_vals(:,1))-pfs(ii,jj))*6.67,peak_vals(:,2),'k*');
                plot(6.67*(motif_onset-pfs(ii,jj)),thisMotif(motif_frames_N == motif_onset),'*r');
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
                ImgSeq_file_name = makeName(makeName('ImgSeq.mat',d.pEvokedFolderTrials{anum,ii}{jj}),getMainDataFolder);
                msf = matfile(ImgSeq_file_name);
                motif_frames = msf.ImgSeq(:,:,motif_frames_N);

                thisROI = getROI(d.animalNumber(an),'ROI_FL');
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
                    ylims = ylim;
                    plot([motif_frames_N(mfn) motif_frames_N(mfn)],[ylims(1) ylims(2)],'r');
                    if motif_frames_N(mfn) == motif_onset
                        subplot 211
                        text(10,10,'Estimated Onset','Color','w');
                    end
%                     if motif_frames_N(mfn) > (pf+26-7) && motif_frames_N(mfn) < (pf+26+7)                    
                    if motif_frames_N(mfn) > (motif_onset-3) && motif_frames_N(mfn) < (motif_onset+3)
                        pause(1);
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
                if roi == 3 && ss == 1
                    baseline_mean = mean(whole_trace((motif_onset-10):motif_onset));
                    baseline_std = std(whole_trace((motif_onset-10):motif_onset));
                else
                    baseline_mean = mean(whole_trace((motif_onset-15):motif_onset));
                    baseline_std = std(whole_trace((motif_onset-15):motif_onset));
                end
                dfval = thisdf(pf)-baseline_mean;
                uvval = thisuv(pfuv);
                pfang = pfuv; % max of Speed
%                 pfang = 4;
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
                fd = [thisuv(pfang); thisuv(pfang-m:pfang-1); thisuv(pfang+1:pfang+6)];
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
%                 trialFolder = makeName(trials{jj},makeName(names{ss}{ii},makeName('pEvoked',makeName(d.animalNumber{an},getMainDataFolder))));
                trialFolder = makeName(d.pEvokedFolderTrials{an,ii}{jj},getMainDataFolder);
                fileName = makeName('ImgSeq.mat',trialFolder);
                maguv = abs(thisuv);
                uvAreaUC = trapz(maguv);
%                 areaval = getArea(fileName,pf+26,mask);
                
            catch
                jj = jj + 1;
                continue;
            end
%             if dfval < amplitude_threshold
            mean_post_onset = mean(whole_trace(motif_onset:(pf+26+15)));
            if mean_post_onset < (baseline_mean + amplitude_threshold*baseline_std)
                jj = jj + 1;
                jjcd = jjcd+1;
%                 allThisdf_discarded{ii}(jjcd,:) = whole_trace((motif_onset-15):(motif_onset+25))-baseline_mean;
                tdfd = [tdfd dfval];
                continue;
            end
            tdf = [tdf dfval];
            tfd = [tfd fd];
            tuv = [tuv 10*abs(uvval)];
            tang = [tang uvval];
            trT = [trT rT];
            tuvAreaUC = [tuvAreaUC uvAreaUC];
            trspT = [trspT (responseFrame-27)*6.67];
            tangdur = [tangdur uvdur];
            jjc = jjc + 1;
            if roi == 3 && ss == 1
                allThisdf{ii}(jjc,:) = whole_trace((motif_onset-10):(motif_onset+65))-baseline_mean;
                allThisuv{ii}(jjc,:) = 10*abs(whole_trace_uv((motif_onset-10):(motif_onset+60)));
            else
                allThisdf{ii}(jjc,:) = whole_trace((motif_onset-15):(motif_onset+65))-baseline_mean;
                allThisuv{ii}(jjc,:) = 10*abs(whole_trace_uv((motif_onset-15):(motif_onset+60)));
            end
            Trials(ii,jjc) = jj; % trials that made it to the end
            riseTime(jjc) = (pf+26-motif_onset)*6.67;
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
        %         allang{kk,ii} = tang;
        allang{kk,ii} = tangdur;
        alldf{kk,ii} = tdf;
        allRiseTime{kk,ii} = riseTime;
        alldfd{kk,ii} = tdfd;
        alluv{kk,ii} = tang;
        allrT{kk,ii} = trT;
        allrspT{kk,ii} = trspT;
        %         allangDur{kk,ii} = tangdur;
    end
end
output.Trials = Trials;
output.allRiseTime = allRiseTime;
output.allThisuv = allThisuv;
output.allThisdf = allThisdf;
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
        
        motif_frames_N = ((pf+26)-pre_frames):((pf+26)+post_frames);
        thisMotif = whole_trace(motif_frames_N);
        baseline_thisMotif = findBleachingTrend(thisMotif);
        intersection_points = InterX([motif_frames_N;thisMotif'],[motif_frames_N;baseline_thisMotif']);
        motif_onset = floor(intersection_points(1,find(intersection_points(1,:) < (pf+26),1,'last')));
        if isempty(motif_onset)
            continue;
        end
        if motif_onset <= 15
            continue;
        end
        baseline_mean = mean(whole_trace((motif_onset-15):motif_onset));
        baseline_std = std(whole_trace((motif_onset-15):motif_onset));
        
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
        Avgalluv{kk,ii} = uvval;
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
try
    output.Avgallrlat = Avgallrlat;
catch
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
temp = makeName(animalNumber{1},getMainDataFolder);
temp = makeName('pEvoked',temp);
folderName = makeName(folderName,temp);
fileName = makeName('stimTemplate.mat',folderName);
try
    temp = load(fileName);
    responseFrame = temp.stimTemplate.frameNumbers(1);
catch
    responseFrame = 28;
end


