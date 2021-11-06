%%
clc
mainDataFolder = getMainDataFolderYFP;
animalNumber =  {'190533'; '190545'};
Fpass1 = 0.5;
Fpass2 = 6;
Fpass1 = round(Fpass1*100);
%%
quit = 0;
FNS = {'ROI_BC'};
for an = 1:length(animalNumber)
    if quit
        break;
    end
    dataFolder = makeName(animalNumber{an},mainDataFolder);
    peDataFolder = makeName('pEvoked',dataFolder);
    ignoreFolders = [];
    ofolders = dir(peDataFolder);
    folders = cleanFolderList(ofolders,ignoreFolders);
    mask = getMask(animalNumber(an),'glut',0,'yfp',1,'ec',1); 
    ss = 1;
    if quit
        break;
    end
    folderName = makeName(folders(ss).name,peDataFolder);
    fileName = sprintf('%s\\stimTemplate_BPF_0%dHz_%dHz.mat',folderName, Fpass1, Fpass2);
    load(fileName); stimTemplate.frameNumbers
    fileName = sprintf('%s\\ImgSeq_BPF_0%dHz_%dHz.mat',folderName, Fpass1, Fpass2);
    temp = load(fileName);
    ImgSeq = temp.ImgSeq;
    ImgSeq = applyMask(ImgSeq,mask);
    avg = mean(getMaskedValues(ImgSeq,mask),2);
    ii = 90;
    mImgSeq = ImgSeq;
    mImgSeq = getMaskedValues(mImgSeq,mask);
    mImgSeq = mean(mImgSeq,2);
    minI = min(mImgSeq)*0;
    maxI = max(mImgSeq)/2;
    multiplier = 1.5;
    while 1
        figure(2);clf;set(gcf,'position',[616 630 515 412]);plot(mImgSeq);hold on;
        plot((ii),mImgSeq(ii),'or');
        figure(1);set(gcf,'position',[60 630 515 412]);clf;
        imagesc(ImgSeq(:,:,ii),[minI maxI]);colormap jet;
        maskVals = getMaskedValues(ImgSeq(:,:,ii),mask);
        threshold = mean(maskVals) + multiplier*std(maskVals);
        BW = im2bw(ImgSeq(:,:,ii),max(0,threshold));
        B = bwboundaries(BW);
        s = regionprops(BW,'centroid');
        centroids = cat(1, s.Centroid);
        hold on;
        for jj = 1:size(centroids,1)
            plot(centroids(jj,1),centroids(jj,2),'.','markersize',15,'color','c');
            text(centroids(jj,1),centroids(jj,2),sprintf('%d',jj));
            plot(B{jj}(:,2),B{jj}(:,1),'k','linewidth',1.5);
        end
        
        hold on;
        axis equal; axis off;
        titleText = sprintf('%s - %s - %d',animalNumber{an},folders(ss).name,ii);
        title(titleText,'interpreter','none');
        figure(3);set(gcf,'position',[616 100 515 412]);clf;
        contour(ImgSeq(:,:,ii));
        [C, Clevels] = contourcs(ImgSeq(:,:,ii));
        set(gca,'YDir','Reverse');
        
        ch = getkey;
        if ch == 29
            ii = ii + 1;
        end
        if ch == 28
            ii = ii - 1;
        end
        if ch == 115 %s
            display('Saving file ...');
            fileName = makeName(sprintf('%s.mat',FNS{ss}),dataFolder);
            cNum = input('Select contour \n');
            ROI = poly2mask(B{cNum}(:,2),B{cNum}(:,1),size(ImgSeq,1),size(ImgSeq,2));
            save(fileName,'ROI');
            display('File Saved');
            break;
        end
        if ch == 27 %Esc
            break;
        end
        if ch == 119 %w
            winopen(folderName);
        end
        if ch == 113 %q
            quit = 1;
            break;
        end
        if ch == 109 %m
            multiplier
            multiplier = input('Multiplier = ? \n');
        end
        if ch == 121 %y
            maxI
            maxI = input('maxI = ? \n');
        end
    end
end
n=0;
