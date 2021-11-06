function saveEvokedTemplate_GFP_BPF(animalNumber,varargin)

GFP = 0;
if nargin>1
    mainDataFolder = getMainDataFolderYFP;
    GFP = 1;
else
    mainDataFolder = getMainDataFolder;
end
filtOptions.Fstop1 = 0.25; filtOptions.Fpass1 = 0.5;
filtOptions.Fpass2 = 6; filtOptions.Fstop2 = 6.5;
filtOptions.dur = 5;
if nargin > 2
    filtOptions = varargin{2};
end
Fpass1 = filtOptions.Fpass1; Fpass2 = filtOptions.Fpass2;
Fpass1 = round(Fpass1*100);
Fpass2 = round(Fpass2);
stim = {'BC'};
quit = 0;
%%
nFTemplate = 10;
Cmap = load('selfmade3.txt');
Cmap = Cmap(:,2:4)/255;
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
    for ss = 1%:length(folders)
        if quit
            break;
        end
        folderName = makeName(folders(ss).name,peDataFolder);
        fileName = sprintf('%s\\ImgSeq_BPF_0%dHz_%dHz.mat',folderName, Fpass1, Fpass2);
        temp = load(fileName);
        ImgSeq = temp.ImgSeq;
        ImgSeq = ImgSeq(:,:,91-26:91-26+100-1);
        ImgSeq = applyMask(ImgSeq,mask);
        
        ii = 25;
        baseFrames = 10;
        f1 = 26-baseFrames+1;
        mImgSeq = ImgSeq(:,:,f1:ii+50-1);
        mImgSeq = getMaskedValues(mImgSeq,mask);
        mImgSeq = mean(mImgSeq,2);
        minI = 0;
        maxI = max(mImgSeq)/2;
        fileName = sprintf('%s\\stimTemplate_BPF_0%dHz_%dHz.mat',folderName, Fpass1, Fpass2);
        if exist(fileName,'file')
            load(fileName)
            stimTemplate;
            continue;
        end
        
        while 1
            figure(2);clf;plot(mImgSeq); hold on; plot((baseFrames)*[1 1],[min(mImgSeq) max(mImgSeq)])
            figure(1);clf;
            imagesc(ImgSeq(:,:,ii),[minI maxI]);
            colormap(Cmap);
            hold on;
            axis equal; axis off;
            titleText = sprintf('%s - %s - %d',animalNumber{an},folders(ss).name,ii);
            title(titleText);
            ch = getkey;
            if ch == 29 % '-->' key for going to next frame
                ii = ii + 1;
            end
            if ch == 28 % '<--' key for going to previous frame
                ii = ii - 1;
            end
            if ch == 115 % 's' for saving template frames and frame numbers
                display('Saving file ...');
                stimTemplate.frames = ImgSeq(:,:,ii:(ii+nFTemplate-1));
                stimTemplate.frameNumbers = ii:(ii+nFTemplate-1);
                stimTemplate.frameNumbers = stimTemplate.frameNumbers + 64;
                fileName = sprintf('%s\\stimTemplate_BPF_0%dHz_%dHz.mat',folderName, Fpass1, Fpass2);
                save(fileName,'stimTemplate');
                display('File Saved');
                break;
            end
            if ch == 27 % 'Esc' key to quit the current animal
                break;
            end
            if ch == 119 % 'w' key to open the folder
                winopen(folderName);
            end
            if ch == 113 % 'q' key to completly exit from this script
                quit = 1;
                break;
            end
        end
        
    end
end
