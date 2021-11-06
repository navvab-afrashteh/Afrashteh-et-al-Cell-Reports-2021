function R1_Spon_SpatioTemp_YFP

mainDataFolder = getMainDataFolderYFP;
% All animals 
animals =  {'190533'; '190545'};

calc = 0;
if calc
    calcST(animals)
end

plotST

function calcST(animals)
mainDataFolder = getMainDataFolderYFP;
for an = 1:2
    anNum = animals{an};
    anFolder = sprintf('%s\\%s',mainDataFolder,anNum);
    sponFiles = dir(anFolder);
    idx = strfind({sponFiles(:).name},'.mat');
    idx = ~(cellfun(@isempty,idx));
    sponFiles = sponFiles(idx);
    mask = imread(sprintf('%s\\Mask.tif',anFolder));
    Mask(:,:,an) = mask;
    for ii = 1%:length(sponFiles)
        sponPath = sprintf('%s\\%s',sponFiles(ii).folder,sponFiles(ii).name);
        load(sponPath)
        
        m(:,:,an,ii) = mean(DF_F0,3);    % spatial mean
        sd(:,:,an,ii) = std(DF_F0,[],3); % spatial standard deviation
        values = getMaskedValues (DF_F0,mask);
        sig(:,an,ii) = mean(values,2);
    end
end
save('R1_Spon_SpatioTemp_YFP.mat','m','sd','sig','Mask')

function plotST
cmap = load('selfmade3.txt')/255; cmap = cmap(:,2:4);
load('R1_Spon_SpatioTemp_YFP.mat','m','sd','sig','Mask')
for an = 1:2
    pdfFileName = sprintf('%s_animal_%d.pdf',mfilename,an);
    pdfFileName = makeName(pdfFileName,getpdffolder)
    if exist(pdfFileName,'file')
        delete(pdfFileName);
    end
    if exist(makeName('temp.pdf',getpdffolder),'file')
        delete(makeName('temp.pdf',getpdffolder));
    end

    mask = Mask(:,:,an);
    nFiles = size(m,4);
    for ii = 1:nFiles        
        mtemp = m(:,:,an,ii).*mask;    % spatial mean
        sdtemp = sd(:,:,an,ii).*mask; % spatial standard deviation
        
        ff = makeFigureWindow__one_axes_only(2,[1 4 2 2],[0.05 0.05 0.9 0.9]);
        imagesc(mtemp,[0,0.1])
        colormap(cmap)
        axis equal
        set(ff.ha,'YDir','reverse','XTick',[],'YTick',[])
        colorbar
        title(sprintf('Mean Spon animal #%d',an))
        save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
        append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));
        
        ff = makeFigureWindow__one_axes_only(2,[1 4 2 2],[0.05 0.05 0.9 0.9]);
        imagesc(sdtemp,[0.4,0.8])
        colormap(cmap)
        axis equal
        set(ff.ha,'YDir','reverse','XTick',[],'YTick',[])
        colorbar
        title(sprintf('SD Spon animal #%d',an))
        save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
        append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));
    end
end