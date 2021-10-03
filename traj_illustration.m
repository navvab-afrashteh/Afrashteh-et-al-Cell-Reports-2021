function traj_illustration()

% Two example ROIs are saved in a file. They are for two consequtive
% frames. For real data, ROIs for each frame is selected as the largest ROI
% containing values greater than 50th percentile.
FolName = 'Trajectory Illustration';
FolPath = makeName(FolName,pwd);
RoiPath = makeName('RoiSet.zip',FolPath);
ROIs = ReadImageJROI(RoiPath);
% for better display ROIs are interpolated with more points to get a
% smoother curve.
N = 300;
roi1 = ROIs{1}.mnCoordinates; roi1(end+1,:) = roi1(1,:);
roi1 = interparc(N,roi1(:,1),roi1(:,2),'spline');
x1c = mean(roi1(:,1)); y1c = mean(roi1(:,2));
x1 = roi1(:,1); y1 = roi1(:,2);
roi2 = ROIs{2}.mnCoordinates; roi2(end+1,:) = roi2(1,:);
roi2 = interparc(N,roi2(:,1),roi2(:,2),'spline');
% finding ddirection of activity flow
x2c = mean(roi2(:,1)); y2c = mean(roi2(:,2));
x2 = roi2(:,1); y2 = roi2(:,2);
% nDir = 4*5;
% flowdir = getFlowDir([roi2(:,2),roi2(:,1)], x2c, y2c,[roi1(:,2),roi1(:,1)], x1c, y1c,nDir);
flowdir = (x2c-x1c) + 1i*(y2c-y1c);
% flowdir = flowdir./abs(flowdir);

% find the intersection of the arrow connecting the centers of mass and the
% second contour 
[xe,ye] = FlowPoint([roi2(:,2),roi2(:,1)],x1c,y1c,flowdir);

% plot ROIs
W = 3; H = 3.3;
hfig = figure(1);
set(hfig,'units','inches','pos',[4,3,W,H]); clf; 
axes('units','inches','pos',[0,0,W,H]);hold on;
hp = plot(x1,y1,x2,y2,'linewidth',2);
plot(x1c,y1c,'*','color',get(hp(1),'color'));
plot(x2c,y2c,'*','color',get(hp(2),'color'));
quiver(x1c,y1c,real(flowdir),imag(flowdir),2,'linewidth',1,'color','m',...
    'MaxHeadSize',0.35)
dimX = [min([x1;x2])*0.65, max([x1;x2])*1.05]; 
dimY = [min([y1;y2])*0.95, max([y1;y2])*1.13]; 
% axis image
xlim(dimX);
ylim(dimY);
set(gca,'ydir','rev')
plot(xe,ye,'kx','linewidth',2)
plot(xe,ye,'ko','linewidth',0.5)
% add texts for explanation
fontsize = 10;
txtStr = sprintf('Activation contour\nof frame #');
text(31,34,[txtStr, '\it n-1'],'fontsize',fontsize,'rotation',0,...
    'horizontalalignment','center')
txtStr = sprintf('Center of\nmass');
text(43,38,[txtStr, '\it n-1'],'fontsize',fontsize,'rotation',0,...
    'horizontalalignment','center')
txtStr = sprintf('Activation contour\nof frame #');
text(25,51,[txtStr, '\it n'],'fontsize',fontsize,'rotation',0,...
    'horizontalalignment','center')
txtStr = sprintf('Center of\nmass');
text(36,53,[txtStr, '\it n'],'fontsize',fontsize,'rotation',0,...
    'horizontalalignment','center')
txtStr = sprintf('Vector connecting\ncenters of masses');
text(29,73,txtStr,'fontsize',fontsize,'rotation',0,...
    'horizontalalignment','center')
txtStr = sprintf('The intersection point of the\nvector and second contour');
text(44,65.5,txtStr,'fontsize',fontsize,'rotation',0,...
    'horizontalalignment','center')
color = get(hfig,'Color');
set(gca,'XColor',color,'YColor',color,'TickDir','out')
set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
%% save as pdf
pdfFilePath = makeName('Trajectory Illustration.pdf',FolPath);
save2pdf(pdfFilePath,hfig);
end

function flowdir = getFlowDir(B1,x1,y1,B2,x2,y2,nDir)
% make sure that first and last points are the same and contour is closed
if ~(B1(1,1)==B1(end,1) && B1(1,2)==B1(end,2))
    B1(end+1,1)= B1(1,1);
    B1(end,2)= B1(1,2);
end
if ~(B2(1,1)==B2(end,1) && B2(1,2)==B2(end,2))
    B2(end+1,1)= B2(1,1);
    B2(end,2)= B2(1,2);
end

alltheta = linspace(0,2*pi,nDir);
alldir = cos(alltheta)+1i*sin(alltheta);
xe1 = zeros(nDir-1,1); ye1 = zeros(nDir-1,1); idx1 = zeros(nDir-1,1);
xe2 = zeros(nDir-1,1); ye2 = zeros(nDir-1,1); idx2 = zeros(nDir-1,1);
for k = 1:nDir-1
    if k==11
        nn=1;
    end
    [xe1(k),ye1(k),idx1(k)] = FlowPoint(B1,x1,y1,alldir(k));
    [xe2(k),ye2(k),idx2(k)] = FlowPoint(B2,x2,y2,alldir(k));
end
[idx1, I1] = sort(idx1);
[idx2, I2] = sort(idx2);

x1dir = zeros(nDir-1,1); y1dir = zeros(nDir-1,1); x2dir = zeros(nDir-1,1); y2dir = zeros(nDir-1,1);
for k = 1:nDir-2
    if ~isnan(idx1(k)) && ~isnan(idx1(k+1))
        x1dir(I1(k)) = mean(B1(idx1(k):idx1(k+1),2));
        y1dir(I1(k)) = mean(B1(idx1(k):idx1(k+1),1));
    end
    if ~isnan(idx2(k)) && ~isnan(idx2(k+1))
        x2dir(I2(k)) = mean(B2(idx2(k):idx2(k+1),2));
        y2dir(I2(k)) = mean(B2(idx2(k):idx2(k+1),1));
    end
end
for k = nDir-1
    if ~isnan(idx1(k)) && ~isnan(idx1(1))
        x1dir(I1(k)) = mean([B1(1:idx1(1),2); B1(idx1(k):end,2)]);
        y1dir(I1(k)) = mean([B1(1:idx1(1),1); B1(idx1(k):end,1)]);
    end
    if ~isnan(idx2(k)) && ~isnan(idx2(1))
        x2dir(I2(k)) = mean([B2(1:idx2(1),2); B2(idx2(k):end,2)]);
        y2dir(I2(k)) = mean([B2(1:idx2(1),1); B2(idx2(k):end,1)]);
    end
end
orig1 = x1 - 1i*y1;
vec1 = x1dir + 1i*y1dir - orig1*0;
orig2 = x2 - 1i*y2;
vec2 = x2dir + 1i*y2dir - orig2*0;
vec = vec1 - vec2;
flowdir = sum(vec);
flowdir = flowdir/abs(flowdir);

% figure(1); clf; hold on;
% plot(B1(:,2),B1(:,1),'g',B2(:,2),B2(:,1),'b','linewidth',2);
% plot(x1,y1,'*g');
% % plot(B2(:,2),B2(:,1),'b','linewidth',2);
% plot(x2,y2,'*b');
% plot(x1dir,y1dir,'r.')
% plot(x2dir,y2dir,'m.')
% quiver(x2,y2,real(flowdir),imag(flowdir),50)
% quiver(x2dir,y2dir,real(vec),imag(vec),5)
% 
% dim1 = 128; dim2 = 128;
% axis image
% xlim([1,dim2]);
% ylim([1,dim1]);
% set(gca,'ydir','rev')

n=0;
end
function [xe,ye,idx] = FlowPoint(B,xe,ye,avguv)
if isempty(B)
    idx = nan;
    xe = nan;
    ye = nan;
    return;
end
if ~(B(1,1)==B(end,1) && B(1,2)==B(end,2))
    B(end+1,1)= B(1,1);
    B(end,2)= B(1,2);
end
m = imag(avguv)/real(avguv);
y = m*(B(:,2)-xe)+ye;
b = B(:,1)-y;
b = b(1:end-1).*b(2:end);
idx = find(b<=0);
b = [B(idx,2), B(idx,1)];
% b = [B(idx,2), y(idx)];
b = b - repmat([xe,ye],size(b,1),1);
b = b * [real(avguv); imag(avguv)];
b = b > 0;
idx = idx(b);
B = B(idx,:);
% B = [y(idx),B(idx,2)];
d = B - repmat([ye,xe],size(B,1),1);
d = sqrt(sum(d.^2,2));
[~,idx1] = max(d);
ye = B(idx1,1);
xe = B(idx1,2);
idx = idx(idx1);
if isempty(idx)
    idx = nan;
    xe = nan;
    ye = nan;
end
end
