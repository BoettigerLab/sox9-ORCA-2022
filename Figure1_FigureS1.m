%% data prep
clear all
analysisFolder2 = 'Z:\Liang-Fu\2021-03-09_h9a7_WT_Sox9_p2\DNA_Expt\Analysis\'

[polys2,maps2,spotData2] = CombineAllFits(analysisFolder2, 'byFOV',true);

ESCFov = [7:21, 38:48];
CNCCFov = [1:5,50:55];

ESCmaps = cat(3,maps2{ESCFov});
CNCCmaps = cat(3,maps2{CNCCFov});

ESCpolys = cat(3,polys2{ESCFov});
CNCCpolys = cat(3,polys2{CNCCFov});

[cMap,nObs] = ContactFrac(ESCmaps,'threshold',250);
figure(1); clf; imagesc(nObs);

[cMap,nObs] = ContactFrac(CNCCmaps,'threshold',250);
figure(2); clf; imagesc(nObs);

badHybes = [16,32,44,53,54,67,81:86];

goodESCmaps = ESCmaps;
goodESCmaps(badHybes,:,:)= NaN;
goodESCmaps(:,badHybes,:)= NaN;

goodCNCCmaps = CNCCmaps;
goodCNCCmaps(badHybes,:,:)= NaN;
goodCNCCmaps(:,badHybes,:)= NaN;

goodESCpolys = ESCpolys;
goodESCpolys (badHybes, :, :) = NaN;

goodCNCCpolys = CNCCpolys;
goodCNCCpolys (badHybes, :, :) = NaN;

ESCrpc = ReadsPerCell(goodESCmaps);
CNCCrpc = ReadsPerCell(goodCNCCmaps);

figure(7); clf; 
subplot (1,2,1);
hist(ESCrpc,20);
title('ESC');

subplot (1,2,2);
hist(CNCCrpc,20);
title('CNCC');

ESChighDetect = ESCrpc>38;
CNCChighDetect = CNCCrpc>38;

%% Distance and contact map
ESCmedianDistance = nanmedian(goodESCmaps(:,:,ESChighDetect),3);
CNCCmedianDistance = nanmedian(goodCNCCmaps(:,:,CNCChighDetect),3);

ESC_Distance_interp = InterpMapNans(ESCmedianDistance);
CNCC_Distance_interp = InterpMapNans(CNCCmedianDistance);

ESC_Contact_interp = InterpMapNans( ContactFrac(goodESCmaps (:,:,ESChighDetect),'threshold',250));
CNCC_Contact_interp = InterpMapNans( ContactFrac(goodCNCCmaps (:,:,CNCChighDetect),'threshold',250));

%% Figure 1E

figure(10); clf; 
% ESC contact frequency 
subplot (1,2,1);
imagesc(ESC_Contact_interp); caxis([.05,.4]); axis([1 80 1 80]);
title(['ESC Contact Interp, n=' num2str(size(goodESCmaps (:,:,ESChighDetect),3))]);
colorbar
GetColorMap ('whiteToRed')
% CNCC contact frequency 
subplot (1,2,2);
imagesc(CNCC_Contact_interp) ; caxis ([.05,.4]); axis([1 80 1 80]);
title(['CNCC Contact Interp, n=' num2str(size(goodCNCCmaps (:,:,CNCChighDetect),3))])
colorbar
GetColorMap ('whiteToRed')

%% Figure 1F

% CNCC-ESC contact frequency subtraction map
figure(31); clf;
imagesc(CNCC_Contact_interp - ESC_Contact_interp); caxis([-.1,.1]); axis([1 80 1 80]);
title(['CNCC - ESC Contact Subtraction Interp']);
colorbar
GetColorMap ('BlueWhiteRed')
set(gcf, 'position',[10,10,900,800]);

% Bootstrap for high confidence map CI>95%
ESCdata3= goodESCmaps(:,:,ESChighDetect);
CNCCdata3 = goodCNCCmaps(:,:,CNCChighDetect);


pCNCCB =zeros(86,86);

for i=1:86;
    for k=1:86;
ESCdata1 =squeeze(ESCdata3(i,k,:));
ESCdata11=~isnan(ESCdata1);
ESCdata111=ESCdata1(ESCdata11);

CNCCdata1 =squeeze(CNCCdata3(i,k,:));
CNCCdata11=~isnan(CNCCdata1);
CNCCdata111=CNCCdata1(CNCCdata11);

if  nansum(ESCdata1)==0
            pCNCCB(i,k)=0
else
    
bootstatESC = bootstrp(1000,@median,ESCdata111);
bootstatCNCC = bootstrp(1000,@median,CNCCdata111);
ESC_CNCC= bootstatCNCC-bootstatESC;

alpha=.05;
CI = prctile(ESC_CNCC,[100*alpha/2,100*(1-alpha/2)]);

%Hypothesis test: Does the confidence interval cover zero?
H = CI(1)>0 | CI(2)<0;

pCNCCB(i,k)=H
end
    end
end;

CNCC1S = CNCC_Contact_interp - ESC_Contact_interp;
CNCC1SS = CNCC1S.*pCNCCB;

% high confidence map. only showing pixels that CI>95%
figure(30); clf;
imagesc(CNCC1SS,'AlphaData',pCNCCB); caxis([-.1,.1]); axis([1 80 1 80]);
title(['CNCC - ESC Contact Subtraction Interp Bootstrape']);
colorbar
GetColorMap ('BlueWhiteRed');
set(gca, 'color',[.85 .85 .85]);
set(gcf, 'position',[10,10,900,800]);

%% Figure 1I
% Radius of Gyration 
CNCCdata = goodCNCCpolys (:, :, CNCChighDetect);
ESCdata = goodESCpolys (:, :, ESChighDetect);

% Radius of Gyration
ESCRg = zeros(1063,1);
for i= 1:1063;
    ESCRg(i)= RadiusOfGyration(ESCdata(11:52,:,i));
end;
 
CNCCRg = zeros(1632,1);
for i= 1:1632;
    CNCCRg(i)= RadiusOfGyration(CNCCdata(11:52,:,i));
end;

% Wilcoxon rank sum test for equal medians
[p,h] = ranksum(CNCCRg,ESCRg);

figure(2); clf;
cdfplot (ESCRg);hold on;
cdfplot (CNCCRg);hold on;
legend('ESC','CNCC');
xlabel ('Radius of Gyration (nm)');
ylabel ('% of obs');
set(gcf, 'position',[10,10,900,800]);

%% Figure S1B

figure(9); clf; 
subplot (1,2,1);
imagesc(ESC_Distance_interp); caxis([100 500]); axis([1 80 1 80]);
title(['ESC Distance Interp, n=' num2str(size(goodESCmaps (:,:,ESChighDetect),3))]);
colorbar
GetColorMap ('redToWhite')

subplot (1,2,2);
imagesc(CNCC_Distance_interp); caxis([100 500]); axis([1 80 1 80]);
title(['CNCC Distance Interp, n=' num2str(size(goodCNCCmaps (:,:,CNCChighDetect),3))]);
colorbar
GetColorMap ('redToWhite')

%% Figure S1C

figure(22); clf;
imagesc (CNCC_Distance_interp-ESC_Distance_interp); caxis([-200,200]); axis([1 80 1 80]);
title(['CNCC-ESC Distance'])
colorbar
GetColorMap ('RedWhiteBlue')
set(gcf, 'position',[10,10,900,800]);

%% Figure S1D

jboxFolder = 'U:\GenomeData\JuiceboxExport\';
exportName = 'H1_hESC_Dekker_Sox9_balanced_5kb'
mapRes = 1E4;
h19_locus = 'chr17:68000000-72000000';
SOX9 = ReadJuiceboxMatrix([jboxFolder,exportName],...  
          'locus',h19_locus,'mapRes',mapRes,'displayRes',mapRes);
figure(110); clf; imagesc(SOX9, [0 25]);
colorbar
GetColorMap ('whiteToRed')
set(gcf, 'position',[10,10,900,800]);

%% Figure S1E

% MDS
ESCMDS = ESC_Distance_interp(2:79,2:79)-diag(diag(ESC_Distance_interp(2:79,2:79)));
ESCcurrPoly = mdscale (ESCMDS,3);

CNCCMDS = CNCC_Distance_interp(2:79,2:79)-diag(diag(CNCC_Distance_interp(2:79,2:79)));
CNCCcurrPoly = mdscale (CNCCMDS,3);

% 3D tube plots

CTCF = [11,49];
EC145 = [12];
EC135=[14];
EC125=[16];
Sox9 = [41];

CNCCCTCF =CNCCcurrPoly(CTCF,:);
ESCCTCF=ESCcurrPoly(CTCF,:);

CNCCSox93D = CNCCcurrPoly(Sox9,:);
ESCSox93D=ESCcurrPoly(Sox9,:);

CNCCEC1253D = CNCCcurrPoly (EC125,:);
ESCEC1253D = ESCcurrPoly (EC125,:);

CNCCEC1353D = CNCCcurrPoly (EC135,:);
ESCEC1353D = ESCcurrPoly (EC135,:);

CNCCEC1453D = CNCCcurrPoly (EC145,:);
ESCEC1453D = ESCcurrPoly (EC145,:);

figure(22);cla; title('ESC');
PlotPolymerTube(ESCcurrPoly,'tubeRadius',6,'showSpheres',false,...
    'colormap','hsv','alpha',.6,'method','pchip','center',false,'lightOn',true); grid off;
colormap(flipud(parula));
set(gca,'color','w'); hold on; 
PlotSpheres(ESCCTCF,'color',[0 0 0],'r',12); hold on;
PlotSpheres(ESCSox93D,'color',[1 0 0],'r',12); hold on;
PlotSpheres(ESCEC1253D,'color',[0 0 1],'r',12); hold on;
PlotSpheres(ESCEC1353D,'color',[0 1 0],'r',12); hold on;
PlotSpheres(ESCEC1453D,'color',[1,0.4,0.1],'r',12); hold on;
set(gcf, 'position',[10,10,900,800]);

figure(21);cla; title('CNCC');

PlotPolymerTube(CNCCcurrPoly,'tubeRadius',6,'showSpheres',false,...
    'colormap','hsv','alpha',.6,'method','pchip','center',false,'lightOn',true); grid off;
colormap(flipud(parula));
set(gca,'color','w'); hold on; 
PlotSpheres(CNCCCTCF,'color',[0 0 0],'r',12); hold on;
PlotSpheres(CNCCSox93D,'color',[1 0 0],'r',12); hold on;
PlotSpheres(CNCCEC1253D,'color',[0 0 1],'r',12); hold on;
PlotSpheres(CNCCEC1353D,'color',[0 1 0],'r',12); hold on;
PlotSpheres(CNCCEC1453D,'color',[1,0.4,0.1],'r',12); hold on;
set(gcf, 'position',[10,10,900,800]);

