% location of ORCA data (update to your local filepaths)
NAS02_Vol3 = 'Z:\';
dataFolder = [NAS02_Vol3,'\Liang-Fu\2021-03-09_h9a7_WT_Sox9_p2\DNA_Expt\Analysis\'];

% location of simulation data (update to your file paths. Code to reproduce these simulations is included, see Readme).
simFolders ={'T:/2022-02-24_LiangFu/StripeSim_UniformShortLived\',...
    'T:/2022-02-24_LiangFu/StripeSim_ReelIn3\'};

ESCFov = [7:21, 38:48];
CNCCFov = [1:5,50:55];
[polysSox9,mapsSox9] = CombineAllFits(dataFolder,'byFOV',true);
badHybes = [16,53,54];


%% Distance and contact map
escMaps = cat(3,mapsSox9{ESCFov});
escContact = ContactFrac(escMaps,'threshold',200);
escDist = nanmedian(escMaps,3);
cncMaps = cat(3,mapsSox9{CNCCFov});
cncContact = ContactFrac(cncMaps,'threshold',200);
cncDist = nanmedian(cncMaps,3);

%% sims data prep
distSims = cell(2,1);
polyStep = 10; % 2
mapStep = 1;
timeStep = 1;
% v2


    repFolders = FindFiles([simFolders{1},'rep*'],'onlyFolders',true,'fullPath',false);
    nRuns = length(repFolders);
    distMap = cell(nRuns,1);
    polyDat = cell(nRuns,1);
     for r=1:nRuns  % this is now a loop over reps. 
        folder = [simFolders{e},repFolders{r}];  % filesep,
        % disp([folder,'  now loading data']);
            if r==1
                disp([folder,'  now loading data for n=',num2str(nRuns)]);
                parsTable = readtable([folder,'simPars.txt']); % all reps have same pars
                ctcfSites = str2num(parsTable{1,5}{1}); %#ok<ST2NM>
                chrChains = parsTable.replicates(1);
                lenPoly = parsTable.monomers(1);
                disp(parsTable(:,1:4));
            end
        
          % %   loads a cell array of 1xN-rpts, each monomers x time
          try
[polyReps1{r},mapReps1{r},lefData{r},loopData{r},imPars{r}] = LoadOpenPolySim(folder,...
                                'polyStep',10,'timeStep',1,'mapStep',1);
catch                
              continue
          end
     end
     

 
for e=1:2
    repFolders = FindFiles([simFolders{e},'rep*'],'onlyFolders',true,'fullPath',false);
    nRuns = length(repFolders);
    distMap = cell(nRuns,1);
    polyDat = cell(nRuns,1);
     for r=1:nRuns  % this is now a loop over reps. 
        folder = [simFolders{e},repFolders{r}];  % filesep,
        % disp([folder,'  now loading data']);
            if r==1
                disp([folder,'  now loading data for n=',num2str(nRuns)]);
                parsTable = readtable([folder,'simPars.txt']); % all reps have same pars
                ctcfSites = str2num(parsTable{1,5}{1}); %#ok<ST2NM>
                chrChains = parsTable.replicates(1);
                lenPoly = parsTable.monomers(1);
                disp(parsTable(:,1:4));
            end
        
          % %   loads a cell array of 1xN-rpts, each monomers x time
          try
            [polyDat{r},distMap{r}]  = LoadPolySimDynamics(folder,...
                                'TADs',ctcfSites,...
                                'lenPoly',lenPoly,...
                                'rptPoly',chrChains,...
                                'polyStep',polyStep,...
                                'mapStep',mapStep,...
                                'timeStep',timeStep,...   %  50
                                'threshold',10,...
                                'showPlots',false,...
                                'contMapFig',0,'distMapFig',0,...
                                'maxBlocks',inf);
                            
            catch                
              continue
          end
     end
    distMaps = cat(1,distMap{:}); % cat chrs (only 1)
    distMaps = cat(3,distMaps{:}); % cat replicates
    polyStk = cat(1,polyDat{:});
    polyStk1{e} = cat(3,polyStk{:});
    maps = cat(3,distMaps);
    distSims{e} = maps;
    dmap = nanmedian(maps(:,:,1:1000),3);

end
%% Figure 5B
t =4;
figure(3); clf; 
subplot(1,3,3); imagesc((cMap1(14:end,14:end))); colorbar; caxis([0,.04]);
subplot(1,3,1); imagesc((cMap2(14:end,14:end))); colorbar; caxis([0,.007]);
aMap = InterpMapNans(ContactFrac(cncMaps(:,:,:),'threshold',200),'badHybes',[16]);
subplot(1,3,2); imagesc((aMap(12:52,12:52))); colorbar;  caxis([0.1,.3]);
GetColorMap('whiteToRed'); 

%% Figure 5C
% ORCA data sort
dM = cncMaps(1:52,1:52,:);
strA = squeeze(dM(42,:,:));
strC = strA < 150;
a = sum(cumsum(strC,1));
[~,id] = sort(a);

[~,first]=max(strC,[],1,'linear');
[rw,cl] = ind2sub(size(strC),first);
[~,id] = sort(rw(1:size(strC,2)));
datSort = strC(:,id);

% 'reel-in' sim sort
strA = squeeze(distSims{2}(81,:,:));
strC = strA <  5;
[~,first]=max(strC,[],1,'linear');
[rw,cl] = ind2sub(size(strC),first);
[~,id] = sort(rw(1:size(strC,2)));
reelSort = strC(:,id);

% 'Multi-loop' sim sort
strA = squeeze(distSims{1}(81,:,:));
strC = strA < 5; 
strCC =sum(strC,1);
[~,first]=max(strC,[],1,'linear');
[rw,cl] = ind2sub(size(strC),first);
[~,id] = sort(rw(1:size(strC,2)));
multiSort = strC(:,id);

figure(1); clf; 
subplot(1,3,1); imagesc(reelSort'); title('Reel-in');  ylim([0,6000]);
subplot(1,3,2); imagesc(datSort'); title('Data');  ylim([0,1000]);
subplot(1,3,3); imagesc(multiSort'); title('Multi-loop');  ylim([0,2700]);
GetColorMap ('whiteToRed');
set(gcf,'color','w','position',[10,10,900,800]);

%% Figure 5D and S5D
% ORCA 15-42 and 35-42
CNCCno42 = isnan(squeeze(CNCCpolys(42,1,:)));
sum(CNCCno42)
CNCC42= squeeze (CNCCmaps(:,42,CNCChighDetect&~CNCCno42));
CNCC42 =  rot90(CNCC42);
CNCC42S = CNCC42<150;
CNCC42A=  CNCC42S(:,:);
[~,pt15] = sort(CNCC42A(:,15),'descend');
[~,pt35] = sort(CNCC42A(:,35),'descend');
CNCC4215 =CNCC42S(pt15,:);
CNCC4235 =CNCC42S(pt35,:);

figure(1);clf;
subplot (1,2,1);
imagesc(CNCC4215);caxis([0 1]);axis([12 50 1 1103]);
title('CNCC 15 Threshold =150');
subplot(1,2,2);
imagesc(CNCC4235);caxis([0 1]);axis([12 50 1 1103]);
title('CNCC 35 Threshold =150');
GetColorMap ('whiteToRed');

% 'reel-in' 21-81 and 67-81
reel = squeeze(distSims{2}(81,:,:));
reel1 = reel <  5;
reel2 = rot90(reel1);
[~,pt21] = sort(reel2(:,21),'descend');
[~,pt67] = sort(reel2(:,67),'descend');
reel21 = reel2(pt21,:);
reel67 = reel2(pt67,:);

figure(77);clf;
subplot (1,2,1);
imagesc(reel21);caxis([0 1]);axis([15 100 1 5000]);
title('reel in 21');
subplot (1,2,2);
imagesc(reel67);caxis([0 1]);axis([15 100 1 5000]);
title('reel in 67');

% 'Multi-loop' 21-81 and 67-81
mul = squeeze(distSims{1}(81,:,:));
mul1 = mul <  5;
mul2 = rot90(mul1);
[~,pt21] = sort(mul2(:,21),'descend');
[~,pt67] = sort(mul2(:,67),'descend');
mul21 = mul2(pt21,:);
mul67 = mul2(pt67,:);

figure(666);clf;
subplot (1,2,1);
imagesc(mul21);caxis([0 1]);axis([15 100 1 3400]);
title('mul 21');
subplot (1,2,2);
imagesc(mul67);caxis([0 1]);axis([15 100 1 3400]);
title('mul 67');
GetColorMap ('whiteToRed');

% sim 67-81
t=8;
a = 67; 
isLoop = squeeze(distSims{1}(81,a,:)<t);
aMap = ContactFrac(distSims{1}(:,:,isLoop),'threshold',t);
figure(2); clf; 
subplot(2,3,3); imagesc(aMap); colorbar; axis([15 100 15 100]); caxis([0,.42]);
isLoop = squeeze(distSims{2}(81,a,:)<t);
aMap = ContactFrac(distSims{2}(:,:,isLoop),'threshold',t);
figure(2); 
subplot(2,3,1); imagesc(aMap); colorbar; axis([15 100 15 100]);caxis([0,.42]);
set(gcf,'color','w')

% sim 21-81
a = 21; 
isLoop = squeeze(distSims{1}(81,a,:)<t);
aMap = ContactFrac(distSims{1}(:,:,isLoop),'threshold',t);
figure(2); subplot(2,3,6); imagesc(aMap); colorbar; axis([15 100 15 100]); caxis([0,.42]);
isLoop = squeeze(distSims{2}(81,a,:)<t);
aMap = ContactFrac(distSims{2}(:,:,isLoop),'threshold',t);
figure(2); subplot(2,3,4); imagesc(aMap); colorbar; axis([15 100 15 100]); caxis([0,.42]);

% ORCA data 35-42
a = 35;  
t = 200;
isLoop = squeeze(cncMaps(42,a,:)<t);
aMap = ContactFrac(cncMaps(1:52,1:52,isLoop),'threshold',t);
figure(2); subplot(2,3,2); imagesc(aMap); colorbar; axis([12 50 12 50]);caxis([0.1,.45]);

% ORCA data 15-42
a = 15;
isLoop = squeeze(cncMaps(42,a,:)<t);
aMap = ContactFrac(cncMaps(1:52,1:52,isLoop),'threshold',t);
figure(2); subplot(2,3,5); imagesc(aMap); colorbar; axis([12 50 12 50]);caxis([0.1,.45]);
GetColorMap('whiteToRed')
set(gcf,'color','w')

%% Figure 5E

mulpoly =polyStk1{1};
reelpoly = polyStk1{2};

% multi-loop
mulDis = zeros(100,3400);

for i= 1:3400;
    centerdata= mulpoly(:,:,i);
    rm = mean(centerdata,1,'omitnan');
    mapp=mulpoly(:,:,i);
    mulDis(:,i)= sqrt((mapp(:,1)-rm(1)).^2 + (mapp(:,2)-rm(2)).^2 + (mapp(:,3)-rm(3)).^2);

end;

mulmedianDis = median(mulDis,2,'omitnan' );

%reel-in
reelDis = zeros(100,9600);

for i= 1:9600;
    centerdata= reelpoly(:,:,i);
    rm = mean(centerdata,1,'omitnan');
    mapp=reelpoly(:,:,i);
    reelDis(:,i)= sqrt((mapp(:,1)-rm(1)).^2 + (mapp(:,2)-rm(2)).^2 + (mapp(:,3)-rm(3)).^2);

end;

reelmedianDis = median(reelDis,2,'omitnan' );

% ORCA data
CNCCdata = goodCNCCpolys (:, :, CNCChighDetect);
CNCCDis = zeros(86,1632);

for i= 1:1632;
    centerdata= CNCCdata(11:52,:,i);
    rm = mean(centerdata,1,'omitnan');
    mapp=CNCCdata(:,:,i);
    CNCCDis(:,i)= sqrt((mapp(:,1)-rm(1)).^2 + (mapp(:,2)-rm(2)).^2 + (mapp(:,3)-rm(3)).^2);

end;

CNCCmedianDis = median(CNCCDis,2,'omitnan' );
CNCCmedianDis2=fillmissing(CNCCmedianDis,'movmean',3);

figure(284); clf;
subplot(3,1,1);
plot(CNCCmedianDis2(11:52));axis([1 42 150 500]);hold on;
ylabel ('Median distance from center/Rg');
title('Median distance to center of SOX9 TAD');
subplot(3,1,2);
plot(mulmedianDis(15:90));axis([1 75 10 20]);hold on;
ylabel ('Median distance from center/Rg');
title('multiloop');
subplot(3,1,3);
plot(reelmedianDis(15:90));axis([1 75 15 30]);hold on;
ylabel ('Median distance from center/Rg');
title('reel in');

%% Figure 5F
CNCCtad= goodCNCCmaps(1:52,1:52,CNCChighDetect);
CNCCtad1 = median(CNCCtad,3, 'omitnan' );
CNCCtad2 = median(CNCCtad1,2, 'omitnan' );
CNCCtad3 = fillmissing(CNCCtad2,'movmean',3);

reelmaps = distSims{2};
reelmap1 = median(reelmaps,3, 'omitnan' );
reelmap2 = median(reelmap1,2, 'omitnan' );

mulmaps = distSims{1};
mulmap1 = median(mulmaps,3, 'omitnan' );
mulmap2 = median(mulmap1,2, 'omitnan' );

figure(354); clf;
subplot(3,1,2); 
plot(CNCCtad3);axis([11 52 300 600]);hold on;
title('Distance to points across SOX9');
subplot(3,1,1);
plot(reelmap2);axis([15 90 30 40]);hold on;
title('reelin');
subplot(3,1,3);
plot(mulmap2); axis([15 90 18 28]);hold on;
title('mulloop');

%% Figure S5C
dM = cncMaps(1:52,1:52,:);
strA = squeeze(dM(42,:,:));
strC = strA < 150; %100,150,200,250
a = sum(cumsum(strC,1));
[~,id] = sort(a);

[~,first]=max(strC,[],1,'linear');
[rw,cl] = ind2sub(size(strC),first);
[~,id] = sort(rw(1:size(strC,2)));
datSort = strC(:,id);
