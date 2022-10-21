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

% take only the SOX9 TAD (11:52) to SOX9 gene (42)
ESCdata2= squeeze (goodESCmaps(42,11:52,ESChighDetect));
ee=~isnan(ESCdata2);
eee=sum(ee,1);
eeee=eee>25;
ESCdata3=ESCdata2(:,eeee);
ESCdata4= fillmissing(ESCdata3,'movmean',5);
ESCdata5=flipud(rot90(ESCdata4));
ESCdata6= sum(~isnan(ESCdata5),2)>41;
ESCdata7 = ESCdata5(ESCdata6,:);
%Y = tsne(ESCdata5,'Perplexity',10);
%figure(22);clf;
%gscatter(Y(:,1),Y(:,2));

CNCCdata2 = squeeze(goodCNCCmaps(42,11:52,CNCChighDetect));
cc=~isnan(CNCCdata2);
ccc=sum(cc,1);
cccc=ccc>25;
CNCCdata3=CNCCdata2(:,cccc);
CNCCdata4= fillmissing(CNCCdata3,'movmean',5);
CNCCdata5=flipud(rot90(CNCCdata4));
CNCCdata6= sum(~isnan(CNCCdata5),2)>41;
CNCCdata7 = CNCCdata5(CNCCdata6,:);
%Y = tsne(CNCCdata5,'Perplexity',10);
%figure(23);clf;
%gscatter(Y(:,1),Y(:,2));

ESC_CNCC=cat(1,ESCdata7,CNCCdata7);

%% Figure 2A

Group3= flipud(rot90(repelem([{'ESC'}, {'CNCC'}],[570 780])));

%TSNE, group by cell type
Y=tsne(ESC_CNCC,'Perplexity',10);
figure(23);clf;
gscatter(Y(:,1),Y(:,2),Group3,[],[],15);
set(gcf, 'position',[10,10,900,800]);

%% Figure 2B
% kmean, 15 clusters
Kidx = kmeans (ESC_CNCC, 15);
% TSNE, by  kmean clusters
figure(31);clf;
gscatter(Y(:,1),Y(:,2),Kidx,[],[],15);
set(gcf, 'position',[10,10,900,800]);

%% Figure 2C

KidxESC=Kidx(1:570);
KidxCNCC=Kidx(571:1350);

ESCclustermap = goodESCmaps(11:52,11:52,ESChighDetect);
ESCclustermap2 = ESCclustermap (:,:,eeee); 
ESCclustermap3 = ESCclustermap2 (:,:,ESCdata6);

CNCCclustermap = goodCNCCmaps(11:52,11:52,CNCChighDetect);
CNCCclustermap2 = CNCCclustermap (:,:,cccc); 
CNCCclustermap3 = CNCCclustermap2 (:,:,CNCCdata6);

ESC_CNCC=cat(1,ESCdata7,CNCCdata7);
ESC_CNCC_Maps=cat(3,ESCclustermap3,CNCCclustermap3);


A = char.empty;
B = char.empty;
C = char.empty;
CN = char.empty;
EN = char.empty;
for i=1:15
    cluster1=ESC_CNCC_Maps(:,:,Kidx==i);
    cluster_Contact_interp_i = InterpMapNans (ContactFrac(cluster1));
    cluster_Distance_interp_i = InterpMapNans (nanmedian(cluster1,3));
A{i} = cluster_Distance_interp_i;
B{i} = cluster_Contact_interp_i;
C{i} =size (cluster1,3);
EN{i} = size(ESCdata7(KidxESC==i),1);
CN{i} = size(CNCCdata7(KidxCNCC==i),1);
ENN{i} = (EN{i}/570)/((EN{i}/570)+(CN{i}/780));
CNN{i} = (CN{i}/780)/((EN{i}/570)+(CN{i}/780));

end

figure(223);clf;
for plotid =1:15
subplot (3,5,plotid);
imagesc(A{plotid}); caxis([100 500]); 
title(strcat('cluster ', num2str(plotid), ' n=', num2str(C{plotid}),' ESC=', num2str(ENN{plotid},'%.2f'),' CNCC=', num2str(CNN{plotid},'%.2f')))
colorbar
GetColorMap ('redToWhite')
end

%% Figure S2A
% whole SOX9 TAD
% re organize data and interpating missing points
ESCdataT= goodESCmaps(11:52,11:52,ESChighDetect);
ESCrpc2 = ReadsPerCell(ESCdataT);
figure(7); clf;
hist(ESCrpc2,20);
title('ESC');
ESChighDetect2= ESCrpc2>25;
sum(ESChighDetect2)

ESCdata3=ESCdataT(:,:,ESChighDetect2);
ESCdata4 =zeros(42,42,803);
for i=1:803;
    K=ESCdata3(:,:,i);
    ESCdata4(:,:,i)=InterpMapNans(K);
end;


ESCdata5=extractUpperToVec(ESCdata4(:,:,:), 1);
ESCdata6= sum(~isnan(ESCdata5),2)>860;
ESCdata7 = ESCdata5(ESCdata6,:);


CNCCdataT= goodCNCCmaps(11:52,11:52,CNCChighDetect);
CNCCrpc2 = ReadsPerCell(CNCCdataT);
figure(7); clf;
hist(CNCCrpc2,20);
title('CNCC');
CNCChighDetect2= CNCCrpc2>25;
sum(CNCChighDetect2)

CNCCdata3=CNCCdataT(:,:,CNCChighDetect2);
CNCCdata4 =zeros(42,42,1193);
for i=1:1193;
    K=CNCCdata3(:,:,i);
    CNCCdata4(:,:,i)=InterpMapNans(K);
end;


CNCCdata5=extractUpperToVec(CNCCdata4(:,:,:), 1);
CNCCdata6= sum(~isnan(CNCCdata5),2)>860;
CNCCdata7 = CNCCdata5(CNCCdata6,:);


ESC_CNCC_all=cat(1,ESCdata7,CNCCdata7);

%Cell type id
Group4= flipud(rot90(repelem([{'ESC'}, {'CNCC'}],[262 411])));

%TSNE
Z=tsne(ESC_CNCC_all,'Perplexity',10);
figure(27);clf;
gscatter(Z(:,1),Z(:,2),Group4,[],[],15);
set(gcf, 'position',[10,10,900,800]);

%% Figure S2C

% kmean, 3 clusters
Kidx3 = kmeans (ESC_CNCC, 3);
% TSNE, by  kmean clusters
figure(31);clf;
gscatter(Y(:,1),Y(:,2),Kidx3,[],[],15);
set(gcf, 'position',[10,10,900,800]);

%% Figure S2D


A = char.empty;
B = char.empty;
C = char.empty;
CN = char.empty;
EN = char.empty;
for i=1:3
    cluster1=ESC_CNCC_Maps(:,:,Kidx3==i);
    cluster_Contact_interp_i = InterpMapNans (ContactFrac(cluster1));
    cluster_Distance_interp_i = InterpMapNans (nanmedian(cluster1,3));
A{i} = cluster_Distance_interp_i;
B{i} = cluster_Contact_interp_i;
C{i} =size (cluster1,3);
EN{i} = size(ESCdata7(KidxESC==i),1);
CN{i} = size(CNCCdata7(KidxCNCC==i),1);
ENN{i} = (EN{i}/570)/((EN{i}/570)+(CN{i}/780));
CNN{i} = (CN{i}/780)/((EN{i}/570)+(CN{i}/780));

end

figure(223);clf;
for plotid =1:3
subplot (3,1,plotid);
imagesc(A{plotid}); caxis([100 500]); 
title(strcat('cluster ', num2str(plotid), ' n=', num2str(C{plotid}),' ESC=', num2str(ENN{plotid},'%.2f'),' CNCC=', num2str(CNN{plotid},'%.2f')))
colorbar
GetColorMap ('redToWhite')
end

%% Figure S2E
figure(139);clf;
for plotid =1:15;
    si(plotid)=(size(ESC_CNCC(Kidx==plotid),1)/size((Kidx),1))*0.6;
    B(plotid)=0.97-sum(si(1:plotid))-0.02*plotid;
ax(plotid)= axes('position', [0.13 B(plotid) 0.775 si(plotid)]);
imagesc(ESC_CNCC(Kidx==plotid,:)); caxis([0 500]); 
set(gca,'YTickLabel',[],'XTickLabel',[]);
title(strcat('ESC CNCC cluster ', num2str(plotid), ' n=', num2str(size(ESC_CNCC(Kidx==plotid),1))))
%ylabel(strcat(' n=', num2str(size(ESC_CNCC(Kidx==plotid),1))));
colorbar;
GetColorMap ('redToWhite');

end
