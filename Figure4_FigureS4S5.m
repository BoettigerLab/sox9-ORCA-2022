%% data prep
clear all
analysisFolder = 'Z:\Liang-Fu\2021-05-23\DNA_Expt\analysis2\'

[polys,maps,spotData] = CombineAllFits(analysisFolder, 'byFOV',true);


A7Fov = [22:38];
B10Fov = [10:20,40:49];
E1Fov = [1:8,51:53];

A7maps = cat(3,maps{A7Fov});
B10Cmaps = cat(3,maps{B10Fov});
E1maps = cat(3,maps{E1Fov});

A7polys = cat(3,polys{A7Fov});
B10polys = cat(3,polys{B10Fov});
E1polys = cat(3,polys{E1Fov});

[cMap,nObs] = ContactFrac(A7maps,'threshold',250);
figure(1); clf; imagesc(nObs);

[cMap,nObs] = ContactFrac(B10Cmaps,'threshold',250);
figure(2); clf; imagesc(nObs);

[cMap,nObs] = ContactFrac(E1maps,'threshold',250);
figure(3); clf; imagesc(nObs);

badHybes = [1:6,16,32,36,40,44,53,54,59:86];

goodA7maps = A7maps;
goodA7maps(badHybes,:,:)= NaN;
goodA7maps(:,badHybes,:)= NaN;

goodB10maps = B10Cmaps;
goodB10maps(badHybes,:,:)= NaN;
goodB10maps(:,badHybes,:)= NaN;

goodE1maps = E1maps;
goodE1maps(badHybes,:,:)= NaN;
goodE1maps(:,badHybes,:)= NaN;


A7rpc = ReadsPerCell(goodA7maps);
B10rpc = ReadsPerCell(goodB10maps);
E1rpc = ReadsPerCell(goodE1maps);

figure(7); clf; 
subplot (1,3,1);
hist(A7rpc,20);
title('A1');

subplot (1,3,2);
hist(B10rpc,20);
title('B10');

subplot (1,3,3);
hist(E1rpc,20);
title('E1');

goodA7polys = A7polys;
goodA7polys (badHybes, :, :) = NaN;

goodB10polys = B10polys;
goodB10polys (badHybes, :, :) = NaN;

goodE1polys = E1polys;
goodE1polys (badHybes, :, :) = NaN;


A7highDetect = A7rpc>20;
B10highDetect = B10rpc>20;
E1highDetect = E1rpc>20;

A7data= goodA7maps(:,:,A7highDetect);
E1data = goodE1maps(:,:,E1highDetect);
B10data= goodB10maps(:,:,B10highDetect);

%% Distance and contact map
A7medianDistance = nanmedian(goodA7maps(:,:,A7highDetect),3);
B10medianDistance = nanmedian(goodB10maps(:,:,B10highDetect),3);
E1medianDistance = nanmedian(goodE1maps(:,:,E1highDetect),3);


A7_Distance_interp = InterpMapNans(A7medianDistance);
B10_Distance_interp = InterpMapNans(B10medianDistance);
E1_Distance_interp = InterpMapNans (E1medianDistance);

A7_Contact_interp = InterpMapNans( ContactFrac(goodA7maps (:,:,A7highDetect),'threshold',250));
B10_Contact_interp = InterpMapNans( ContactFrac(goodB10maps (:,:,B10highDetect),'threshold',250));
E1_Contact_interp = InterpMapNans( ContactFrac(goodE1maps (:,:,E1highDetect),'threshold',250));

%% Figure 4C

figure(114);clf;
imagesc(A7_Contact_interp); caxis ([.05,.4]); axis([8 57 8 57]);
title(['A7 Contact (WT) Interp, n=' num2str(size(goodA7maps (:,:,A7highDetect),3))])
colorbar
GetColorMap ('whiteToRed')
set(gcf, 'position',[10,10,900,800]);

figure(114);clf;
imagesc(E1_Contact_interp); caxis ([.05,.4]); axis([8 57 8 57]);
title(['E1 Contact (SSE1.35 delta CTCF) Interp, n=' num2str(size(goodE1maps (:,:,E1highDetect),3))])
colorbar
GetColorMap ('whiteToRed')
set(gcf, 'position',[10,10,900,800]);

%% Figure 4D

pE1B =zeros(86,86);

for i=1:86;
    for k=1:86;
A7data1 =squeeze(A7data(i,k,:));
A7data11=~isnan(A7data1);
A7data111=A7data1(A7data11);

E1data1 =squeeze(E1data(i,k,:));
E1data11=~isnan(E1data1);
E1data111=E1data1(E1data11);

if  nansum(A7data1)==0
            pE1B(i,k)=0
else
    
bootstatA7 = bootstrp(1000,@median,A7data111);
bootstatE1 = bootstrp(1000,@median,E1data111);
E1_A7= bootstatE1-bootstatA7;

alpha=.05;
CI = prctile(E1_A7,[100*alpha/2,100*(1-alpha/2)]);

%Hypothesis test: Does the confidence interval cover zero?
H = CI(1)>0 | CI(2)<0;

pE1B(i,k)=H
end
    end
end;

E1S = E1_Contact_interp - A7_Contact_interp;
E1SS = E1S.*pE1B;

figure(42); clf;
imagesc(E1_Contact_interp - A7_Contact_interp); caxis([-.1,.1]); axis([8 57 8 57]);
title(['E1 - A7 Contact Subtraction Interp']);
colorbar
GetColorMap ('BlueWhiteRed')
set(gcf, 'position',[10,10,900,800]);
%exportgraphics(gcf,'E1_A7_Contact_Subtraction_Interp.pdf')


figure(41); clf;
imagesc(E1SS,'alphadata', pE1B); caxis([-.1,.1]); axis([8 57 8 57]);
title(['E1 - A7 Contact Subtraction Interp Bootstrape']);
colorbar
GetColorMap ('BlueWhiteRed');
set(gca, 'color',[.85 .85 .85]);
set(gcf, 'position',[10,10,900,800]);
%exportgraphics(gcf,'E1_A7_Contact_Subtraction_Interp_Bootstrape_v2.pdf')

%% Figure 4E
figure(114);clf;
imagesc(B10_Contact_interp); caxis ([.05,.4]); axis([8 57 8 57]);
title(['B10 Contact (SSEProm delta CTCF) Interp, n=' num2str(size(goodB10maps (:,:,B10highDetect),3))])
colorbar
GetColorMap ('whiteToRed')
set(gcf, 'position',[10,10,900,800]);

%% Figure 4F


pB10B =zeros(86,86);

for i=1:86;
    for k=1:86;
A7data1 =squeeze(A7data(i,k,:));
A7data11=~isnan(A7data1);
A7data111=A7data1(A7data11);

B10data1 =squeeze(B10data(i,k,:));
B10data11=~isnan(B10data1);
B10data111=B10data1(B10data11);

if  nansum(A7data1)==0
            pB10B(i,k)=0
else
    
bootstatA7 = bootstrp(1000,@median,A7data111);
bootstatB10 = bootstrp(1000,@median,B10data111);
B10_A7= bootstatB10-bootstatA7;

alpha=.05;
CI = prctile(B10_A7,[100*alpha/2,100*(1-alpha/2)]);

%Hypothesis test: Does the confidence interval cover zero?
H = CI(1)>0 | CI(2)<0;

pB10B(i,k)=H
end
    end
end;

B10S = B10_Contact_interp - A7_Contact_interp;
B10SS = B10S.*pB10B;

figure(11); clf;
imagesc(B10_Contact_interp - A7_Contact_interp); caxis([-.1,.1]); axis([8 57 8 57]);
title(['B10 - A7 Contact Subtraction Interp']);
colorbar
GetColorMap ('BlueWhiteRed')
set(gcf, 'position',[10,10,900,800]);
%exportgraphics(gcf,'B10_A7_Contact_Subtraction_Interp.pdf')


figure(10); clf;
imagesc(B10SS,'alphadata',pB10B); caxis([-.1,.1]); axis([8 57 8 57]);
title(['B10 - A7 Contact Subtraction Interp Bootstrap']);
colorbar
GetColorMap ('BlueWhiteRed');
set(gca, 'color',[.85 .85 .85]);
set(gcf, 'position',[10,10,900,800]);
%exportgraphics(gcf,'B10_A7_Contact_Subtraction_Interp_Bootstrap_v2.pdf')

%% Figure 4H

A7Rg = zeros(1537,1);
for i= 1:1537;
    A7Rg(i)= RadiusOfGyration(A7data(11:52,:,i));
end;
 
B10Rg = zeros(1474,1);
for i= 1:1474;
    B10Rg(i)= RadiusOfGyration(B10data(11:52,:,i));
end;

E1Rg = zeros(1450,1);
for i= 1:1450;
    E1Rg(i)= RadiusOfGyration(E1data(11:52,:,i));
end;

[pE1,hE1] = ranksum(E1Rg,A7Rg);
[pB10,hB10] = ranksum(B10Rg,A7Rg);

figure(3); clf;
cdfplot (A7Rg);hold on;
cdfplot (E1Rg);hold on;
cdfplot (B10Rg);hold on;
legend('A7','E1','B10');
xlabel ('Radius of Gyration (nm)');
ylabel ('% of obs');
xlim ([0 800]);
set(gcf, 'position',[10,10,900,800]);
%exportgraphics(gcf,'A7_E1_B10_CDF_x800_v2.pdf');

%% Figure 4I

A7data = goodA7polys (:, :, A7highDetect);
B10data = goodB10polys (:,:, B10highDetect);
E1data =goodE1polys (:,:, E1highDetect);

delData = cell(3,1);
delData{1} = A7data;
delData{2} = E1data;
delData{3} = B10data;

delDis = cell(3,1);
for d=1:3;
data = delData{d};
ncells = size(data,3);
Dis = zeros (86, ncells);
for i = 1:ncells;
    centerdata = data (11:52,:,i);
    rm =mean(centerdata,1,'omitnan');
    mapp=data(:,:,i);
    Dis(:,i)= sqrt((mapp(:,1)-rm(1)).^2 + (mapp(:,2)-rm(2)).^2 + (mapp(:,3)-rm(3)).^2);
end
medianDis = median(Dis,2,'omitnan');
medianDis =fillmissing(medianDis,'movmean',3);
delDis{d} = medianDis;
end

figure(11); clf;
plot(delDis{1}(11:52));axis([1 42 100 500]);hold on;
plot(delDis{2}(11:52));axis([1 42 100 500]);hold on;
plot(delDis{3}(11:52));axis([1 42 100 500]);hold on;
legend('A7','E1','B10');
xlabel ('Barcodes');
ylabel ('Median distance from center/Rg');
title('Median distance to center of SOX9 TAD');

%% Figure 4J

delMaps = cell(3,1);
delMaps{1} = goodA7maps (11:52, 11:52, A7highDetect);
delMaps{3} = goodB10maps (11:52,11:52, B10highDetect);
delMaps{2} = goodE1maps (11:52,11:52, E1highDetect);

delmean = cell(3,1);
for e=1:3;
tad = delMaps{e};
tad1= median(tad,3,'omitnan');
tad2 = mean(tad1,2,'omitnan');
tad3= fillmissing (tad2, 'movmean',3);

delmean{e} = tad3;
end


figure(12); clf;
plot(delmean{1});axis([1 42 300 600]);hold on;
plot(delmean{2});axis([1 42 300 600]);hold on;
plot(delmean{3});axis([1 42 300 600]);hold on;
legend('A7','E1','B10');
%xlabel ('Barcodes');
%ylabel ('Median distance from center/Rg');
title('Mean distance to other points');

%% Figure S5A-D

figure(112);clf;
imagesc(A7_Contact_interp); caxis ([.05,.4]); axis([1 80 1 80]);;
title(['A7 Contact (WT) Interp, n=' num2str(size(goodA7maps (:,:,A7highDetect),3))])
colorbar
GetColorMap ('whiteToRed')
set(gcf, 'position',[10,10,900,800]);
exportgraphics(gcf,'A7_Interp_All80.pdf')

figure(115);clf;
imagesc(E1_Contact_interp); caxis ([.05,.4]); axis([1 80 1 80]);;
title(['E1 Contact (SSE1.35 delta CTCF) Interp, n=' num2str(size(goodE1maps (:,:,E1highDetect),3))])
colorbar
GetColorMap ('whiteToRed')
set(gcf, 'position',[10,10,900,800]);
exportgraphics(gcf,'E1_Interp_All80.pdf')

figure(42); clf;
imagesc(E1_Contact_interp - A7_Contact_interp); caxis([-.1,.1]); axis([1 80 1 80]);
title(['E1 - A7 Contact Subtraction Interp']);
colorbar
GetColorMap ('BlueWhiteRed')
set(gcf, 'position',[10,10,900,800]);
exportgraphics(gcf,'E1_A7_Contact_Subtraction_Interp_All80.pdf')

figure(114);clf;
imagesc(B10_Contact_interp); caxis ([.05,.4]); axis([1 80 1 80]);
title(['B10 Contact (SSEProm delta CTCF) Interp, n=' num2str(size(goodB10maps (:,:,B10highDetect),3))])
colorbar
GetColorMap ('whiteToRed')
set(gcf, 'position',[10,10,900,800]);
exportgraphics(gcf,'B10_Interp_All80.pdf')

figure(15); clf;
imagesc(B10_Contact_interp - A7_Contact_interp); caxis([-.1,.1]); axis([1 80 1 80]);
title(['B10 - A7 Contact Subtraction Interp']);
colorbar
GetColorMap ('BlueWhiteRed')
set(gcf, 'position',[10,10,900,800]);
exportgraphics(gcf,'B10_A7_Contact_Subtraction_Interp_All80.pdf')











