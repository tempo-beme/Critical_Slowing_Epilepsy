% this demo will produce some of the figures from Maturana etal 2020. Note
% that some of the figures will not be reproduced exactly due to only
% including a subset of the total data

load('demo_data1.mat')
d = DataACFWidth(:,:,7); %change this to DataVariance, DataEnergy, DataACF to see other plots
%channel used to generate ACFW plot in Figure 2 is 16
channel = 16; %change between 1 and 16 to see different channels
dch = d(:,channel);

figure(1); clf
set(gcf, 'OuterPosition', [100 100 1400 800]);
plot(x_days, dch,'Color',[1 1 1]*0.7);hold on;

%compute the slow and fast signals
causal = 0;
if(causal)
    dataSig_slow = movmean(dch, [720*2 0]);
    dataSig_fast = movmean(dch, [20 0] ) ;
else
    dataSig_slow = movmean(dch, [720*2]);
    dataSig_fast = movmean(dch, [20] ) ;
end

%plot(x_days, dataSig_slow,'Color', [.83, .82, .78]);hold on;
plot(x_days,dataSig_fast,'Color', [.2,.3,.49],'LineWidth',2);
plot(x_days, dataSig_slow,'Color',[.85, .33, .1] ,'LineWidth',2);
scatter(sz_days, 25*ones(1,length(sz_days)),150, 'vr', 'filled');
xticks(260:10:300); yticks(0:10:40); box off;
xlim([260 300]); 
ylim([0 31])
xticks(260:10:300); 

SzDayID = nan(1,length(sz_days));    
for sz=1:length(sz_days)
   [mn id] = min(abs(x_days-sz_days(sz)));
   SzDayID(sz) = id;    
end

horz_wind = 2;

% hilbert transform
hbsig1 = hilbert(dataSig_fast-mean(dataSig_fast));
hasig1 = angle(hbsig1);
hbsig2 = hilbert(dataSig_slow-mean(dataSig_slow));
hasig2 = angle(hbsig2);

% find the signal phase relative to sizure times
sz_phases_c = hbsig1(SzDayID-horz_wind); sz_phases_c = sz_phases_c./abs(sz_phases_c);
SI1 = sum(sz_phases_c); SI1 = abs(SI1)/length(sz_phases_c);
sz_phase1 = angle( sz_phases_c); 

sz_phases_c = hbsig2(SzDayID-horz_wind); sz_phases_c = sz_phases_c./abs(sz_phases_c);
SI2 = sum(sz_phases_c); SI2 = abs(SI2)/length(sz_phases_c);
sz_phase2 = angle(sz_phases_c); 
disp(['SI fast: ' num2str(SI1) ', SI slow: ' num2str(SI2)])

%cmpute the polar histograms
anglebins = -pi:2*pi/20:pi;
figure(1);clf;
c = histogram(sz_phase1,anglebins);
BinCounts1 = c.BinCounts;
c = histogram(sz_phase2,anglebins);
BinCounts2 = c.BinCounts;
mx1 = max([max(BinCounts1) max(BinCounts2)]);
BinCounts1 = BinCounts1/mx1;  
BinCounts2 = BinCounts2/mx1;    
c = histogram(hasig1,anglebins);
BinCounts3 = c.BinCounts;BinCounts3=BinCounts3/max(BinCounts3);
c = histogram(hasig2,anglebins);
BinCounts4 = c.BinCounts;BinCounts4=BinCounts4/max(BinCounts4);

% plot polar histograms
figure(1);clf;
pax = polarhistogram('BinCounts',BinCounts3, 'BinEdges',anglebins,'EdgeColor','b',...
   'FaceColor','none','DisplayStyle', 'stairs','LineWidth',1,'FaceAlpha',0.7);hold on;
pax = polarhistogram('BinCounts',BinCounts4, 'BinEdges',anglebins,'EdgeColor','r',...
       'FaceColor','none','DisplayStyle', 'stairs','LineWidth',1,'FaceAlpha',0.7);hold on;
pax = polarhistogram('BinCounts',BinCounts1, 'BinEdges',anglebins,'EdgeColor','k',...
   'FaceColor','b','LineStyle', '-','LineWidth',2,'FaceAlpha',0.7);hold on;
pax = polarhistogram('BinCounts',BinCounts2, 'BinEdges',anglebins,'EdgeColor','k',...
    'FaceColor','r','LineStyle', '-','LineWidth',2,'FaceAlpha',0.7);hold on;
pax.Parent.ThetaTickLabel = [];
pax.Parent.RTickLabel = [];
pax.Parent.GridAlpha = 0.3;
pax.Parent.RLim = [0 1.05];
set(gcf, 'OuterPosition', [200 100 800 800]);


%% will reproduce figure 1D ACFW
load('demo_data2.mat')
dt = 1/Fs;

% filter seizure signal
dfilt = designfilt('lowpassiir','FilterOrder',6,'HalfPowerFrequency',170,...
                'SampleRate',Fs); 
do = Data;
b = find(isnan(Data));
d = Data - nanmean(Data); d(b) = nanstd(d)*randn(1,length(b));
Data = filtfilt(dfilt,d);
Data(b) = nan;
ttsz = dt:dt:dt*length(Data);ttsz = ttsz+(xmin*60); ttsz = ttsz/60;

figure(1);
subplot(3,1,2);
plot(ttsz, do,'Color',[1 1 1]*.7);hold on
plot(ttsz, Data,'k');hold on;                
a = find(ttsz > 0 & ttsz < szdur/(60));
plot(ttsz(a), Data(a), 'r');
xlim([xmin xmax]);
box off;
if(xmin<-10)
    xticks(xmin:10:xmax)
    xticklabels({});           
else
    xticks(xmin:xmax);xticklabels({});
end
set(gca,'FontName', 'Arial', 'FontSize', 26); 
drawnow;

%compute ACFW and var
acfw = []; vars = [];idx = 1; step = 20;
win = 5*round(Fs);
for ii = win+1:step:length(Data)
    i1 = ii - win; i2 = ii;
    d = Data(i1:i2);
    try
        [r, lag] = autocorr(d, floor(Fs));
        [width, height] = fwhm(lag,r,1);
        acfw(idx) = width; 
    catch
    end
    vars(idx) = var(d);
    idx = idx + 1;
end
Mlen = 5;
macf = movmean(acfw, [Mlen 0]);
mvars = movmean(vars, [Mlen 0]);
tacf = 1:length(macf); tacf = (step*dt)*tacf; tacf = tacf/60; tacf = tacf + xmin;
a = find(tacf > 0 & tacf < szdur/(60));
subplot(3,1,1)
plot(tacf, macf, 'k');hold on;
plot(tacf(a), macf(a), 'r');hold on;   
xlim([xmin xmax]); ylim([0 max(macf)+1]);
xticks(xmin:10:xmax)
xticklabels({});set(gca,'FontName', 'Arial', 'FontSize', 26); box off

subplot(3,1,3)
plot(tacf, mvars, 'k');hold on;
plot(tacf(a), mvars(a), 'r');hold on;               
set(gca,'FontName', 'Arial', 'FontSize', 26); 
set(gca,'YScale','log')
ylim([100 max(mvars)+1000]);xlim([xmin xmax]);
box off
