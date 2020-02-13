

%% Mean SI across patients
filename = 'PtSI.mat';
load(filename)
PtSI(3,:) =[];

filenamephases = ['PT_MEAN_PHASE_ANGLE.mat'];
load(filenamephases);
MEAN_SI_PHASE(3,:) =[];

col1 = [0.2 0.3 0.49];
col2 = [0.8 0. 0.2];
col3 = [0.0 0.75 0.75];
COL = [col1;col2;col3];

figure(1);clf;
set(gcf, 'OuterPosition', [100 100 800 800]);
ha1 = axes('Position',[0.1 0.1 0.4 0.85]);
ha2 = axes('Position',[0.55 0.1 0.4 0.85]);
% short cycle
axes(ha1);
for i=4:6
    SI = PtSI(:,i);
    mx = mean(SI);
    stdx = std(SI);
    
    PHASES = MEAN_SI_PHASE(:,i);
    mph = mean(PHASES);
    sph = std(PHASES);
    
    ii = i-3;
    xx = ii*ones(1,length(SI)) + randn(1,length(SI))*0.00;
    scatter3(xx, SI, PHASES, 20, 'MarkerFaceColor',COL(ii,:), 'MarkerEdgeColor',COL(ii,:));
    hold on;
    plot3([ii ii], mx+[-stdx stdx], [mph mph],'Color',COL(ii,:), 'LineWidth', 2)
    plot3([ii ii], [mx mx], mph+[-sph sph],'Color',COL(ii,:), 'LineWidth', 2)
    scatter3(ii, mx,mph,150, 'sk', 'MarkerFaceColor',COL(ii,:), 'MarkerEdgeColor',COL(ii,:));
end
xlim([0.5 3.5]);
ylim([0 1]);
xticks(1:3); xticklabels({})
zlim([-pi pi]);zticks([-pi 0 pi]);zticklabels({})
yticks(0:0.5:1); yticklabels({})
set(gca,'FontName', 'Arial', 'FontSize', 26);
view(-65,10)
grid on;

% long cycles

for ii=1:3
    SI = PtSI(:,ii);
    mx = mean(SI);
    stdx = std(SI);
    
    PHASES = MEAN_SI_PHASE(:,ii);
    mph = mean(PHASES);
    sph = std(PHASES);
    
    axes(ha2);
%     xx = i*ones(1,length(SI)) + randn(1,length(SI))*0.05;
%     scatter(xx, SI, 20, 'MarkerFaceColor',COL(i,:), 'MarkerEdgeColor',COL(i,:)); hold on;
%     plot([i i], mx+[-stdx stdx],'Color',COL(i,:), 'LineWidth', 2)
%     scatter(i, mx,150, 'sk', 'MarkerFaceColor',COL(i,:), 'MarkerEdgeColor',COL(i,:));
    xx = ii*ones(1,length(SI)) + randn(1,length(SI))*0.00;
    scatter3(xx, SI, PHASES, 20, 'MarkerFaceColor',COL(ii,:), 'MarkerEdgeColor',COL(ii,:));
    hold on;
    plot3([ii ii], mx+[-stdx stdx], [mph mph],'Color',COL(ii,:), 'LineWidth', 2)
    plot3([ii ii], [mx mx], mph+[-sph sph],'Color',COL(ii,:), 'LineWidth', 2)
    scatter3(ii, mx,mph,150, 'sk', 'MarkerFaceColor',COL(ii,:), 'MarkerEdgeColor',COL(ii,:));
    
end
xlim([0.5 3.5]);
ylim([0 1]);
xticks(1:3); xticklabels({})
yticks(0:0.5:1); yticklabels({})
zlim([-pi pi]);zticks([-pi 0 pi]);zticklabels({})
set(gca,'FontName', 'Arial', 'FontSize', 26);
view(-65,10)
grid on;




figure(2);clf;
set(gcf, 'OuterPosition', [100 100 400 800]);
ha1 = axes('Position',[0.1 0.1 0.4 0.85]);
ha2 = axes('Position',[0.55 0.1 0.4 0.85]);
% autocorr
ii = [1, 4];
x1 = 1;
for i=ii
    if(x1 ==1)
        axes(ha1)
    else
        axes(ha2)
    end
    c1 = 1;
    SI = PtSI(:,i);
    mx = mean(SI);
    stdx = std(SI);
    PHASES = MEAN_SI_PHASE(:,i);
    mph = mean(PHASES);
    sph = std(PHASES);
    xx = x1*ones(1,length(SI)) + randn(1,length(SI))*0.00;
    scatter3(xx, SI, PHASES, 20, 'MarkerFaceColor',COL(c1,:), 'MarkerEdgeColor',COL(c1,:));
    hold on;
    plot3([x1 x1], mx+[-stdx stdx], [mph mph],'Color',COL(c1,:), 'LineWidth', 2)
    plot3([x1 x1], [mx mx], mph+[-sph sph],'Color',COL(c1,:), 'LineWidth', 2)
    scatter3(x1, mx,mph,150, 'sk', 'MarkerFaceColor',COL(c1,:), 'MarkerEdgeColor',COL(c1,:));
    view(90,0)
    ylim([0 1]);
    yticks(0:0.5:1); yticklabels({})
    zlim([-pi pi]);zticks([-pi 0 pi]);zticklabels({})
    set(gca,'FontName', 'Arial', 'FontSize', 26);
    grid on;
    
    x1 = x1+1;
end

figure(3);clf;
set(gcf, 'OuterPosition', [100 100 400 800]);
ha1 = axes('Position',[0.1 0.1 0.4 0.85]);
ha2 = axes('Position',[0.55 0.1 0.4 0.85]);
% autocorr
ii = [2, 5];
x1 = 1;
for i=ii
    if(x1 ==1)
        axes(ha1)
    else
        axes(ha2)
    end
    c1 = 2;
    SI = PtSI(:,i);
    mx = mean(SI);
    stdx = std(SI);
    PHASES = MEAN_SI_PHASE(:,i);
    mph = mean(PHASES);
    sph = std(PHASES);
    xx = x1*ones(1,length(SI)) + randn(1,length(SI))*0.00;
    scatter3(xx, SI, PHASES, 20, 'MarkerFaceColor',COL(c1,:), 'MarkerEdgeColor',COL(c1,:));
    hold on;
    plot3([x1 x1], mx+[-stdx stdx], [mph mph],'Color',COL(c1,:), 'LineWidth', 2)
    plot3([x1 x1], [mx mx], mph+[-sph sph],'Color',COL(c1,:), 'LineWidth', 2)
    scatter3(x1, mx,mph,150, 'sk', 'MarkerFaceColor',COL(c1,:), 'MarkerEdgeColor',COL(c1,:));
    view(90,0)
    ylim([0 1]);
    yticks(0:0.5:1); yticklabels({})
    zlim([-pi pi]);zticks([-pi 0 pi]);zticklabels({})
    set(gca,'FontName', 'Arial', 'FontSize', 26);
    grid on;
    
    x1 = x1+1;
end

figure(4);clf;
set(gcf, 'OuterPosition', [100 100 400 800]);
ha1 = axes('Position',[0.1 0.1 0.4 0.85]);
ha2 = axes('Position',[0.55 0.1 0.4 0.85]);
% autocorr
ii = [3, 6];
x1 = 1;
for i=ii
    if(x1 ==1)
        axes(ha1)
    else
        axes(ha2)
    end
    c1 = 3;
    SI = PtSI(:,i);
    mx = mean(SI);
    stdx = std(SI);
    PHASES = MEAN_SI_PHASE(:,i);
    mph = mean(PHASES);
    sph = std(PHASES);
    xx = x1*ones(1,length(SI)) + randn(1,length(SI))*0.00;
    scatter3(xx, SI, PHASES, 20, 'MarkerFaceColor',COL(c1,:), 'MarkerEdgeColor',COL(c1,:));
    hold on;
    plot3([x1 x1], mx+[-stdx stdx], [mph mph],'Color',COL(c1,:), 'LineWidth', 2)
    plot3([x1 x1], [mx mx], mph+[-sph sph],'Color',COL(c1,:), 'LineWidth', 2)
    scatter3(x1, mx,mph,150, 'sk', 'MarkerFaceColor',COL(c1,:), 'MarkerEdgeColor',COL(c1,:));
    view(90,0)
    ylim([0 1]);
    yticks(0:0.5:1); yticklabels({})
    zlim([-pi pi]);zticks([-pi 0 pi]);zticklabels({})
    set(gca,'FontName', 'Arial', 'FontSize', 26);
    grid on;
    
    x1 = x1+1;
end


%% corr coeff
MC1 = [];
MC2 = [];
MC3 = [];
filename = 'meanCorrCoeff.mat';
load(filename);

[a,b,stats] = anova2([MC1(:) MC2(:) MC3(:) ])
d3 = multcompare(stats)

col1 = [0.2 0.3 0.49];
col2 = [0.8 0. 0.2];
col3 = [0.0 0.75 0.75];
COL = [col1;col2;col3];
m1 = mean(MC1); s1 = std(MC1);
m2 = mean(MC2); s2 = std(MC2);
m3 = mean(MC3); s3 = std(MC3);

close all
figure(1);clf;
set(gcf, 'OuterPosition', [100 100 400 800]);
plot([1 1], m1+[s1 -s1],'Color', col1,'LineWidth',3);hold on;
scatter(1, m1,300,'sk','MarkerFaceColor', col1, 'MarkerEdgeColor', 'none');
scatter(ones(1,length(MC1))+randn(1,length(MC1))*0.05, MC1,...
    'MarkerFaceColor', col1, 'MarkerEdgeColor', 'none','MarkerFaceAlpha', 0.5);

plot(2*[1 1], m2+[s2 -s2],'Color', col2,'LineWidth',3);hold on;
scatter(2*1, m2,300,'sk','MarkerFaceColor', col2, 'MarkerEdgeColor', 'none');
scatter(2*ones(1,length(MC1))+randn(1,length(MC1))*0.05, MC2,...
    'MarkerFaceColor', col2, 'MarkerEdgeColor', 'none','MarkerFaceAlpha', 0.5);

plot(3*[1 1], m3+[s3 -s3],'Color', col3,'LineWidth',3);hold on;
scatter(3*1, m3,300,'sk','MarkerFaceColor', col3, 'MarkerEdgeColor', 'none');
scatter(3*ones(1,length(MC1))+randn(1,length(MC1))*0.05, MC3,...
    'MarkerFaceColor', col3, 'MarkerEdgeColor', 'none','MarkerFaceAlpha', 0.5);
xlim([0.5 3.5]);
ylim([0 1]);
xticks(1:3); xticklabels({})
yticks(0:0.25:1); %yticklabels({})
set(gca,'FontName', 'Arial', 'FontSize', 26);
box off;
grid on;


%% plot overall forecasting results (figure 5 B, C, D)

filename = 'RiskTransitions.mat';
load(filename);
filename = 'PtSzLevel2.mat';
load(filename);
RandSzLevel(3,:)=[];%remove patient 3
RandTimeRisk(3,:)=[];
RandSzLevel(4,:)=[];%remove patient 5 due to too few seizures
RandTimeRisk(4,:)=[];


% compare stats for iterative
figure(1);clf
% Seizures in high risk
P1 = [Forecast_acor_var_spikes(:,1),...
    Forecast_iterative(:,1) RandSzLevel(:,1)];
mean3 = mean(P1(:,1));std3 = std(P1(:,1));
mean4 = mean(P1(:,2));std4 = std(P1(:,2));
mean5 = mean(P1(:,3));std5 = std(P1(:,3));
plot([1 1]*0.9, mean3+[std3 -std3],'k');hold on
scatter(0.9, mean3,300,'dk', 'filled');
plot([1 1]*1, mean4+[std4 -std4],'k');hold on
scatter(1, mean4,300,'ok', 'filled');
plot([1 1]*1.1, mean5+[std5 -std5],'k');hold on
scatter(1.1, mean5,300,'sk', 'filled');
[a,b,stats] = anova2(P1);
figure(2);clf
c1 = multcompare(stats);


% seizure in low risk
P2 = [Forecast_acor_var_spikes(:,3),...
    Forecast_iterative(:,3) RandSzLevel(:,3)];
figure(1);
mean3 = mean(P2(:,1));std3 = std(P2(:,1));
mean4 = mean(P2(:,2));std4 = std(P2(:,2));
mean5 = mean(P2(:,3));std5 = std(P2(:,3));
plot([1 1]*1.9, mean3+[std3 -std3],'k');hold on
scatter(1.9, mean3,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 'd');
plot([1 1]*2, mean4+[std4 -std4],'k');hold on
scatter(2, mean4,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
plot([1 1]*2.1, mean5+[std5 -std5],'k');hold on
scatter(2.1, mean5,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 's');
box off; xlim([0.5 2.5]); ylim([0 1])
set(gcf, 'OuterPosition', [200 100 800 800]);
set(gca,'FontName', 'Arial', 'FontSize', 32);
yticks(0:.25:1); xticks([]);xticklabels({})
figure(2);clf
[a,b,stats] = anova2(P2);
c1 = multcompare(stats);


%time in high risk level
figure(4);clf
% Seizures in high risk
P3 = [ Forecast_acor_var_spikes(:,4),...
    Forecast_iterative(:,4) RandTimeRisk(:,1)];

mean3 = mean(P3(:,1));std3 = std(P3(:,1));
mean4 = mean(P3(:,2));std4 = std(P3(:,2));
mean5 = mean(P3(:,3));std5 = std(P3(:,3));
plot([1 1]*0.9, mean3+[std3 -std3],'k');hold on
scatter(0.9, mean3,300,'dk', 'filled');
plot([1 1]*1, mean4+[std4 -std4],'k');hold on
scatter(1, mean4,300,'ok', 'filled');
plot([1 1]*1.1, mean5+[std5 -std5],'k');hold on
scatter(1.1, mean5,300,'sk', 'filled');
[a,b,stats] = anova2(P3);
figure(2);clf
c1 = multcompare(stats);


% time in low risk
P4 = [ Forecast_acor_var_spikes(:,6),...
    Forecast_iterative(:,6) RandTimeRisk(:,3)];
figure(4);
set(gcf, 'OuterPosition', [200 100 800 800]);
set(gca,'FontName', 'Arial', 'FontSize', 32);
mean3 = mean(P4(:,1));std3 = std(P4(:,1));
mean4 = mean(P4(:,2));std4 = std(P4(:,2));
mean5 = mean(P4(:,3));std5 = std(P4(:,3));
plot([1 1]*1.9, mean3+[std3 -std3],'k');hold on
scatter(1.9, mean3,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 'd');
plot([1 1]*2, mean4+[std4 -std4],'k');hold on
scatter(2, mean4,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
plot([1 1]*2.1, mean5+[std5 -std5],'k');hold on
scatter(2.1, mean5,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 's');
box off;
xlim([0.5 2.5]); ylim([0 1])
yticks(0:.25:1);xticks([]);xticklabels({})
figure(2);clf
[a,b,stats] = anova2(P4);
c1 = multcompare(stats);

Pratio1 = P1(:,1).*P4(:,1);%ratio of seizure in high risk, time in low risk
Pratio2 = P1(:,2).*P4(:,2);%ratio of seizure in high risk, time in low risk
Pratio3 = P1(:,3).*P4(:,3);%ratio of seizure in high risk, time in low risk
figure(6);clf;
mean1 = mean(Pratio1);std1 = std(Pratio1);
mean2 = mean(Pratio2);std2 = std(Pratio2);
mean3 = mean(Pratio3);std3 = std(Pratio3);
plot([1 1]*0.9, mean1+[std1 -std1],'k');hold on
scatter(0.9, mean1,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 'd');
plot([1 1]*1, mean2+[std2 -std2],'k');hold on
scatter(1, mean2,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
plot([1 1]*1.1, mean3+[std3 -std3],'k');hold on
scatter(1.1, mean3,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 's');
box off;
xlim([0.85 1.15]); ylim([0 1])
yticks(0:.25:1);xticks([]); xticklabels({})
set(gcf, 'OuterPosition', [200 100 800 800]);
set(gca,'FontName', 'Arial', 'FontSize', 32);
figure(2);clf
[a,b,stats] = anova2([Pratio1(:) Pratio2(:) Pratio3(:)]);
c1 = multcompare(stats);


%% forecasting with ACF & Var, compared to spikes and ACF+VAR+SPIKES
%supp figure 21

filename = 'PtSzLevel2.mat';
load(filename);

F1 =[]; %Acor+var
F2 = []; %Spikes only
F3 = []; %combined spikes+crit slowing
for iPt = [1,2,4:15]
    pt = sprintf('%02.0f', iPt);
    load(['Forecast/Pt_' pt '/Forecast.mat']);
    F1(iPt,:) = ForeCastSummaryACF;
    F2(iPt,:) = ForeCastSummarySpikesOnly;
    F3(iPt,:) = ForeCastSummaryWithSpikes;
end

figure(1);clf
% Seizures in high risk
P1 = [Forecast_acor_var(:,1) Forecast_spikes(:,1),...
    Forecast_acor_var_spikes(:,1)];


mean3 = mean(P1(:,1));std3 = std(P1(:,1));
mean4 = mean(P1(:,2));std4 = std(P1(:,2));
mean5 = mean(P1(:,3));std5 = std(P1(:,3));
plot([1 1]*0.9, mean3+[std3 -std3],'k');hold on
scatter(0.9, mean3,300,'dk', 'filled');
plot([1 1]*1, mean4+[std4 -std4],'k');hold on
scatter(1, mean4,300,'ok', 'filled');
plot([1 1]*1.1, mean5+[std5 -std5],'k');hold on
scatter(1.1, mean5,300,'sk', 'filled');
[a,b,stats] = anova2(P1);
figure(2);clf
c1 = multcompare(stats);


% seizure in low risk
P2 = [Forecast_acor_var(:,3) Forecast_spikes(:,3),...
    Forecast_acor_var_spikes(:,3)];

figure(1);
mean3 = mean(P2(:,1));std3 = std(P2(:,1));
mean4 = mean(P2(:,2));std4 = std(P2(:,2));
mean5 = mean(P2(:,3));std5 = std(P2(:,3));
plot([1 1]*1.9, mean3+[std3 -std3],'k');hold on
scatter(1.9, mean3,300,'MarkerFaceColor','w', 'MarkerEdgeColor', 'k', 'Marker', 'd');
plot([1 1]*2, mean4+[std4 -std4],'k');hold on
scatter(2, mean4,300,'MarkerFaceColor','w', 'MarkerEdgeColor', 'k', 'Marker', 'o');
plot([1 1]*2.1, mean5+[std5 -std5],'k');hold on
scatter(2.1, mean5,300,'MarkerFaceColor','w', 'MarkerEdgeColor', 'k', 'Marker', 's');
box off; xlim([0.5 2.5]); ylim([0 1])
set(gcf, 'OuterPosition', [200 100 800 800]);
set(gca,'FontName', 'Arial', 'FontSize', 32);
yticks(0:.25:1); xticks([]);xticklabels({})
figure(2);clf
[a,b,stats] = anova2(P2);
c1 = multcompare(stats);



%time in high risk level
figure(4);clf
% Seizures in high risk
P3 = [Forecast_acor_var(:,4) Forecast_spikes(:,4),...
    Forecast_acor_var_spikes(:,4)];

mean3 = mean(P3(:,1));std3 = std(P3(:,1));
mean4 = mean(P3(:,2));std4 = std(P3(:,2));
mean5 = mean(P3(:,3));std5 = std(P3(:,3));
plot([1 1]*0.9, mean3+[std3 -std3],'k');hold on
scatter(0.9, mean3,300,'dk', 'filled');
plot([1 1]*1, mean4+[std4 -std4],'k');hold on
scatter(1, mean4,300,'ok', 'filled');
plot([1 1]*1.1, mean5+[std5 -std5],'k');hold on
scatter(1.1, mean5,300,'sk', 'filled');
[a,b,stats] = anova2(P3);
figure(2);clf
c1 = multcompare(stats);


% time in low risk
P4 = [Forecast_acor_var(:,6) Forecast_spikes(:,6),...
    Forecast_acor_var_spikes(:,6)];
figure(4);
set(gcf, 'OuterPosition', [200 100 800 800]);
set(gca,'FontName', 'Arial', 'FontSize', 32);
mean3 = mean(P4(:,1));std3 = std(P4(:,1));
mean4 = mean(P4(:,2));std4 = std(P4(:,2));
mean5 = mean(P4(:,3));std5 = std(P4(:,3));
plot([1 1]*1.9, mean3+[std3 -std3],'k');hold on
scatter(1.9, mean3,300,'MarkerFaceColor','w', 'MarkerEdgeColor', 'k', 'Marker', 'd');
plot([1 1]*2, mean4+[std4 -std4],'k');hold on
scatter(2, mean4,300,'MarkerFaceColor','w', 'MarkerEdgeColor', 'k', 'Marker', 'o');
plot([1 1]*2.1, mean5+[std5 -std5],'k');hold on
scatter(2.1, mean5,300,'MarkerFaceColor','w', 'MarkerEdgeColor', 'k', 'Marker', 's');
box off;
xlim([0.5 2.5]); ylim([0 1])
yticks(0:.25:1);xticks([]);xticklabels({})
figure(2);clf
[a,b,stats] = anova2(P4);
c1 = multcompare(stats);

Pratio1 = P1(:,1).*P4(:,1);%ratio of seizure in high risk, time in low risk
Pratio2 = P1(:,2).*P4(:,2);%ratio of seizure in high risk, time in low risk
Pratio3 = P1(:,3).*P4(:,3);%ratio of seizure in high risk, time in low risk
figure(6);clf;
% G1 = [ones(10,1);2*ones(10,1);3*ones(10,1)];
mean1 = mean(Pratio1);std1 = std(Pratio1);
mean2 = mean(Pratio2);std2 = std(Pratio2);
mean3 = mean(Pratio3);std3 = std(Pratio3);
plot([1 1]*0.9, mean1+[std1 -std1],'k');hold on
scatter(0.9, mean1,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 'd');
plot([1 1]*1, mean2+[std2 -std2],'k');hold on
scatter(1, mean2,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
plot([1 1]*1.1, mean3+[std3 -std3],'k');hold on
scatter(1.1, mean3,300,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'Marker', 's');
box off;
xlim([0.85 1.15]); ylim([0 1])
yticks(0:.25:1);xticks([]); xticklabels({})
set(gcf, 'OuterPosition', [200 100 800 800]);
set(gca,'FontName', 'Arial', 'FontSize', 32);
figure(2);clf
[a,b,stats] = anova2([Pratio1(:) Pratio2(:) Pratio3(:)]);
c1 = multcompare(stats);


%% mean duration of slow and fast cycles, Supp figure 20A
load(['CycleDurations.mat'])

col1 = [0.2 0.3 0.49];
col2 = [0.8 0. 0.2];
col3 = [0.0 0.75 0.75];
COL = [col1; col2; col3];
figure(2);clf
rw = 1;
xi = [1.9 2 2.1 0.9 1 1.1];
for iPt = 1:15
    if(iPt ==3)
        continue; 
    end
    Mean = CycleDurations.Mean(rw,:);
    STD = CycleDurations.STD(rw,:);
   
    
    for f = 1:6
        figure(2);
        r = rem(f,3);if(r==0), r = 3; end
        scatter(xi(f)+randn*0.05, Mean(f), 'MarkerFaceColor', COL(r,:),'MarkerEdgeColor','none');hold on
    end
    
    set(gcf, 'OuterPosition', [200 100 1000 1000]);
    rw = rw + 1;
end
set(gcf, 'OuterPosition', [200 100 1000 1000]);

Long = [CycleDurations.Mean(:,1:3)];
Short = [CycleDurations.Mean(:,4:6)];
M1 = mean(Long(:))
S1 = std(Long(:))
M2 = mean(Short(:))
S2 = std(Short(:))
figure(2); 
scatter(2, M1, 300,'+k','LineWidth',3);
plot(2*[1 1],M1+[-S1 S1],'k','LineWidth',3);
scatter(1, M2, 300,'+k','LineWidth',3);
plot(1*[1 1],M2+[-S2 S2],'k','LineWidth',3);
xticks([1 2]); xticklabels([]);
yticks([0.5 1 10]); yticklabels({'0.5','1','10'})
set(gca, 'YScale', 'Log');
set(gcf, 'OuterPosition', [200 100 800 800]);
set(gca,'FontName', 'Arial', 'FontSize', 26);
grid on
