% Code to run the forecast using a combination of autocorrelation width,
% variance and spike rates. 

function CSD_Forecasting(iPt)
    % add the ieeg.org matlab library
    addpath(genpath('ieeg-matlab-1.13.2'));

    %ieeg.org login details
    login = '';
    pword = '.bin';
    
    SAVE = 0; % save data
    SAVEFIG = 0; %save figures
    
    Patient{iPt} = '';
    pt = sprintf('%02.0f', iPt)
    curPt = Patient{iPt};
    patient = IEEGSession(['session_name'],login,pword);
    Fs = patient.data.sampleRate;
    dt = 1/Fs;
    
    % load seizure timing information
    szinfo = load(['Portal Annots/' curPt '_Annots']);
    SzDur = szinfo.SzDur;
    SzIndices = szinfo.SzIndices;
    SzType = szinfo.SzType;
    % use Type 1 and 2 seizures
    SzIndL = find( (SzType==1 | SzType==2) );
    SzDur = SzDur(SzIndL);
    SzIndices = SzIndices(SzIndL);

    %load original data. this is just to tell us where nans are
    fileListName = ['CSD_data/Pt_' pt '/DataTSCompiled.mat'];
    dO = load(['CSD_data/Pt_' pt '/DataTimeSeries.mat']);
    dO = dO.DataTimeSeries.DataACFWidth;   
    % load compiled data. analysis computed on this data
    d = load(fileListName);
    DataACFWidth = d.DataACFWidth;
    DataVariance = d.DataVariance; 
    t = (d.T*dt);%seconds
    x_days = t/(60*60*24); 
    
    szid = SzIndices;
    sz_days = szid*dt/(60*60*24);
    
    % compute seizure times relative to downsampled time x_days
    SzDayID = [];    
    for sz=1:length(sz_days)
       [mn id] = min(abs(x_days-sz_days(sz)));
       SzDayID(sz) = id;    
    end
    
    % phase data will be binned into anglebins
    anglebins = -pi:2*pi/20:pi;
    da = diff(anglebins);
    anglebins2mid = anglebins(1:end-1)+da;
    
    % load spike rates
    filename = ['CSD_data/Rates_' pt '.mat'];
    Rates = load(filename); Rates = Rates.Rates;

    % find best electrode and compute probability of seizure relative to
    % signals
    [elecACF1, elecVAR1,esp1, hbACFW1, hbVAR1, hbSP1, PACF_sp1,...
        PVAR_sp1, PSP_sp1] = forecast1(DataACFWidth(:,:,7),DataVariance(:,:,7),Rates); %slow cycle
    [elecACF2, elecVAR2,esp2, hbACFW2, hbVAR2,hbSP2, PACF_sp2,...
        PVAR_sp2, PSP_sp2] = forecast2(DataACFWidth(:,:,7),DataVariance(:,:,7),Rates); % fast cycle

    % compute seizure probability iteratively. use electrodes taken from
    % forecsast above
    IterativeForecast(elecACF1, elecVAR1, elecACF2, elecVAR2, esp1, esp2)
    
    
    P_slow = PACF_sp1(:)*PVAR_sp1(:)';
    P_fast = PACF_sp2(:)*PVAR_sp2(:)';

    P_Hor = 1; % prediction horizon
    %compute prediction from slow cycles
    %ACF
    [PcSzHigh1,PcSzLow1,PcTimeLow1, PcTimeHigh1,...
                pTimeACF1] = TimeProbability(hbACFW1, PACF_sp1, 0);
    % Var
    [PcSzHigh2,PcSzLow2,PcTimeLow2, PcTimeHigh2,...
                pTimeVar1] = TimeProbability(hbVAR1, PVAR_sp1, 0);
    Ptime_slow = pTimeACF1.*pTimeVar1;

    [xth,PcSzHighSlow,PcSzLowSlow,PcTimeLowSlow, PcTimeHighSlow]  = ProbTimeEval(Ptime_slow, SzDayID-P_Hor);
   
    %Spikes
    [PcSzHigh3,PcSzLow3,PcTimeLow3, PcTimeHigh3,...
                pTimeSP1] = TimeProbability(hbSP1, PSP_sp1, 0);
    [xth,PcSzHighSlow,PcSzLowSlow,PcTimeLowSlow, PcTimeHighSlow]  = ProbTimeEval(pTimeSP1, SzDayID-P_Hor);
   
    %compute prediction from fast cycles

    [PcSzHigh1,PcSzLow1,PcTimeLow1, PcTimeHigh1,...
                pTimeACF2] = TimeProbability(hbACFW2, PACF_sp2, 0);
    [PcSzHigh2,PcSzLow2,PcTimeLow2, PcTimeHigh2,...
                pTimeVar2] = TimeProbability(hbVAR2, PVAR_sp2, 0);
    Ptime_fast = pTimeACF2.*pTimeVar2;

    [xth,PcSzHighFast,PcSzLowFast,PcTimeLowFast, PcTimeHighFast]  = ProbTimeEval(Ptime_fast, SzDayID-P_Hor);
   
    [PcSzHigh1,PcSzLow1,PcTimeLow1, PcTimeHigh1,...
                pTimeSP2] = TimeProbability(hbSP2, PSP_sp2, 0);
    [xth,PcSzHighFast,PcSzLowFast,PcTimeLowFast, PcTimeHighFast]  = ProbTimeEval(pTimeSP2, SzDayID-P_Hor);
    

    % slow and fast spikes
    P_SPIKES = PSP_sp1(:)*PSP_sp2(:)';
    figure(7);clf
    contourf(anglebins2mid, anglebins2mid, P_SPIKES, 30,'EdgeColor','none');
    xticks([-pi 0 pi]);yticks([-pi 0 pi]);xticklabels([]);yticklabels([]);
    set(gcf, 'OuterPosition', [200 100 800 800]); 
    caxis([0 max(max(P_fast))])
    colormap(hot), hold on;

    % combine slow and fast cycles ACF+Var
    P_combinedACF = Ptime_slow.*Ptime_fast;
    [xth,~,~,~, ~]  = ProbTimeEval(P_combinedACF, SzDayID-P_Hor);


    [PcSzHighComb,PcSzMedComb,PcSzLowComb,...
                PcTimeHighComb, PcTimeMedComb, PcTimeLowComb] = RiskRatio(Prisk1);

    % combine cycles with spikes
    P_spikes_combined = pTimeSP1.*pTimeSP2;
    P_ACF_VAR_spikes = P_combinedACF.*P_spikes_combined;
    [xth,~,~,~, ~]  = ProbTimeEval(P_ACF_VAR_spikes, SzDayID-P_Hor);
    Prisk_ACF_VAR_spikes = zeros(1,length(P_ACF_VAR_spikes));
    Prisk_ACF_VAR_spikes(P_ACF_VAR_spikes<xth(1)) = 1;Prisk_ACF_VAR_spikes(P_ACF_VAR_spikes>=xth(2)) = 3;
    Prisk_ACF_VAR_spikes(P_ACF_VAR_spikes>=xth(1) & P_ACF_VAR_spikes<xth(2)) = 2;

    [PcSzHighComb,PcSzMedComb,PcSzLowComb,...
                PcTimeHighComb, PcTimeMedComb, PcTimeLowComb] = RiskRatio(Prisk_ACF_VAR_spikes);

    ForeCastSummaryWithSpikes = [PcSzHighComb,PcSzMedComb,...
        PcSzLowComb,PcTimeHighComb, PcTimeMedComb, PcTimeLowComb];

    str = ['Combined (ACF+VAR+SPK): ' num2str(PcSzHighComb) ' seizures in high, ' num2str(PcSzLowComb) ' seizures in low, ',...
        num2str(PcTimeHighComb) ' time in high, ' num2str(PcTimeLowComb) ' time in low, '];
    disp(str);


    % run a random markov model
    MarkovModel(Prisk_ACF_VAR_spikes)


%-----------------------------------------------------------------------
    function MarkovModel(P_risk)
        
        % transition probability
        filename = ['CSD_data/RiskTransitions.mat'];
        
        if(exist(filename,'file'))
            d = load(filename);
            RandSzLevel = d.RandSzLevel;
            RandTimeRisk = d.RandTimeRisk;
        else
            RandSzLevel = [];
        end
        risk = P_risk;

        a = find(risk == 1);
        N0 = 0; N1=0; N2 = 0; N=0;
        %transition probability when at 1
        for ii=1:length(a)
            if(a(ii)+1>length(risk)), continue, end
            if(risk(a(ii)+1) == 1)
                N0 = N0+1;
            elseif(risk(a(ii)+1) == 2)
                N1 = N1+1;
            elseif(risk(a(ii)+1) == 3)
                N2 = N2 + 1;
            end
            N = N+1;
        end
        RISKTRANS(1,:) = [N0/N N1/N N2/N];
        
        %transition probability when at 2
        N0 = 0; N1=0; N2 = 0; N=0;
        a = find(risk ==2);
        for ii=1:length(a)
            if((a(ii)+1)>length(risk)), continue, end
            if(risk(a(ii)+1) == 2)
                N0 = N0+1;
            elseif(risk(a(ii)+1) == 1)
                N1 = N1+1;
            elseif(risk(a(ii)+1) == 3)
                N2 = N2 + 1;
            end
            N = N+1;
        end
        RISKTRANS(2,:) = [N0/N N1/N N2/N];
        
        %transition probability when at 3
        N0 = 0; N1=0; N2 = 0; N=0;
        a = find(risk ==3);
        for ii=1:length(a)
            if((a(ii)+1)>length(risk)), continue, end
            if(risk(a(ii)+1) == 3)
                N0 = N0+1;
            elseif(risk(a(ii)+1) == 2)
                N1 = N1+1;
            elseif(risk(a(ii)+1) == 1)
                N2 = N2 + 1;
            end
            N = N+1;
        end
        RISKTRANS(3,:) = [N0/N N1/N N2/N];
        assignin('base', 'RISKTRANS', RISKTRANS)

        randomRisk = ones(1,length(risk));
        for ii=2:length(randomRisk)
            p=[];
            if(randomRisk(ii-1) == 1)
                p1 = RISKTRANS(1,1); %probability of no change
                p2 = RISKTRANS(1,2); %probability of change to 2
                p3 = RISKTRANS(1,3); %probability of change to 3
                p(1) = p1; p(2) = p1+p2; p(3) = 1;
            elseif(randomRisk(ii-1) == 2)
                p1 = RISKTRANS(2,1); %probability of no change
                p2 = RISKTRANS(2,2); %probability of change to 1
                p3 = RISKTRANS(2,3); %probability of change to 3
                p(1) = p1; p(2) = p1+p2; p(3) = 1;
            elseif(randomRisk(ii-1) == 3)
                p1 = RISKTRANS(3,1); %probability of no change
                p2 = RISKTRANS(3,2); %probability of change to 2
                p3 = RISKTRANS(3,3); %probability of change to 1
                p(1) = p1; p(2) = p1+p2; p(3) = 1;
            end
            
            rr = rand;
            plevel = randomRisk(ii-1);
            if(rr<p(1)) %no change
                randomRisk(ii) = plevel;
            elseif(rr>=p(1) && rr<p(2))
                if(plevel == 1)
                    randomRisk(ii) = 2;
                elseif(plevel == 2)
                    randomRisk(ii) = 1;
                else
                    randomRisk(ii) = 2;
                end
            else
                if(plevel == 1)
                    randomRisk(ii) = 3;
                elseif(plevel == 2)
                    randomRisk(ii) = 3;
                else
                    randomRisk(ii) = 1;
                end
            end
        end
        
        assignin('base', 'randomRisk', randomRisk)
        szin1 = length(find(randomRisk(SzDayID)==1));
        szin2 = length(find(randomRisk(SzDayID)==2));
        szin3 = length(find(randomRisk(SzDayID)==3));
        RandSzLevel(iPt,:) = [szin3 szin2 szin1]/length(SzDayID);

        szin1 = length(find(randomRisk()==1));
        szin2 = length(find(randomRisk()==2));
        szin3 = length(find(randomRisk()==3));
        RandTimeRisk(iPt,:) = [szin3 szin2 szin1]/length(randomRisk);

        save(filename, 'RandSzLevel','RandTimeRisk');
    end


%-----------------------------------------------------------------------
    function IterativeForecast(ch1, ch2, ch3, ch4, ch5, ch6)
        
        dataSig1 = (DataACFWidth(:, ch1,7));%slow cycle
        dataSig2 = (DataVariance(:, ch2,7));%slow cycle
        dataSig3 = (DataACFWidth(:, ch3,7));%fast cycle
        dataSig4 = (DataVariance(:, ch4,7));%fast cycle
        dataSig5 = (Rates(ch5,:));%long
        dataSig6 = (Rates(ch6,:));%short
        
        
        Mlen =20;
        CausalSigACFSlow = movmean(dataSig1,[720*2 0]) - movmean(dataSig1, [720*30 0]);
        CausalSigVARSlow = movmean(dataSig2,[720*2 0]) - movmean(dataSig2, [720*30 0]);
        CausalSigACFFast = movmean(dataSig3,[Mlen 0]) - movmean(dataSig3,[720*2 0]);
        CausalSigVARFast = movmean(dataSig4,[Mlen 0]) - movmean(dataSig4,[720*2 0]);
        CausalSigSPSlow = movmean(dataSig5,[720*2 0]) - movmean(dataSig5, [720*30 0]);
        CausalSigSPFast = movmean(dataSig6,[Mlen 0]) - movmean(dataSig6,[720*2 0]);
        
        d1 = CausalSigACFSlow;
        d2 = CausalSigVARSlow;
        d3 = CausalSigACFFast;
        d4 = CausalSigVARFast;
        d5 = CausalSigSPSlow;
        d6 = CausalSigSPFast;
        dh1 = hilbert(d1 - mean(d1));dha1 = angle(dh1);
        dh2 = hilbert(d2 - mean(d2));dha2 = angle(dh2);
        dh3 = hilbert(d3 - mean(d3));dha3 = angle(dh3);
        dh4 = hilbert(d4 - mean(d4));dha4 = angle(dh4);
        dh5 = hilbert(d5 - mean(d5));dha5 = angle(dh5);
        dh6 = hilbert(d6 - mean(d6));dha6 = angle(dh6);
        
        
        sz_phases_c = dh1(SzDayID-1); sz_phases_c = sz_phases_c./abs(sz_phases_c);
        SI1 = sum(sz_phases_c); SI1 = abs(SI1)/length(sz_phases_c);
        sz_phases_c = dh2(SzDayID-1); sz_phases_c = sz_phases_c./abs(sz_phases_c);
        SI2 = sum(sz_phases_c); SI2 = abs(SI2)/length(sz_phases_c);
        sz_phases_c = dh3(SzDayID-1); sz_phases_c = sz_phases_c./abs(sz_phases_c);
        SI3 = sum(sz_phases_c); SI3 = abs(SI3)/length(sz_phases_c);
        sz_phases_c = dh4(SzDayID-1); sz_phases_c = sz_phases_c./abs(sz_phases_c);
        SI4= sum(sz_phases_c); SI4= abs(SI4)/length(sz_phases_c);
        sz_phases_c = dh5(SzDayID-1); sz_phases_c = sz_phases_c./abs(sz_phases_c);
        SI5 = sum(sz_phases_c); SI5 = abs(SI5)/length(sz_phases_c);
        sz_phases_c = dh6(SzDayID-1); sz_phases_c = sz_phases_c./abs(sz_phases_c);
        SI6= sum(sz_phases_c); SI6= abs(SI6)/length(sz_phases_c);
            
        SzStart = 10;
        
        P_sp1 = ProbDist(dha1(1:SzDayID(SzStart)), dha1(SzDayID(1:SzStart)-1));
        P_sp2 = ProbDist(dha2(1:SzDayID(SzStart)), dha2(SzDayID(1:SzStart)-1));
        P_sp3 = ProbDist(dha3(1:SzDayID(SzStart)), dha3(SzDayID(1:SzStart)-1));
        P_sp4 = ProbDist(dha4(1:SzDayID(SzStart)), dha4(SzDayID(1:SzStart)-1));
        P_sp5 = ProbDist(dha5(1:SzDayID(SzStart)), dha5(SzDayID(1:SzStart)-1));
        P_sp6 = ProbDist(dha6(1:SzDayID(SzStart)), dha6(SzDayID(1:SzStart)-1));
        [~,~,~, ~,P_time1] = TimeProbability(dh1(1:SzDayID(SzStart)), P_sp1, 0);
        [~,~,~, ~,P_time2] = TimeProbability(dh2(1:SzDayID(SzStart)), P_sp2, 0);
        [~,~,~, ~,P_time3] = TimeProbability(dh3(1:SzDayID(SzStart)), P_sp3, 0);
        [~,~,~, ~,P_time4] = TimeProbability(dh4(1:SzDayID(SzStart)), P_sp4, 0);
        [~,~,~, ~,P_time5] = TimeProbability(dh5(1:SzDayID(SzStart)), P_sp5, 0);
        [~,~,~, ~,P_time6] = TimeProbability(dh6(1:SzDayID(SzStart)), P_sp6, 0);
        
        %assumption of independence: combine probabilities through
        %multiplication
        Ptime_slowI = P_time1.*P_time2;
        Ptime_fastI = P_time3.*P_time4;
        Ptime_spikesI = P_time5.*P_time6;
        Ptime_combI_ = Ptime_slowI.*Ptime_fastI;
        Ptime_combI_ = Ptime_combI_.*Ptime_spikesI;
        Ptime_combI_ = Ptime_combI_/max(Ptime_combI_);
        
        %find optimal thresholds to begin with
        [xth,~,~,~, ~]  = ProbTimeEval(Ptime_combI_, SzDayID(1:SzStart)-1);
        if(xth(1) == 1e4), xth(1) = max(Ptime_combI_)/30,xth(2) = max(Ptime_combI_)/20; end
       
        % divide into three risk states
        PriskI = zeros(1,SzDayID(end)-1);
        PriskI(Ptime_combI_<xth(1)) = 1;PriskI(Ptime_combI_>=xth(2)) = 3;
        PriskI(Ptime_combI_>=xth(1) & Ptime_combI_<xth(2)) = 2;
                
        PRiskSz = zeros(1, length(SzDayID));
        PRiskSz(1:SzStart) = PriskI(SzDayID(1:SzStart) - 1); %sample prior to the seizures
        Ptime_combI = zeros(1, SzDayID(end));
        Ptime_combI(1:SzDayID(SzStart)) = Ptime_combI_;
        
        Pthresh = zeros(2, SzDayID(end));
        Pthresh(1,1:SzDayID(SzStart)) = xth(1);
        Pthresh(2,1:SzDayID(SzStart)) = xth(2);
        
        dhasz1 = nan(1,length(SzDayID));dhasz1(1:SzStart) = dha1(SzDayID(1:SzStart)-1);
        dhasz2 = nan(1,length(SzDayID));dhasz2(1:SzStart) = dha2(SzDayID(1:SzStart)-1);
        dhasz3 = nan(1,length(SzDayID));dhasz3(1:SzStart) = dha3(SzDayID(1:SzStart)-1);
        dhasz4 = nan(1,length(SzDayID));dhasz4(1:SzStart) = dha4(SzDayID(1:SzStart)-1);
        dhasz5 = nan(1,length(SzDayID));dhasz5(1:SzStart) = dha5(SzDayID(1:SzStart)-1);
        dhasz6 = nan(1,length(SzDayID));dhasz6(1:SzStart) = dha6(SzDayID(1:SzStart)-1);
        PLOT = 0;
        
        % now compute sizure probability iteratively
        for sz = SzStart+1:length(SzDayID)
            sz
            
            dataSig1 = (DataACFWidth(1:SzDayID(sz), ch1,7));
            dataSig2 = (DataVariance(1:SzDayID(sz), ch2,7));
            dataSig3 = (DataACFWidth(1:SzDayID(sz), ch3,7));
            dataSig4 = (DataVariance(1:SzDayID(sz), ch4,7));            
            dataSig5 = (Rates(ch5,1:SzDayID(sz)));%long
            dataSig6 = (Rates(ch6,1:SzDayID(sz)));%short
            

            CausalSigACFSlow = movmean(dataSig1,[720*2 0]) - movmean(dataSig1, [720*30 0]);
            CausalSigVARSlow = movmean(dataSig2,[720*2 0]) - movmean(dataSig2, [720*30 0]);
            CausalSigACFFast = movmean(dataSig3,[Mlen 0]) - movmean(dataSig3,[720*2 0]);
            CausalSigVARFast = movmean(dataSig4,[Mlen 0]) - movmean(dataSig4,[720*2 0]);
            CausalSigSPSlow = movmean(dataSig5,[720*2 0]) - movmean(dataSig5, [720*30 0]);
            CausalSigSPFast = movmean(dataSig6,[Mlen 0]) - movmean(dataSig6,[720*2 0]);
        
            d1 = CausalSigACFSlow;
            d2 = CausalSigVARSlow;
            d3 = CausalSigACFFast;
            d4 = CausalSigVARFast;
            d5 = CausalSigSPSlow;
            d6 = CausalSigSPFast;
            
            dh1 = hilbert(d1 - mean(d1));dha1 = angle(dh1);
            dh2 = hilbert(d2 - mean(d2));dha2 = angle(dh2);
            dh3 = hilbert(d3 - mean(d3));dha3 = angle(dh3);
            dh4 = hilbert(d4 - mean(d4));dha4 = angle(dh4);
            dh5 = hilbert(d5 - mean(d5));dha5 = angle(dh5);
            dh6 = hilbert(d6 - mean(d6));dha6 = angle(dh6);

            if(PLOT)
                figure(1);clf;
                subplot(4,1,1)
                plot(x_days(1:SzDayID(sz)), abs(dh1), 'b'); hold on;
                scatter(x_days(SzDayID(1:sz)), ones(1,sz)*pi);
                subplot(4,1,2)
                plot(x_days(1:SzDayID(sz)), abs(dh2), 'b'); hold on;
                scatter(x_days(SzDayID(1:sz)), ones(1,sz)*pi);
                subplot(4,1,3)
                plot(x_days(1:SzDayID(sz)), abs(dh3), 'b'); hold on;     
                scatter(x_days(SzDayID(1:sz)), ones(1,sz)*pi);
                subplot(4,1,4)
                plot(x_days(1:SzDayID(sz)), abs(dh4), 'b'); hold on;
                scatter(x_days(SzDayID(1:sz)), ones(1,sz)*pi);

                figure(2);clf;
                subplot(4,1,1)
                plot(x_days(1:SzDayID(sz)), dha1,'k');hold on;
                scatter(x_days(SzDayID(1:sz)), ones(1,sz)*pi);
                subplot(4,1,2)
                plot(x_days(1:SzDayID(sz)), dha2,'k');hold on;
                scatter(x_days(SzDayID(1:sz)), ones(1,sz)*pi);
                subplot(4,1,3)
                plot(x_days(1:SzDayID(sz)), dha3,'k');hold on;
                scatter(x_days(SzDayID(1:sz)), ones(1,sz)*pi);
                subplot(4,1,4)
                plot(x_days(1:SzDayID(sz)), dha4,'k');hold on;
                scatter(x_days(SzDayID(1:sz)), ones(1,sz)*pi);
            end
            
            dhasz1(1:sz) = dha1(SzDayID(1:sz)-1);
        	dhasz2(1:sz) = dha2(SzDayID(1:sz)-1);
        	dhasz3(1:sz) = dha3(SzDayID(1:sz)-1);
            dhasz4(1:sz) = dha4(SzDayID(1:sz)-1);
        	dhasz5(1:sz) = dha5(SzDayID(1:sz)-1);
        	dhasz6(1:sz) = dha6(SzDayID(1:sz)-1);
            
            P_sp1 = ProbDist(dha1(1:SzDayID(sz)), dhasz1);
            P_sp2 = ProbDist(dha2(1:SzDayID(sz)), dhasz2);
            P_sp3 = ProbDist(dha3(1:SzDayID(sz)), dhasz3);
            P_sp4 = ProbDist(dha4(1:SzDayID(sz)), dhasz4);
            P_sp5 = ProbDist(dha5(1:SzDayID(sz)), dhasz5);
            P_sp6 = ProbDist(dha6(1:SzDayID(sz)), dhasz6);
            
            [~,~,~, ~,P_time1] = TimeProbability(dh1(1:SzDayID(sz)), P_sp1, 0);
            [~,~,~, ~,P_time2] = TimeProbability(dh2(1:SzDayID(sz)), P_sp2, 0);
            [~,~,~, ~,P_time3] = TimeProbability(dh3(1:SzDayID(sz)), P_sp3, 0);
            [~,~,~, ~,P_time4] = TimeProbability(dh4(1:SzDayID(sz)), P_sp4, 0);
            [~,~,~, ~,P_time5] = TimeProbability(dh5(1:SzDayID(sz)), P_sp5, 0);
            [~,~,~, ~,P_time6] = TimeProbability(dh6(1:SzDayID(sz)), P_sp6, 0);
            
            %assumption of independence: combine probabilities through
            %multiplication
            Ptime_slowI = P_time1.*P_time2;
            Ptime_fastI = P_time3.*P_time4; 
            Ptime_spikesI = P_time5.*P_time6;
            Ptime_combI_ = Ptime_slowI.*Ptime_fastI;
            Ptime_combI_ = Ptime_combI_.*Ptime_spikesI;
            %Ptime_combI_ = Ptime_combI_/max(Ptime_combI_);
            

            prisk = zeros(1,length(Ptime_combI_));
            prisk(Ptime_combI_<xth(1)) = 1;prisk(Ptime_combI_>=xth(2)) = 3;
            prisk(Ptime_combI_>=xth(1) & Ptime_combI_<xth(2)) = 2;
            
            Pthresh(1,SzDayID(sz-1)+1:SzDayID(sz)) = xth(1);
            Pthresh(2,SzDayID(sz-1)+1:SzDayID(sz)) = xth(2);
            assignin('base', 'Pthresh', Pthresh)
            
            Ptime_combI(SzDayID(sz-1)+1:SzDayID(sz)) = Ptime_combI_(SzDayID(sz-1)+1:end);  
            PriskI(SzDayID(sz-1)+1:SzDayID(sz)) = prisk(SzDayID(sz-1)+1:end);            
            PRiskSz(sz) = prisk(end-1); %sample prior to the seizures
            MeanRisk = mean(nonzeros(PRiskSz))
            
            if(PLOT)
                figure(3);clf;
                subplot(2,1,1);
                pplot = Ptime_combI(1:SzDayID(sz)); pplot(pplot==0) = 1e-20;
                plot(x_days(1:SzDayID(sz)), log(pplot), 'k');hold on;
                scatter(x_days(SzDayID(1:sz)), 0*ones(1,sz),150, 'vr','filled');
                grid on
                subplot(2,1,2);
                plot(x_days(1:SzDayID(sz)), PriskI(1:SzDayID(sz)));hold on;
                scatter(x_days(SzDayID(1:sz)), 3.3*ones(1,sz),150, 'vr','filled');
                grid on
                drawnow
                pause
            end
            
            %update thresholds
            if(sz ~= length(SzDayID))
                [xth_,~,~,~, ~]  = ProbTimeEval(Ptime_combI_, SzDayID(1:sz));
                xth_(1) = xth_(1)/2;
                if(xth_(1) ~= 10000 && xth_(2) ~= 10000)
                    xth = xth_;
                end
            end
        end
        
        PR = PriskI(SzDayID -1);
        figure(2);clf;        
        plot(1:length(SzDayID), movmean(PR,5),'k','LineWidth',4)
        yticks(1:3); xticks([0 length(SzDayID)]); box off;
        yticklabels({});xticklabels({});
        xlim([0 length(SzDayID)+1]); ylim([0.8 3.2]);
        set(gca,'FontName', 'Arial', 'FontSize', 32);
        set(gcf, 'OuterPosition', [100 100 800 800]);
        set(gca, 'YGrid', 'on');
        
                
        figure(1);clf;
        subplot(2,1,1);
        plot(x_days(1:SzDayID(end)), Ptime_combI, 'k'); hold on;
        scatter(x_days(SzDayID), max(Ptime_combI)*ones(1,length(SzDayID)),150, 'vr','filled');
        plot([x_days(1) x_days(end)], [xth(1) xth(1)],'r')
        plot([x_days(1) x_days(end)], [xth(2) xth(2)],'--r')
        grid on
        subplot(2,1,2);
        plot(x_days(1:SzDayID(end)), PriskI);hold on;
        scatter(x_days(SzDayID), 3.3*ones(1,length(SzDayID)),150, 'vr','filled');
        grid on
        drawnow
        
        PR1 = PR(PR==1);
        PR2 = PR(PR==2);
        PR3 = PR(PR==3);
        
        PrSzHigh = length(PR3)/length(SzDayID)
        PrSzMed = length(PR2)/length(SzDayID)
        PrSzLow = length(PR1)/length(SzDayID)
        
        PrTimeHigh = length(find(PriskI==3))/length(PriskI)
        PrTimeMed = length(find(PriskI==2))/length(PriskI)
        PrTimeLow = length(find(PriskI==1))/length(PriskI)
        
        IterativeSummary = [PrSzHigh PrSzLow PrTimeHigh PrTimeLow];
        
        
    end
%-----------------------------------------------------------------------
    function [PcSzHigh1,PcSzLow1,PcTimeLow1, PcTimeHigh1,...
            P_time1] = TimeProbability(hbSig1, P_sp1, PLOT)
        
        %compute optimal risk threshold and time spent in each risk
        %category
        
        P_time1=zeros(1,length(hbSig1));
        phase_time1 = angle(hbSig1);

        for tt =1:length(hbSig1)        
            a1 = phase_time1(tt);
            if(isnan(a1) ), continue; end
            [n1, x_1] = histc(a1,anglebins);
            i1 = find(n1);
            p1 = P_sp1(i1);
            P_time1(tt) = p1;
        end
        
        % return optimal risk thresholds (x) and time othe stats
        [x,PcSzHigh1,PcSzLow1,PcTimeLow1, PcTimeHigh1]  = ProbTimeEval(P_time1, SzDayID);
        ACF_TH1 = x(1);
        ACF_TH2 = x(2);

        if(PLOT)
            figure(1); clf;
            subplot(2,1,1);
            plot(x_days, P_time1);hold on;
            scatter(x_days(SzDayID), max(P_time1)*ones(1,length(SzDayID)),150, 'vr','filled');
            plot([x_days(1) x_days(end)], [ACF_TH1 ACF_TH1],'r')
            plot([x_days(1) x_days(end)], [ACF_TH2 ACF_TH2],'--r')
            grid on
        end

    end
%-----------------------------------------------------------------------
    function [elecACF, elecVAR,elecSpikes, hbsig1, hbsig2,hbsig3,...
            P_sp1, P_sp2, P_sp3] = forecast2(sig1, sig2,sig3)  % short period, causal
        
        SI = [];        
        for ch=1:16
            dataSig1 = sig1(:,ch);
            dataSig2 = sig2(:,ch);
            dataSigSpikes = sig3(ch,:);
            
            %subtract causal 2-day moving average
            Dd1 = dataSig1 - movmean(dataSig1,[720*2 0]); 
            Dd2 = dataSig2 - movmean(dataSig2,[720*2 0]);
            Dd3 = dataSigSpikes - movmean(dataSigSpikes,[720*2 0]);
            
            %compute hilbert transform            
            hbsig1 = hilbert(Dd1-mean(Dd1));
            hasig1 = angle(hbsig1);
            hbsig2 = hilbert(Dd2-mean(Dd2));
            hasig2 = angle(hbsig2);
            hbsig3 = hilbert(Dd3-mean(Dd3));
            hasig3 = angle(hbsig3);
           
            %compute phases relative to seizure times
            sz_phases1 = hasig1(SzDayID-1);
            sz_phases_c = hbsig1(SzDayID-1); sz_phases_c = sz_phases_c./abs(sz_phases_c);
            SI1 = sum(sz_phases_c); SI1 = abs(SI1)/length(sz_phases_c);
            sig_phase1 = hasig1; 
            
            sz_phases2 = hasig2(SzDayID-1);
            sz_phases_c = hbsig2(SzDayID-1);sz_phases_c = sz_phases_c./abs(sz_phases_c);
            SI2 = sum(sz_phases_c); SI2 = abs(SI2)/length(sz_phases_c);
            sig_phase2 = hasig2; 
            
            sz_phases3 = hasig3(SzDayID-1);
            sz_phases_c = hbsig3(SzDayID-1);sz_phases_c = sz_phases_c./abs(sz_phases_c);
            SI3 = sum(sz_phases_c); SI3 = abs(SI3)/length(sz_phases_c);
            sig_phase3 = hasig3;
            
            SI(ch, 1) = SI1;
            SI(ch, 2) = SI2;
            SI(ch, 3) = SI3;
        end
        
        SI_shortACF(iPt, :) = SI(:,1)';
        SI_shortVAR(iPt, :) = SI(:,2)';
        SI_shortSP(iPt, :) = SI(:,3)';
        
        [SI1, elecACF] = max(SI(:,1));
        [SI2, elecVAR] = max(SI(:,2));
        [SI3, elecSpikes] = max(SI(:,3));
        
        dataSig1 = sig1(:,elecACF);
        dataSig2 = sig2(:,elecVAR);
        dataSigSpikes = sig3(elecSpikes,:);
        
        %recompute hilbert transform for the signal with highest SI
        Dd1 = dataSig1 - movmean(dataSig1,[720*2 0]);            
        Dd2 = dataSig2 - movmean(dataSig2,[720*2 0]);
        Dd3 = dataSigSpikes - movmean(dataSigSpikes,[720*2 0]);
        hbsig1 = hilbert(Dd1-mean(Dd1));
        hasig1 = angle(hbsig1);
        hbsig2 = hilbert(Dd2-mean(Dd2));
        hasig2 = angle(hbsig2);
        hbsig3 = hilbert(Dd3-mean(Dd3));
        hasig3 = angle(hbsig3);
        sz_phases1 = hasig1(SzDayID-1);
        sz_phases2 = hasig2(SzDayID-1);
        sz_phases3 = hasig3(SzDayID-1);
        
        % compute the mean cycle time
        anan = isnan(dO(:,elecACF,7));
        C1 = MeanCycleTime(hbsig1,2,anan);
        anan = isnan(dO(:,elecVAR,7));
        C2 = MeanCycleTime(hbsig2,2,anan);
        C3 = MeanCycleTime(hbsig3,2,anan);

        A1 = hbsig1(SzDayID-1);A1 = angle(sum(A1)/length(SzDayID));
        A2 = hbsig2(SzDayID-1);A2 = angle(sum(A2)/length(SzDayID));
        A3 = hbsig3(SzDayID-1);A3 = angle(sum(A3)/length(SzDayID));
        
        P_sp1 = ProbDist(hasig1,sz_phases1);
        P_sp2 = ProbDist(hasig2,sz_phases2);
        P_sp3 = ProbDist(hasig3,sz_phases3);
        
        %compute synchronisation indices
        SI1A = hbsig1./abs(hbsig1);
        SI1A = sum(SI1A);
        SI1A = abs(SI1A)/length(hbsig1);
        SI2A = hbsig2./abs(hbsig2);
        SI2A = sum(SI2A); 
        SI2A = abs(SI2A)/length(hbsig2);
        SI3A = hbsig3./abs(hbsig3);
        SI3A = sum(SI3A); 
        SI3A = abs(SI3A)/length(hbsig3);
        
        PtSI(iPt,4:6) = [SI1, SI2, SI3];
        
        % plot polar diagrams
        col = [0.2 0.3 0.49];
        PolarPlot(sz_phases1, hasig1, SI1, col, 2,4,SI1A,C1)
        col = [0.8 0. 0.2];
        PolarPlot(sz_phases2, hasig2, SI2, col, 2,5,SI2A,C2)
        col = [0.0 0.75 0.75];
        PolarPlot(sz_phases3, hasig3, SI3, col, 2,6,SI3A,C3)
       
        MEAN_SI_PHASE(iPt,4:6) = [A1 A2 A3];
    end
%-----------------------------------------------------------------------
    function [elecACF, elecVAR,elecSpikes, hbsig1, hbsig2, hbsig3,...
            P_sp1, P_sp2, P_sp3] = forecast1(sig1, sig2,sig3) % long period, causal
        
        SI = [];
        for ch=1:16
            dataSig1 = sig1(:,ch); % autocorr width
            dataSig2 = sig2(:,ch); % variance
            dataSig3 = sig3(ch,:); % spike rate
            
            % compute the 2-day moving average: causal
            Dd1 = movmean(dataSig1,[720*2 0]); 
            Dd2 = movmean(dataSig2,[720*2 0]); 
            Dd3 = movmean(dataSig3,[720*2 0]); 
            
            % compute the hilbert transform
            hbsig1 = hilbert(Dd1-mean(Dd1));
            hasig1 = angle(hbsig1);
            hbsig2 = hilbert(Dd2-mean(Dd2));
            hasig2 = angle(hbsig2);
            hbsig3 = hilbert(Dd3-mean(Dd3));
            hasig3 = angle(hbsig3);
            
            % find the signal phase relative to sample prior to sizures
            sz_phases_c = hbsig1(SzDayID-1); sz_phases_c = sz_phases_c./abs(sz_phases_c);
            SI1 = sum(sz_phases_c); SI1 = abs(SI1)/length(sz_phases_c);
            sig_phase1 = hasig1; 
            
            sz_phases_c = hbsig2(SzDayID-1);sz_phases_c = sz_phases_c./abs(sz_phases_c);
            SI2 = sum(sz_phases_c); SI2 = abs(SI2)/length(sz_phases_c);
            sig_phase2 = hasig2; 
            
            sz_phases_c = hbsig3(SzDayID-1);sz_phases_c = sz_phases_c./abs(sz_phases_c);
            SI3 = sum(sz_phases_c); SI3 = abs(SI3)/length(sz_phases_c);
            sig_phase3 = hasig3; 
            
            SI(ch, 1) = SI1;
            SI(ch, 2) = SI2;
            SI(ch, 3) = SI3;
        end
        
        SI_longACF(iPt,:) = SI(:,1);
        SI_longVAR(iPt,:) = SI(:,2);
        SI_longSP(iPt, :) = SI(:,3)';
        [SI1, elecACF] = max(SI(:,1));
        [SI2, elecVAR] = max(SI(:,2));
        [SI3, elecSpikes] = max(SI(:,3));
        
        % recompute phases for the electrode with highest SI
        dataSig1 = sig1(:,elecACF);
        dataSig2 = sig2(:,elecVAR);
        dataSigSpikes = sig3(elecSpikes,:);
        
        Dd1 = movmean(dataSig1,[720*2 0]);
        Dd2 = movmean(dataSig2,[720*2 0]);
        Dd3 = movmean(dataSigSpikes,[720*2 0]);
            
        hbsig1 = hilbert(Dd1-mean(Dd1));
        hasig1 = angle(hbsig1);
        hbsig2 = hilbert(Dd2-mean(Dd2));
        hasig2 = angle(hbsig2);
        hbsig3 = hilbert(Dd3-mean(Dd3));
        hasig3 = angle(hbsig3);
        sz_phases1 = hasig1(SzDayID-1);
        sz_phases2 = hasig2(SzDayID-1);
        sz_phases3 = hasig3(SzDayID-1);
        
        % compute the mean cycle time
        anan = isnan(dO(:,elecACF,7));
        C1 = MeanCycleTime(hbsig1,1,anan);assignin('base', 'C1', C1)
        anan = isnan(dO(:,elecVAR,7));
        C2 = MeanCycleTime(hbsig2,1,anan);assignin('base', 'C2', C2)
        C3 = MeanCycleTime(hbsig3,1,anan);assignin('base', 'C3', C3)
        
        A1 = hbsig1(SzDayID-1);A1 = angle(sum(A1)/length(SzDayID));
        A2 = hbsig2(SzDayID-1);A2 = angle(sum(A2)/length(SzDayID));
        A3 = hbsig3(SzDayID-1);A3 = angle(sum(A3)/length(SzDayID));
        
        % compute the probability distribution of seizure relative to phase
        P_sp1 = ProbDist(hasig1,sz_phases1);
        P_sp2 = ProbDist(hasig2,sz_phases2);
        P_sp3 = ProbDist(hasig3,sz_phases3);
        
        %compute synchronisation indices
        SI1A = hbsig1./abs(hbsig1); 
        SI1A = sum(SI1A); 
        SI1A = abs(SI1A)/length(hbsig1);
        SI2A = hbsig2./abs(hbsig2);
        SI2A = sum(SI2A); 
        SI2A = abs(SI2A)/length(hbsig2);
        SI3A = hbsig3./abs(hbsig3);
        SI3A = sum(SI3A); 
        SI3A = abs(SI3A)/length(hbsig3);
        
        % plot polar diagrams            
        col = [0.2 0.3 0.49];
        PolarPlot(sz_phases1, hasig1, SI1, col, 2,1,SI1A,C1)
        col = [0.8 0. 0.2];
        PolarPlot(sz_phases2, hasig2, SI2, col, 2,2,SI2A,C2)        
        col = [0.0 0.75 0.75];
        PolarPlot(sz_phases3, hasig3, SI3, col, 2,3,SI3A,C3)
        
        PtSI(iPt,1:3) = [SI1, SI2, SI3];
        MEAN_SI_PHASE(iPt,1:3) = [A1 A2 A3];
        
        
    end
 
%-----------------------------------------------------------------------
    function C = MeanCycleTime(hbsig, T,anan)
        C = [];
        figure(15);clf;
        hasig = angle(hbsig);
        plot(x_days, hasig); hold on;
        hd = [0; abs(diff(hasig(:)))];
        hd(anan) = 0;
        
        if(T == 1) %long rhythms
            [pk, loc] = findpeaks(hd, 'MinPeakDist', 720*2,'MinPeakHeight', 4);
        elseif(T == 2)%short rhythms
            [pk, loc] = findpeaks(hd, 'MinPeakDist', 100,'MinPeakHeight', 4);
        end
        scatter(x_days(loc), hasig(loc));
        C(1) = mean(diff(x_days(loc)));
        C(2) = std(diff(x_days(loc)));
        plot(x_days, abs(hbsig),'k');
        
        
    end

%-----------------------------------------------------------------------
    function P_sp = ProbDist(hasig, sz_phases)
      
        L = length(SzDayID);
        for i=2:length(anglebins)
            a1 = anglebins(i-1);
            a2 = anglebins(i);
            A = find(hasig>= a1 & hasig<a2);
            asz = find(sz_phases>= a1 & sz_phases<a2);
            
            % ACF: prob phase given seizure
            p_phase(i-1) = length(A)/length(hasig);
            P_ps(i-1) = length(asz)/L;
            % prob of seizure=1 given phase
            P_sp(i-1) =  length(asz)/length(A);
            
        end
    end
%-----------------------------------------------------------------------
    function PolarPlot(sz_phase, sig_phase, SI, col, fignum, spn, pval, Cyc)
       
        print_pval=1;
        if(nargin==5)
            figure(fignum);
            print_pval =0; Cyc = [0 0];
        elseif(nargin ==6)
            figure(fignum); subplot(2,3,spn); cla;
            print_pval =0;Cyc = [0 0];
        elseif(nargin >= 7)
            figure(fignum); subplot(2,3,spn); cla;
        end  
        
        sigphaseplot=1;
        if(isempty(sig_phase))
            sigphaseplot = 0;
        end
        
        figure(fignum+10);clf
        c = histogram(sz_phase,anglebins);
        BinCounts1 = c.BinCounts;
        %mx1 = max((BinCounts1));
        %BinCounts1 = BinCounts1/mx1;   
        if(sigphaseplot)
            c = histogram(sig_phase,anglebins);
            BinCounts3 = c.BinCounts;
            BinCounts3=BinCounts3/max(BinCounts3);
        end
        %BinCounts1 = BinCounts1./BinCounts3; 
        BinCounts1 = BinCounts1/max(BinCounts1);
        
        close(gcf)
        
              
        
        %figure(fignum);clf;
   
        if(sigphaseplot)
        pax = polarhistogram('BinCounts',BinCounts3, 'BinEdges',anglebins,'EdgeColor',col,...
           'FaceColor','none','DisplayStyle', 'stairs','LineWidth',2,'FaceAlpha',0.7);hold on;
        end 
        pax = polarhistogram('BinCounts',BinCounts1, 'BinEdges',anglebins,'EdgeColor','k',...
           'FaceColor',col,'LineStyle', '-','LineWidth',2,'FaceAlpha',0.7);hold on;
        pax.Parent.ThetaTickLabel = [];
        pax.Parent.RTickLabel = [];
        pax.Parent.GridAlpha = 0.3;
        pax.Parent.RLim = [0 1.05];
        si = sprintf('%01.2f', SI);
        cyc1 = sprintf('%01.2f', Cyc(1));
        cyc2 = sprintf('%01.2f', Cyc(2));
        cyc_time = [ cyc1 ' ' char(177) ' ' cyc2];
        
        if(print_pval)            
%             title(['SI: ' si newline '{\it p} value: ' num2str(pval)],'fontsize',18)
            si2 = sprintf('%01.2f', pval)
            %title(['SI_1: ' si newline 'SI_2: ' si2],'fontsize',18)
            coltxt = [num2str(col(1)) ' ' num2str(col(2)) ' ' num2str(col(3))];
            title(['{\fontsize{18} {\color[rgb]{' coltxt '} SI_1:' si '}' newline,...
                '{\color[rgb]{0 0 0} SI_2:' si2 '} }'])
            text(0.2, -0.1, cyc_time, 'Units', 'normalized','FontSize',14)
        else
            title(['SI: ' si],'fontsize',18) 
        end
        set(gcf, 'OuterPosition', [200 100 800 800]);
        %figname = ['Ch:' num2str(ch) ', SI: ' num2str(SI)];
        set(gcf, 'OuterPosition', [200 100 800 800]);
    end
   
%-----------------------------------------------------------------------
    function [PcSzHigh,PcSzMed,PcSzLow,...
            PcTimeHigh, PcTimeMed, PcTimeLow] = RiskRatio(Prisk)
        
        PR = Prisk(SzDayID-P_Hor);
        PR1 = PR(PR==1);
        PR2 = PR(PR==2);
        PR3 = PR(PR==3);
        PcSzHigh = length(PR3)/length(SzDayID);
        PcSzMed = length(PR2)/length(SzDayID);
        PcSzLow = length(PR1)/length(SzDayID);
        PcTimeHigh = length(find(Prisk==3))/length(Prisk);
        PcTimeMed = length(find(Prisk==2))/length(Prisk);
        PcTimeLow = length(find(Prisk==1))/length(Prisk);
    end
end