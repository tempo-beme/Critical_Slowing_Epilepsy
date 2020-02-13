% function to download data from iEEG.org and compute the signal features; 
% autcorrelation, signal energy, variance
% iPt is the patient number


function CSD_compute_features(iPt) 

% add the ieeg.org matlab library
addpath(genpath('ieeg-matlab-1.13.2'));

%ieeg.org login details
login = '';
pword = '.bin';


PLOT = 0; %1 to plot data while downloading, 0 to not plot
SAVEDATA = 1; %1 to save data, 0 to skip saving

% load information
% Patients identifier
Patient{iPt} = 'insert_identifier';
curPt = Patient{iPt};
patient = IEEGSession(['session_name'],login,pword);
Fs = patient.data.sampleRate;
dt = 1/Fs;

% load seizure timing information
szinfo = load(['Portal Annots/' curPt '_Annots']);
Ptinfo.SzDur = szinfo.SzDur;
Ptinfo.SzIndices = szinfo.SzIndices;
Ptinfo.SzType = szinfo.SzType;
% only use Type1 and 2 seizures
SzIndL = find( (Ptinfo.SzType==1 | Ptinfo.SzType==2) );
Ptinfo.SzIndL = Ptinfo.SzIndices(SzIndL);

pt = sprintf('%02.0f', iPt);
Block = floor(30*60*Fs);

%filters
d_all(1).dfilt = designfilt('lowpassiir','FilterOrder',12,'HalfPowerFrequency',8,'SampleRate',Fs); 
d_all(2).dfilt = designfilt('bandpassiir','FilterOrder',12,'HalfPowerFrequency1',8,'HalfPowerFrequency2',12,'SampleRate',Fs); 
d_all(3).dfilt = designfilt('bandpassiir','FilterOrder',12,'HalfPowerFrequency1',16,'HalfPowerFrequency2',24,'SampleRate',Fs);  
d_all(4).dfilt = designfilt('bandpassiir','FilterOrder',12,'HalfPowerFrequency1',24,'HalfPowerFrequency2',45,'SampleRate',Fs);  
d_all(5).dfilt = designfilt('bandpassiir','FilterOrder',12,'HalfPowerFrequency1',45,'HalfPowerFrequency2',80,'SampleRate',Fs);  
d_all(6).dfilt = designfilt('bandpassiir','FilterOrder',12,'HalfPowerFrequency1',80,'HalfPowerFrequency2',170,'SampleRate',Fs);  
%cutoff at 170Hz to cut interference when patients recharged their devices
d_all(7).dfilt = designfilt('lowpassiir','FilterOrder',12,'HalfPowerFrequency',160,'SampleRate',Fs); 


download_data();
    
   %%-------------------------------------------
    function download_data() 

        folder = ['CSD_data/Pt_' pt];
        ptfolder = folder;
        ptfile = [ptfolder '/DataTimeSeries.mat'];
        buffer = 3000;
        
        if(exist(ptfile, 'file'))
            %if file already exist, continue
            d = load(ptfile);
            DataTimeSeries = d.DataTimeSeries;
            BlockStart = DataTimeSeries.lastBlock;
            BlockEnd = BlockStart + Block -1;
            dataTime = DataTimeSeries.T;
            DataEnergy = DataTimeSeries.DataEnergy;
            DataACF = DataTimeSeries.DataACF;
            DataACFWidth = DataTimeSeries.DataACFWidth;
            DataVariance = DataTimeSeries.DataVariance;
            n = DataTimeSeries.n;
            lim = n+500;
        else
            % else begin from the beginning
            BlockStart = buffer+1;
            BlockEnd = BlockStart + Block -1;
            step = floor(Fs*60*2);
            dataTime = buffer+1+floor(Fs):step:imax;
            DataEnergy = nan(length(dataTime),16,7);
            DataACF = nan(length(dataTime),16,7);
            DataACFWidth = nan(length(dataTime),16,7);
            DataVariance = nan(length(dataTime),16,7);
            DataTimeSeries.T = dataTime;
            lim = 500;
        end
                

        while(BlockStart < imax)
            
            % get data in chunks of size Block
            idblock1 = BlockStart;
            idblock2 = BlockEnd;
            
            try
                % connect to iEEG and get data
                Data = patient.data.getvalues(idblock1-buffer:idblock2+buffer, 1:16);
            catch
                % if failed, try again after some time
                pause(120);
                try
                    Data = patient.data.getvalues(idblock1-buffer:idblock2+buffer, 1:16);
                catch
                    pause(300);
                    Data = patient.data.getvalues(idblock1-buffer:idblock2+buffer, 1:16);
                end
            end
            
            a = find(isnan(Data(:,1)));
            if(length(a) == length(Data(:,1)))
                BlockStart = BlockEnd + 1;
                BlockEnd = BlockStart + Block;
                disp('data full of nan');
                continue; 
            end
            
            % filter data
            FiltData = [];
            % only using low pass filter
            dfilt = d_all(7).dfilt;
            for ch =1:16
                d1 = Data(:,ch);
                anan = find(isnan(d1));
                bnan = find(~isnan(d1));
                d = d1;
                d(anan) = [];
                if(length(d)>40)%filtfilt restriction
                    df = filtfilt(dfilt, d);
                    d1 = nan(1,length(d1));
                    d1(bnan) = df;
                    d1 = d1(:,buffer+1:end-buffer);
                    FiltData(ff).chan(ch).data = d1(:);
                else
                    d1 = d1(buffer+1:end-buffer);
                    FiltData(ff).chan(ch).data = nan(length(d1),1);
                end
            end
            
            blockIds = idblock1:idblock2;
            if(PLOT)
                figure(1);clf
                plot(dt*(blockIds-BlockStart), FiltData(7).chan(1).data); hold on;
            end
                        
            a = find(ismember(blockIds,dataTime));
            
            for ii=1:length(a)
                dtime = blockIds(a(ii));
                n = find(dataTime == dtime);
                disp([num2str(n) ' of ' num2str(length(dataTime))]);
                
                id1 = a(ii)- floor(Fs);
                id2 = a(ii);
                if(id1<1), continue;end
                
                
                for ch=1:16
                    d = FiltData(7).chan(ch).data;
                    if(length(d) == 0), continue; end
                    d = d(id1:id2);
                    anan = find(isnan(d));
                    if(~isempty(anan))
                        %disp('d has nans')
                        continue; 
                    end
                    
                    d1 = DataEnergy(n,ch,ff);
                    d2 = DataACF(n,ch,ff);
                    d3 = DataACFWidth(n,ch,ff);
                    d4 = DataVariance(n,ch,ff);

                    % Compute features
                    Energy = sum(d.*d);
                    DataEnergy(n,ch,ff) = Energy;
                    
                    [r, lag] = autocorr(d, floor(Fs));
                    [width, height] = fwhm(lag,r,1);
                    ACF1 = r(2);
                    
                    %store data
                    DataACF(n,ch,ff) = ACF1;
                    DataACFWidth(n,ch,ff) = width;
                    sigvar = var(d);
                    DataVariance(n,ch,ff) = sigvar;

                    if(PLOT)
                        figure(1)
                        scatter(dt*(blockIds(id1)-BlockStart), 0.1, 'r');
                        figure(2);clf
                        tt = dt:dt:dt*length(d);
                        subplot(2,1,1)
                        plot(tt, d);
                        drawnow
                        [r, lag] = autocorr(d, floor(Fs));
                        [width, height] = fwhm(lag,r,1);
                        figure(2);
                        LAG = [-lag(:); lag(:)];LAG = unique(sort(LAG));
                        R = [fliplr(r(:)') r(2:end)'];
                        subplot(2,1,2)
                        plot(LAG, R); hold on;
                        x1 = ceil(length(LAG)/2);
                        X1 = -width; X2 = width;
                        plot([X1 X2], [0.5 0.5]);
                        pause
                    end
                end
                
            end
            
            
            if(SAVEDATA && n > lim)
                DataTimeSeries.lastBlock = BlockStart;
                DataTimeSeries.n=n;
                DataTimeSeries.DataEnergy = DataEnergy;
                DataTimeSeries.DataACF = DataACF;
                DataTimeSeries.DataACFWidth = DataACFWidth;
                DataTimeSeries.DataVariance= DataVariance;
                
                if(~exist(ptfolder, 'dir')), mkdir(ptfolder),end
                save(ptfile, 'DataTimeSeries');
				ptfile_backup = [ptfolder '/DataTimeSeries_backup.mat'];
				save(ptfile_backup, 'DataTimeSeries');
				
                lim = n+500;
            end
            BlockStart = BlockEnd + 1;
            BlockEnd = BlockStart + Block;
        end
        
        if(SAVEDATA)
            DataTimeSeries.lastBlock = BlockStart;
            DataTimeSeries.n=n;
            DataTimeSeries.DataEnergy = DataEnergy;
            DataTimeSeries.DataACF = DataACF;
            DataTimeSeries.DataACFWidth = DataACFWidth;
            DataTimeSeries.DataVariance= DataVariance;
            
            if(~exist(ptfolder, 'dir')), mkdir(ptfolder),end
            save(ptfile, 'DataTimeSeries');
            
            % compile data
            [DataEnergy,DataACF,DataACFWidth,DataVariance, T] =Data_compile(DataTimeSeries);
            dFs=720; %a sample every 2 mins
            
            % compute seizure phase relationships
            SzPhaseRelationship(DataACFWidth, DataACF, DataEnergy,...
                DataVariance,T,SzIndices,dFs,dt,pt);


            disp('sz phase relationships complete')
        end
    end

%%-------------------------------------------
    function [DataEnergy,DataACF,DataACFWidth,DataVariance, T, dFs] = Data_compile(DataTimeSeries)
        
        OutputFileName = ['Pow_data/Pt_' pt '/DataTSCompiled.mat'];
        T = DataTimeSeries.T;
        DataEnergy = DataTimeSeries.DataEnergy;
        DataACF = DataTimeSeries.DataACF;
        DataACFWidth = DataTimeSeries.DataACFWidth;
        DataVariance = DataTimeSeries.DataVariance;

        Acorr=nan(size(DataACF));
        Energy = nan(size(DataEnergy));
        ACFWidth = nan(size(DataACFWidth));
        Variance = nan(size(DataVariance));


        L = length(T);
        chan = 1:16;
        %fill gaps with noise
        for ff=7:7  % only using low pass filter
            for ch=chan
                d1 = DataEnergy(:, ch, ff);
                d2 = DataACF(:, ch, ff);
                d3 = DataACFWidth(:, ch, ff);
                d4 = DataVariance(:, ch, ff);

                a1 = find(isnan(d1)); 
                a2 = find(isnan(d2));
                a3 = find(isnan(d3));
                a4 = find(isnan(d4));

                std1 = nanstd(d1);
                std2 = nanstd(d2);
                std3 = nanstd(d3);
                std4 = nanstd(d4);
                
                m1 = nanmean(d1);
                m2 = nanmean(d2);
                m3 = nanmean(d3);
                m4 = nanmean(d4);

                %Energy
                gaps = find(diff(a1) > 1);
                gap_end = a1(gaps);
                gap_start = [a1(1); a1(gaps+1)];            
                d1out = d1;
                for ii=1:length(gap_end)
                    g1 = gap_start(ii);
                    g2 = gap_end(ii);

                    if(g1==g2 && g1 ~= 1 && g2 ~=length(d1out))
                        p1 = d1out(g1-1);
                        p2 = d1out(g2+1);
                        d1out(g1) = std1*rand+ m1;
                    elseif(g1==1)
                        p2 = d1out(g2+1);
                        tnew = g1:g2;
                        d1out(g1:g2) = std1*rand(1,length(tnew))+m1;
                    elseif(g2 == length(d1out))
                        p1 = d1out(g1-1);
                        tnew = g1:g2;
                        d1out(g1:g2) = std1*rand(1,length(tnew))+m1;
                    else
                        p1 = d1out(g1-1);
                        p2 = d1out(g2+1);
                        t = [g1-1 g2+1];
                        tnew = g1-1:g2+1;
                        xnew = interp1(t,[p1, p2],tnew); 
                        d1out(g1:g2) = std1*rand(1,length(tnew)-2) +  m1;
                    end                
                end
                a = find(isnan(d1out));d1out(a) = std1*rand(1,length(a)) + m1;
                
                %ACF
                gaps = find(diff(a2) > 1);
                gap_end = a2(gaps);
                gap_start = [a2(1); a2(gaps+1)];            
                d2out = d2;
                for ii=1:length(gap_end)
                    g1 = gap_start(ii);
                    g2 = gap_end(ii);

                    if(g1==g2)
                        p1 = d2out(g1-1);
                        p2 = d2out(g2+1);
                        d2out(g1) = std2*rand+ m2;
                    elseif(g1==1)
                        p2 = d2out(g2+1);
                        tnew = g1:g2;
                        d2out(g1:g2) = std2*rand(1,length(tnew))+m2;
                    elseif(g2 == length(d2out))
                        p1 = d2out(g1-1);
                        tnew = g1:g2;
                        d2out(g1:g2) = std2*rand(1,length(tnew))+m2;
                    else
                        p1 = d2out(g1-1);
                        p2 = d2out(g2+1);
                        t = [g1-1 g2+1];
                        tnew = g1-1:g2+1;
                        xnew = interp1(t,[p1, p2],tnew); 
                        d2out(g1:g2) = std2*rand(1,length(tnew)-2) +  m2;
                    end                
                end
                a = find(isnan(d2out));d2out(a) = std2*rand(1,length(a))+m2;
                
                %ACF width
                gaps = find(diff(a3) > 1);
                gap_end = a3(gaps);
                gap_start = [a3(1); a3(gaps+1)];            
                d3out = d3;
                for ii=1:length(gap_end)
                    g1 = gap_start(ii);
                    g2 = gap_end(ii);

                    if(g1==g2)
                        p1 = d3out(g1-1);
                        p2 = d3out(g2+1);
                        d3out(g1) = std3*rand +m3;
                    elseif(g1==1)
                        p2 = d3out(g2+1);
                        tnew = g1:g2;
                        d3out(g1:g2) = std3*rand(1,length(tnew))+ m3;
                    elseif(g2 == length(d3out))
                        p1 = d3out(g1-1);
                        tnew = g1:g2;
                        d3out(g1:g2) = std3*rand(1,length(tnew))+m3;
                    else
                        p1 = d3out(g1-1);
                        p2 = d3out(g2+1);
                        t = [g1-1 g2+1];
                        tnew = g1-1:g2+1;
                        xnew = interp1(t,[p1, p2],tnew); 
                        d3out(g1:g2) = std3*rand(1,length(tnew)-2) +  m3;
                    end                
                end
                a = find(isnan(d3out));d3out(a) = std3*rand(1,length(a))+m3;
                
                %Variance
                gaps = find(diff(a4) > 1);
                gap_end = a4(gaps);
                if(isempty(a4)), continue; end
                gap_start = [a4(1); a4(gaps+1)];            
                d4out = d4;
                for ii=1:length(gap_end)
                    g1 = gap_start(ii);
                    g2 = gap_end(ii);

                    if(g1==g2)
                        p1 = d4out(g1-1);
                        p2 = d4out(g2+1);
                        d4out(g1) = std4*rand + m4;
                    elseif(g1==1)
                        p2 = d4out(g2+1);
                        tnew = g1:g2;
                        d4out(g1:g2) = std4*rand(1,length(tnew))+m4;
                    elseif(g2 == length(d4out))
                        p1 = d4out(g1-1);
                        tnew = g1:g2;
                        d4out(g1:g2) = std4*rand(1,length(tnew))+m4;
                    else
                        p1 = d4out(g1-1);
                        p2 = d4out(g2+1);
                        t = [g1-1 g2+1];
                        tnew = g1-1:g2+1;
                        xnew = interp1(t,[p1, p2],tnew); 
                        d4out(g1:g2) = std4*rand(1,length(tnew)-2) +  m4;
                    end                
                end
                a = find(isnan(d4out));d4out(a) = std4*rand(1,length(a))+m4;

                Energy(:,ch,ff) = d1out;
                Acorr(:,ch,ff) = d2out;
                ACFWidth(:,ch,ff) = d3out;
                Variance(:,ch,ff) = d4out;

            end
        end

        DataEnergy = Energy;
        DataACF = Acorr;
        DataACFWidth = ACFWidth;
        DataVariance = Variance;


        save(OutputFileName, 'DataEnergy','DataACF','DataACFWidth','DataVariance', 'T');


        disp('Compile complete');

    end



end