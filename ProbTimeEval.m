function [x,PercSzInHigh,PercSzInLow,PercTimeInLow,PercTimeInHigh]  = ProbTimeEval(Ptime, SzInd)
    x(1) = 10000;
    x(2) = 10000;

    if(SzInd(end) > length(Ptime))
        a = find(SzInd<=length(Ptime));
        SzInd = SzInd(a);
    end
    PercSzInLow=[];
    PercSzInMed=[];
    PercSzInHigh=[];
    PercTimeInLow=[];
    PercTimeInMed=[];
    PercTimeInHigh=[];
    
  
    
    u = sort(unique(Ptime));
    if(u(1) == 0), u(1) = []; end
    amin = floor(log10(min(u)));
    if(isempty(amin)), amin = -10; end
    th1 = u - min(u)/10;
    
    if(length(th1)>50)
        nl = floor(length(th1)/50);
        th1 = th1(1:nl:end);
    else
        th1 = logspace(amin, 0, 50);
    end
    
    th = th1;
    
    %brute force search through thresholds
    F = zeros(length(th), length(th));
    for ii=1:length(th)
        for jj=1:length(th)
            if(th(jj)<=th(ii))
                continue;
            end
            f = ftime([th(ii) th(jj)]);
            F(ii, jj) = f;
        end
    end

    if(max(F(:)) == 0)
        PercSzInLow=1;
        PercSzInMed=0;
        PercSzInHigh=0;
        PercTimeInLow=0;
        PercTimeInMed=0;
        PercTimeInHigh=1;
        x(1) = 10000;%max(Ptime)/3;
        x(2) = 10000;%2*max(Ptime)/3;
    else

        [maxValue, linearIndexesOfMaxes] = max(F(:));
        [rw cl] = find(F == maxValue);
        x(1) = th(rw(1));
        x(2) = th(cl(1));

        f = ftime([x(1) x(2)]);
    end
    
    % fmincon search for optimal thresholds
%     x0 = [0 max(Ptime)];
%     A = [1 -1];
%     b = 0;
%     lb = [0 0];
%     ub = [max(Ptime) max(Ptime)];
%     %x = fmincon(@ftime, x0);
%     x = fmincon(@ftime, x0, [], [],[],[],lb,ub);
    
    %% function to maximise
    function F = ftime(xin)
        
        TH1 = xin(1);
        TH2 = xin(2);
        a = (Ptime<TH1);
        szLow = a(SzInd);
        numSzInLow = sum(szLow);
        PercSzInLow = numSzInLow/length(SzInd);
        PercTimeInLow = sum(a)/length(Ptime);
        
        a = (Ptime>=TH2);
        szHigh = a(SzInd);
        numSzInHigh = sum(szHigh);
        PercSzInHigh = numSzInHigh/length(SzInd);
        PercTimeInHigh = sum(a)/length(Ptime);
        
        numSzInMed = length(SzInd) - (numSzInHigh+numSzInLow);
        PercSzInMed = numSzInMed/length(SzInd);
        PercTimeInMed = 1 - (PercTimeInHigh + PercTimeInLow);
        
        if(PercTimeInLow < PercTimeInMed || PercTimeInMed < PercTimeInHigh)
            F = 0.5*(PercTimeInLow) * (PercSzInHigh) * (1 - PercTimeInHigh);
        elseif(numSzInLow > numSzInMed || numSzInMed > numSzInHigh)
            F = 0.5*(PercTimeInLow) * (PercSzInHigh) * (1 - PercTimeInHigh);
        else
            F = (PercTimeInLow) * (PercSzInHigh) * (1 - PercTimeInHigh);
        end
        
        
    end

end