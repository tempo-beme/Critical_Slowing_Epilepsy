%% Critical slowing example
%create noisy signal
dt = 0.0005;
Tsamp = 1/dt; 
rng(1); %repeatable
y = randn(1,Tsamp);
Y1 = y; Y2 = y;
x = dt:dt:dt*Tsamp;

[pks, a] = findpeaks(abs(y), 'MinPeakHeight', 0.5);
tau1 = 1e-3; %system time constant when close to K1
tau2 = 10e-3; %system time constant when close to K2
xx = 0:dt:0.03; %signal lasting 10 ms
sig1 = exp(-xx./tau1);
sig2 = exp(-xx./tau2);
figure(1);clf; 
xx2 = 0:dt:dt*(-1+length(sig1)*2);
noise_std = 0.03;
noisysig1 = noise_std*randn(1,length(sig1)*2);
id1 = floor(length(sig1)/2); id2 = id1 + length(sig1)-1;
noisysig1(id1:id2) = noisysig1(id1:id2) - sig1; 
noisysig2 = noise_std*randn(1,length(sig2)*2);
id1 = floor(length(sig2)/2); id2 = id1 + length(sig2)-1;
noisysig2(id1:id2) = noisysig2(id1:id2) - sig2; 
plot(xx2,noisysig1);hold on;
plot(xx2, noisysig2);
for ii=1:length(a)
    id1 = a(ii);
    id2 = a(ii) + length(sig1) - 1;
    if(id2 > length(y)), continue; end
    if(y(a(ii))<0)
        Y1(id1:id2) = Y1(id1:id2) - sig1;
        Y2(id1:id2) = Y2(id1:id2) - sig2;
    else
        Y1(id1:id2) = Y1(id1:id2) + sig1;
        Y2(id1:id2) = Y2(id1:id2) + sig2;
    end
end

figure(2); clf
plot(x, Y1); hold on; plot(x, Y2+20)
var(Y1)
var(Y2)

[r1, lag1] = autocorr(Y1, length(Y1)-10);r1 = r1(:);
[width1, height1] = fwhm(lag1,r1,1);
[r2, lag2] = autocorr(Y2, length(Y2)-10);r2 = r2(:);
[width2, height2] = fwhm(lag2,r2,1);
LAG1 = [-lag1(:); lag1(:)];LAG1 = unique(sort(LAG1));
R1 = [flipud(r1(:)); r1(2:end)];
LAG2 = [-lag2(:); lag2(:)];LAG2 = unique(sort(LAG2));
R2 = [flipud(r2(:)); r2(2:end)];

figure(3); clf;
plot(LAG1, R1); hold on;
plot(LAG2, R2)
xlim([-200 200])



%%  bifrucation diagram Figure 6D

k = linspace(-1000, 1000, 300); %sort(unique(dataSig)); %autocorr
r = linspace(-1000, 1000, 300);%sort(unique(sigSpikes)); %spikerate


figure(4);clf
c3 = -1;
c2 = 1;
c1 = 0.1;
c0 = 1;
for kk=1:length(k)
    XS = [];
    XU = [];
    X =[];
    
    for rr=1:length(r)
        % dv/dt = k+rx-x^3
        x = roots([c3, c2, c1*(r(rr)-0), c0*(k(kk)-0)]);
        for jj=1:3
            if(isreal(x(jj)))
                X(end+1,:) = [r(rr) x(jj)];
            end
            if(isreal(x(jj)))
               a = x(jj);
               XS(end+1,:) = [r(rr) x(jj)];
               
            end
        end
        
    end
    
    if(~isempty(XS))
    	scatter3((k(kk))*ones(1,length(XS(:,1))),XS(:,1), XS(:,2),1,...
            'MarkerFaceColor',.7*[1 1 1],'MarkerEdgeColor','none');hold on;
    end
    if(kk == 149)
        
        plot3((k(kk))*ones(1,length(XS(:,1))),XS(:,1), XS(:,2), 'r.');hold on;
    end
    ax = gca;
    ax.YDir = 'reverse';
    %scatter(XS(:,1), XS(:,2), 'k')
    drawnow
end

XU= []; XS = [];
rr=244;
for kk=1:length(k)  
    x = roots([c3, c2, c1*(r(rr)-0), c0*(k(kk)-0)]);
    for jj=1:3
        if(isreal(x(jj)))
            X(end+1,:) = [k(kk) x(jj)];
        end
        if(isreal(x(jj)))
           a = x(jj);
           if( -3*a^2 + 0.1 > 1)  
               XU(end+1,:) = [k(kk) x(jj)];
           else
               XS(end+1,:) = [k(kk) x(jj)];
           end
        end
    end

end
scatter3(XS(:,1),(r(rr))*ones(1,length(XS(:,1))), XS(:,2), 'k.');hold on;  

xlabel('slow parameter (k)');
ylabel('fast parameter (r)')