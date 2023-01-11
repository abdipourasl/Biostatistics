%% Initialize parameters
Fs=100;                                  % sampling rate
t_length=75000;                              % data length (4 s)
TW=1:1:t_length;
TW_p=round(TW*Fs);
n_run=20;                                % number of used runs
sti_f=[12 8.57 5.45];             % stimulus frequencies 10, 9, 8, 6 Hz
n_sti=length(sti_f);                     % number of stimulus frequencies
n_correct2=zeros(23,length(TW));
n_correct3=zeros(23,length(TW));
n_correct4=zeros(23,length(TW));
n_correct5=zeros(23,length(TW));

load('mobile_scalp_dataset.mat') ;
%%
TEI2ave=0;
TEI3ave=0;
TEI4ave=0;
TEI5ave=0;

TEI2ave1=0;
TEI3ave1=0;
TEI4ave1=0;
TEI5ave1=0;

TEI2ave2=0;
TEI3ave2=0;
TEI4ave2=0;
TEI5ave2=0;



%%
ts = 0.01;
fs = 1/ts;
t = [0:ts:750];

% ETOTAL = (norm(data)^2) ;
% plot(data) ;
% title('Raw Data') ;
% figure ()
% xlabel('Time(T)');
% ylabel('Amplitude');
% grid on ;




% x1 = data
% na=length(x1) ;
% n=nextpow2(na) ;
% ffta=fft(x1,2^n);
% ampffta=abs(ffta(1:length(ffta)/2));
% f=1/2*linspace(0,1,length(ampffta)) ;
% plot(f,ampffta) ;
% xlabel('Frequency(Hz)');
% ylabel('Amplitude');
% title('FFT Of Raw Data');
% grid on;
% figure ();



% pdata = pwelch(data);
% plot(pdata) ;
% xlabel('Normalized Frequency');
% ylabel('Power/frequency');
% title('Welch power Spectral');
% grid on;
% figure()

%c)
% pxx = periodogram(data)
% plot(pxx) ;
% xlabel('Normalized Frequency');
% ylabel('Power/frequency');
% title('Priodogram power Spectral');
% grid on;
% figure()




for j=1:23
    for m=1:46
pTHETA(j).session2(m) = bandpower(Scalp(j).session2(m,:),Fs,[4 8]);

pALPHA(j).session2(m) = bandpower(Scalp(j).session2(m,:),Fs,[8 12]);

pBETA(j).session2(m) = bandpower(Scalp(j).session2(m,:),Fs,[12 30]);

pTHETA(j).session3(m) = bandpower(Scalp(j).session3(m,:),Fs,[4 8]);

pALPHA(j).session3(m) = bandpower(Scalp(j).session3(m,:),Fs,[8 12]);

pBETA(j).session3(m) = bandpower(Scalp(j).session3(m,:),Fs,[12 30]);

pTHETA(j).session4(m) = bandpower(Scalp(j).session4(m,:),Fs,[4 8]);

pALPHA(j).session4(m) = bandpower(Scalp(j).session4(m,:),Fs,[8 12]);

pBETA(j).session4(m) = bandpower(Scalp(j).session4(m,:),Fs,[12 30]);






TEI(j).session2(m) = pBETA(j).session2(m)/(pALPHA(j).session2(m) + pTHETA(j).session2(m));
TEI(j).session3(m) = pBETA(j).session3(m)/(pALPHA(j).session3(m) + pTHETA(j).session3(m));
TEI(j).session4(m) = pBETA(j).session4(m)/(pALPHA(j).session4(m) + pTHETA(j).session4(m));



    end 
    %c is for channels
    


end

for j=1:15
    for m=1:46
        
pTHETA(j).session5(m) = bandpower(Scalp(j).session5(m,:),Fs,[4 8]);

pALPHA(j).session5(m) = bandpower(Scalp(j).session5(m,:),Fs,[8 12]);

pBETA(j).session5(m) = bandpower(Scalp(j).session5(m,:),Fs,[12 30]);
TEI(j).session5(m) = pBETA(j).session5(m)/(pALPHA(j).session5(m) + pTHETA(j).session5(m));
    end
    
end



for m=1:46
pTHETA(18).session5(m) = bandpower(Scalp(j).session5(m,:),Fs,[4 8]);

pALPHA(18).session5(m) = bandpower(Scalp(j).session5(m,:),Fs,[8 12]);

pBETA(18).session5(m) = bandpower(Scalp(j).session5(m,:),Fs,[12 30]);
TEI(18).session5(m) = pBETA(j).session5(m)/(pALPHA(j).session5(m) + pTHETA(j).session5(m));
end


%%
clear all
close all
clc
load('C:/Users/ASUS/Desktop/BSc Thesis/BSc project/matlab.mat');
skew = skewness(tei(:,:));

for i=1:4
    if i<4
        figure();
        histogram(tei(:,i));
        h(i) = lillietest(tei(:,i));
        s(i) = jbtest(tei(:,i));
        st(i) = std(tei(:,i));
        m(i) = mean(tei(:,i));


    else
        figure();
        histogram(tei([1:15,18],i));
        h(i) = lillietest(tei([1:15,18],i));
        s(i) = jbtest(tei([1:15,18],i));
        st(i) = std(tei([1:15,18],i));
        m(i) = mean(tei([1:15,18],i));

    end
end
skew = skewness(tei(:,:));
figure();
histogram(tei);

figure();
hold on
plot(skew,'*-','linewidth',1.5) ; 
title('Skewness of data');
xlabel('Session(speed)');
ylabel('Skewness');
grid on ;

%% test

Spooled = ((22*(st(3)^2)+15*(st(4)^2))/37)^0.5;
T = (m(3)-m(4))/(Spooled*(((1/23)^0.5+(1/16)^0.5))^0.5)

[p,tbl,stats] = anova1(tei);