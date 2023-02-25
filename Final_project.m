clear
close all
clc
tic
L=zeros(9,20);
LH=zeros(9,20);
H=zeros(9,20);
S=zeros(9,20);

opw_ttest1 =H;
opw_ttest2=H;
opw_wilcoxon2=H;
opw_wilcoxon1=H;
oP_ttest1=H;
oP_ttest2=H;
oP_wilcoxon1=H;
oP_wilcoxon2=H;
opw_LH1=H;
opw_LH2=H;
opw_H1=H;
opw_H2=H;
opw_L1=H;
opw_L2=H;
oP_LH1=H;
oP_LH2=H;
oP_H1=H;
oP_H2=H;
oP_L1=H;
oP_L2=H;
%%
for q=1:1000
    
    pw_ttest1 =S;
pw_ttest2=S;
pw_wilcoxon2=S;
pw_wilcoxon1=S;
P_ttest1=S;
P_ttest2=S;
P_wilcoxon1=S;
P_wilcoxon2=S;
pw_LH1=S;
pw_LH2=S;
pw_H1=S;
pw_H2=S;
pw_L1=S;
pw_L2=S;
P_LH1=S;
P_LH2=S;
P_H1=S;
P_H2=S;
P_L1=S;
P_L2=S;
   
    
    
    
    
for n=5:5:20
    % n: Number of samples
    % 1:uniform , 2:norm , 3:lap , 4: T5, 5:T6 , 6: T7 , 7:T8 , 8:T9 , 9:T10
    clear samp
    samp(1,:) = unifrnd(-1,1,1,n);
    samp(2,:) = normrnd(0,1,1,n);
    samp(3,:) = randl(1, n);
    samp(4,:) = random('T',5,1,n);
    samp(5,:) = random('T',6,1,n);
    samp(6,:) = random('T',7,1,n);
    samp(7,:) = random('T',8,1,n);
    samp(8,:) = random('T',9,1,n);
    samp(9,:) = random('T',10,1,n);
    
    for j=1:9


        %Calculate Kurtosis
        k(j) = kurtosis(samp(j,:));
        
        if n==5
            high = 2.8767;
            low = 1.2780;

        elseif n==10
            high = 3.9409;
            low = 1.5638;

        elseif n==15
            high = 4.1186;
            low = 1.7218;

        else
            high = 4.1493;
            low = 1.8308;

        end
        
        
        if k(j)>high || k(j)<low
            LH(j,n)=LH(j,n)+1;
            P_LH2(j,n) = signrank(samp(j,:));
            P_LH1(j,n) = signrank(samp(j,:),0,'tail','left');
            pw_LH1(j,n) = sampsizepwr('Z',[0 1],0.3,[],n,'tail','left');
            pw_LH2(j,n) = sampsizepwr('Z',[0 1],0.3,[],n);
            if k(j)>high
                H(j,n)=H(j,n)+1;
                P_H2(j,n) = signrank(samp(j,:));
                P_H1(j,n) = signrank(samp(j,:),0,'tail','left');
                pw_H1(j,n) = sampsizepwr('Z',[0 1],0.3,[],n,'tail','left');
                pw_H2(j,n) = sampsizepwr('Z',[0 1],0.3,[],n);
            else
                P_H2(j,n) = ttest(samp(j,:));
                P_H1(j,n) = ttest(samp(j,:),0,'tail','left');
                pw_H1(j,n) = sampsizepwr('t',[0 1],0.3,[],n,'tail','left');
                pw_H2(j,n) = sampsizepwr('t',[0 1],0.3,[],n);
            end
            if k(j)<low
                L(j,n)=L(j,n)+1;
                P_L2(j,n) = signrank(samp(j,:));
                P_L1(j,n) = signrank(samp(j,:),0,'tail','left');
                pw_L1(j,n) = sampsizepwr('Z',[0 1],0.3,[],n,'tail','left');
                pw_L2(j,n) = sampsizepwr('Z',[0 1],0.3,[],n);
            else
                P_L2(j,n) = ttest(samp(j,:));
                P_L1(j,n) = ttest(samp(j,:),0,'tail','left');
                pw_L1(j,n) = sampsizepwr('t',[0 1],0.3,[],n,'tail','left');
                pw_L2(j,n) = sampsizepwr('t',[0 1],0.3,[],n);
            end
        else
            P_LH2(j,n) = ttest(samp(j,:));
            P_LH1(j,n) = ttest(samp(j,:),0,'tail','left');
            pw_LH1(j,n) = sampsizepwr('t',[0 1],0.3,[],n,'tail','left');
            pw_LH2(j,n) = sampsizepwr('t',[0 1],0.3,[],n);
        end
        
        
        
        
        P_wilcoxon2(j,n) = signrank(samp(j,:));
        [h P_ttest2(j,n)] = ttest(samp(j,:));
        P_wilcoxon1(j,n) = signrank(samp(j,:),0,'tail','left');
        [h P_ttest1(j,n)] = ttest(samp(j,:),0,'tail','left');
        
        
        
        
        pw_ttest1(j,n) = sampsizepwr('t',[0 1],0.3,[],n,'tail','left');
        pw_ttest2(j,n) = sampsizepwr('t',[0 1],0.3,[],n);
        pw_wilcoxon1(j,n) = sampsizepwr('Z',[0 1],0.3,[],n,'tail','left');
        pw_wilcoxon2(j,n) = sampsizepwr('Z',[0 1],0.3,[],n);
        
    end
end

opw_ttest1=opw_ttest1+pw_ttest1;
opw_ttest2=opw_ttest2+pw_ttest2;
opw_wilcoxon2=opw_wilcoxon2+pw_wilcoxon2;
opw_wilcoxon1=opw_wilcoxon1+pw_wilcoxon1;
oP_ttest1=oP_ttest1+P_ttest1;
oP_ttest2=oP_ttest2+P_ttest2;
oP_wilcoxon1=oP_wilcoxon1+P_wilcoxon1;
oP_wilcoxon2=oP_wilcoxon2+P_wilcoxon2;
opw_LH1=opw_LH1+pw_LH1;
opw_LH2=opw_LH2+pw_LH2;
opw_H1=opw_H1+pw_H1;
opw_H2=opw_H2+pw_H2;
opw_L1=opw_L1+pw_L1;
opw_L2=opw_L2+pw_L2;
oP_LH1=oP_LH1+P_LH1;
oP_LH2=oP_LH2+P_LH2;
oP_H1=oP_H1+P_H1;
oP_H2=oP_H2+P_H2;
oP_L1=oP_L1+P_L1;
oP_L2=oP_L2+P_L2;

end
%%
opw_ttest1=opw_ttest1/1000;
opw_ttest2=opw_ttest2/1000;
opw_wilcoxon2=opw_wilcoxon2/1000;
opw_wilcoxon1=opw_wilcoxon2/1000;
oP_ttest1=oP_ttest1/1000;
oP_ttest2=oP_ttest2/1000;
oP_wilcoxon1=oP_wilcoxon1/1000;
oP_wilcoxon2=oP_wilcoxon2/1000;
opw_LH1=opw_LH1./LH;
opw_LH2=opw_LH2./LH;
opw_H1=opw_H1./H;
opw_H2=opw_H2./H;
opw_L1=opw_L1./L;
opw_L2=opw_L2./L;
oP_LH1=oP_LH1./LH;
oP_LH2=oP_LH2./LH;
oP_H1=oP_H1./H;
oP_H2=oP_H2./H;
oP_L1=oP_L1./L;
oP_L2=oP_L2./L;

toc