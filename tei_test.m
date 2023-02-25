close all
% clear
%load tei (calculated on midterm)

high = 4.1186;
low = 1.7218;
for i=1:4
   k(i) = kurtosis(tei(1:15,i)); 
    signrank(tei(1:15,i))
end
p1=signrank(tei(1:15,1));
p2=signrank(tei(1:15,2));
[h t1]=ttest(tei(1:15,1))
[h t2]=ttest(tei(1:15,2))