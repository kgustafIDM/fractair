% find the connectivityjapan (see Ferreira et al.) 
% ArrivalAirport == japanorigin; DepartureAirport == japandestination

[uniO numuO] = count_unique(japanorigin);
[uniD numuD] = count_unique(japandestination);

connectivityjapan = zeros(size(uniO,1),1);

for k=1:size(uniO,1)
  
    noonu = strcmp(japanorigin,uniO(k));
    fnoo = find(noonu==1);
    nar = count_unique(japandestination(fnoo));
    connectivityjapan(k) = size(nar,1);
end


[uniO numuO] = count_unique(chinaorigin);
[uniD numuD] = count_unique(chinadestination);

connectivitychina = zeros(size(uniO,1),1);

for k=1:size(uniO,1)
  
    noonu = strcmp(chinaorigin,uniO(k));
    fnoo = find(noonu==1);
    nar = count_unique(chinadestination(fnoo));
    connectivitychina(k) = size(nar,1);
end


[uniO numuO] = count_unique(indiaorigin);
[uniD numuD] = count_unique(indiadestination);

connectivityindia = zeros(size(uniO,1),1);

for k=1:size(uniO,1)
  
    noonu = strcmp(indiaorigin,uniO(k));
    fnoo = find(noonu==1);
    nar = count_unique(indiadestination(fnoo));
    connectivityindia(k) = size(nar,1);
end


[uniO numuO] = count_unique(usaorigin);
[uniD numuD] = count_unique(usadestination);

connectivityusa = zeros(size(uniO,1),1);

for k=1:size(uniO,1)
  
    noonu = strcmp(usaorigin,uniO(k));
    fnoo = find(noonu==1);
    nar = count_unique(usadestination(fnoo));
    connectivityusa(k) = size(nar,1);
end


[uniO numuO] = count_unique(germanyorigin);
[uniD numuD] = count_unique(germanydestination);

connectivitygermany = zeros(size(uniO,1),1);

for k=1:size(uniO,1)
  
    noonu = strcmp(germanyorigin,uniO(k));
    fnoo = find(noonu==1);
    nar = count_unique(germanydestination(fnoo));
    connectivitygermany(k) = size(nar,1);
end


[uniO numuO] = count_unique(franceorigin);
[uniD numuD] = count_unique(francedestination);

connectivityfrance = zeros(size(uniO,1),1);

for k=1:size(uniO,1)
  
    noonu = strcmp(franceorigin,uniO(k));
    fnoo = find(noonu==1);
    nar = count_unique(francedestination(fnoo));
    connectivityfrance(k) = size(nar,1);
end


% [uniO numuO] = count_unique(nigeriaorigin);
% [uniD numuD] = count_unique(nigeriadestination);
% 
% connectivitynigeria = zeros(size(uniO,1),1);
% 
% for k=1:size(uniO,1)
%   
%     noonu = strcmp(nigeriaorigin,uniO(k));
%     fnoo = find(noonu==1);
%     nar = count_unique(nigeriadestination(fnoo));
%     connectivitynigeria(k) = size(nar,1);
% end


[uniO numuO] = count_unique(aussieorigin);
[uniD numuD] = count_unique(aussiedestination);

connectivityaustralia = zeros(size(uniO,1),1);

for k=1:size(uniO,1)
  
    noonu = strcmp(aussieorigin,uniO(k));
    fnoo = find(noonu==1);
    nar = count_unique(aussiedestination(fnoo));
    connectivityaustralia(k) = size(nar,1);
end


[uniO numuO] = count_unique(brazilorigin);
[uniD numuD] = count_unique(brazildestination);

connectivitybrazil = zeros(size(uniO,1),1);

for k=1:size(uniO,1)
  
    noonu = strcmp(brazilorigin,uniO(k));
    fnoo = find(noonu==1);
    nar = count_unique(brazildestination(fnoo));
    connectivitybrazil(k) = size(nar,1);
end


[uniO numuO] = count_unique(canadaorigin);
[uniD numuD] = count_unique(canadadestination);

connectivitycanada = zeros(size(uniO,1),1);

for k=1:size(uniO,1)
  
    noonu = strcmp(canadaorigin,uniO(k));
    fnoo = find(noonu==1);
    nar = count_unique(canadadestination(fnoo));
    connectivitycanada(k) = size(nar,1);
end


[histcjapan,bincjapan] = hist(connectivityjapan,15);
figure; loglog(bincjapan,histcjapan,'r','DisplayName','japan');
xlabel('connectivity'); ylabel('histogram'); title('japan');

[histcchina,bincchina] = hist(connectivitychina,15);
figure; loglog(bincchina,histcchina,'r','DisplayName','china');
xlabel('connectivity'); ylabel('histogram'); title('china');

[histcusa,bincusa] = hist(connectivityusa,15);
figure; loglog(bincusa,histcusa,'r','DisplayName','usa');
xlabel('connectivity'); ylabel('histogram'); title('usa');

[histccanada,binccanada] = hist(connectivitycanada,15);
figure; loglog(binccanada,histccanada,'r','DisplayName','canada');
xlabel('connectivity'); ylabel('histogram'); title('canada');

[histcaustralia,bincaustralia] = hist(connectivityaustralia,15);
figure; loglog(bincaustralia,histcaustralia,'r','DisplayName','australia');
xlabel('connectivity'); ylabel('histogram'); title('australia');

[histcbrazil,bincbrazil] = hist(connectivitybrazil,15);
figure; loglog(bincbrazil,histcbrazil,'r','DisplayName','brazil');
xlabel('connectivity'); ylabel('histogram'); title('brazil');

[histcgermany,bincgermany] = hist(connectivitygermany,15);
figure; loglog(bincgermany,histcgermany,'r','DisplayName','germany');
xlabel('connectivity'); ylabel('histogram'); title('germany');

[histcfrance,bincfrance] = hist(connectivityfrance,15);
figure; loglog(bincfrance,histcfrance,'r','DisplayName','france');
xlabel('connectivity'); ylabel('histogram'); title('france');

[histcindia,bincindia] = hist(connectivityindia,15);
figure; loglog(bincindia,histcindia,'r','DisplayName','india');
xlabel('connectivity'); ylabel('histogram'); title('india');

mcnigeria = 0;
mcjapan = mean(connectivityjapan);
mcusa = mean(connectivityusa);
mcindia = mean(connectivityindia);
mccanada = mean(connectivitycanada);
mcfrance = mean(connectivityfrance);
mcgermany = mean(connectivitygermany);
mcbrazil = mean(connectivitybrazil);
mcchina = mean(connectivitychina);
mcaussie = mean(connectivityaussie);
meancntvty = [mcaussie mcbrazil mccanada mcchina mcfrance mcgermany mcindia mcnigeria mcjapan mcusa];

tcfrance = 0;
tcusa = 1.9;
tcjapan = 1.7;
tcchina = 2.1;
tccanada = 2.4;
tcaussie = 2.5;
tcindia = 2.53;
tcbrazil = 1.3;
tcgermany = 0.5;
tcnigeria = 0;
tailcntvty = [tcaussie tcbrazil tccanada tcchina tcfrance tcgermany tcindia tcnigeria tcjapan tcusa];
