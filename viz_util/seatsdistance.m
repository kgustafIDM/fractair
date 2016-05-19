efmenu install
addpath ~/Research/Naef/paperHiSeq/matlabscripts/

load India/data/indiadata.mat
% [histIndiaw binsIndiaw] = hist(indiamiles.*indiaseats,100);
% figure; loglog(binsIndiaw,histIndiaw,'k')
[histIndiaww binsIndiaww] = hist(indiamiles.*indiaseats.*indiadays./7,100);
figure; loglog(binsIndiaww,histIndiaww./sum(histIndiaww),'k'); title('India flight'); 
xlabel('seat-miles/week'); ylabel('frequency');
indiasmpw_exp = 1.55;
indiagdp_ppp = 4.96e12;
indiapop = 1.273e9; 
indiasurface = 3.29e6; % km^2
indiaroads = 3.32e6; % km
indiacars = 8;

load Australia/data/aussiedata.mat
% [histAustralw binsAustralw] = hist(aussiemiles.*aussieseats,100);
% figure; loglog(binsAustralw,histAustralw,'k')
[histAustralww binsAustralww] = hist(aussiemiles.*aussieseats.*aussiedays./7,100);
figure; loglog(binsAustralww,histAustralww./sum(histAustralww),'k'); title('Australia flight'); 
xlabel('seat-miles/week'); ylabel('frequency');
australiasmpw_exp = 2.2; % note some flattening of tail
australiagdp_ppp = 1e12;
aussiepop = 23.28e6; 
aussiesurface = 7.7e6; % km^2
aussieroads = 0.813e6;
aussiecars = 545.44;

load Germany/data/germanydata.mat
% [histGermanyw binsGermanyw] = hist(germanymiles.*germanyseats,100);
% figure; loglog(binsGermanyw,histGermanyw,'k')
[histGermanyww binsGermanyww] = hist(germanymiles.*germanyseats.*germanydays./7,100);
figure; loglog(binsGermanyww,histGermanyww./sum(histGermanyww),'k'); title('Germany flight'); 
xlabel('seat-miles/week'); ylabel('frequency');
germanysmpw_exp = 1.46; % very unlike france
germanygdp_ppp = 3.23e12;
germanypop = 80.52e6; 
germanysurface = 0.357e6; % km^2
germanyroads = 0.644e6;
germanycars = 556.1;

load China/data/chinadata.mat
% [histChinaw binsChinaw] = hist(chinamiles.*chinaseats,100);
% figure; loglog(binsChinaw,histChinaw,'k');
[histChinaww binsChinaww] = hist(chinamiles.*chinaseats.*chinadays./7,100);
figure; loglog(binsChinaww,histChinaww./sum(histChinaww),'k'); title('China flight'); 
xlabel('seat-miles/week'); ylabel('frequency');
chinasmpw_exp = 3.3; % steep and consistent
chinagdp_ppp = 13.4e12;
chinapop = 1.35e9; 
chinasurface = 9.7e6; % km^2
chinaroads = 3.86e6;
chinacars = 22.47; % per 1000 people

load USA/data/usadata.mat
% [histUSAw binsUSAw] = hist(usamiles.*usaseats,100);
% figure; loglog(binsUSAw,histUSAw,'k');
[histUSAww binsUSAww] = hist(usamiles.*usaseats.*usadays./7,100);
figure; loglog(binsUSAww,histUSAww./sum(histUSAww),'k'); title('USA flight'); 
xlabel('seat-miles/week'); ylabel('frequency');
usasmpw_exp = (1.73+5.53)/2; % clear break in slope
usagdp_ppp = 16.7e12;
usapop = 317.2e6; 
usasurface = 9.82e6 ; % km^2
usaroads = 6.5e6;
usacars = 450.67;

load Brazil/data/brazildata.mat
% [histBrazilw binsBrazilw] = hist(brazilmiles.*brazilseats,1000);
% figure; loglog(binsBrazilw,histBrazilw,'k')
[histBrazilww binsBrazilww] = hist(brazilmiles.*brazilseats.*brazildays./7,100);
figure; loglog(binsBrazilww,histBrazilww./sum(histBrazilww),'k'); title('Brazil flight'); 
xlabel('seat-miles/week'); ylabel('frequency');
brazilsmpw_exp = 1.1;
brazilgdp_ppp = 2.42e12;
brazilpop = 201e6; 
brazilsurface = 8.52e6 ; % km^2
brazilroads = 1.75e6;
brazilcars = 158.1;

load Canada/data/canadadata.mat
% [histCanadaw binsCanadaw] = hist(canadamiles.*canadaseats,1000);
% figure; loglog(binsCanadaw,histCanadaw,'k')
[histCanadaww binsCanadaww] = hist(canadamiles.*canadaseats.*canadadays./7,100);
figure; loglog(binsCanadaww,histCanadaww./sum(histCanadaww),'k'); title('Canada flight'); 
xlabel('seat-miles/week'); ylabel('frequency');
canadasmpw_exp = 1.48; % very flat curve
canadagdp_ppp = 1.52e12;
canadapop = 35.1e6; 
canadasurface = 9.98e6; % km^2
canadaroads = 1.04e6;
canadacars = 371.98;

load France/data/francedata.mat
% [histFrancew binsFrancew] = hist(francemiles.*franceseats,1000);
% figure; loglog(binsFrancew,histFrancew,'k')
[histFranceww binsFranceww] = hist(francemiles.*franceseats.*francedays./7,100);
figure; loglog(binsFranceww,histFranceww./sum(histFranceww),'k'); title('France flight'); 
xlabel('seat-miles/week'); ylabel('frequency');
francesmpw_exp = 2; % note the major flattening of long tail
francegdp_ppp = 2.27e12;
francepop = 63.46e6; 
francesurface = 0.64e6; % km^2
franceroads = 1.027e6; % source http://chartsbin.com/view/2770
francecars = 497.51; % source http://chartsbin.com/view/1113

load Japan/data/japandata.mat
% [histJapanw binsJapanw] = hist(japanmiles.*japanseats,1000);
% figure; loglog(binsJapanw,histJapanw,'k')
[histJapanww binsJapanww] = hist(japanmiles.*japanseats.*japandays./7,100);
figure; loglog(binsJapanww,histJapanww./sum(histJapanww),'k'); title('Japan flight'); 
xlabel('seat-miles/week'); ylabel('frequency');
japansmpw_exp = 1.68; % note the major flattening of long tail
japangdp_ppp = 4.7e12;
japanpop = 127e6; 
japansurface = 0.378e6; % km^2
japanroads = 1.20e6;
japancars = 324.56;

load Nigeria/data/nigeriadata.mat
% [histNigeriaw binsNigeriaw] = hist(nigeriamiles.*nigeriaseats,1000);
% figure; loglog(binsNigeriaw,histNigeriaw,'k')
[histNigeriaww binsNigeriaww] = hist(nigeriamiles.*nigeriaseats.*nigeriadays./7,100);
figure; loglog(binsNigeriaww,histNigeriaww./sum(histNigeriaww),'k'); title('Nigeria flight'); 
xlabel('seat-miles/week'); ylabel('frequency');
nigeriasmpw_exp = 1; % note the major flattening of long tail
nigeriagdp_ppp = 0.4785e12;
nigeriapop = 174.5e6; 
nigeriasurface = 0.923e6; % km^2
nigeriaroads = 0.193e6;
nigeriacars = 30.81;

smpw = [australiasmpw_exp brazilsmpw_exp canadasmpw_exp chinasmpw_exp francesmpw_exp germanysmpw_exp indiasmpw_exp nigeriasmpw_exp japansmpw_exp usasmpw_exp];
populations = [aussiepop brazilpop canadapop chinapop francepop germanypop indiapop japanpop nigeriapop usapop];
surfacearea = [aussiesurface brazilsurface canadasurface chinasurface francesurface germanysurface indiasurface japansurface nigeriasurface usasurface];
densitypeople = populations./surfacearea;
gdp_ppp = [australiagdp_ppp brazilgdp_ppp canadagdp_ppp chinagdp_ppp francegdp_ppp germanygdp_ppp indiagdp_ppp japangdp_ppp nigeriagdp_ppp usagdp_ppp];
roadskm = [aussieroads brazilroads canadaroads chinaroads franceroads germanyroads indiaroads japanroads nigeriaroads usaroads];
carsper1000 = [aussiecars brazilcars canadacars chinacars francecars germanycars indiacars japancars nigeriacars usacars];
roadspersurf = roadskm./surfacearea;