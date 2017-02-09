
clear;
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultTextFontSize',16);
fsize=16

d=200;                   % GW source distance Mpc
zGW=0.0455
%d=75;                   % GW source distance Mpc
%zGW=0.017;
oneday=86400; % sec
dcm=d*3.08e24; % cm
% r SDSS
freq_zero_r=4.8471e14;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%off-axis

%R-band flux zeropint in erg/cm^2/s/Hz
flusso_zero_R=2941.*1.e-23;


% fields={'tempooff','flussooff'};
% catm2=readcatalog('sgrbEjets1e50n1thetajet0d2opticaltheta0d40.txt',fields,',');
raw = load('sgrbEjets1e50n1thetajet0d2opticaltheta0d40.txt');
catm2.tempooff = raw(:,1);
catm2.flussooff = raw(:,2);
magooff=-2.5*log10([catm2.flussooff]/(flusso_zero_R/1.0e-26))+5.*log10(d*1.0e6*3.08d18/1.0d28);
tempooff2 = catm2.tempooff*(1+zGW)/(1+0.56);

% fields={'temposhortoff','flussoshortoff'};
% catm3=readcatalog('sgrbEjets1e50n1e-3thetajet0d2opticaltheta0d40.txt',fields,',');
raw = load('sgrbEjets1e50n1e-3thetajet0d2opticaltheta0d40.txt');
catm3.temposhortoff = raw(:,1);
catm3.flussoshortoff = raw(:,2);
magooff_short=-2.5*log10([catm3.flussoshortoff]/(flusso_zero_R/1.0e-26))+5.*log10((d*1.0e6*3.08d18)/1.0d28);
tempooff3 = catm3.temposhortoff*(1+zGW)/(1+0.56);

%% afterglow models
time=logspace(2,8,50); % sec
time/oneday
time_stripe=[time(1)  time(1) time(end) time(end) time(1)];
shb_hi=SHBafterglow(d,time,0.0);
shb_lo=SHBafterglow(d,time,8.0);
upleft=shb_hi(1); bottomleft=shb_lo(1);
upright=shb_hi(end); bottomright=shb_lo(end);
shb_stripe=[bottomleft upleft upright bottomright bottomleft];


%% make plots

%figure(1)
%set(gcf, 'PaperSize',[9 6])
%set(gcf, 'PaperPosition', [0 0 9 6])
clf
semilogx(time/oneday,shb_hi,'b', 'linewidth', 5)
hold on
semilogx(time/oneday,shb_lo,'b', 'linewidth', 5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%;;;;;;;;;;;;;;;;;KILONOVA MODEL;;;;;;;;;;;;;;;;;;;;
%Now overplotting the kilonova model by Metzger et al. 2010 stored
% in the file named BB_R_1e-2 that he provided us

% fields={'tempokilo','lumkilo','empty'};
% kilo1=readcatalog('BBR1e-2_new.dat',fields,',');
raw = load('BBR1e-2_new.dat');
kilo1.tempokilo = raw(:,1);
kilo1.lumkilo = raw(:,2);
kilo1.empty = raw(:,3);

%To convert from luminosity into R-band magnitudes, I use the
%following values for the flux zeropint and central frequency

%R-band flux zeropint in erg/cm^2/s/Hz
flusso_zero_R=2941.*1.e-23;

flusso_zero_R_mJy=flusso_zero_R*1.e26;

%R-band central frequency in Hz
freq_zero=4.3e14;

%calculating corresponding observed magnitudes
magkilo=-2.5*log10([kilo1(:).lumkilo]/(4.*pi*(d*1.d6*3.08*1.d18)^2.0))+2.5*log10(flusso_zero_R*freq_zero);

%Now overplotting the kilonova model by Metzger et al. 2010 stored
% in the file Fe_R_1e-2.dat that he provided us

% fields={'tempokilo','lumkilo','empty'};
% kilo2=readcatalog('Fe_R_1e-2_new.dat',fields,',');
raw = load('Fe_R_1e-2_new.dat');
kilo2.tempokilo = raw(:,1);
kilo2.lumkilo = raw(:,2);
kilo2.empty = raw(:,3);

magkilosmall=-2.5*log10([kilo2(:).lumkilo]/(4.*pi*(d*1.d6*3.08*1.d18)^2.0))+2.5*log10(flusso_zero_R*freq_zero);

%Now overplotting the kilonova model by Piran et al. 2012 stored
% in the file BH14_NS14.txt that they provided us. Here I make the approximation that all the bolometric luminosity goes in R band

% fields={'tempokilo','lumkilo','empty'};
% kilo2p=readcatalog('BH10_NS14.txt',fields,',');
raw = load('BH10_NS14.txt');
kilo2p.tempokilo = raw(:,1);
kilo2p.lumkilo = raw(:,2);
kilo2p.empty = raw(:,3);

magkilo_bhns=-2.5*log10([kilo2p.lumkilo]/(4.*pi*(d*1.d6*3.08*1.d18)^2.0))+2.5*log10(flusso_zero_R*freq_zero);

%Now the kilonova model by Barnes & Kasen stored
% in the file lv_l.txt that they provided us. 

% fields={'tempokilo','pip','empty','magkilo'};
% kilobarnes=readcatalog('lv_l.dat',fields,' ');
raw = load('lv_l.dat');
kilobarnes.tempokilo = raw(:,1);
kilobarnes.magkilo = raw(:,2);

magkilobarnes=[kilobarnes(:).magkilo]+5.0*log10(d*1.d5);

% kilonova model by Kasen, Fernandez & Metzger 2015 (r) stored
% in the file Kasen15_t0.dat (day, nu L_nu / 1e40  erg/s) 
% (BNS --> BH)

raw = load('Kasen15_t0_R.dat');
tempokilobr = raw(:,1)*(1+zGW);
lumkilobr = raw(:,2)*1.e+40;
magkilobarnesr=-2.5*log10(lumkilobr/(4.*pi*(dcm)^2.0))+2.5*log10(flusso_zero_R*freq_zero_r);



% kilonova model by Kasen, Fernandez & Metzger 2015 (r--) stored
% in the file Kasen15_tinf.dat (day, nu L_nu / 1e40  erg/s) 
% (BNS --> Hypermassive NS)

raw = load('Kasen15_tinf_R.dat');
tempokilobr2 = raw(:,1)*(1+zGW);
lumkilobr2 = raw(:,2)*1.e+40;
magkilobarnesr2=-2.5*log10(lumkilobr2/(4.*pi*(dcm)^2.0))+2.5*log10(flusso_zero_R*freq_zero_r);


% kilonova model by Kasen, Fernandez & Metzger 2015 (r) stored
% in the file Kasen15_t100.dat (day, nu L_nu / 1e40  erg/s) 
% (BNS --> Hypermassive NS for 100 ms --> BH)

%raw = load('Kasen15_t100_R.dat');
%tempokilobr3 = raw(:,1)*(1+zgw);
%lumkilobr3 = raw(:,2)*1.e+40;
%magkilobarnesr3=-2.5*log10(lumkilobr3/(4.*pi*(dcm)^2.0))+2.5*log10(flusso_zero_R*freq_zero_r);

%;;;;;;;;;;;;;;;;;OFF-axis Afterglow MODEL;;;;;;;;;;;;;;;;;;;;

% Distance of GRB off-axis
%130603B
dd=5.907e27;
zdd=0.3564

%raw = csvread('sgrbEjets1e50n1e-3thetajet0d2opticaltheta0d40.txt',49);
%raw = csvread('/Users/giulia/local_sw/boxfit/boxfitoutput/Results/130603B/lc1.txt',82);
raw = load('130603B_4thj_rSDSS.txt');
tempooffr1 = raw(:,2)*(1+zGW)/(1+zdd);
flussooffr1 = raw(:,4);
%magoff1=-2.5*log10(flussooff1/flusso_zero_R_mJy)+5.*log10(dcm/1.0d28);
magoffr1=-2.5*log10(flussooffr1/flusso_zero_R_mJy)+5.*log10(dcm/dd);

%raw = csvread('/Users/giulia/local_sw/boxfit/boxfitoutput/Results/130603B/lc2.txt',82);
raw = load('130603B_2thj_r.txt');
tempooffr2 = raw(:,2)*(1+zGW)/(1+zdd);
flussooffr2 = raw(:,4);
%magoff2=-2.5*log10(flussooff2/flusso_zero_R_mJy)+5.*log10(dcm/1.0d28);
magoffr2=-2.5*log10(flussooffr2/flusso_zero_R_mJy)+5.*log10(dcm/dd);



%% make plots


%semilogx([kilo1(:).tempokilo],magkilo,'r-', 'linewidth', 8)
%semilogx([kilo2p(:).tempokilo],magkilo_bhns,'r-', 'linewidth', 8)
semilogx([kilobarnes(:).tempokilo],magkilobarnes,'r-', 'linewidth', 8)
%semilogx([tempooff2],magooff,'k--', 'linewidth', 10)
%semilogx([tempooff3],[magooff_short],'k--', 'linewidth', 3)
semilogx([tempokilobr2],[magkilobarnesr2],'r--', 'linewidth', 8);
semilogx([tempokilobr],[magkilobarnesr],'r--', 'linewidth', 5);
semilogx(tempooffr1,magoffr1,'g--', 'linewidth', 5);
semilogx(tempooffr2,magoffr2,'g--', 'linewidth', 8);


hold off
%set(gca,'xtick',10.^(-1:1:2))
set(gca,'YDir','reverse')
xlabel('Observed Time (T-T_0) [days]');
ylabel('R mag');
%title(sprintf('afterglow light curves (source distance d=%d Mpc)',d));


%axis([1000.*9/oneday 10000000*0.4/oneday 5 30]); %questo
%axis([96.*9/oneday 10000000*0.4/oneday 15 28]); %mich
axis([1000.*12/86400 100000000*0.1/86400 8 30]);
%axis([1000.*12/86400 100000000*0.1/86400 8 26]);
%print('-color','-dpng','OpticalbEWASS200.png',sprintf('-FHelvetica:%d',fsize));
%print("-color","-deps","Opticalbw.eps",sprintf("-FHelvetica:%d",fsize));
% saveas(gcf,'kilonova2_200MPc.pdf');
% saveas(gcf,'kilonova2_200MPc.png');
