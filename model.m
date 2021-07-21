function rf=model(x)
% close all
% clear all
% clc

%Design of heat exchanger to minimize fouling, minimize area/size, maximize
%U By :  Totok R. Biyanto, Ph.D

%variabel optimisasi

% do=0.0254;
% nb=8;
% ds=0.7;

ds=x(1);
do=x(2);
nb=x(3);

%tube side
lbi=0.6;
lbo=0.761;
di=do-(2*0.00277);
ltp=1.25*do;
lbb=(12+0.005*ds)/1000;
dotl=ds-lbb;
dctl=dotl-do;
nt=(0.3008*(dctl^2))/(ltp^2);
at=((22/7)*(di^2)*nt)/4;
gt=((35.3675*1.38/4)*4)/at;
ret=(gt*di)/0.000516;
prt=0.57*0.000516/0.00001975;
ht=(0.023*(ret^0.8)*(prt^0.4)*(0.00001975/di)*1.008700185)*3600*1.163;

%shell side
lta=4.5;
lbc=(lta/(nb+1));
lti=((nb-1)*lbc)+lbi+lbo;
tetads=2.158476633;
tetactl=2*(acos((ds/dctl)*(1-(2*26.4/100))));
sm=lbc*((lbb+(dctl/ltp)*(ltp-do)));
fw=(tetactl/(2*(22/7)))-((sin(tetactl))/((2*(22/7))));
fc=1-(2*fw);
ntcc=(ds/ltp)*(1-(2*26.4/100));
sb=lbc*(ds-dotl+0);
fsbp=(sb/sm);
lsb=(3.1+(0.004*ds))/1000;
ssb=(22/7)*ds*(lsb/2)*(((2*(22/7))-(2*tetads))/(2*(22/7)));
stb=((22/7)/4)*(((do+0.00079)^2)-(do^2))*nt*(1-fw);
jc=0.55+(0.72*fc);
rs=ssb/(ssb+stb);
rlm=(ssb+stb)/sm;
jl=(0.44*(1-rs))+((1-(0.44*(1-rs)))*exp(-2.2*rlm));
gs=(9.4/2)/sm;
res=(gs*do)/0.00001855;
jb=exp((-1.35*fsbp*(1-(2*rs))));
jr=1;
libintang=lbi/lbc;
lobintang=lbo/lbc;
js=((nb-1)+(libintang^(1-0.6))+(lobintang^(1-0.6)))/((nb-1)+(libintang-1)+(lobintang-1));
prs=(0.00001855*0.655)/0.0000107361111111;
ji=0.236*(res^(-0.346));
hi=(ji*0.655*gs*0.97/(prs^(2/3)))*(3600*1.163);
hs=hi*jc*jl*jb*js*jr;

%pressure drop
fs=exp(0.576-(0.19*log(res)));
ps=((fs*(gs^2)*ds*(nb+1))/(9.0119*do*0.8))*0.00001019716213*14.2233;
ft=((1.58*log(ret))-3.28)^(-2);
pt=(((2*ft*5*4/di)+(2*4))*872.6500*(1.3^2))*0.00001019716213*14.2233;

%fouling
alfa=277.8;
ea=48;
gamma=(4.17*(10^(-13)));
r=0.008314;

drft=(alfa*(ret^(-0.8))*(prt^(-1/3))*(5.53685*(10^(-5)))-(gamma*(ret^0.8)));
rftu=(drft*300);

drfs=(alfa*(res^(-0.8))*(prs^(-1/3))*(7.51479*(10^(-6)))-(gamma*(res^0.8)));
rfsh=(drfs*300);

rf=(rftu+rfsh);

cond=(do*(log(do/di)))/(2*20.8);
uf=1/((do/(di*ht))+((do*rftu)/di)+cond+rfsh+(1/hs));

%heat duty
ao=((22/7)*do*lti*nt*2);
lmtdcorr=59.67995442;
q=(uf*lmtdcorr*ao)/1000000;
























