% Figs4_6_7.m

% This matlab script is used to produce 
% Figures 4,6 and 7 in "Investigating Heterogeneity of Intracellular
% Calcium Dynamics in Anterior Pituitary Lactotrophs Using a Combined
% Modelling/Experimental Approach", by M. Tomaiuolo, R. Bertram,
% A. E. Gonzalez-Iglesias, and J. Tabak, in J. Neuroendocrinology,
% 22:1279-1289, 2010.

% the cell mode of Matlab is used
% on Matlab R2007a

%%
% Figure 4

clear all;

ton=0;
toff=39999;
fc=0.01;
fer=0.01;
nu=30;
alpha=4.5e-06;
pl=0.0002;
dpip=0.0004;
vhold=-60;
vm=-20; 
sm=12;
vca=25;

allpar=[];
allres=[];

t=0:10:100000;
rest=t;

gca=0;
pip=10;
kc=0.15;
ks=0.3;

for i=1:21,
    
    cer_init = 39 + 9.1*(i-1);
    
    ip3=pip*dpip;
    minf = 1/(1+exp((vm-vhold)/sm));
    ica = gca*minf*(vhold-vca);
    jin=-alpha*ica;

    params=[cer_init,ip3,kc,ks];
    lambda3=fer*nu*(ks+ip3+pl)*jin./(kc+ks+ip3+pl);
    lambda2=-fer*nu*kc.*(pl+ip3)./(kc+ks+ip3+pl);
    lambda1=lambda3./lambda2+cer_init;
    
    peak=(jin+(pl+ip3).*cer_init)./(kc+ks+ip3+pl);
    decay=abs(lambda2);
    ka=1./(kc+ks+ip3+pl);
    area = cer_init/(fer*nu*kc);
   
    results=[peak,decay,area];
    allpar=[allpar,params];
    allres=[allres,results];
end

allpar=reshape(allpar,4,i);
allpar=allpar';
allres=reshape(allres,3,i);
allres=allres';

cer_s = [allpar(:,1) allres];
for j=1:size(cer_s,2),
    cer_sen(:,j)=(cer_s(:,j)-cer_s(11,j))/cer_s(11,j);
end

allpar=[];
allres=[];

t=0:10:100000;
rest=t;

gca=0;
kc=0.15;
ks=0.3;
cer_init = 132;

for i=1:21,
    
    pip = 3 + 0.7*(i-1);
    
    ip3=pip*dpip;
    minf = 1/(1+exp((vm-vhold)/sm));
    ica = gca*minf*(vhold-vca);
    jin=-alpha*ica;

    params=[cer_init,ip3,kc,ks];
    lambda3=fer*nu*(ks+ip3+pl)*jin./(kc+ks+ip3+pl);
    lambda2=-fer*nu*kc.*(pl+ip3)./(kc+ks+ip3+pl);
    lambda1=lambda3./lambda2+cer_init;
    
    peak=(jin+(pl+ip3).*cer_init)./(kc+ks+ip3+pl);
    decay=abs(lambda2);
    ka=1./(kc+ks+ip3+pl);
    area = cer_init/(fer*nu*kc);
   
    results=[peak,decay,area];
    allpar=[allpar,params];
    allres=[allres,results];
end

allpar=reshape(allpar,4,i);
allpar=allpar';
allres=reshape(allres,3,i);
allres=allres';

ip3_s = [allpar(:,2) allres];
for j=1:size(ip3_s,2),
    ip3_sen(:,j)=(ip3_s(:,j)-ip3_s(11,j))/ip3_s(11,j);
end

allpar=[];
allres=[];

t=0:10:100000;
rest=t;

gca=0;
ks=0.3;
cer_init = 132;
pip = 10;

for i=1:21,
    
    kc = 0.045 + 0.0105*(i-1);
    
    ip3=pip*dpip;
    minf = 1/(1+exp((vm-vhold)/sm));
    ica = gca*minf*(vhold-vca);
    jin=-alpha*ica;

    params=[cer_init,ip3,kc,ks];
    lambda3=fer*nu*(ks+ip3+pl)*jin./(kc+ks+ip3+pl);
    lambda2=-fer*nu*kc.*(pl+ip3)./(kc+ks+ip3+pl);
    lambda1=lambda3./lambda2+cer_init;
    
    peak=(jin+(pl+ip3).*cer_init)./(kc+ks+ip3+pl);
    decay=abs(lambda2);
    ka=1./(kc+ks+ip3+pl);
    area = cer_init/(fer*nu*kc);
   
    results=[peak,decay,area];
    allpar=[allpar,params];
    allres=[allres,results];
end

allpar=reshape(allpar,4,i);
allpar=allpar';
allres=reshape(allres,3,i);
allres=allres';

kc_s = [allpar(:,3) allres];
for j=1:size(kc_s,2),
    kc_sen(:,j)=(kc_s(:,j)-kc_s(11,j))/kc_s(11,j);
end

allpar=[];
allres=[];

t=0:10:100000;
rest=t;

gca=0;
kc=0.15;
cer_init = 132;
pip = 10;

for i=1:21,
    
    ks = 0.09 + 0.021*(i-1);
    
    ip3=pip*dpip;
    minf = 1/(1+exp((vm-vhold)/sm));
    ica = gca*minf*(vhold-vca);
    jin=-alpha*ica;

    params=[cer_init,ip3,kc,ks];
    lambda3=fer*nu*(ks+ip3+pl)*jin./(kc+ks+ip3+pl);
    lambda2=-fer*nu*kc.*(pl+ip3)./(kc+ks+ip3+pl);
    lambda1=lambda3./lambda2+cer_init;
    
    peak=(jin+(pl+ip3).*cer_init)./(kc+ks+ip3+pl);
    decay=abs(lambda2);
    ka=1./(kc+ks+ip3+pl);
    area = cer_init/(fer*nu*kc);
   
    results=[peak,decay,area];
    allpar=[allpar,params];
    allres=[allres,results];
end

allpar=reshape(allpar,4,i);
allpar=allpar';
allres=reshape(allres,3,i);
allres=allres';

ks_s = [allpar(:,4) allres];
for j=1:size(ks_s,2),
    ks_sen(:,j)=(ks_s(:,j)-ks_s(11,j))/ks_s(11,j);
end

subplot(2,2,1);
title('A) Cer(0)','FontName','Arial','FontSize',16);
hold all;
plot(cer_sen(:,1),cer_sen(:,2),'*k');
plot(cer_sen(:,1),cer_sen(:,3),'ok');
plot(cer_sen(:,1),cer_sen(:,4),'dk');
axis([-0.8 0.8 -1 2.5]);
xlabel('relative change parameter (\alpha)','FontName','Arial','FontSize',14);
ylabel('relative change output','FontName','Arial','FontSize',14);
subplot(2,2,2);
title('B) p_{ip3}','FontName','Arial','FontSize',16);
hold all;
plot(ip3_sen(:,1),ip3_sen(:,2),'*k');
plot(ip3_sen(:,1),ip3_sen(:,3),'ok');
plot(ip3_sen(:,1),ip3_sen(:,4),'dk');
axis([-0.8 0.8 -1 2.5]);
xlabel('relative change parameter (\alpha)','FontName','Arial','FontSize',14);
ylabel('relative change output','FontName','Arial','FontSize',14);
subplot(2,2,3);
title('C) k_{serca}','FontName','Arial','FontSize',16);
hold all;
plot(ks_sen(:,1),ks_sen(:,2),'*k');
plot(ks_sen(:,1),ks_sen(:,3),'ok');
plot(ks_sen(:,1),ks_sen(:,4),'dk');
axis([-0.8 0.8 -1 2.5]);
xlabel('relative change parameter (\alpha)','FontName','Arial','FontSize',14);
ylabel('relative change output','FontName','Arial','FontSize',14);
subplot(2,2,4);
hold all;
title('D) k_{pmca}','FontName','Arial','FontSize',16);
plot(kc_sen(:,1),kc_sen(:,2),'*k');
plot(kc_sen(:,1),kc_sen(:,3),'ok');
plot(kc_sen(:,1),kc_sen(:,4),'dk');
axis([-0.8 0.8 -1 2.5]);
xlabel('relative change parameter (\alpha)','FontName','Arial','FontSize',14);
ylabel('relative change output','FontName','Arial','FontSize',14);

%%
% figure 6

rcer_s=[];
rip3=[];
rkc=[];
rks=[];

for i=1:size(ks_s,2),
    if(std(cer_s(:,i)) > 0.00001)
        rcer_s(:,i)=(cer_s(:,i)-mean(cer_s(:,i)))/std(cer_s(:,i));
    else
        rcer_s(:,i)=(cer_s(:,i)-mean(cer_s(:,i)));
    end
    
    if(std(ip3_s(:,i)) > 0.00001)
        rip3(:,i)=(ip3_s(:,i)-mean(ip3_s(:,i)))/std(ip3_s(:,i));
    else
        rip3(:,i)=(ip3_s(:,i)-mean(ip3_s(:,i)));
    end
    
    if(std(kc_s(:,i)) > 0.00001)
        rkc(:,i)=(kc_s(:,i)-mean(kc_s(:,i)))/std(kc_s(:,i));
    else
        rkc(:,i)=(kc_s(:,i)-mean(kc_s(:,i)));
    end
    
    if(std(ks_s(:,i)) > 0.00001)
        rks(:,i)=(ks_s(:,i)-mean(ks_s(:,i)))/std(ks_s(:,i));
    else
        rks(:,i)=(ks_s(:,i)-mean(ks_s(:,i)));
    end
end

subplot(3,1,1);
plot(rkc(:,3),rkc(:,2),'sk',...
    rks(:,3),rks(:,2),'ok',...
    rip3(:,3),rip3(:,2),'dk',...
    rcer_s(:,3),rcer_s(:,2),'*k');
axis([-2 2 -2 2]);
title('A','FontName','Arial','FontSize',16); 
xlabel('decay','FontName','Arial','FontSize',14);
ylabel('peak','FontName','Arial','FontSize',14);
subplot(3,1,2);
plot(rkc(:,4),rkc(:,2),'sk',...
    rks(:,4),rks(:,2),'ok',...
    rip3(:,4),rip3(:,2),'dk',...
    rcer_s(:,4),rcer_s(:,2),'*k');
axis([-2 2 -2 2]);
title('B','FontName','Arial','FontSize',16);
%title('area vs peak','FontName','Arial','FontSize',16); 
xlabel('area','FontName','Arial','FontSize',14);
ylabel('peak','FontName','Arial','FontSize',14);
subplot(3,1,3);
plot(rkc(:,4),rkc(:,3),'sk',...
    rks(:,4),rks(:,3),'ok',...
    rip3(:,4),rip3(:,3),'dk',...
    rcer_s(:,4),rcer_s(:,3),'*k');
axis([-2 2 -2 2]);
title('C','FontName','Arial','FontSize',16);
%title('area vs decay','FontName','Arial','FontSize',16); 
xlabel('area','FontName','Arial','FontSize',14);
ylabel('decay','FontName','Arial','FontSize',14);

%%
% figure 7

clear all;

ton=0;
toff=39999;
fc=0.01;
fer=0.01;
nu=30;
alpha=4.5e-06;
pl=0.0002;
dpip=0.0004;
vhold=-60;
vm=-20; 
sm=12;
vca=25;
allpar=[];
allres=[];

t=0:10:100000;
rest=t;

gca=0;
pip=10;
kc=0.15;
ks=0.3;

rand('state',1);
for i=1:200,

    cer_init = 66 + 132*rand;
    pip = 5 + 10*rand;
    kc = 0.075 + 0.15*rand;
    ks = 0.15 + 0.3*rand;
    
    ip3=pip*dpip;
    minf = 1/(1+exp((vm-vhold)/sm));
    ica = gca*minf*(vhold-vca);
    jin=-alpha*ica;

    params=[cer_init,ip3,kc,ks];
    lambda3=fer*nu*(ks+ip3+pl)*jin./(kc+ks+ip3+pl);
    lambda2=-fer*nu*kc.*(pl+ip3)./(kc+ks+ip3+pl);
    lambda1=lambda3./lambda2+cer_init;
    
    peak=(jin+(pl+ip3).*cer_init)./(kc+ks+ip3+pl);
    decay=abs(lambda2);
    ka=1./(kc+ks+ip3+pl);
    area = cer_init/(fer*nu*kc);
   
    results=[peak,decay,area];
    allpar=[allpar,params];
    allres=[allres,results];
end

allpar=reshape(allpar,4,i);
allpar=allpar';
allres50=reshape(allres,3,i);
allres50=allres50';

[rpd50,ppd50]=corrcoef(allres50(:,1),allres50(:,2));
[rpa50,ppa50]=corrcoef(allres50(:,1),allres50(:,3));
[rda50,pda50]=corrcoef(allres50(:,2),allres50(:,3));

allpar=[];
allres=[];

for i=1:200,

    kc = 0.075 + 0.15*rand;
    cer_init = 99 + 66*rand;
    pip = 7.5 + 5*rand;
    ks = 0.225 + 0.15*rand;
    
    ip3=pip*dpip;
    minf = 1/(1+exp((vm-vhold)/sm));
    ica = gca*minf*(vhold-vca);
    jin=-alpha*ica;

    params=[cer_init,ip3,kc,ks];
    lambda3=fer*nu*(ks+ip3+pl)*jin./(kc+ks+ip3+pl);
    lambda2=-fer*nu*kc.*(pl+ip3)./(kc+ks+ip3+pl);
    lambda1=lambda3./lambda2+cer_init;
    
    peak=(jin+(pl+ip3).*cer_init)./(kc+ks+ip3+pl);
    decay=abs(lambda2);
    ka=1./(kc+ks+ip3+pl);
    area = cer_init/(fer*nu*kc);
   
    results=[peak,decay,area];
    allpar=[allpar,params];
    allres=[allres,results];
end

allpar=reshape(allpar,4,i);
allpar=allpar';
allres25=reshape(allres,3,i);
allres25=allres25';

[rpd25,ppd25]=corrcoef(allres25(:,1),allres25(:,2));
[rpa25,ppa25]=corrcoef(allres25(:,1),allres25(:,3));
[rda25,pda25]=corrcoef(allres25(:,2),allres25(:,3));

subplot(2,3,1);
plot(allres50(:,2),allres50(:,1),'ok','MarkerFaceColor','Black');
title('A','FontName','Arial','FontSize',16);
xlabel('decay','FontName','Arial','FontSize',14);
ylabel('peak','FontName','Arial','FontSize',14);
subplot(2,3,2);
plot(allres50(:,3),allres50(:,1),'ok','MarkerFaceColor','Black');
title('B','FontName','Arial','FontSize',16);
xlabel('area','FontName','Arial','FontSize',14);
ylabel('peak','FontName','Arial','FontSize',14);
subplot(2,3,3);
plot(allres50(:,3),allres50(:,2),'ok','MarkerFaceColor','Black');
title('C','FontName','Arial','FontSize',16);
xlabel('area','FontName','Arial','FontSize',14);
ylabel('decay','FontName','Arial','FontSize',14);
subplot(2,3,4);
plot(allres25(:,2),allres25(:,1),'ok','MarkerFaceColor','Black');
title('D','FontName','Arial','FontSize',16);
xlabel('decay','FontName','Arial','FontSize',14);
ylabel('peak','FontName','Arial','FontSize',14);
subplot(2,3,5);
plot(allres25(:,3),allres25(:,1),'ok','MarkerFaceColor','Black');
title('E','FontName','Arial','FontSize',16);
xlabel('area','FontName','Arial','FontSize',14);
ylabel('peak','FontName','Arial','FontSize',14);
subplot(2,3,6);
plot(allres25(:,3),allres25(:,2),'ok','MarkerFaceColor','Black');
title('F','FontName','Arial','FontSize',16);
xlabel('area','FontName','Arial','FontSize',14);
ylabel('decay','FontName','Arial','FontSize',14);

