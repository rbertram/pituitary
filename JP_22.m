% JP_22.m
%
% This computer programs is for the pituitary corticotroph model used in the
% publication "Chronic Stress Facilitates Bursting Electrical Activity in 
% Pituitary Corticotrophs", Journal of Physiology, 600:313-332, 2022. 
% Authors are Peter Duncan, Mehran Fazli, Nicola Romano, Paul Le Tissier, 
% Richard Bertram, and Michael Shipston.

%% parameters
cm = 7.0; % membrane capacitance (pf)
vca = 60.0; % calcium Nernst potential (mV)
vk = -70.0; % potassium Nernst potential (mV)
vns = -20.0; % half-activation voltage of non-selective channels (mV)
vm = -20.0; % half-activation voltage of L-type calcium channels (mV)
vn = -5.0; % half-activation voltage of delayed rectifier potassium channels (mV)
vir = -50.0; % half-activation voltage of inward rectifier potassium channels (mV)
vl = -50.0; % half-activation voltage of leak currents (mV)
vs = -20.0; % half-activation voltage of BK-STREX channels (mV)
vz = -5.0; % half-activation voltage of BK-ZERO channels (mV)
sm = 12.0; % slope factor for activation function of L-type calcium channels
sn = 10.0; % slope factor for activation function of delayed rectifier potassium channels
sbs = 2; % slope factor for activation function of BK-STREX channels
sbz = 2; % slope factor for activation function of BK-ZERO channels
kik = 0.4; % activation coefficient of intermediate conductance potassium channel
sir = -1; % slope factor for activation function of 
taun = 30; % time scale of delayed rectifier potassium channel (ms)
taubf = 1000; % time scale of BK-far channels opening (ms)
taubn = 5; % time scale of BK-near channels opening (ms)
tauoc = 5; % time scale of all BK channels closing (ms)
kc = 0.12; % the Calcium pump rate
alpha = 0.0015; % coefficient to convert current to concentration
fc = 0.005; % the fraction of cytosolic Calcium that is free
gl = 0.2; % leak conductance (ns)
gkdr = 6.5; % delayed rectifier potassium channel maximum conductance (ns)
gir=0.93; % inward rectifier potassium channel maximum conductance (ns)
gik=0.5; % intermediate conductance potassium channel maximum conductance (ns)
aff=10; % binding affinity ratio of paxilline to close channels compare to open channel
gbk=0.3; % single BK channel conductance (ns)

p_=[1 2.1 0.12 0.2 0.2 5 20 5 10 2];
CRH        = p_(1); % 1 means in the presence of CRH and 0 means in basal condition
gca        = p_(2); % L-type calcium channel maximum conductance (ns)
gns        = p_(3); % Non Selective channel maximum conductance (ns)
betaz      = p_(4); % fraction of colocalised BK-ZERO channels with L-type Calcium channel
betas      = p_(5); % fraction of colocalised BK-STREX channels with L-type Calcium channel
Nz_min      = p_(6); % total number of BK-ZERO channel in basal condition
Nz_max      = p_(7); % total number of BK-ZERO channel in the presence of CRH
Ns_min      = p_(8); % total number of BK-STREX channel in the presence of CRH
Ns_max      = p_(9); % total number of BK-STREX channel in basal condition

%% BK channel status
Ns = (CRH*(Ns_min)+(1-CRH)*(Ns_max)); % total number of BK-STREX
Nz = (CRH*(Nz_max)+(1-CRH)*(Nz_min)); % total number of BK-ZERO
ofa        = p_(10); % 0 means all  of them are blocked and 1 means 3 (CRH) or 2 (basal) are unblocked and 2 means without paxilline 
Nzn=round(betaz*Nz); % total number of BK-ZERO-near
Nsn=round(betas*Ns); % total number of BK-STREX-near
Nzf=Nz-Nzn; % total number of BK-STREX-far
Nsf=Nz-Nsn; % total number of BK-ZERO-far
%% variables and time step
sec=20; % total time (sec)
dt=0.05; % solver time steps (ms)
tsec=0:dt*(1000^(-1)):sec;
msec=sec*1000;
N=msec*(1/dt);

v=zeros(1,N+1); % membrane voltage
Osn=zeros(1,N+1); % number of open BK-STREX-near channel
Ozn=zeros(1,N+1); % number of open BK-ZERO-near channel
n=zeros(1,N+1); % fraction of open delayed rectifier potassium channel
c=zeros(1,N+1); % calcium concentration
Osf=zeros(1,N+1); % number of open BK-STREX-far channel
Ozf=zeros(1,N+1); % number of open BK-ZERO-far channel

% initial condition
x0 = [-60 0 0 0.1 0.1 0 0];


v(1)          = x0(1);
Osn(1)        = x0(2);
Ozn(1)        = x0(3);
n(1)          = x0(4);
c(1)          = x0(5);
Osf(1)        = x0(6);
Ozf(1)        = x0(7);


if CRH==1 % in presence of CRH number o 
   nbc=3; % number of unblocked BK channel in the presence of CRH and paxilline
else
   nbc=2; % number of unblocked BK channel in the basal condition and presence of paxilline
end
%% Euler method ODE solver
for t=1:N
    
    % infinity functions
    zinf = (1+exp((vz-v(t))*sbz^(-1)))^(-1);
    sinf = (1+exp(-(v(t)-vs)*sbs^(-1)))^(-1);
    minf = (1+exp(-(v(t)-vm)*sm^(-1)))^(-1);
    ninf = (1+exp(-(v(t)-vn)*sn^(-1)))^(-1);
    kirinf = (1+exp(-(v(t)-vir)*sir^(-1)))^(-1);
    ikinf = (kik^2+c(t)^2)^(-1)*c(t)^2;
    
    % currents
    ica = (v(t)-vca)*minf*gca;
    ikdr = n(t)*(v(t)-vk)*gkdr;
    ibf = (v(t)-vk)*Osf(t)*gbk + (v(t)-vk)*Ozf(t)*gbk;
    ibn = (v(t)-vk)*Ozn(t)*gbk + (v(t)-vk)*Osn(t)*gbk;
    il = gl*(v(t)-vl);
    ins = (v(t)-vns)*gns;
    ikir = kirinf*gir*(v(t)-vk);
    iik = ikinf*(v(t)-vk)*gik;
    
    % Deterministic variables
    v(t+1) = v(t)+dt*(-cm^(-1)*(iik+il+ikdr+ins+ibf+ikir+ibn+ica));
    n(t+1) = n(t)+dt*(-(n(t)-ninf)*taun^(-1));
    c(t+1) = c(t)+dt*(-fc*(alpha*ica+c(t)*kc));
    
    % Stochastic variables
    if ofa~=0 
        Ozn(t+1) = Ozn(t)-binornd(Ozn(t),dt*(1-zinf)/tauoc)+binornd(Nzn-Ozn(t),dt*zinf/taubn); % computing number of open BK-ZERO-near channel
        Ozf(t+1) = Ozf(t)-binornd(Ozf(t),dt*(1-zinf)/tauoc)+binornd(Nzf-Ozf(t),dt*zinf/taubf); % computing number of open BK-ZERO-far channel
        Osn(t+1) = Osn(t)-binornd(Osn(t),dt*(1-sinf)/tauoc)+binornd(Nsn-Osn(t),dt*sinf/taubn); % computing number of open BK-STREX-near channel
        Osf(t+1) = Osf(t)-binornd(Osf(t),dt*(1-sinf)/tauoc)+binornd(Nsf-Osf(t),dt*sinf/taubf); % computing number of open BK-STREX-near channel
        if(ofa==1) % Presence of paxilline 
            pax_zn=0;
            pax_zf=0;
            pax_sn=0;
            pax_sf=0;
            for k=1:nbc % algorithm for determining type of unblocked BK channel
                chantype=randi(Nz+Ns);
                zchanprox=randi(Nz);
                schanprox=randi(Ns);
                if(chantype<=Nz)  
                    if(zchanprox<=Nzn)
                        zochanprox=randi(aff*Ozn(t+1)+Nzn-Ozn(t+1));
                        if (zochanprox<=aff*Ozn(t+1))&&(pax_zn<Nzn)
                            pax_zn=1+pax_zn;
                        end
                    else
                        zochanprox=randi(aff*Ozf(t+1)+Nzf-Ozf(t+1));
                        if (zochanprox<=aff*Ozf(t+1))&&(pax_zf<Nzf)
                            pax_zf=1+pax_zf;
                        end
                    end
                else
                    if(schanprox<=Nsn)
                        sochanprox=randi(aff*Osn(t+1)+Nsn-Osn(t+1));
                        if (sochanprox<=aff*Osn(t+1))&&(pax_sn<Nsn)
                            pax_sn=1+pax_sn;
                        end
                    else
                        sochanprox=randi(aff*Osf(t+1)+Nsf-Osf(t+1));
                        if (sochanprox<=aff*Osf(t+1))&&(pax_sf<Nsf)
                            pax_sf=1+pax_sf;
                        end
                    end
                end
            end
            Ozn(t+1) = pax_zn; 
            Ozf(t+1) = pax_zf;
            Osn(t+1) = pax_sn; 
            Osf(t+1) = pax_sf;
        end
    end
end
