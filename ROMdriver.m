%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% 7/30/2017                                    %
% Author: Muhammad Mohebujjaman                %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




 
%endTimestep = 166;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all
clc

%load EFRROMtestSV35kN10_1P_vis_801  

%load EFRROMtestSV35kN10_1P_vis

load EFRROMtestSV35kN10_11P_vis
load points5d.txt;
load weights5d.txt;

nn = 1;
filter = 'version-2';

sample_points = points5d;

m=801;
c = 1.0;
L = 1e-2;
dt = 0.002;

%delta = 0.0069;
delta = 0.0118;
chi = 0.002;

endTime = 10; % over a single period

N=4;

    MassROM = MassROM(1:N,1:N);
    StiffROM = StiffROM(1:N,1:N);
    StiffROMsin = StiffROMsin(1:N,1:N);
    StiffROMsin2 = StiffROMsin2(1:N,1:N);
    StiffROMcos = StiffROMcos(1:N,1:N);
    StiffROMcos2 = StiffROMcos2(1:N,1:N);
    TriLinROM2 = TriLinROM2(1:N,1:N,1:N);

    NLlift = NLlift(1:N,1:N);
    NLdrag = NLdrag(1:N,1:N);

    vdmass = vdmass(1:N);
    vlmass = vlmass(1:N);

    vdstiff = vdstiff(1:N);
    vdstiffsin = vdstiffsin(1:N);
    vdstiffsin2 = vdstiffsin2(1:N);
    vdstiffcos = vdstiffcos(1:N);
    vdstiffcos2 = vdstiffcos2(1:N);


    vlstiff = vlstiff(1:N);
    vlstiffsin = vlstiffsin(1:N);
    vlstiffsin2 = vlstiffsin2(1:N);
    vlstiffcos = vlstiffcos(1:N);
    vlstiffcos2 = vlstiffcos2(1:N);


for i=1:2
    sqrt_xi(i)=sqrt(sqrt(pi)*L)*exp(-(i*pi*L)^2/8.0);
end

for IC=1:m
    IC
    ran_var = sample_points(IC,:);

    nu1 = (c+sqrt(sqrt(pi)*L/2)*ran_var(1))/1000.0;

    nu2 = sqrt_xi(1)*ran_var(2)/1000.0;
    nu3 = sqrt_xi(1)*ran_var(3)/1000.0;
    nu4 = sqrt_xi(2)*ran_var(4)/1000.0;
    nu5 = sqrt_xi(2)*ran_var(5)/1000.0;

    EFRfilename = sprintf('vissnapT_7_%d.txt',IC);
    vvv = load(EFRfilename);
   
    
    % L2 project the initial condition (held in GlobalV) into ROM basis - coeff
    % vector is put into "velInit"

    RHS = zeros(N,1);
    for i=1:N
        RHS(i) = vvv' * (MassMatrix * PhiR(:,i) );
    end
    A = MassROM;
    RHS = RHS - A(:,1)*1;
    A(1,:)=0;
    A(:,1)=0;
    A(1,1)=1;
    RHS(1)=1;
    velInit = A  \ RHS;

    velPrevPrev=velInit;
    velPrevPrev_Init = velInit;

    

    T=endTime;
    numTimeSteps =  round(T/dt);
 
    mm = nn;
    EFR_ROM_vis_time_evolution;
    
    mm = nn + 1;
    velPrevPrev = velPrevPrev_Init;
    EFR_ROM_vis_time_evolution;
    
    mm = nn + 2;
    velPrevPrev = velPrevPrev_Init;
    EFR_ROM_vis_time_evolution;
    
    mm = nn + 3;
    velPrevPrev = velPrevPrev_Init;
    EFR_ROM_vis_time_evolution;
    
    velPrevPrev = velPrevPrev_Init;
    G_ROM_vis_time_evolution;
  

end

weights = weights5d;
       
numfiles = m;
EFRdata = cell(1,numfiles);
DNSdata = cell(1,numfiles);

for k = 1:m
    
    HODF = sprintf('visHODF1_%d.txt',k);
    EFRHODF1{k} = importdata(HODF);
    
    HODF = sprintf('visHODF2_%d.txt',k);
    EFRHODF2{k} = importdata(HODF);
    
    HODF = sprintf('visHODF3_%d.txt',k);
    EFRHODF3{k} = importdata(HODF);
    
    HODF = sprintf('visHODF4_%d.txt',k);
    EFRHODF4{k} = importdata(HODF);
    
    GROMfilename = sprintf('GROMQIvis%d.txt',k);
    GROMdata{k} = importdata(GROMfilename);
    
    DNSfilename = sprintf('visQI_%d.txt',k);
    DNSdata{k} = importdata(DNSfilename);
end

HODF1 = 0.0;
HODF2 = 0.0;
HODF3 = 0.0;
HODF4 = 0.0;

DNS_data = 0.0;
GROM_data = 0.0;

fact=1/(2*sqrt(3))^5;

for k=1:m
     HODF1 = HODF1 + fact* weights(k)*EFRHODF1{k}(:,2:end);
     HODF2 = HODF2 + fact* weights(k)*EFRHODF2{k}(:,2:end);
     HODF3 = HODF3 + fact* weights(k)*EFRHODF3{k}(:,2:end);
     HODF4 = HODF4 + fact* weights(k)*EFRHODF4{k}(:,2:end);
     
     GROM_data = GROM_data + fact* weights(k)*GROMdata{k}(:,2:end);
     DNS_data = DNS_data + fact* weights(k)*DNSdata{k}(:,2:end);
end

% figure;
% plot(EFRdata{1}(1:end,1),EFR_data(1:end,1),'r',EFRdata{1}(1:end,1),EFRHODF_data(1:end,1),'b-.',GROMdata{1}(:,1),GROM_data(:,1),'k',DNSdata{1}(3500:end,1)-7,DNS_data(3500:end,4),'b','LineWidth',2);
% xlabel('t','FontSize',20)
% ylabel('Lift','FontSize',20)
% title(['r=' num2str(N)],'FontSize',20)
% K=legend('EFR-ROM','HODF-ROM','G-ROM','DNS')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

figure;
plot(EFRHODF1{1}(:,1),HODF1(:,3),'r',EFRHODF1{1}(:,1),HODF2(:,3),'k-.',EFRHODF1{1}(:,1),HODF3(:,3),'r-.',EFRHODF1{1}(:,1),HODF4(:,3),'g',EFRHODF1{1}(:,1),GROM_data(:,3),'k',DNSdata{1}(3500:end,1)-7,DNS_data(3500:end,1),'b','LineWidth',2);
xlabel('t','FontSize',20)
ylabel('Energy','FontSize',20)
title(['r=' num2str(N), '  ' filter '  '],'FontSize',20)
K=legend('EFR-ROM','HODF-ROM2','HODF-ROM3','HODF-ROM4','G-ROM','DNS')
set(K,'Interpreter','Latex');
set(gca,'FontSize',20)
axis tight

 
% figure;
% plot(EFRdata{1}(:,1),EFR_data(:,2),'r',EFRdata{1}(1:end,1),EFRHODF_data(1:end,2),'b-.',GROMdata{1}(:,1),GROM_data(:,2),'k',DNSdata{1}(3500:end,1)-7,DNS_data(3500:end,3),'b','LineWidth',2);
% xlabel('t','FontSize',20)
% ylabel('Drag','FontSize',20)
% title(['r=' num2str(N)],'FontSize',20)
% K=legend('EFR-ROM','HODF-ROM','G-ROM','DNS')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight
%  
% 
% figure;
% plot(EFRdata{1}(:,1),EFR_data(:,3),'r',EFRdata{1}(1:end,1),EFRHODF_data(1:end,3),'k-.',GROMdata{1}(:,1),GROM_data(:,3),'k',DNSdata{1}(3500:end,1)-7,DNS_data(3500:end,1),'b','LineWidth',2);
% xlabel('t','FontSize',20)
% ylabel('Energy','FontSize',20)
% title(['r=' num2str(N)],'FontSize',20)
% K=legend('EFR-ROM','HODF-ROM','G-ROM','DNS')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight


% figure;
% plot(EFRdata{1}(:,1),EFR_data(:,3),'r',EFRdata{1}(1:end,1),EFRHODF_data(1:end,2),'b-.',GROMdata{1}(:,1),GROM_data(:,3),'k',DNSdata{1}(3500:end,1)-7,DNS_data(3500:end,1),'b',0,0.5,'w',0,0.6,'w','LineWidth',2);
% xlabel('t','FontSize',20)
% ylabel('Energy','FontSize',20)
% title(['r=' num2str(N)],'FontSize',20)
% K=legend('EFR-ROM','HODF-ROM','G-ROM','DNS')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

% figure
% plot(dataTableEFR_vis(1:end,1),dataTableEFR_vis(1:end,4),'r','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Energy','FontSize',20)
% title(['r=' num2str(N)],'FontSize',20)
% K = legend('EFR-ROM')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight
% 
% figure
% plot(dataTableEFR_vis(1:end,1),dataTableEFR_vis(1:end,2),'r','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Lift','FontSize',20)
% title(['r=' num2str(N)],'FontSize',20)
% K = legend('EFR-ROM')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight
% 
% figure
% plot(dataTableEFR_vis(1:end,1),dataTableEFR_vis(1:end,3),'r','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Drag','FontSize',20)
% title(['r=' num2str(N)],'FontSize',20)
% K = legend('EFR-ROM')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

