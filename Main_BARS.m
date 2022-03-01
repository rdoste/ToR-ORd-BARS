clear
close all
load ('sf_cal');   %load the current multipliers from the populations of models
% Setting parameters
DataRep = DataReporter;

index=35;
sf=sf_cal(35,:);

% Simulation and output extraction
 TimebyBeats=1; %Specify simulated time by number of beats (1) or interval duration (0)
              
    Sinterval = [1000]; %BCL
    t= [150]; 
    if TimebyBeats==1
         nBeats(1)=350;
    end
    
   %Starting Values
       %
        X0_Signaling= [0.00685533455220118,0.0184630401160325,0.000731426797266862,0.00746268345094940,0.0191020788696145,0.00115141961261304,0.000607348898749425,0.000639038753581265,0.000419992815346110,0.344257659177271,9.62262190945345,0.474028267773051,0.0148474493496437,0.203015985725400,0.00944459882156118,1.76916170568022e-10,8.36801924219387e-10,5.01719476559362e-11,0.0898875954193269,0.00272422687231002,0.225041998219388,0.0322381220102320,0.192803876209063,0.205457169881723,0.174050238959555,0.817148066343192,0.567236181979005,0.249911884364108,0.0646981828366729,0.0664977613514265,0.489057632872032,0.362107101574369,0.126950531297531,0.0233954992478800,0.0128401592747216,0.00629647926927854,4.29166115483698e-05,0.00917030613498568,0.0123536190564101,0.000664158274421826,0.000765842738691197,0.666165471397222,0.673477756497978,0.236980176272067,0.124710628782511,0.00404925913372347,0.0589106047787742,0.0274407333314977,6.32124571143896e-10,0.00159025206300466,0.00209267447556971,0.000502422412564797,0.0110248493074472,8.04005829146876e-11,0.000364313646402573,0.000705654325530757,0.000341340679127998];
        X0{1}=[getStartingState('Torord_endo_dynCl_BARS_1000') X0_Signaling ]; 
        X0{2}=[getStartingState('Torord_endo_dynCl_BARS_1000') X0_Signaling ]; 
        X0{3}=[getStartingState('Torord_endo_dynCl_HCM_BARS_1000') X0_Signaling ]; 
        X0{4}=[getStartingState('Torord_endo_dynCl_HCM_BARS_1000') X0_Signaling ]; 
 
 
  
%% BARS 

%ISO values
ISO=zeros(1,nBeats);
initISO=round(20); %start of ISO application
ISO(initISO:end)=0.1; %Iso values for each beat  (uM/l)
ISO(initISO:200)=linspace(0,max(ISO),200-initISO+1);%Iso values for each beat  (uM/l)  -->Lineal variation of the ISO values over time

%Apply PKA phosphorylation cascade from first non-negativa ISO value
NOISO= sum(cumprod(ISO==0,1),2); %detect first beats without BARS
runSignalingPathway=ones(1,nBeats);  
runSignalingPathway(1:NOISO)=0;  %0--> NO BARS   1->BARS

%%
%APD simulations for each interval
tic
  for j = 1:length(Sinterval)
     
    % param is the default model parametrization here
    param.bcl = Sinterval(j);
    param.model = @model_Torord_dynCl_BARS;
    param.ISO=ISO;
    param.runSignalingPathway=runSignalingPathway;
    param.verbose = 1; % It is possible to use verbose output even when parfor is used for simulations (the numbers from threads are mixed together, but they still given an idea of the progress).
    param.cellType=0;
    params(1:4) = param;
    params(1).ICaL_Multiplier=sf(1);  params(1).INa_Multiplier=sf(2); params(1).Ito_Multiplier=sf(3); params(1).INaL_Multiplier=sf(4);
                params(1).IKr_Multiplier=sf(5); params(1).IKs_Multiplier=sf(6); params(1).IK1_Multiplier=sf(7);  params(1).INaCa_Multiplier=sf(8); 
                params(1).INaK_Multiplier=sf(9); params(1).Jrel_Multiplier=sf(10); params(1).Jup_Multiplier=sf(11);

    params(2)=params(1);
    params(2).model=@model_Torord_dynCl_BARS;
    params(3)=params(1);
    params(3).model=@model_Torord_dynCl_BARS_HCM;
    params(4)=params(1);
    params(4).model=@model_Torord_dynCl_BARS_HCM;
 
    params(1).runSignalingPathway=zeros(size(runSignalingPathway)); % NO BARS !!!!!!!!!; 
    params(3).runSignalingPathway=zeros(size(runSignalingPathway)); % NO BARS !!!!!!!!!
 
    options = [];
          
     if  TimebyBeats==1
            Beats=nBeats(j);
     else 
            Beats = round(t(j)/(Sinterval(j)/1000));% number of beats
     end
    ignoreFirst = 0;
   
    %simulations
    parfor i = 1:length(params)       
        [time{j,i}, X{j,i}] = modelRunner_BARS(X0{i}, options, params(i),  Beats, ignoreFirst);
        currents{j,i} = getCurrentsStructure_BARS(time{j,i}, X{j,i}, params(i),  Beats-1);  %last beat
        currents_complete{j,i} = getCurrentsStructure_BARS(time{j,i}, X{j,i}, params(i),  0);%all beats
    end
     
    %change of initial conditions
    %         for i=1:length(params)          
    %                      %X0{i}=X{j,i}{end,1}(end,:);   %asign state after simulation
    %         end
  end        
toc  

   figure(1)
 cm2=colormap(lines(4));
 
    title('AP')
    hold on
    plot(time{1,1}{end,1}(1:end)-min(time{1,1}{end,1}),currents{1,1}.V,'color',cm2(1,:),'linewidth',2);
    plot(time{1,2}{end,1}(1:end)-min(time{1,2}{end,1}),currents{1,2}.V,'color',cm2(2,:),'linewidth',2);
    plot(time{1,3}{end,1}(1:end)-min(time{1,3}{end,1}),currents{1,3}.V,'color',cm2(3,:),'linewidth',2);
    plot(time{1,4}{end,1}(1:end)-min(time{1,4}{end,1}),currents{1,4}.V,'color',cm2(4,:),'linewidth',2);
    legend('Control','Control+BARS','HCM','HCM+BARS')
    ylabel('Membrane potential (mV)'), xlabel('Time  (ms)')

 figure(2)
    title('Ca2+')
    hold on
    plot(time{1,1}{end,1}(1:end)-min(time{1,1}{end,1}),currents{1,1}.Cai,'color',cm2(1,:),'linewidth',2);
    plot(time{1,2}{end,1}(1:end)-min(time{1,2}{end,1}),currents{1,2}.Cai,'color',cm2(2,:),'linewidth',2);
    plot(time{1,3}{end,1}(1:end)-min(time{1,3}{end,1}),currents{1,3}.Cai,'color',cm2(3,:),'linewidth',2);
    plot(time{1,4}{end,1}(1:end)-min(time{1,4}{end,1}),currents{1,4}.Cai,'color',cm2(4,:),'linewidth',2);
    legend('Control','Control+BARS','HCM','HCM+BARS')
    ylabel('Ca^{2+}  (mM)'), xlabel('Time  (ms)')
  
  figure(3) % PKA phosphorylation 
   for ind=1:2
       subplot(1,2,ind)
                 hold on
                  plot(currents_complete{1,ind}.time,currents_complete{1,ind}.fICaLP,'linewidth',2);
                  plot(currents_complete{1,ind}.time,currents_complete{1,ind}.fIKsP,'g','linewidth',2);
                  plot(currents_complete{1,ind}.time,currents_complete{1,ind}.fPLBP,'r','linewidth',2);
                  plot(currents_complete{1,ind}.time,currents_complete{1,ind}.fTnIP,'c','linewidth',2);
                  plot(currents_complete{1,ind}.time,currents_complete{1,ind}.fINaP,'y','linewidth',2);
                  plot(currents_complete{1,ind}.time,currents_complete{1,ind}.fINaKP,'m','linewidth',2);
                  plot(currents_complete{1,ind}.time,currents_complete{1,ind}.fIKurP,'Color',[204 204 204]/255,'linewidth',2);
                  plot(currents_complete{1,ind}.time,currents_complete{1,ind}.fRyRP,'k','linewidth',2);
   end
  
  title('Phosphorylation')
  legend('fICaL','fIKsP','fPLBP','fTnIP','fINaP','fINaKP','fRyRP','fIKurP')
  ylabel('f'), xlabel('time  ms')
  
  

 