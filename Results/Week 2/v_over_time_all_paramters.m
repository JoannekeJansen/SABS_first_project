close all
clear all

parameters ={
%% Fractional currents
'Fjunc'
'Fjunc_CaL'

%% Currents
'GNa'
'GNaB'
'IbarNaK'
'KmKo'
'KmNaip'
'Q10KmNai'
'Q10NaK'
'gkp'
'pNaK'
'epi'
'GClB'
'GClCa'
'KdClCa'
'Q10CaL'
'pCa'
'pK'
'pNak'
'IbarNCX'
'Kdact'
'KmCai'
'KmCao'
'KmNai'
'KmNao'
'Q10NCX'
'ksat'
'nu'
'IbarSLCaP'
'KmPCa'
'Q10SLCaP'
'GCaB'
'Kmf'
'Kmr'
'MaxSR'
'MinSR'
'Q10SRCaP'
'Vmax_SRCaP'
'ec50SR'
'hillSRCaP'
'kiCa'
'kim'
'koCa'
'kom'
'ks'

%% Ion concentrations
'Nao'
'Ko'
'Cao'
'Cli'
'Clo'
'Mgi'
%Ki?? state variable instead of parameter?

%% Buffering
%'Bmax_Naj'
%'Bmax_Nasl'
%'koff_na'
%'kon_na'
%'Bmax_CaM'
%'Bmax_SR'
%'Bmax_TnChigh'
%'Bmax_TnClow'
%'Bmax_myosin'
%'koff_cam'
%'koff_myoca'
%'koff_myomg'
%'koff_sr'
%'koff_tnchca'
%'koff_tnchmg'
%'koff_tncl'
%'kon_cam'
%'kon_myoca'
%'kon_myomg'
%'kon_sr'
%'kon_tnchca'
%'kon_tnchmg'
%'kon_tncl'
%'Bmax_SLhighj0'
%'Bmax_SLhighsl0'
%'Bmax_SLlowj0'
%'Bmax_SLlowsl0'
%'koff_slh'
%'koff_sll'
%'kon_slh'
%'kon_sll'
%'Bmax_Csqn0'
%'koff_csqn'
%'kon_csqn'

%% Environmental parameters:
%'cellLength'
%'cellRadius'
%'distJuncSL'
%'distSLcyto'
%'junctionLength'
%'junctionRadius'
%'DnaJuncSL'
%'DnaSLcyto'
%'DcaJuncSL'
%'DcaSLcyto'
%'J_ca_juncsl'
%'J_ca_slmyo'
%'J_na_juncsl'
%'J_na_slmyo'

%% Physical constants:
%'Cmem'
%'Frdy'
%'R'
%'Temp'

%% Stim constants:
%'stim_amplitude'
%'stim_duration'
%'stim_period'
%'stim_start'
};

C=load('values_GNA_23.txt');

for i = (1:1:length(parameters))
    a =char(strcat('values_', cellstr(parameters(i)), '90.txt'));
    b= char(strcat('values_', cellstr(parameters(i)), '110.txt'));
    A=load(a);
    B=load(b);
    figure
    subplot(1,2,1)
    hold on
    plot(C(:,1),A(:,3), 'Linewidth', 3)
    plot(C(:,1),C(:,3), 'Linewidth', 3)
    plot(C(:,1),B(:,3), 'Linewidth', 3)
    xlim([1 499])
    xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
    ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex') 
    set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
    h=legend('$90\%$','$100\%$','$110\%$');
    set(h,'FontSize',20,'Interpreter', 'Latex')
    title(parameters(i))
    subplot(1,2,2)
    hold on
    plot(C(:,1),A(:,2), 'Linewidth', 3)
    plot(C(:,1),C(:,2), 'Linewidth', 3)
    plot(C(:,1),B(:,2), 'Linewidth', 3)
    xlim([1 499])
    xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
    ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex') 
    set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
    h=legend('$90\%$','$100\%$','$110\%$');
    set(h,'FontSize',20,'Interpreter', 'Latex')
    title(parameters(i))
end

shg