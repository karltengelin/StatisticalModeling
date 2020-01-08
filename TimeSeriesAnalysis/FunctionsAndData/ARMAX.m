load proj18.mat

nvdi = ElGeneina.nvdi;
nvdi = (nvdi)/(255/2)-1; %gör om punkterna till att passa intervallet [-1 1]
nvdi = log(nvdi);       %Normalfördelar den
nvdi = nvdi-mean(nvdi); %Gör om till zero mean

modeldata = nvdi(1:floor(length(nvdi)*0.7)); % skapar vår modeldata och valdata
valdata = nvdi(length(nvdi)*0.7+1:end);
%% regndatan
rain_inp = zsave(793:end) - mean(zsave(793:end));
rain_inp = rain_inp(1:length(modeldata));
plot(rain_inp,'r--')
% raininp = sqrt(raininp);
% raininp = raininp - mean(raininp);
% 
% rain_m = raininp(1:length(raininp)*0.8);
% rain_v = raininp(length(raininp)*0.8+1:end);

myplotter(rain_inp,100)
%% Vi vill ta bort periodicitet på 36 - s for "season"
Arain36 = [1 zeros(1,35) -0.6]; %(-0.6 funkade bra annars)
rain_inps = filter(Arain36,1,rain_inp);
myplotter(rain_inps,100)
%% vi vill nu skapa "rainw" dvs prewhitea regnet, vi tar reda på hur regnet modeleras genom att kolla på myplotter(rain_inps)
% vi vill också hitta d,r och s som används i main script
Arainw = [1 0 0];
Crainw = [1 zeros(1,36)];

initmodelrain = idpoly(Arainw,[],Crainw);
initmodelrain.structure.C.Free = [1 zeros(1,35) 1];
%initmodelrain.structure.A.Free = [1 1 0 1];

ARrain_model = pem(rain_inp',initmodelrain);
present(ARrain_model)
%%
rainw = filter(ARrain_model.A,ARrain_model.C,rain_inp);
Arainw = ARrain_model.A;
Crainw = ARrain_model.C;
myplotter(rainw);
whitenessTest(rainw)
%%
nvdiw_modeldata = filter(conv(ARrain_model.A,Arain36),ARrain_model.C,modeldata);                    %filter(conv(ARrain_model.A,Arain36),ARrain_model.C,nvdi);

figure
M= 40; 
stem(-M:M,crosscorr(rainw,nvdiw_modeldata,M)); 
title('Cross correlation function'), xlabel('Lag')
hold on
plot(-M:M, 2/sqrt(length(rainw))*ones(1,2*M+1),'-') %confidensintervallen 2 genom roten ur N
plot(-M:M, -2/sqrt(length(rainw))*ones(1,2*M+1),'-') 
hold off
% detta ger d=3, r=2, s=1 OBS detta behöver inte vara helt korrekt


%%
% nvdiw_modeldata = filter(ARrain_model.A,ARrain_model.C,modeldata);                    %filter(conv(ARrain_model.A,Arain36),ARrain_model.C,nvdi);
% 
% figure
% M= 40; 
% stem(-M:M,crosscorr(rainw,nvdiw_modeldata,M)); 
% title('Cross correlation function'), xlabel('Lag')
% hold on
% plot(-M:M, 2/sqrt(length(rainw))*ones(1,2*M+1),'-') %confidensintervallen 2 genom roten ur N
% plot(-M:M, -2/sqrt(length(rainw))*ones(1,2*M+1),'-') 
% hold off
% detta ger d=3, r=2, s=1 OBS detta behöver inte vara helt korrekt