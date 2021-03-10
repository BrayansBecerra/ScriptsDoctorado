% Modelo sintetico para evaluaciones de las señales de STFT, WT y HHT
 fs=4;
 %%t=0:1/fs:300-(1/fs); %Para el modelo HRV
%t=0:1/fs:300;
 %  f1=0.05;
%  f2=0.3;
 
%  m=1+cos(2*pi*f1*t);
%  c=cos(2*pi*f2*t);
%  y=m.*c;
%  
%  
%  figure
%  subplot(3,1,1)
%  plot(t,m);
%  subplot(3,1,2)
%  plot(t,c,'r')
%  subplot(3,1,3)
%  plot(t,y,'black')
 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%  Modelo para HRV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% t1=0:1/fs:150-(1/fs);
%  t2=150:1/fs:300-(1/fs);
% % 
%  f1=0.05; %f1
%  f2=0.28; %f2
%  f3=0.09; %f3
%  f4=0.18; %f4
% % 
%  Y1=1.2228+0.15*cos(2*pi*f1*t1)+0.2*cos(2*pi*f2*t1);
%  Y2=0.921+0.05*cos(2*pi*f4*t2)+0.1*cos(2*pi*f3*t2);
% % 
% % 
%  yyLF1=0.15*cos(2*pi*f1*t1);
% % yyLF2=0.1*cos(2*pi*f3*t2);
% % yyHF1=0.2*cos(2*pi*f2*t1);
%  yyHF2=0.05*cos(2*pi*f4*t2);
% % 
% % norm(x,2)^2/numel(x)
% % 
% % l2normLF = (norm(yyLF1,2)^2/numel(yyLF1))*(1000)^2
% % l2normHF = (norm(yyHF2,2)^2/numel(yyHF2))*(1000)^2
% 
% PotenciaLF1=bandpower(yyLF1)*(1000)^2
% % PotenciaLF2=bandpower(yyLF2,fs,[0.04 0.15])*(1000)^2
% % PotenciaHF1=bandpower(yyHF1,fs,[0.15 0.4])*(1000)^2
% PotenciaHF2=bandpower(yyHF2)*(1000)^2
% % 
% % 
%  y=cat(1,Y1',Y2');
%  y=y';
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% plot(t1,Y1,'blue');
% hold on
% plot(t2,Y2,'red');
% 
% figure
% plot(t,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  Señal real  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fin=length(t);
load('custom_ascii_ecg_data_hrv.mat');
t=cell2mat(Res.HRV.Data.T_RRi)';
y=cell2mat(Res.HRV.Data.RRi)';

x=y-mean(y);
x=detrend(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Modelo de cosenoidales estaticos %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f1=0.1;
% f2=0.2;
% 
% Y1=0.3*cos(2*pi*f1*t);
% Y2=0.2*cos(2*pi*f2*t);
% y=1+(0.3*cos(2*pi*f1*t)+0.2*cos(2*pi*f2*t)); 
%  
% PotenciaLF=bandpower(Y1)*(1000)^2
% PotenciaHF=bandpower(Y2)*(1000)^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Modelo de amplitud lineal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f1=0.05;
% f2=0.25;
% 
% A=0.4;
% B=0.42;
% m=0.01;
% 
% %y=(cos(2*pi*t*f1)./(A*(1+(m*t))))+((B*(1+(m*t))).*cos(2*pi*t*f2)); 
% Y1=(cos(2*pi*t*f1)./(A*(1+(m*t))));
% Y2=((B*(1+(m*t))).*cos(2*pi*t*f2));
% y=(cos(2*pi*t*f1)./(A*(1+(m*t))))+((B*(1+(m*t))).*cos(2*pi*t*f2));
% 
% PotenciaLF=(bandpower(Y1))*(1000)^2 
% PotenciaHF=(bandpower(Y2))*(1000)^2 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Modelo chirp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f0=0.04;
% f1=0.12;
% 
% f2=0.25;
% f3=0.18;
% 
% t1=300;
% p=1;
% phi=0;
% 
% beta   = (f1-f0).*(t1.^(-p));
% beta2   = (f3-f2).*(t1.^(-p));
% 
% Y1=cos(2*pi * ( beta./(1+p).*(t.^(1+p)) + f0.*t + phi/360));
% Y2=cos(2*pi * ( beta2./(1+p).*(t.^(1+p)) + f2.*t + phi/360));
% y = cos(2*pi * ( beta./(1+p).*(t.^(1+p)) + f0.*t + phi/360))+ cos(2*pi * ( beta2./(1+p).*(t.^(1+p)) + f2.*t + phi/360));
% 
% PotenciaLF=(bandpower(Y1))*(1000)^2 
% PotenciaHF=(bandpower(Y2))*(1000)^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Modelo compuesto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B=0.42;
% m=0.01;
% f1=0.1;
% 
% f2=0.25;
% f3=0.18;
% p=1;
% phi=0;
% t1=300;
% beta2 = (f3-f2).*(t1.^(-p));
% 
% 
%  
% %y = ((B*(1+(m*t)))).*cos(2*pi * ( beta2./(1+p).*(t.^(1+p)) + f2.*t + phi/360));
% Y1=((B*(1+(m*t))).*cos(2*pi*t*f1));
% Y2=2*cos(2*pi * ( beta2./(1+p).*(t.^(1+p)) + f2.*t + phi/360));
% y = ((B*(1+(m*t))).*cos(2*pi*t*f1)) + 2*cos(2*pi * ( beta2./(1+p).*(t.^(1+p)) + f2.*t + phi/360));
% 
% PotenciaLF=(bandpower(Y1))*(1000)^2 
% PotenciaHF=(bandpower(Y2))*(1000)^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=y;

y2=y-mean(y);

figure
plot(t,y,'*-blue');
hold on
plot(t,y2,'r');
title('Tacograma real interpolado - paciente sano (cambio de posición)');
legend('Tacograma crudo','Tacograma con remocion de media');
%ylim([0 2]);

segmentLength=length(y2);
[p,f]=pwelch(y,hann(segmentLength),segmentLength/2,2400,fs);
[p2,f2]=pwelch(y2,hann(segmentLength),segmentLength/2,2400,fs);

figure
%plot(f,p);
%hold on
plot(f2,p2);
xlim([0 0.5])

%x=h;

%%%%%%%%%%%%%%%%%%%% PARAMETROS PARA PROCESAMIENTO CON STFT %%%%%%%%%%%%%%% 
%fs=4;
duracion=100; %% TAMAÑO DE LA VENTANA DE ANALISIS EN SEGUNDOS
ventana=duracion*fs;

  a=x(1:ventana/2);
  aa=flip(a);
  b=x(end-(ventana/2):end);
  c=flip(b);

xX=cat(1,aa',x',c');
yy=detrend(xX);%yy=señal acondicionada, puntos de la ventana,traslape,nfft

gridSegundos=1; %% DESPLAZAMIENTO DEL GRID EN SEGUNDOS
puntoscorri=(gridSegundos*4); %+((gridSegundos/2)*4)
grid=ventana-puntoscorri;
nfft=2*ventana; %% NUMERO DE PUNTOS DE LA STFT

%%%%%%%%%%%%% APLICACION DEL ESPECTROGRAMA PARA EL TACOGRAMA AJUSTADO %%%%
[S,F,T,P]=spectrogram(yy,hann(ventana),grid,nfft,fs,'yaxis');
% % % 
figure
surf((T-T(1)),F,P,'FaceColor', 'interp',...
     'EdgeColor','none',...
     'FaceLighting','phong');
 axis tight;
 ylabel('Frecuencia (Hz)');
xlabel('Tiempo (s)');
title('Representacion tiempo-frecuencia de la HRV con STFT');
 colorbar
 colormap('jet');
 ylim([0 0.4]);
 caxis([0 max(max(P(:,:)))/20]);
 %zlim([0 0.3]);
 % % % axis ij;
% % % view(-45,35);
view(0,90);

VLF = [0.003,0.04]; %VLF = [0.003,0.041]; valor original
LF  = [0.041,0.15]; 
HF  = [0.151, 0.4];

iVLF = (F>=VLF(1)) & (F<VLF(2));
 iLF = (F>=LF(1)) &  (F<=LF(2));
 iHF = (F>=HF(1)) &  (F<=HF(2));
    %Find peaks
      %VLF Peak
      tmpF=F(iVLF);
      tmppsd=P(iVLF,1);
      [pks,ipks] = zipeaks(tmppsd);
      if ~isempty(pks)
        [tmpMax i]=max(pks);        
        peakVLF=tmpF(ipks(i));
      else
        [tmpMax i]=max(tmppsd);
        peakVLF=tmpF(i);
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %LF Peak
      tmpF=F(iLF);
      tmppsd=P(iLF,1);
      
      [pks,ipks] = zipeaks(tmppsd);
      if ~isempty(pks)
        [tmpMax i]=max(pks);
        peakLF=tmpF(ipks(i));
      else
        [tmpMax i]=max(tmppsd);
        peakLF=tmpF(i);
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
                 %HF Peak
      tmpF=F(iHF);
      tmppsd=P(iHF,1);
      
      [pks,ipks] = zipeaks(tmppsd);
      if ~isempty(pks)
        [tmpMax i]=max(pks);
        peakHF=tmpF(ipks(i));
      else
        [tmpMax i]=max(tmppsd);
        peakHF=tmpF(i);
      end
      
            % calculate raw areas (power under curve), within the freq bands (ms^2)
    aVLF=trapz(F(iVLF),P(iVLF,:))*(1000)^2
     aLF=trapz(F(iLF),P(iLF,:))*(1000)^2
     aHF=trapz(F(iHF),P(iHF,:))*(1000)^2
%     aTotal=aVLF+aLF+aHF;      
      

% % figure
% % plot(P(:,1));
% % hold on
% % plot(P(:,2),'r');
% % hold on
% % plot(P(:,3),'green');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ylim([0 0.5])
%xlim([0 15])

% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%%%%%%%%%%%% PARAMETROS PARA PROCESAMIENTO PARA CON CWT %%%%%%%%%%%%

[Sw,Fw,Pw]=cwt(h,fs);
sw=abs(Sw);


figure
surf(t/60,Fw,(sw), 'FaceColor', 'interp',...
  'EdgeColor','none',...
  'FaceLighting','phong');
           
axis tight;
view(0,90);
ylim([0 0.4]);
%zlabel('PSD (s^2/Hz)');
ylabel('Frecuencia (Hz)');
xlabel('Tiempo (s)');
title('Representacion tiempo-frecuencia de la HRV con CWT');
% % % %zlim([0 0.08]);
% % % title('Representacion tiempo-frecuencia con CWT');
colormap('jet');
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cfsVLF,Fw] = cwt(h,fs);
[cfsLF,Fw] = cwt(h,fs);
[cfsHF,Fw] = cwt(h,fs);
%[cfs2,f] = cwt(y,fs);
%T1 = 0;  %T2 = t(end,1);
F1 = 0.003;   F2 = 0.04;
F3=0.15; F4=0.4; 
%cfs(f < F1 & f > F2, t> T1 & t < T2) = 0;


cfsVLF(Fw < F1, t>0) = 0;
cfsVLF(Fw>F2, t>0)=0;

cfsLF(Fw < F2, t>0) = 0;
cfsLF(Fw>F3, t>0)=0;

cfsHF(Fw < F3, t>0) = 0;
cfsHF(Fw>F4, t>0)=0;

xrecVLF = icwt(cfsVLF);
xrecLF = icwt(cfsLF);
xrecHF = icwt(cfsHF);

PotVLFw=bandpower(xrecVLF)*(1000)^2
PotLFw=bandpower(xrecLF)*(1000)^2
PotHFw=bandpower(xrecHF)*(1000)^2
% % figure
% % cwt(xrecVLF,fs);
% % colormap('jet');
% % 
% % figure
% % cwt(xrecLF,fs);
% % colormap('jet');
% % 
% % figure
% % cwt(xrecHF,fs);
% % colormap('jet');

% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % 
datafile=y2;
 %==============Run EEMD=================================  
%give noise level
 noiselevel=0;
%give Number of ensemble 
 Nensemble=1;
%do eemd  and show result
 EEMDIMF=eemd(datafile,noiselevel,Nensemble);
%plot the EEMD result 
%      figure 
%      anoiselevel=num2str(noiselevel);
%      aNensemble=num2str(Nensemble);
%      strips(EEMDIMF);
%      title(['EEMD result    ,noise level=',anoiselevel,' Number of ensemble=',aNensemble]);
%      clear noiselevel Nensemble
 disp('===== EEMD calculation complete! =====')    
 
%==============Run Hilbert-Energy Spectrum for EEMD =======================   
%give samplerate 
  samplerate=4;
%give frequency-axis resolution for hilbert-spectrum
  freqsol=600;%200; %120;% 200 %1080;%medida original    %para ver bandas es suficiente en freqsol =300
%give time-axis resolution for hilbert-spectrum
  timesol=150;%50;% 100 %2160; %%original 500 %%medida orginal de 333     %para ver en tiempo es suficiente en timesol=600
%give frequency-axis maximun value
  hsp_fre1=2;

% dealing with EEMD result 
     au=size(EEMDIMF);
     nIMF=au(2)-2;
     nPT=au(1)-1;
     totalt=(nPT+1)/samplerate;
     Xlow=1/totalt;
     Xhig=samplerate/2;
     
%  figure
%  for i=1:nIMF
%  subplot(5,2,i);
%  plot(EEMDIMF(:,i))
%  title('IMFs')
%  end
 
%Calculate Hilbert-Energy-Spectrum(NNSPE.m)
%EEMDIMF(1:nPT,2:nIMF+1)
[nte,tae,fae]=nnspe(EEMDIMF(1:nPT,2:3), 0, totalt, freqsol, timesol, 0.00001, hsp_fre1,0,totalt);     
%[nte,tae,fae]=nnspe(EEMDIMF(1:nPT,[2,3:nIMF-4]), 0, totalt, freqsol, timesol, 0.00001, hsp_fre1,0,totalt);
     %nnspe(EEMDIMF(1:nPT,2:nIMF-1), 0, totalt, freqsol, timesol, 0.00001,%hsp_fre1,0,totalt);
 
 %Smooth the pic ,use an Gaussian Filter "q",let nt been filtered twice,HESP spectrum will be good for reading
      q=fspecial('gaussian', 7, 0.9);
      nse=filter2(q, nte);
      nsue=filter2(q, nse);
  
% fig = figure;nsue.^.5
 figure;  
 subplot(2,1,1)
 plot(t/60,y2,'red');
 subplot(2,1,2)
 surf(tae/60,fae,nsue.^.5,'FaceColor', 'interp',...
     'EdgeColor','none',...
     'FaceLighting','phong');
colormap('jet');
% % ylim([0 0.5])
% % axis tight;
%axis([0 totalt/60 0 0.4 0 0.5]); %0.5
%caxis([0 0.08])
%view(100, 36);
view([0 90]);
ylim([0 0.4]);
xlim([0 5]);
xlabel('Tiempo (s)');
title('Representacion tiempo-frecuencia de la HRV con HHT');
zlabel('PSD (s^2/Hz)');
ylabel('Frecuencia (Hz)');
colorbar
%grid off

PotLFhht=bandpower(EEMDIMF(:,3))*(1000)^2
PotHFhht=bandpower(EEMDIMF(:,2))*(1000)^2
%PotLF1hht=bandpower(EEMDIMF(:,4))*(1000)^2
