clear all
clc

%nuevoCorre=load ('RRrem_p12s1.txt'); %RRremST13s1mod;
%nuevo(:,2)=detrend(nuevoCorre(:,2));
%nuevo(:,1)=nuevoCorre(:,1);
nuevo=load ('RRrem_ST1s1.txt');  %%corregir linea de tiempo del archivo 12
Horas=input('Ingrese el tiempo de terapia(horas): ')
Minutos=input('Ingrese el delta de tiempo (minutos): ')
Soloterapia=nuevo(Minutos*60*4:(Minutos*60*4)+(Horas*60)*60*4,:);

figure
plot(nuevo(:,1)/60,nuevo(:,2));
hold on
plot(Soloterapia(:,1)/60,Soloterapia(:,2));
xlabel('tiempo (minutos)')
ylabel('Intervalo RR (s)')
title('Tacograma real sin tendencia')
legend('Registro de tacograma completo','Registros de tacograma durante la sesión')

tamano=length(Soloterapia(:,1));
trem=Soloterapia(:,1);
x=Soloterapia(:,2); %(Aqui ya estan dado en segundos)
%X=x*1000;
%xCorr=ada_f(x);
%xCorr=detrend(x);
%xCorr=xCorr-mean(xCorr);
xCorr=x;
%figure
%plot(nuevo(:,1));
%hold on
%plot(nuevoST(:,1))


figure

plot(trem/60,xCorr);
xlabel('tiempo (minutos)')
ylabel('Intervalo RR (s)')
title('Tacograma solo durante la terapia')

%%%%%%%%%%%%%%%%%%%% PARAMETROS PARA PROCESAMIENTO CON STFT %%%%%%%%%%%%%%% 
fs=4;
duracion=600; %% TAMAÑO DE LA VENTANA DE ANALISIS EN SEGUNDOS
ventana=duracion*fs;

  a=x(1:ventana/2);
  aa=flip(a);
  b=x(end-(ventana/2):end);
  c=flip(b);

xX=cat(1,aa,xCorr,c);
yy=detrend(xX);%yy=señal acondicionada, puntos de la ventana,traslape,nfft

gridSegundos=360; %% DESPLAZAMIENTO DEL GRID EN SEGUNDOS
puntoscorri=(gridSegundos*4); %+((gridSegundos/2)*4)
grid=ventana-puntoscorri;
nfft=2*ventana; %% NUMERO DE PUNTOS DE LA STFT

%%%%%%%%%%%%% APLICACION DEL ESPECTROGRAMA PARA EL TACOGRAMA AJUSTADO %%%%
[S,F,T,P]=spectrogram(yy,hann(ventana),grid,nfft,fs,'yaxis');
% % % 
figure

surf((T-T(1))/60,F(5:end,1),P(5:end,:),'FaceColor', 'interp',...
     'EdgeColor','none',...
     'FaceLighting','phong');
 
% surf((T-T(1))/60,F(5:end,1),P(5:end,:),'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud'); %,...
     %'EdgeColor','none',...
     %'FaceLighting','phong');
%surf((T/60)-2.5,F(1:end,1),P(1:end,:));
axis tight;
% % % 
%axis ij;
%view(-45,35);
view(0,90);
ylim([0 0.4]);
% % % title('Representacion tiempo-frecuencia con STFT');
colormap('jet');
colorbar

%label={'apple'};
%lcolorbar(label,'fontweight','bold');
% % % %caxis([0 max(max(P(:,:)))/10]);
%caxis([0 max(max(P(5:end,:)))/10]);
caxis([0 0.044]);
%labels = {'ms^2/Hz'};
%lcolorbar('s^2/Hz','fontweight','bold');
xlabel('time (minutes)','fontweight','bold');
ylabel('Frequency (Hz)','fontweight','bold');
title('Spectrogram  HRV with STFT','fontweight','bold');
%print(gcf,'-dtiff',strcat('C:\PROYECTOS\Fibromialgia\Datos2aParte\Figuras\',nombre,'RR24hrs.tif'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%print(gcf,'-dtiff',strcat('C:\PROYECTOS\STFThemo\p1s1Hemo.tif'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%axis ij   %inversion de ejes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VLF = [0.003,0.04]; %VLF = [0.003,0.041]; valor original
LF  = [0.04,0.15]; 
HF  = [0.15, 0.4];

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
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
            % calculate raw areas (power under curve), within the freq bands (ms^2)
    aVLF=trapz(F(iVLF),P(iVLF,:))*(1000)^2;
     aLF=trapz(F(iLF),P(iLF,:))*(1000)^2;
     aHF=trapz(F(iHF),P(iHF,:))*(1000)^2;
%     aTotal=aVLF+aLF+aHF;
        
    %fprintf('\n Power VLF(ms^2): %5.5f\n\r',aVLF);
    figure
    plot(((T-T(1))/60),aVLF,'-*');  
    hold on
    plot(((T-T(1))/60),aLF,'-o');
    hold on
    plot(((T-T(1))/60),aHF,'-^');
    legend('VLF','LF','HF')
    xlabel('Tiempo (minutos)')
    ylabel('Potencia (ms^2)');
    
    %PromedioVLF=sum(aVLF)/length(T);
    %PromedioLF=sum(aLF)/length(T);
    %PromedioHF=sum(aHF)/length(T);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aPotTotal=aVLF+aLF+aHF;
aLFHF=(aLF./aHF);

PotNormLF=(aLF./(aPotTotal-aVLF))*100;
PotNormHF=(aHF./(aPotTotal-aVLF))*100;
%aLFHFnorm=PotNormLF./PotNormHF;
%promedioLFHF=sum(aLFHF)/length(T);
figure
subplot(3,1,1)
plot((T-T(1))/60,PotNormLF)
legend('Potencia normalizada LF')
xlabel('tiempo (minutos)')
ylabel('Potencia (u.n.)')

subplot(3,1,2)
plot((T-T(1))/60,PotNormHF)
legend('Potencia normalizada HF')
xlabel('tiempo (minutos)')
ylabel('Potencia (u.n.)')

subplot(3,1,3)
plot((T-T(1))/60,aLFHF)
hold on
legend('LF/HF')
xlabel('tiempo (minutos)')
ylabel('LF/HF')

%subplot(4,1,4)
%plot((T-T(1))/60,aLFHFnorm)
%legend('LF/HF')
%xlabel('tiempo (minutos)')
%ylabel('LF/HF normalizado')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VectoresPotencia=[aVLF; aLF; aHF; aPotTotal; aLFHF]';
VectoresNormalizados=[PotNormLF; PotNormHF]';


%save PotNorma-p20s1.txt VectoresNormalizados -ascii;

%Potencias=[PromedioVLF; PromedioLF; PromedioHF; promedioLFHF];

% Auto=aVLF';
% save PruebaIndependencia.txt Auto -ascii
% 
% 
% YYY=xcorr(Auto,Auto); % matlab built in function
% len=length(YYY);
% index=len/2;
% h=YYY(index:1:end); %extacting one side of the result
% figure;
% plot((1:length(h))/6,(h./max(h)),'o');
% title('MATLAB xcorr function OUTPUT'); % output plot as per matlab function
% [m1,n1]=max(h) % max value in the correlation and its index
