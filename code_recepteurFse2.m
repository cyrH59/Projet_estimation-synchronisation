%%% FSE = 2
clear;
close all;
clc;

%% Chargement du fichier contenant le signal reçu
load 'signal_recu.mat';
load bitSynchro.mat;

%% La modulation utilisée pour encoder les données
M = 4;                                                                      % Ordre de la modulation
theMod = comm.PSKModulator(M, pi/4,'BitInput',0,'SymbolMapping','Gray');    % Le modulateur
theDemod = comm.PSKDemodulator(M, pi/4, 'DecisionMethod','HARD DECISION');  % Le démodulateur

%% Traitement du signal par le récepteur
% paramètres
Ns=65536;
Fse=2;
% Filtre adapté
sps=2;  
span=16;
rolloff=0.5;
b = rcosdesign(rolloff, span, sps,"sqrt");
signal=conv(signal_recu,b);
rnavantsynchtf=signal;
rn=signal(1:5:end);

%%% Synchronisation temporelle

N=floor(length(rn)/Fse);
h1=[0.5 -0.5 -0.5 0.5];
h2=[-0.5 -0.5 1.5 -0.5];
h3=[0 1 0 0];
terme1=filter(h1,1,rn);
terme2=filter(h2,1,rn);
terme3=filter(h3,1,rn);
ltau=0:0.1:1.9;
Ftau = zeros(1,length(ltau));
for ii = 1:length(ltau)
    tau = ltau(ii);
    dec1=(1:N-2)*Fse+floor(tau)+1;
    dec2=(1:N-2)*Fse+Fse+floor(tau)+1;
    dec3=(1:N-2)*Fse+Fse/2+floor(tau)+1;
    ftau = tau - floor(tau);
    y1=ftau^2*terme1(dec1)+ftau*terme2(dec1)+terme3(dec1);
    y2=ftau^2*terme1(dec2)+ftau*terme2(dec2)+terme3(dec2);
    y3=ftau^2*terme1(dec3)+ftau*terme2(dec3)+terme3(dec3);
    Ftau(ii)=mean(real((y1-y2).*conj(y3)));

end

tau = 1.7;
dec1=(1:N-2)*Fse+floor(tau)+1;
dec2=(1:N-2)*Fse+Fse+floor(tau)+1;
dec3=(1:N-2)*Fse+Fse/2+floor(tau)+1;
ftau = tau - floor(tau);
y1=ftau^2*terme1(dec1)+ftau*terme2(dec1)+terme3(dec1);





%% Detection maximum de corrélation (Synchronisation trame)

Spilote=theMod(bi2de(bitSynchro'));                                         % Séquence pilote
rn = rn(tau:Fse:end);                                                       % Sous echantillonnage 
N=150;                                                                      % Nombre d'echantillons à appliquer sur notre algo
[dt, ck]=preambuledetect(Spilote,rn,N);                                     % Recherche du début de l'information
PiloteEstimee=rn(dt:1:dt+length(Spilote)-1);                                % Séquence pilote de la trame reçu
rn = rn(dt+length(Spilote)+8:1:end);                                        % Application de la synchro Trame

%% Implémentation decallage fréquenciel :
rsanssynchrof=rn;
% Synchro frequence
rnbis=rn.^4; %
scatterplot(rn);
title("rn avant")
phiinit=0; %Phase initiale
alpha=3; 
beta=0.5;

phi=zeros(1,length(rn ));
phi(1)=phiinit;
for k=2:length(phi)
    etmp=imag(rnbis(k)*exp(-1i*phi(k-1)));
    b1=[alpha+beta -alpha];
    a1=[1 -1];
    if(k==2)
        [vn, zf1]=filter(b1,a1,etmp);
        b2=[0 1];
        a2=[1 -1];
        [phitmp, zf2]=filter(b2,a2,vn);
        phi(k)=phitmp/2;
    end
    if (k~=2)
        [vn, zf1]=filter(b1,a1,etmp,zf1);
        b2=[0 1];
        a2=[1 -1];
        [phitmp, zf2]=filter(b2,a2,vn,zf2);
        phi(k)=phitmp/2;
    end
    rn(k)=rn(k)*exp(-1i*phi(k)/4);


end
%% Ambiguïté de phase
Ambiguite=sum((PiloteEstimee).*conj(Spilote)./abs(Spilote).^2);             % Calcul de l'ambiguïté de phase
Ambiguite=pi*round(angle(Ambiguite)/pi);                                    % Arrondissement à 0 ou à pi (phase la plus proche)
rn=rn*exp(-1i*Ambiguite);                                                   % Correction de l'ambiguïté de phase





%% Décodage de source avec synchro f
demodulation=pskdemod(rn(1:Ns), 4,0,'gray');
hatB=de2bi(demodulation,2);
hatB=hatB';
hatB=hatB(:)';

hatMatBitImg = reshape(hatB(:),[],8);
matImg = bi2de(hatMatBitImg);
Img = reshape(matImg,128,128);
figure;
imshow(uint8(Img))
title('image après toutes les synchronisations');

%% Décodage de source sans synchro f
demodulation=pskdemod(rsanssynchrof(1:Ns), 4,pi/4,'gray');
hatB=de2bi(demodulation,2);
hatB=hatB';
hatB=hatB(:)';

hatMatBitImg = reshape(hatB(:),[],8);
matImg = bi2de(hatMatBitImg);
Img = reshape(matImg,128,128);

figure;
imshow(uint8(Img))
title('image sans synchronisation fréquentiel');


%% Affichage

figure;
plot(phi);
title("Evolution de phi")

scatterplot(rn);
title("rn après")


figure;
subplot(1,3,3)
plot(rn,'g.');
title("Synchronisation temporelle et fréquentielle (Fse =2)")
subplot(1,3,1)
plot(rnavantsynchtf,'r.');
title("Sans synchronisation (Fse=2)")
subplot(1,3,2)
plot(rsanssynchrof,'y.');
title("Synchronisation temporelle (Fse=2)")
figure;
plot(Fse-ltau,Ftau,'r') %dérivée seconde négative
title("F("+char([0xD835 0xDF0F])+")");
xlabel(char([0xD835 0xDF0F])); 

figure; 
plot(abs(y1))


