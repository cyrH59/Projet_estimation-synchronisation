%%% FSE = 10
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
Fse = 10;
Ns=65536; 
% Filtre adapté
sps=10;   
span=16;
rolloff=0.5;
b = rcosdesign(rolloff, span, sps,"sqrt");
rn=conv(signal_recu,b);

%% Synchronisation temporelle (recherche tau)

rnavantsynchtf=rn;
N=floor(length(rn)/Fse); 
Ftau=zeros(1,Fse);                                                          % Contient l'ensemble des valeurs de F(tau) 
for tau=1:Fse
    for m=1:N-1
        Ftau(tau)=Ftau(tau)-(1/Ns)*((real((rn(m*Fse-Fse+tau)-rn(m*Fse+tau))*conj(rn(m*Fse-(Fse/2)+tau)))));
    end
end
tau=ZeroPassage(Ftau);                                                      % Determination automatique de tau à l'aide de la courbe passant par 0 (cf F(tau))


%% Detection maximum de corrélation (Synchronisation trame)

Spilote=theMod(bi2de(bitSynchro'));                                         % Séquence pilote
rn = rn(tau:Fse:end);                                                       % Sous echantillonnage 
N=150;                                                                      % Nombre d'echantillons à appliquer sur notre algo
[dt, ck]=preambuledetect(Spilote,rn,N);                                     % Recherche du début de l'information
PiloteEstimee=rn(dt:1:dt+length(Spilote)-1);                                % Séquence pilote de la trame reçu
rn = rn(dt+length(Spilote)+8:1:end);                                        % Application de la synchro Trame


%% Implémentation Synchro phase :
rsanssynchrof=rn;
% Boucle à quadrature
rnbis=rn.^4; 
phiinit=0;                                                                  % Phase initiale 
alpha=1;
beta=0.1;
phi=zeros(1,length(rn ));                                                   % Phase
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
% Ambiguïté de phase
Ambiguite=sum((PiloteEstimee).*conj(Spilote)./abs(Spilote).^2);             % Calcul de l'ambiguïté de phase
Ambiguite=pi*round(angle(Ambiguite)/pi);                                    % Arrondissement à 0 ou à pi (phase la plus proche)
rn=rn*exp(-1i*Ambiguite);                                                   % Correction de l'ambiguïté de phase



%% Décodage de source sans synchro f et sans synchro t
demodulation=pskdemod(rnavantsynchtf(1:Ns), 4,0,'gray');
hatB=de2bi(demodulation,2);
hatB=hatB';
hatB=hatB(:)';

hatMatBitImg = reshape(hatB(:),[],8);
matImg = bi2de(hatMatBitImg);
Imgas2 = reshape(matImg,128,128);


%% Décodage de source avec synchro f
demodulation=pskdemod(rn(1:Ns), 4,0,'gray');
hatB=de2bi(demodulation,2);
hatB=hatB';
hatB=hatB(:)';

hatMatBitImg = reshape(hatB(:),[],8);
matImg = bi2de(hatMatBitImg);
Imgas = reshape(matImg,128,128);


%% Décodage de source sans synchro f
demodulation=pskdemod(rsanssynchrof(1:Ns), 4,pi/4,'gray');
hatB=de2bi(demodulation,2);
hatB=hatB';
hatB=hatB(:)';

hatMatBitImg = reshape(hatB(:),[],8);
matImg = bi2de(hatMatBitImg);
Imgss = reshape(matImg,128,128);




%% Affichage des differents résultats 


% Affichage image : 
figure;
subplot(1,3,2)
imshow(uint8(Imgss))
title('Synchronisation temporelle');
subplot(1,3,1)
imshow(uint8(Imgas2))
title('Aucune Synchronisation');
subplot(1,3,3)
imshow(uint8(Imgas))
title('Synchronisation temporelle et fréquentielle')

figure;
plot(phi);
title("Evolution de la Phase")
xlabel("n°Iteration")
ylabel("Phase en °")


% Affichage constellation 
figure;
subplot(1,3,3)
plot(rn,'g.');
title("Synchronisation temporelle et fréquentielle")
subplot(1,3,1)
plot(rnavantsynchtf,'r.');
title("Sans synchronisation")
subplot(1,3,2)
plot(rsanssynchrof,'y.');
title("Synchronisation temporelle")

% Affichage F(tau) :
figure; 
plot(Ftau) %dérivée seconde positive 
title("F("+char([0xD835 0xDF0F])+")");
xlabel(char([0xD835 0xDF0F])); 

% Affichage produit de convolution 
figure();
plot(1:N,abs(ck),'r')
title("Produit de convolution");
xlabel("Indice");

