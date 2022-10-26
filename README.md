# Projet_estimation-synchronisation

À travers ce TP, nous cherchons à mettre en œuvre nos connaissances en traitement du signal
afin de correctement restituer une information modulée par une QPSK reçue via une chaîne de
communication numérique classique.
Une synchronisation temporelle, paquet et en phase a été effectué durant ce TP pour obtenir un 
signal correct. 

Cette implémentation a été effectué pour deux Fse différentes ce qui a permis d’utiliser
plusieurs méthodes notamment au niveau de la synchronisation temporelle

Synchronisation en phase : Une boucle en quadrature a été utilisé.
Synchronisation temporelle : dérivée de l'erreur quadratique moyenne, filtre de Farrow.
Synchronisation paquet : utilisation d'une séquence pilote pour détecter le début de la trame. 
Utilisation d'une intercorrélation pour détecter le début des données.
