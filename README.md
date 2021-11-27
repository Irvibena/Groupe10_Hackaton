# Groupe10 Hackaton AMI2B 2021

Nous détaillons ici dans une première partie le mode opératoire pour lancer le workflow. Dans une seconde partie, nous donnons les instructions pour se connecter à notre VM où les résultats du workflow sont déjà présents.

Pour consulter le workflow, voir les fichiers [workflow.nf](workflow.nf) et la configuration [nextflow.config](nextflow.config).

## Lancer le workflow _de novo_

### MobaXTerm

Lancer **MobaXTerm**.  
Sous **Tools**, cliquer sur **MobaKeyGen**. Générer une clef publique et une clef privée et les sauvegarder.  

### IFB Cloud

Puis il faut créer une machine virtuelle sur [IFB Cloud](https://biosphere.france-bioinformatique.fr/).  
Se connecter en sprécifiant qu'on étudie à Université Paris Saclay.  
Dans son profil, sélectionner **Paramètres** et coller la clef publique dans **PubKey**.  
Sous [RAINBio](https://biosphere.france-bioinformatique.fr/catalogue/) > [BioPipes](https://biosphere.france-bioinformatique.fr/catalogue/appliance/119/), cliquer sur **Lancer** puis sur **Déploiement Avancé**. Donner un nom à la VM, choisir **ReproHack2021** comme groupe, **ifb-core-cloud** comme cloud, et choisir le gabarit **ifb.m4.4xlarge**.  
Attendre que la VM se lance, puis ouvrir la page de la VM et copier l'adresse IP.  

### MobaXTerm

De retour sous MobaXTerm, cliquer sur **Session** > **SSH**.  
Coller l'adresse IP de la VM dans **Remote Host**.  
Cocher **Specify Username** et écrire `ubuntu`.  
Sous **Advanced SSH settings**, cocher **Use private key** et aller chercher le fichier contenant la clef privée.  
Cliquer sur OK.

### GitHub

Une fois connecté à la VM depuis MobaXTerm, éxécuter `ssh-keygen -t rsa` puis `cat ~/.ssh/id_rsa.pub`. Copier la clef SSH.  
Dans son profil [GitHub](https://github.com/), sous **Settings** > **SSH and GPG keys**, cliquer sur **New SSH key** et ajouter la clef SSH copié précédemment.  
Sur la VM (via MobaXTerm), se placer dans le bon dossier, _i.e._ exécuter `cd data/mydatalocal/`.  
Enfin, exécuter `git clone git@github.com:Irvibena/Groupe10_Hackaton.git`.  
Taper `cd Groupe10_Hackaton/` pour se placer dans le répertoire.

### Le workflow

Pour pouvoir utiliser le workflow, taper `conda activate`, puis `conda install nextflow` et `conda update nextflow`.  
Enfin, lancer l'éxécution du workflow avec `nextflow run workflow.nf`.
