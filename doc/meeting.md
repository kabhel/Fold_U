# Meeting 1 : lundi 1/10

## Répartition des rôles :

* **Manager** ~ Hélène
 - Organisation des meetings
 - Gestion de l'avancement du projet
 - Gestion du temps

* **Expert code/github** ~ Gabriel/Arnold
 - Répartition du code entre contributeurs
 - Gestion des branches dans GitHub
 - Vérification du code
 - Choix des classes / fichiers

* **Expert biblio** ~ Flora/Arnold
 - Recherche d'articles intéressants

* **Expert rapport/présentation** ~ Tom
 - Bon niveau en anglais / correction anglais
 - Gestion du rapport en latex (overleaf)
 - Gestion présentation diapos
  
## Commandes GitHub

#### Cloner le repository
```shell
git clone https://github.com/meetU-MasterStudents/2018---2019-Equipe-4.git
```

#### Ajout / Modification fichiers
```shell
cd repository
git add nom_du_fichier
git commit -m "commentaire"
git push
```

#### "Tirer" le repository
```shell
cd repository
git pull
```

### Create a new repository  
```shell
git init
git add README.md
git commit -m "first commit"
git remote add origin https://github.com/...
git push -u origin master
```
### Astuce : Sauvegarde des log de GitHub pendant 8h

Dans le ~/.gitconfig :
```
[user]
	email = mail@gmail.com
	name = pseudo
[core]
	pager = less
[url "https://pseudo@github.com"]
	insteadOf = https://github.com
[credential]
helper = cache --timeout=28800
```

---

# Meeting 2 : jeudi 11 octobre

## Ajout des datas
data test : Agglutini

## Pseudo-code et répartition du code
- génération des classes ~ Hélène
- parsing_folrec ~ Hélène
- parsing_metafold ~ Tom
- parsing_DOPE ~ Arnold
- calc_threading_score
  - calc_dist : Gabriel
  - transf_dist_E : FLora


## Git pull request
1. Création d'une branche :
Une branche = une feature
```shell
git checkout -b dev_feature
```
2. Merge
Le merge se fait une fois le pull request accepté par tout le monde

