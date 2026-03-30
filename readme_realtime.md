Install Xlib:
sudo apt-get install libx11-dev

DISPLAY=:0.0 export DISPLAY

To compile:
g++ main_rt.cpp -pthread -O2 -I/usr/include/X11/lib -L/usr/lib/x86_64-linux-gnu/ -lX11

To run:
./a.out disk texture (png) sensitivity data (txt) point spread function (kernel for convolution) (png) scene parameters (txt) rendering parameters (txt) output file name (without the .png extension)

Note: the number of measurements as well as the step size is hardcoded for the sensitivity data file.

Example with the repo files:
./a.out ./data/adisk.png ./data/sensitivity.txt ./data/psf.png ./scenes/scene4.txt ./scenes/qualityH_a.txt res

The scene file contains in order:
bh.a scn.camera.x scn.camera.y scn.camera.z disk.R_max disk.betamax disk.Tmax disk.texture_rep rdr.R_inf

bh.a: Kerr parameter of the black hole (between 0 and 1)
scn.camera.x scn.camera.y scn.camera.z: camera position (approximate due to Cartesian to BL conversion approximation)
disk.R_max: maximum radius of the dust disk
disk.betamax: velocity of the disk closest to the black hole (in natural units)
disk.TMax: temperature of the disk closest to the black hole
disk.texture_rep: number of repetitions of the disk texture along its length to avoid pixelation
rdr.R_inf: distance beyond which the ray is considered at infinity and the simulation stops; must be greater than the camera distance and the disk’s maximum radius

The rendering parameters file contains in order:
rdr.width rdr.height rdr.ChunkSize rdr.stepmax rdr.stepmin rdr.delta rdr.SamplesPerPixels

rdr.width: render width (horizontal size)
rdr.height: render height (vertical size)
rdr.chunkSize: number of pixels per block processed in parallel (must divide the total number of pixels)
rdr.stepmax: maximum step size
rdr.stepmin: minimum step size
rdr.delta: angular spacing (in pixels) between rays within the same beam
rdr.SamplesPerPixels: samples per pixel (one sample = a beam of 5 rays)
1: no anti-aliasing
16: very high quality


Installer Xlib: sudo apt-get install libx11-dev

DISPLAY=:0.0
export DISPLAY

Pour compiler:  g++ main_rt.cpp -pthread -O2 -I/usr/include/X11/lib -L/usr/lib/x86_64-linux-gnu/ -lX11
Pour executer:  ./a.out texture du disque(png) données de sensibilité(txt) fonction d'étalement du point (noyau pour la convolution)(png) parametres de la scene(txt) parametres de rendu(txt) nom du fichier de sortie(sans l'exension .png)

Remarque: le nombre de mesures ainsi que le pas est hardcodé pour le fichier des données de sensibilité


exemple avec les fichiers du repo: 
./a.out ./data/adisk.png ./data/sensitivity.txt ./data/psf.png ./scenes/scene4.txt ./scenes/qualityH_a.txt res 



le fichier scene contient dans l'ordre:
bh.a scn.camera.x scn.camera.y scn.camera.z disk.R_max disk.betamax disk.Tmax disk.texture_rep rdr.R_inf 

bh.a: Parametre de Kerr du trou noir (entre 0 et 1)
scn.camera.x scn.camera.y scn.camera.z position de la camera (approx car passage cartesienne à BL approx)
disk.R_max: Rayon maximal du disque de poussière
disk.betamax: Vitesse du disque au plus proche du trou noir (dans le système d'unité naturel)
disk.TMax: Temperature du disque au plus proche du trou noir
disk.texture_rep: Nombre de répétition de la texture du dique en longeur pour ne pas qu'elle soit pixelisée
rdr.R_inf: Distance à partir de laquelle on considere etre a l'infini et on arrete la simulation du rayon, doit être plus grand que la distance de la camera et le rayon maximal du disque



le fichier parametre de rendu contient dans l'ordre:
rdr.width rdr.height rdr.ChunkSize rdr.stepmax rdr.stepmin rdr.delta rdr.SamplesPerPixels
      
rdr.width: Taille du rendu (horizontale)
rdr.height: Taille du rendu (verticale)
rdr.chunkSize: Nombre de pixels par blocs traités en parallele (doit diviser le nombre de pixel total)
rdr.stepmax: Pas max
rdr.stepmin: Pas min
rdr.delta: Ecart angulaire (en pixel) entre les rayons d'un même faisceau
rdr.SamplesPerPixels: echantillons par pixels (un echantillon = un faisceau de 5 rayons) 1: sans antialisaing, 16: pour une très bonne qualité

