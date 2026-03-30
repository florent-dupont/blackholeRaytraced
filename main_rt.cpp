#include <stdint.h>
#include <iostream>
#include <math.h>
#include <tgmath.h>
#include <string.h>
#include <fstream>

#include <chrono>
using namespace std::chrono;

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#define STB_IMAGE_IMPLEMENTATION
#include "./libs/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./libs/stb_image_write.h"

#include "./libs/par_for.h"

//Un des deux ou aucun des deux
#define DRAWSTARS 0 //Si 1 etoiles sur le fond celeste, sinon rien
#define DRAWGRID 0  //Si 1 grille sur le fond celeste

//Un des deux ou aucun des deux
#define ADISK_NORMAL 1 //Si 1, texture + temperature, sinon rien
#define ADISK_GRID 0   //Si 1, grille, sinon rien
//1 0 0 0 ne marche pas
#include "fct.h"
#include "raytracing.h"

#include <stdio.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include "graph.h"

Window root, win;
GC graphcontext;
Display *disp;

//XEvent ev;

unsigned long *imageToDisplay;
int *reconstructionMatrix;

#define CHANNEL_NUM 3

#define MAXITER 500000

static const int Neq = 6; //Nombre de variables pour l'integration numerique

struct Rendu rdr;
struct Scene scn;
struct Blackhole bh;
struct Disk disk;

int adisk_width, adisk_height, adisk_bpp;

int kernel_width, kernel_height, kernel_bpp;

uint8_t *adisk = NULL;
uint8_t *kernel = NULL;

static const int SpectrumSampleSize = 76; //nombre d'échantillons
double deltaWaveLength = 4.;              //(En nm, espacement des échantillons)
double wavelengthSamples[SpectrumSampleSize], wavelengthSamples5[SpectrumSampleSize], sensitivitySamplesR[SpectrumSampleSize], sensitivitySamplesG[SpectrumSampleSize], sensitivitySamplesB[SpectrumSampleSize];

void readSensitivityData(char *filename, double *wavelengthSamples, double *wavelengthSamples5, double *sensitivitySamplesR, double *sensitivitySamplesG, double *sensitivitySamplesB)
{
    int cnt = 0;
    ifstream source;       // build a read-Stream
    source.open(filename); // open data

    if (source)
    {
        for (std::string line; std::getline(source, line);) //read stream line by line
        {
            std::istringstream in(line); //make a stream for the line itself

            in >> wavelengthSamples[cnt] >> sensitivitySamplesR[cnt] >> sensitivitySamplesG[cnt] >> sensitivitySamplesB[cnt]; //now read the whitespace-separated doubles

            wavelengthSamples5[cnt] = pow(wavelengthSamples[cnt], 5);
            cnt++;
        }
        cout << "Données de sensibilité trouvées" << endl;
    }
    else
    {
        cout << "Données de sensibilité non trouvées" << endl;
        exit(0);
    }
}

void getBodyColor(double *rgbR, double *rgbG, double *rgbB, double temperature, double brightness)
{
    double I;

    for (int l = 0; l < SpectrumSampleSize; ++l)
    {

        I = (9.e14 / wavelengthSamples5[l]) / exp(1.43913e7 / (wavelengthSamples[l] * temperature) - 1.) * brightness; //6e14 pour renormaliser les composante des pixels

        I *= deltaWaveLength;
        *rgbR += sensitivitySamplesR[l] * I;
        *rgbG += sensitivitySamplesG[l] * I;
        *rgbB += sensitivitySamplesB[l] * I;
    }
}

//Affiche une grille sur le disque
void getDiskColorGrid(double phi, double r, double *rgbR, double *rgbG, double *rgbB)
{
    phi = (phi + M_PI) / (2. * M_PI);

    bool a = int((100 * phi)) % 2; //100 alternances de couleur par tour

    bool b = ((r - disk.R_min) / (disk.R_max - disk.R_min) > .5); //Séparer le disque en 2 radialement

    if (a ^ b) //Ou exclusif
    {
        *rgbR = 255.;
        *rgbG = 0.;
    }
    else
    {
        *rgbR = 0.;
        *rgbG = 255.;
    }
}

//Passage aux coordonnees de Boyer Lindquist
//Pour l'instant juste passage au coordonnées sphériques (plus simple)
//Utilisé uniquement pour trouver la pos initiale de la camera, pas important
void cartesianToBl(double x, double y, double z, double *r, double *theta, double *phi)
{

    double r2 = x * x + y * y + z * z;
    *r = sqrt(r2);

    *phi = atan2(y, x);
    *theta = acos(z / (*r));
}

//Passage des coordonnées de Boyer Lindquist au coordonnées cartésienne (exacte)
void blToCartesian(double r, double theta, double phi, double *x, double *y, double *z)
{

    double sintheta = sin(theta);
    double costheta = cos(theta);
    double cosphi = cos(phi);
    double sinphi = sin(phi);
    double temp = sintheta * sqrt(r * r + bh.a2);

    *x = temp * cosphi;
    *y = temp * sinphi;
    *z = r * costheta;
}

/* Fonction utilisé pour l'intégration */
void geodesic(double L, double kappa, double *y, double *dydx)
{
    double r, theta, pr, ptheta;

    r = y[0];
    theta = y[1];
    pr = y[4];
    ptheta = y[5];

    double r2 = r * r;
    double twor = 2.0 * r;

    double sintheta, costheta;
    sintheta = sin(theta);
    costheta = cos(theta);

    double costheta2 = costheta * costheta;
    double sintheta2 = sintheta * sintheta;

    double sigma = r2 + bh.a2 * costheta2;
    double delta = r2 - twor + bh.a2;
    double sd = sigma * delta;
    double siginv = 1.0 / sigma;
    double bot = 1.0 / sd;

    /* Prevent problems with the axis */
    if (sintheta < 1e-8)
    {

        sintheta = 1e-8;
        sintheta2 = 1e-16;
    }

    dydx[0] = -pr * delta * siginv;
    dydx[1] = -ptheta * siginv;
    dydx[2] = -(twor * bh.a + (sigma - twor) * L / sintheta2) * bot;
    dydx[3] = -(1.0 + (twor * (r2 + bh.a2) - twor * bh.a * L) * bot);
    dydx[4] = -(((r - 1.0) * (-kappa) + twor * (r2 + bh.a2) - 2.0 * bh.a * L) * bot - 2.0 * pr * pr * (r - 1.0) * siginv);
    dydx[5] = -sintheta * costheta * (L * L / (sintheta2 * sintheta2) - bh.a2) * siginv;
}

/* Conditions initiales pour un rayon */
void initial(double r0, double theta0, double *L, double *kappa, double *y0, double *ydot0, double x, double y)
{
    y0[0] = r0;
    y0[1] = theta0;
    y0[2] = 0;
    y0[3] = 0;
    y0[4] = cos(y) * cos(x);
    y0[5] = sin(y) / r0;

    double sintheta, costheta;
    sintheta = sin(theta0);
    costheta = cos(theta0);
    double costheta2 = costheta * costheta;
    double sintheta2 = sintheta * sintheta;

    double rdot0 = y0[4];
    double thetadot0 = y0[5];

    double r2 = r0 * r0;
    double sigma = r2 + bh.a2 * costheta2;
    double delta = r2 - 2.0 * r0 + bh.a2;
    double s1 = sigma - 2.0 * r0;

    y0[4] = rdot0 * sigma / delta;
    y0[5] = thetadot0 * sigma;

    ydot0[0] = rdot0;
    ydot0[1] = thetadot0;
    ydot0[2] = cos(y) * sin(x) / (r0 * sin(theta0));

    double phidot0 = ydot0[2];
    double energy2 = s1 * (rdot0 * rdot0 / delta + thetadot0 * thetadot0) + delta * sintheta2 * phidot0 * phidot0;

    double energy = sqrt(energy2);

    /* Energie de 1 */
    y0[4] = y0[4] / energy;
    y0[5] = y0[5] / energy;

    /* Angular Momentum with E = 1 */
    *L = ((sigma * delta * phidot0 - 2.0 * bh.a * r0 * energy) * sintheta2 / s1) / energy;

    *kappa = y0[5] * y0[5] + bh.a2 * sintheta2 + (*L) * (*L) / sintheta2;

    /* Hack - make sure everything is normalized correctly by a call to geodesic */

    geodesic(*L, *kappa, y0, ydot0);
}

//Obtenir la direction du rayon à partir des coord de BL (y) et de leurs dérivées (dydx)
//Retourne un vecteur normé
void getDirection(double *y, double *dydx, double *xp, double *yp, double *zp)
{

    double costheta = cos(y[1]);
    double sintheta = sin(y[1]);
    double cosphi = cos(y[2]);
    double sinphi = sin(y[2]);
    double r2 = y[0] * y[0];
    double R2 = r2 + bh.a2;
    double R = sqrt(r2 + bh.a2);

    *xp = R * dydx[1] * cosphi * costheta - R * dydx[2] * sinphi * sintheta + sintheta * cosphi * y[0] * dydx[0] / R;
    *yp = R * dydx[1] * sinphi * costheta + R * dydx[2] * cosphi * sintheta + sintheta * sinphi * y[0] * dydx[0] / R;
    *zp = -y[0] * sintheta * dydx[1] + dydx[0] * costheta;
    normalise(xp, yp, zp);
}

//Simulation complete d'un rayon (trajectoire+collision)
void sim(double r0, double theta0, double phi0, double xpixel, double ypixel, double *xp, double *yp, double *zp, double *pixel_transpr, double *pixel_transpg, double *pixel_transpb, double *canalAlpha, bool *ReachedInfinity, int *nbCollision, bool mode)
{

    int N = 0; //Nombre d'itérations
    int k = 0; //nombre de collision

    double oldtheta; //pour detecter le passage à travers le plan z=0 (pour tracer le disque)

    bool zSignChange, diskDistance, diskCollision;

    //Euler

    //double y[Neq];
    //double dydx[Neq];

    //RK4
    double y[Neq];
    double ak[Neq];
    double dydx1[Neq], dydx2[Neq], dydx3[Neq], dydx4[Neq];
    double ytemp[Neq];

    double L, kappa;

    initial(r0, theta0, &L, &kappa, y, ak, xpixel, ypixel);

    double currentStep = rdr.step(y[0], y[1], mode);

    int l = 0; //Juste pour effectuer des boucles sans redeclarer de var

    while ((N < MAXITER) && (rdr.R_min < y[0]) && (y[0] < rdr.R_inf))
    {

        N += 1;

        oldtheta = y[1];

        currentStep = rdr.step(y[0], y[1], mode);

        //Euler
        /*for (int l = 0; l < Neq; l++)
	    {
		    double hdydx = currentStep * dydx[l];
		    y[l] = y[l] + hdydx;
	    }

	    geodesic(L,kappa,y, dydx);*/

        //RK4
        geodesic(L, kappa, y, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx1[l] = ak[l];
            ytemp[l] = y[l] + .5 * currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx2[l] = ak[l];
            ytemp[l] = y[l] + .5 * currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx3[l] = ak[l];
            ytemp[l] = y[l] + currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx4[l] = ak[l];
            dydx1[l] = currentStep / 6. * (dydx1[l] + 2. * dydx2[l] + 2. * dydx3[l] + dydx4[l]); //Recuperer une bonne approx de la derivée dans dydx1
            y[l] = y[l] + dydx1[l];
        }

        zSignChange = (oldtheta > M_PI_2) != (y[1] > M_PI_2);                    //on traverse le plan z=0
        diskDistance = (y[0] < disk.R_max_margin) && (y[0] > disk.R_min_margin); //On l'a traversé la ou est le disque
        diskCollision = zSignChange && diskDistance;

        if (diskCollision)
        {
            double xpos, ypos, zpos;
            blToCartesian(y[0], y[1], y[2], &xpos, &ypos, &zpos); //Obtenir la position en coord cartesiennes
            getDirection(y, dydx1, xp, yp, zp);                   //Obtenir la direction (en coord cartesiennes)

            //Trouver le point de colision (en allant tout droit entre les deux point au dessus et en dessous du disque (valable si le pas est assez petit))
            double lambda = -zpos / (*zp);
            double coll_x = xpos + lambda * (*xp);
            double coll_y = ypos + lambda * (*yp);
            double coll_z = zpos + lambda * (*zp);
            double r = norm(coll_x, coll_y, coll_z);

            diskCollision = (r < disk.R_max) && (r > disk.R_min); //#reverification plus précise
            if (diskCollision)
            {
                if (k < rdr.maxtransparency)
                {

                    double phi = atan2(coll_y / r, coll_x / r) + scn.camera.phi; //Coordonné du point d'impact (en sphérique)

#if ADISK_NORMAL == 1

                    double diskspeed_x, diskspeed_y, diskspeed_z;

                    //rotation prograde du disque

                    sphericalToCartesian(M_PI / 2., phi + M_PI / 2. - scn.camera.phi, &diskspeed_x, &diskspeed_y, &diskspeed_z); //Direction de la vitesse en ce point, sens trigo,mvt circulaire
                    //Prendre en compte le phi camera
                    double beta = disk.RotationSpeed(r);                                                         //Obtenir la norme de la vitesse des poussières en ce point
                    double costhetadoppler = -((*xp) * diskspeed_x + (*yp) * diskspeed_y + (*zp) * diskspeed_z); //Obtenir le cosinus de l'angle entre le rayon et la vitesse des particules (Les vecteurs sont normés)
                    //- à cause du sens des rayons

                    //Temperature du disque à cet endroit
                    double Temp = disk.Temp(r);

                    //Décalage en fréquence égal à décalage en temperature, valable aussi pour l'intensité

                    Temp = Temp * (1 + beta * costhetadoppler) / sqrt(1 - beta * beta); //Effet doppler relativiste

                    //calcul du redshift gravitationnel
                    double gtt = 1. - 1. / r;
                    double gtphi = bh.a / r;
                    double gphiphi = -(r * r + bh.a2 + bh.a2 / r);
                    double omega = beta / r;                                               //Vitesse angulaire
                    Temp = Temp * sqrt(gtt + 2 * gtphi * omega + gphiphi * omega * omega); //Redshift gravitationnel

                    //Recuperer la texture du disque à cet endroit, qu'on utilise uniquemenet pour obtenir la transparence du disque
                    int cy = int((r - disk.R_min) * (adisk_height - 1) / (disk.R_max - disk.R_min));
                    int cx = int((adisk_width - 1) * mod(disk.texture_rep * phi - M_PI, 2. * M_PI) / (2. * M_PI));

                    int loc = (cy * adisk_width + cx) * CHANNEL_NUM;

                    pixel_transpr[k] = 0.;
                    pixel_transpg[k] = 0.;
                    pixel_transpb[k] = 0.;

                    getBodyColor(&pixel_transpr[k], &pixel_transpg[k], &pixel_transpb[k], Temp, 1.); //Obtenir la couleur

                    canalAlpha[k] = (double)adisk[loc] / 255.; //Choix: on utilise la composante verte (le +1) pour reconstituer la transparence

#endif

#if ADISK_GRID == 1

                    getDiskColorGrid(phi, r, &pixel_transpr[k], &pixel_transpg[k], &pixel_transpb[k]); //Motif de grille sur le disque et calcul de la transparence
                    canalAlpha[0] = 1.;
#endif

                    k++;
                }
            }
        }
    }

    *ReachedInfinity = (y[0] >= rdr.R_inf);
    *nbCollision = k;

    if (*ReachedInfinity)
    {
        //Pour visualiser quels rayons ont atteints R_inf
        //*nbCollision=1;
        //pixel_transpr[0]=255;
        //canalAlpha[0]=.5;

        getDirection(y, dydx1, xp, yp, zp);
        normalise(xp, yp, zp); //Ne pas effectuer inutilement cette opération
    }
}

//Simulation d'un rayon (trajectoire uniquement)
void sim_opt(double r0, double theta0, double phi0, double xpixel, double ypixel, double *xp, double *yp, double *zp, bool *ReachedInfinity, bool mode)
{

    int N = 0;

    //Euler
    //double y[Neq];
    //double dydx[Neq];

    //RK4
    double y[Neq], ak[Neq];
    double dydx1[Neq], dydx2[Neq], dydx3[Neq], dydx4[Neq];
    double ytemp[Neq];

    double L, kappa;

    initial(r0, theta0, &L, &kappa, y, ak, xpixel, ypixel); //CI

    double currentStep = rdr.step(y[0], y[1], mode);

    int l;

    while ((N < MAXITER) && (rdr.R_min < y[0]) && (y[0] < rdr.R_inf))
    {

        N += 1;

        currentStep = rdr.step(y[0], y[1], mode);

        //Euler
        /*for (int l = 0; l < Neq; l++)
	    {
		    double hdydx = currentStep * dydx[l];
		    y[l] = y[l] + hdydx;
	    }
	    geodesic(L,kappa,y, dydx);*/

        //RK4

        geodesic(L, kappa, y, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx1[l] = ak[l];
            ytemp[l] = y[l] + .5 * currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx2[l] = ak[l];
            ytemp[l] = y[l] + .5 * currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx3[l] = ak[l];
            ytemp[l] = y[l] + currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx4[l] = ak[l];
            dydx1[l] = currentStep / 6. * (dydx1[l] + 2. * dydx2[l] + 2. * dydx3[l] + dydx4[l]); //Recuperer une bonne aprrox de la derivée
            y[l] = y[l] + dydx1[l];
        }
    }

    *ReachedInfinity = (y[0] >= rdr.R_inf);

    if (*ReachedInfinity)
    { //Ne pas effectuer inutilement cette opération
        getDirection(y, dydx1, xp, yp, zp);
        normalise(xp, yp, zp);
    }
}

//Simulation d'un faisceu de rayons (uniquement dans le cas ou on dessine des étoiles)
void sim_bundle(int i0, int j0, double x, double y, double z, double *pixelr, double *pixelg, double *pixelb)
{ //arguments: i0 et j0 position du pixel, x,y,z position initiale du rayon dans l'espace

    double r0, theta0, phi0;
    cartesianToBl(x, y, z, &r0, &theta0, &phi0); //Approximativement (compliqué de passer de cartesien à BL), mais pas important

    double xpcentre0, ypcentre0, zpcentre0;
    double xpcentref, ypcentref, zpcentref;
    double xp0, yp0, zp0;
    double xpf, ypf, zpf;

    bool ReachedInfinity = true;
    int nbCollision = 0;

    //double mindist2=1000.;
    //double mindist2_0=1000.;

    /*double Distances2[4];
    double Distances2_0[4];
    double dilatation=1.;*/

    double pixel_transpr[rdr.maxtransparency], pixel_transpg[rdr.maxtransparency], pixel_transpb[rdr.maxtransparency];
    double canalAlpha[rdr.maxtransparency];

    double range = 0.075 * 20. / (rdr.width - 1.0); //FOV

    //for (int v = 0; v < rdr.SamplesPerPixels; v++)
    //{
    double tmpR = 0;
    double tmpG = 0;
    double tmpB = 0;

    double i = i0 + doubleRand(1.) - 0.5;
    double j = j0 + doubleRand(1.) - 0.5;

    double maxdist2 = 0.;
    double maxdist2_0 = 0.;

    double xpixel = -(i - (rdr.width + 1.0) / 2) * range;
    double ypixel = -(j - (rdr.height + 1.0) / 2) * range;

    bool mode = (abs(xpixel) < 1e-2); //Si le rayon se dirige vers les poles

    sim(r0, theta0, phi0, xpixel, ypixel, &xpcentref, &ypcentref, &zpcentref, pixel_transpr, pixel_transpg, pixel_transpb, canalAlpha, &ReachedInfinity, &nbCollision, mode); //Simulation du rayon principal

#if DRAWSTARS == 1

    int Nray = 0;

    while (Nray < 4 && ReachedInfinity)
    {

        if (Nray == 0)
        {
            xpixel = -(i + rdr.delta - (rdr.width + 1.0) / 2) * range;
            ypixel = -(j - (rdr.height + 1.0) / 2) * range;
        }
        if (Nray == 1)
        {
            xpixel = -(i - rdr.delta - (rdr.width + 1.0) / 2) * range;
            ypixel = -(j - (rdr.height + 1.0) / 2) * range;
        }
        if (Nray == 2)
        {
            xpixel = -(i - (rdr.width + 1.0) / 2) * range;
            ypixel = -(j + rdr.delta - (rdr.height + 1.0) / 2) * range;
        }
        if (Nray == 3)
        {
            xpixel = -(i - (rdr.width + 1.0) / 2) * range;
            ypixel = -(j - rdr.delta - (rdr.height + 1.0) / 2) * range;
        }

        Nray++;

        xpf = xp0;
        ypf = yp0;
        zpf = zp0;

        mode = (abs(xpixel) < 1e-2);

        sim_opt(r0, theta0, phi0, xpixel, ypixel, &xpf, &ypf, &zpf, &ReachedInfinity, mode); //Simulation des rayons secondaires

        if (ReachedInfinity)
        {
            //Calcule de la taille de la zone d'impact des differents rayons (zone en terme de direction)
            //On trouve le rayon maxdist du cercle centré sur la direction du rayon principal, qui contient les autres rayons (tjrs en terme de direction)
            //Calcul en coord cartesienne pour éviter les discontinuités liées au modulo 2PI

            //double dx0=xpcentre0-xp0;
            //double dy0=ypcentre0-yp0;
            //double dz0=zpcentre0-zp0;

            //double dist2_0=sqrnorm(dx0,dy0,dz0);

            double dx = xpcentref - xpf;
            double dy = ypcentref - ypf;
            double dz = zpcentref - zpf;

            double dist2 = sqrnorm(dx, dy, dz);

            //if (maxdist2_0<dist2_0){ maxdist2_0=dist2_0; }
            if (maxdist2 < dist2)
            {
                maxdist2 = dist2;
            }
        }
    }

    if (ReachedInfinity)
    { //Si tous les rayons se sont échappés, on regarde quelles étoiles sont dans le cercle formé par le faisceau sur le ciel

        double sqrdistanceToStar;

        double theta, phi;
        cartesianToSpherical(xpcentref, ypcentref, zpcentref, &theta, &phi);

        int loc = scn.compute_hash(theta, phi + scn.camera.phi);

        for (int t = 0; t < scn.hashindex[loc]; ++t)
        {
            int index = scn.hashtable[loc][t];

            sqrdistanceToStar = sqrnorm(scn.starx[index] - xpcentref, scn.stary[index] - ypcentref, scn.starz[index] - zpcentref);
            if (sqrdistanceToStar < maxdist2)
            {
                getBodyColor(&tmpR, &tmpG, &tmpB, scn.starTemp[index], scn.starBrightness[index] * 1.); //(Récupere la couleur et effectue l'addition composante par composante)
            }
        }

        double f = 5e-8 / maxdist2; //Angle solide initial environ constant égal à 5e-8
        tmpR *= f;
        tmpG *= f;
        tmpB *= f;
    }
#endif
#if DRAWGRID == 0

    getFinalColor(pixelr, pixelg, pixelb, tmpR, tmpG, tmpB, pixel_transpr, pixel_transpg, pixel_transpb, canalAlpha, nbCollision); // Calcul de transparence

#else
    getFinalColorGrid(xpcentref, ypcentref, zpcentref, scn.camera.phi, pixelr, pixelg, pixelb, pixel_transpr, pixel_transpg, pixel_transpb, nbCollision, ReachedInfinity); // Calcul de la transparence, et calcul de la couleur de la grille
#endif
    //}
    //}

    /* *pixelr *= rdr.invSamplesPerPixels;
    *pixelg *= rdr.invSamplesPerPixels;
    *pixelb *= rdr.invSamplesPerPixels;*/
}

int exp(int a, int b)
{
    int res = 1;
    for (int i = 0; i < b; ++i)
    {
        res *= a;
    }
    return res;
}

void computeChunk(int h, double *image, int *index)
{
    int index0 = h * rdr.chunkSize;
    int index1 = index0 + rdr.chunkSize;

    for (int l = index0; l < index1; ++l)
    {
        double pixelR = 0;
        double pixelG = 0;
        double pixelB = 0;

        int i = index[l] % rdr.width;
        int j = index[l] / rdr.width;

        sim_bundle(i, j, scn.camera.x, scn.camera.y, scn.camera.z, &pixelR, &pixelG, &pixelB); //i horiz, j vert

        int loc = index[l] * CHANNEL_NUM;

        image[loc] += pixelR * rdr.invSamplesPerPixels;
        image[loc + 1] += pixelG * rdr.invSamplesPerPixels;
        image[loc + 2] += pixelB * rdr.invSamplesPerPixels;
    }
}

void displayChunks(int h, double *image, int *index, int pass)
{
    int index0 = h * rdr.chunkSize;
    int index1 = index0 + rdr.chunkSize;

    double f = rdr.SamplesPerPixels / pass;

    if (pass == 1)
    {

        for (int l = index0; l < index1; ++l)
        {

            int i = index[l] % rdr.width;
            int j = index[l] / rdr.width;

            int loc = index[l] * CHANNEL_NUM;

            unsigned long rgb = createRGB(int(image[loc] * f), int(image[loc + 1] * f), int(image[loc + 2] * f));

            int div = exp(2, reconstructionMatrix[index[l]]);

            static const int bigCubeSize = 32; //Doit être une puissance de 2

            int smallCubeSize = bigCubeSize / div;

            int w0 = bigCubeSize * int(i / bigCubeSize);
            int h0 = bigCubeSize * int(j / bigCubeSize);

            w0 = w0 + int((i - w0) / smallCubeSize) * smallCubeSize;
            h0 = h0 + int((j - h0) / smallCubeSize) * smallCubeSize;

            int w1 = w0 + smallCubeSize;
            int h1 = h0 + smallCubeSize;

            if (w1 > rdr.width)
            {
                w1 = rdr.width;
            }
            if (h1 > rdr.height)
            {
                h1 = rdr.height;
            }

            XSetForeground(disp, graphcontext, rgb);
            XFillRectangle(disp, win, graphcontext, w0, h0, w1 - w0, h1 - h0);

            for (int w = w0; w < w1; w++)
            {
                int location = h * rdr.width;
                for (int h = h0; h < h1; h++)
                {
                    int location = h * rdr.width + w;

                    reconstructionMatrix[location] += 1;
                }
            }
        }
    }
    else
    {
        for (int l = index0; l < index1; ++l)
        {

            int i = index[l] % rdr.width;
            int j = index[l] / rdr.width;

            int loc = index[l] * CHANNEL_NUM;

            unsigned long rgb = createRGB(int(image[loc] * f), int(image[loc + 1] * f), int(image[loc + 2] * f));
            drawpixel(disp, win, graphcontext, i, j, rgb);
        }
    }
    XFlush(disp);
}

void estimateComputationTime(double *image, int *index)
{

    cout << "Estimation du temps de calcul... " << endl;
    auto startTest = high_resolution_clock::now();

    computeChunk(0, image, index);
    auto stopTest = high_resolution_clock::now();
    auto testTime = duration_cast<milliseconds>(stopTest - startTest);

    unsigned nb_threads_hint = thread::hardware_concurrency();
    unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);
    cout << "Temps estimé: " << int(testTime.count() * rdr.TotalChunknumber / (nb_threads * 1000)) << " secondes" << endl;
}

void convolution_2D(double *image, double *kernel, double *resultat, int width, int height, int kwidth, int kheight)
{

    // find center position of kernel (half of kernel size)
    int kCenterX = kwidth / 2;
    int kCenterY = kheight / 2;

    for (int i = 0; i < height; ++i) // rows
    {
        for (int j = 0; j < width; ++j) // columns
        {
            int loc1 = (i * width + j) * CHANNEL_NUM;

            for (int m = 0; m < kheight; ++m) // kernel rows
            {
                for (int n = 0; n < kwidth; ++n) // kernel columns
                {

                    int ii = i + (m - kCenterY);
                    int jj = j + (n - kCenterX);

                    if (ii >= 0 && ii < height && jj >= 0 && jj < width)
                    {
                        int loc2 = (ii * width + jj) * CHANNEL_NUM;
                        int loc3 = (m * kwidth + n) * CHANNEL_NUM;

                        resultat[loc1] += image[loc2] * kernel[loc3];
                        resultat[loc1 + 1] += image[loc2 + 1] * kernel[loc3 + 1];
                        resultat[loc1 + 2] += image[loc2 + 2] * kernel[loc3 + 2];
                    }
                }
            }
        }
    }
}

void postprocess(char *name, double *image)
{
    cout << "Postprocessing..." << endl;

    int Npixels = kernel_width * kernel_height;

    int size = Npixels * CHANNEL_NUM;

    double *kernel_double = new double[kernel_width * kernel_height * CHANNEL_NUM];

    double sum = 0;

    for (int k = 0; k < size; ++k)
    {
        sum += kernel[k];
    }
    sum = sum / (3); //Diviser par le nombre de composantes

    for (int k = 0; k < size; ++k)
    {
        kernel_double[k] = (double)kernel[k] / sum; //Normaliser
    }

    double *image_convolve = new double[rdr.Npixels * CHANNEL_NUM];

    int Ncomp = CHANNEL_NUM * rdr.Npixels;
    for (int l = 0; l < Ncomp; ++l)
    {
        image_convolve[l] = 0.; //Mettre l'image à zero
    }

    convolution_2D(image, kernel_double, image_convolve, rdr.width, rdr.height, kernel_width, kernel_height);

    char nametmp[50];
    int n;
    n = sprintf(nametmp, "%s.png", name);
    save(nametmp, image_convolve, rdr.height, rdr.width, CHANNEL_NUM);
}

void render(char *name, double *image)
{

    int *index = new int[rdr.Npixels];

    for (int l = 0; l < rdr.Npixels; ++l)
    {
        reconstructionMatrix[l] = 0;
        index[l] = l;
    }

    int Ncomp = CHANNEL_NUM * rdr.Npixels;
    for (int l = 0; l < Ncomp; ++l)
    {
        image[l] = 0.; //Mettre l'image à zero
    }

    cout << "Mélange des pixels" << endl;
    for (int i = rdr.Npixels - 1; i != 0; i--)
    {
        int j = rand() % i;
        swap(index[i], index[j]);
    }

    estimateComputationTime(image, index);

    //auto start = high_resolution_clock::now();

    //Select what events the window will listen to
    //XSelectInput(disp, win, KeyPressMask | PointerMotionMask);
    //XEvent ev;
    int quit = 0;

    int c = 0;
    int Chunkcomputed = 0;

    int x, y, prevx, prevy, winx, winy;
    long unsigned int child;
    unsigned int mask;

    int pass = 1;

    int bloc = 8;

    while (!quit)
    {

        /*XNextEvent( disp, &ev );
        switch( ev.type ) {
            case MotionNotify:
                printf("x %d y %d\n", ev.xmotion.x, ev.xmotion.y );
                break;
        }*/

        XQueryPointer(disp, win, &root, &child, &x, &y, &winx, &winy, &mask);

        if (x != prevx || y != prevy)
        {

            prevx = x;
            prevy = y;

            cartesianToBl(scn.camera.x, scn.camera.y, scn.camera.z, &scn.camera.r, &scn.camera.theta, &scn.camera.phi);

            scn.camera.theta = M_PI * y / rdr.height;
            scn.camera.phi = 2 * M_PI * x / rdr.width;

            sphericalToCartesianNotNormalized(scn.camera.r, scn.camera.theta, scn.camera.phi, &scn.camera.x, &scn.camera.y, &scn.camera.z);

            pass = 1;
            Chunkcomputed = 0;
            bloc = 8;

            for (int l = 0; l < rdr.Npixels; ++l)
            {
                reconstructionMatrix[l] = 0;
            }

            for (int l = 0; l < Ncomp; ++l)
            {
                image[l] = 0.; //Mettre l'image à zero
            }
        }

        if (Chunkcomputed < rdr.TotalChunknumber)
        {

            int numberofblocs = bloc;
            int rest = rdr.TotalChunknumber - Chunkcomputed;
            if (rest < bloc)
            {
                numberofblocs = rest;
            }

            pl::async_par_for(0, numberofblocs, [&](unsigned h) {
                //for (int h = 0; h < numberofblocs; ++h)
                //{

                double pc = 100. * (double)Chunkcomputed / ((double)rdr.TotalChunknumber);
                int currentChunk = Chunkcomputed + h;
                printf("\rFile %s - Pass %d - Chunk %04d/%04d started - %2.2f pourcents - r: %02.2f, theta: %.2f,phi: %.2f", name, pass, currentChunk, rdr.TotalChunknumber, pc, scn.camera.r, scn.camera.theta, scn.camera.phi);

                fflush(stdout);
                computeChunk(currentChunk, image, index);
            });
            //}

            for (int v = 0; v < numberofblocs; ++v)
            {
                displayChunks(Chunkcomputed + v, image, index, pass);
            }

            Chunkcomputed += numberofblocs;
        }
        else if (pass < rdr.SamplesPerPixels)
        {
            pass++;
            Chunkcomputed = 0;
            bloc = 4 * 32;//Plus gors bloc, moins de rafraichissements
        }
        else
        {
            printf("\33[2K\rFin");
        }
    }

    /*if (ev.type == KeyPress)
        {
            printf("KeyPress: %x\n", ev.xkey.keycode);
            

            
            if (ev.xkey.keycode == 0x09)
                quit = 1;
        }
        if (ev.type == Expose)
        {
        }*/


    /*auto stop = high_resolution_clock::now();
    auto totalTime = duration_cast<seconds>(stop - start);

    printf("\33[2K\rTemps de calcul : %ld secondes\n", totalTime.count());

    char nametmp[50];
    int n;
    n = sprintf(nametmp, "%s_raw.png", name);
    save(nametmp, image, rdr.height, rdr.width, CHANNEL_NUM);*/
}

void readParams(char *file)
{
    ifstream monFlux;
    monFlux.open(file);

    if (monFlux)
    {
        monFlux >> rdr.width >> rdr.height >> rdr.chunkSize >> rdr.stepmax >> rdr.stepmin0 >> rdr.delta >> rdr.SamplesPerPixels;
        cout << "Paramètres de rendu trouvés" << endl;
    }
    else
    {
        cout << "Paramètres de rendu non trouvés" << endl;
        exit(0);
    }
}

void readScene(char *file)
{
    ifstream monFlux;
    monFlux.open(file);

    if (monFlux)
    {
        monFlux >> bh.a >> scn.camera.x >> scn.camera.y >> scn.camera.z >> disk.R_max >> disk.betamax >> disk.TMax >> disk.texture_rep >> rdr.R_inf;
        cout << "Paramètres de la scène trouvés" << endl;
    }
    else
    {
        cout << "Paramètres de la scène non trouvés" << endl;
        exit(0);
    }
}

int main(int argc, char *argv[])
{
    if (argc == 7)
    {
        adisk = stbi_load(argv[1], &adisk_width, &adisk_height, &adisk_bpp, 3);
        if (adisk == NULL)
        {
            cout << "Texture du disque non trouvée" << endl;
            exit(0);
        }
        else
        {
            cout << "Texture du disque trouvée" << endl;
        }
        readSensitivityData(argv[2], wavelengthSamples, wavelengthSamples5, sensitivitySamplesR, sensitivitySamplesG, sensitivitySamplesB);

        kernel = stbi_load(argv[3], &kernel_width, &kernel_height, &kernel_bpp, 3);
        if (kernel == NULL)
        {
            cout << "PSF non trouvée" << endl;
            exit(0);
        }
        else
        {
            cout << "PSF trouvée" << endl;
        }

        readScene(argv[4]);
        readParams(argv[5]);

        /*
            rdr.height: Taille du rendu
            rdr.width: Taille du rendu
            rdr.R_inf: Distance à partir de laquelle on considere etre a l'infini 
            rdr.chunkSize: nombre de pixels par blocs traités en parallele

            rdr.stepmax: Pas max
            rdr.stepmin: Pas min

            rdr.delta: Ecart angulaire (en pixel) entre les rayons d'un même faisceau

            bh.a: Parametre de Kerr du trou noir (entre 0 et 1)

            disk.R_max: Rayon maximal du disque de poussière
            disk.betamax: Vitesse du disque au plus proche du trou noir (c'est la qu'est la vitesse max pour un profil de vitesse en r^-1/2)
            disk.TMax: Temperature du disque au plus proche du trou noir

            disk.texture_rep: Nombre de répétition de la texture du dique en longeur pour ne pas qu'elle soit pixelisée

            scn.FOV: Pas utilisé
            */

        scn.precalc();

        bh.precalc(); //A executer avant inner orbit

        disk.R_min = bh.inner_orbit(); //Trouver l'orbite la plus proche du trou noir encore stable (prograde)

        rdr.precalc(bh.a2); //quelques calculs pour avoir les carrés de certaines qtité et le fov en radian etc..
        disk.precalc();

        scn.generateSky();

        cout << "Prêt" << endl;

        double *image = new double[rdr.Npixels * CHANNEL_NUM]; //Tableau qui va contenir les données

        imageToDisplay = new unsigned long[rdr.Npixels];
        reconstructionMatrix = new int[rdr.Npixels];


        //initScreen(disp, rdr.width,rdr.height,&win, &root,&graphcontext);
        //cout<<win<<","<<root<<","<<graphcontext<<endl;

        //Open Display
        disp = XOpenDisplay(getenv("DISPLAY"));
        if (disp == NULL)
        {
            printf("Couldn't open display.\n");
            //return -1;
        }

        //Create
        int const x = 0, y = 0, border_width = 1;
        int sc = DefaultScreen(disp);
        root = DefaultRootWindow(disp);
        win = XCreateSimpleWindow(disp, root, x, y, rdr.width, rdr.height, border_width, BlackPixel(disp, sc), WhitePixel(disp, sc));
        XMapWindow(disp, win);           //Make window visible
        XStoreName(disp, win, "Render"); // Set window title

        //Prepare the window for drawing
        graphcontext = XCreateGC(disp, root, 0, NULL);

        XClearWindow(disp, win);

        //postprocess(argv[6], image);

        cout << "Ok" << endl;

        render(argv[6], image);

        stbi_image_free(adisk);

        XFreeGC(disp, graphcontext);
        XDestroyWindow(disp, win);
        XCloseDisplay(disp);

        cout << "Fin" << endl;
    }
    else
    {
        cout << "Pas assez d'arguments" << endl;
    }

    return 0;
}
