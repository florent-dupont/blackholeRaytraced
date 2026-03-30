#include <stdint.h>
#include <iostream>
#include <math.h>
#include <tgmath.h>





void normalisepixel(double *rgbR, double *rgbG, double *rgbB)
{
    double maxcomp = *rgbR;
    if (*rgbG > maxcomp)
    {
        maxcomp = *rgbG;
    }
    if (*rgbB > maxcomp)
    {
        maxcomp = *rgbB;
    }

    if (maxcomp > 255.)
    {
        *rgbR = 255. * (*rgbR) / maxcomp;
        *rgbG = 255. * (*rgbG) / maxcomp;
        *rgbB = 255. * (*rgbB) / maxcomp;
    }
}


void save(char *name, double *image, int height, int width, int channel)
{

    int Npixels = height * width;

    uint8_t *image_byte = new uint8_t[Npixels * channel]; //Image finale

    int loc = 0;

    for (int k = 0; k < Npixels; ++k)
    {
        double r=image[loc];
        double g=image[loc+1];
        double b=image[loc+2];

        normalisepixel(&r,&g,&b);

        image_byte[loc] = char(r);
        loc++;

        image_byte[loc] = char(g);
        loc++;

        image_byte[loc] = char(b);
        loc++;
    }

    stbi_write_png(name, width, height, channel, image_byte, width * channel);
}

double doubleRand(double MaxVal)
{
    return ((double)rand() / (double)RAND_MAX) * MaxVal;
}

void cartesianToSpherical(double x, double y, double z, double *theta, double *phi)
{
    *theta = acos(z);
    *phi = atan2(y, x);
}

void sphericalToCartesian(double theta, double phi, double *x, double *y, double *z)
{
    double tmp = sin(theta);
    *x = tmp * cos(phi);
    *y = tmp * sin(phi);
    *z = cos(theta);
}

void sphericalToCartesianNotNormalized(double r, double theta, double phi, double *x, double *y, double *z)
{
    double tmp = r*sin(theta);
    *x = tmp * cos(phi);
    *y = tmp * sin(phi);
    *z = r*cos(theta);
}

double mod(double x, double y)
{
    int resultat = floor(x / y);
    return x - (resultat * y);
}

double sqrnorm(double x, double y, double z)
{
    return x * x + y * y + z * z;
}
double norm(double x, double y, double z)
{
    return sqrt(x * x + y * y + z * z);
}

void normalise(double *x, double *y, double *z)
{
    double tmp = sqrt(sqrnorm(*x, *y, *z));
    *x /= tmp;
    *y /= tmp;
    *z /= tmp;
}

double avg(double *array, int size)
{
    double res = 0;
    for (int i = 0; i < size; ++i)
    {
        res += array[i];
    }
    res = res / (double)size;
    return res;
}


//Utilisé pour ajouter à l'image du fond celest le disque de poussière en superposant dans l'ordre les masques
void getFinalColor(double *out_rgbR, double *out_rgbG, double *out_rgbB,double in_rgbR, double in_rgbG, double in_rgbB, double *pixel_transpr, double *pixel_transpg, double *pixel_transpb, double *canalAlpha, int nbCollision)
{

    //normalisepixel(rgbR, rgbG, rgbB);
 

    for (int l = nbCollision - 1; l >= 0; --l)
    {
        normalisepixel(&pixel_transpr[l], &pixel_transpg[l], &pixel_transpb[l]);

        in_rgbR = canalAlpha[l] * pixel_transpr[l] + (1 - canalAlpha[l]) * in_rgbR;
        in_rgbG = canalAlpha[l] * pixel_transpg[l] + (1 - canalAlpha[l]) * in_rgbG;
        in_rgbB = canalAlpha[l] * pixel_transpb[l] + (1 - canalAlpha[l]) * in_rgbB;
    }
    *out_rgbR+=in_rgbR;
    *out_rgbG+=in_rgbG;
    *out_rgbB+=in_rgbB;
}

void getFinalColorGrid(double xpcentref, double ypcentref, double zpcentref, double cameraPhi,double *rgbR, double *rgbG, double *rgbB, double *pixel_transpr, double *pixel_transpg, double *pixel_transpb, int nbCollision, bool ReachedInfinity)
{
    double theta, phi;
    cartesianToSpherical(xpcentref, ypcentref, zpcentref, &theta, &phi);
    phi = (phi + cameraPhi + M_PI) / (2. * M_PI);
    theta /= M_PI;

    bool a = int((100 * theta)) % 2; //100 changements par tour
    bool b = int((100 * phi)) % 2;   //100 changements par tour

#if ADISK_GRID == 1

    if (nbCollision > 0)
    {
        *rgbR += pixel_transpr[0];
        *rgbG += pixel_transpg[0];
        *rgbB += pixel_transpb[0];
    }
    else if (ReachedInfinity)
    {
        *rgbR += b * 200;
        *rgbG += 0.;
        *rgbB += a * 200;
    }
#else

    if (ReachedInfinity)
    {
        *rgbR += b * 200;
        *rgbG += 0.;
        *rgbB += a * 200;
    }
#endif
}
