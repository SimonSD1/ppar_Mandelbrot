/*
 * Sorbonne Université  --  PPAR
 * Computing the Mandelbrot set, sequential version
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h> /* compile with -lm */
#include <sys/time.h>
#include <arpa/inet.h> /* htonl */
#include <mpi.h>

#include "rasterfile.h"

char info[] = "\
Usage:\n\
      ./mandel dimx dimy xmin ymin xmax ymax depth\n\
\n\
      dimx,dimy : dimensions of the image to generate\n\
      xmin,ymin,xmax,ymax : domain to be computed, in the complex plane\n\
      depth : maximal number of iterations\n\
\n\
Some examples of execution:\n\
      ./mandel 3840 3840 0.35 0.355 0.353 0.358 200\n\
      ./mandel 3840 3840 -0.736 -0.184 -0.735 -0.183 500\n\
      ./mandel 3840 3840 -1.48478 0.00006 -1.48440 0.00044 100\n\
      ./mandel 3840 3840 -1.5 -0.1 -1.3 0.1 10000\n\
";

double wallclock_time()
{
    struct timeval tmp_time;
    gettimeofday(&tmp_time, NULL);
    return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6);
}

unsigned char cos_component(int i, double freq)
{
    double iD = i;
    iD = cos(iD / 255.0 * 2 * M_PI * freq);
    iD += 1;
    iD *= 128;
    return iD;
}

/**
 *  Save the data array in rasterfile format
 */
void save_rasterfile(char *name, int largeur, int hauteur, unsigned char *p)
{
    FILE *fd = fopen(name, "w");
    if (fd == NULL)
    {
        perror("Error while opening output file. ");
        exit(1);
    }

    struct rasterfile file;
    file.ras_magic = htonl(RAS_MAGIC);
    file.ras_width = htonl(largeur);            /* width of the image, in pixels */
    file.ras_height = htonl(hauteur);           /* height of the image, in pixels */
    file.ras_depth = htonl(8);                  /* depth of each pixel (1, 8 or 24 ) */
    file.ras_length = htonl(largeur * hauteur); /* size of the image in nb of bytes */
    file.ras_type = htonl(RT_STANDARD);         /* file type */
    file.ras_maptype = htonl(RMT_EQUAL_RGB);
    file.ras_maplength = htonl(256 * 3);
    fwrite(&file, sizeof(struct rasterfile), 1, fd);

    /* Color palette: red component */
    for (int i = 255; i >= 0; i--)
    {
        unsigned char o = cos_component(i, 13.0);
        fwrite(&o, sizeof(unsigned char), 1, fd);
    }

    /* Color palette: green component */
    for (int i = 255; i >= 0; i--)
    {
        unsigned char o = cos_component(i, 5.0);
        fwrite(&o, sizeof(unsigned char), 1, fd);
    }

    /* Color palette: blue component */
    for (int i = 255; i >= 0; i--)
    {
        unsigned char o = cos_component(i + 10, 7.0);
        fwrite(&o, sizeof(unsigned char), 1, fd);
    }

    fwrite(p, largeur * hauteur, sizeof(unsigned char), fd);
    fclose(fd);
}

/**
 * Given the coordinates of a point $c = a + ib$ in the complex plane,
 * the function returns the corresponding color which estimates the
 * distance from the Mandelbrot set to this point.
 * Consider the complex sequence defined by:
 * \begin{align}
 *     z_0     &= 0 \\
 *     z_{n+1} &= z_n^2 + c
 *   \end{align}
 * The number of iterations that this sequence needs before diverging
 * is the number $n$ for which $|z_n| > 2$.
 * This number is brought back to a value between 0 and 255, thus
 * corresponding to a color in the color palette.
 */
unsigned char xy2color(double a, double b, int depth)
{
    double x = 0;
    double y = 0;
    for (int i = 0; i < depth; i++)
    {
        /* save the former value of x (which will be erased) */
        double temp = x;
        /* new values for x and y */
        double x2 = x * x;
        double y2 = y * y;
        x = x2 - y2 + a;
        y = 2 * temp * y + b;
        if (x2 + y2 > 4.0)
            return (i % 255); /* diverges */
    }
    return 255; /* did not diverge in depth steps */
}

/*
 * Main: in each point of the grid, apply xy2color
 */
int main(int argc, char **argv)
{

    int my_rank; /* rank of the process */
    int p;       /* number of processes */

    MPI_Status status;
    /* Initialisation */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (argc == 1)
    {
        fprintf(stderr, "%s\n", info);
        return 1;
    }

    /* Default values for the fractal */
    double xmin = -2; /* Domain of computation, in the complex plane */
    double ymin = -2;
    double xmax = 2;
    double ymax = 2;
    int w = 3840; /* dimensions of the image (4K HTDV!) */
    int h = 3840;
    int depth = 10000; /* depth of iteration */

    /* Retrieving parameters */
    if (argc > 1)
        w = atoi(argv[1]);
    if (argc > 2)
        h = atoi(argv[2]);
    if (argc > 3)
        xmin = atof(argv[3]);
    if (argc > 4)
        ymin = atof(argv[4]);
    if (argc > 5)
        xmax = atof(argv[5]);
    if (argc > 6)
        ymax = atof(argv[6]);
    if (argc > 7)
        depth = atoi(argv[7]);

    /* Computing steps for increments */
    double xinc = (xmax - xmin) / (w - 1);
    double yinc = (ymax - ymin) / (h - 1);

    /* show parameters, for verification */
    // fprintf(stderr, "Domain: [%g,%g] x [%g,%g]\n", xmin, ymin, xmax, ymax);
    // fprintf(stderr, "Increment: %g, %g\n", xinc, yinc);
    // fprintf(stderr, "depth: %d\n", depth);
    // fprintf(stderr, "Image size: %d x %d\n", w, h);

    /* Allocate memory for the output array */
    unsigned char *image = malloc(w * h);
    if (image == NULL)
    {
        perror("Error while allocating memory for array: ");
        exit(1);
    }

    int nbComputed = 0;

    int mon_rang = my_rank;

    int nbTiles = 20;

    int ligne_par_process = h/nbTiles;

    //  numero de travail -1
    int noMoreWork = 0;

    // la premiere case contient le numero du travail
    unsigned char *localBuff = malloc(ligne_par_process*w * sizeof(unsigned char) + 1);
    unsigned char *localBuffRecepetion = malloc(ligne_par_process*w * sizeof(unsigned char) + 1);

    // tag utilisés
    // demande = 0
    // resultat = 1
    // noMoreWork =3

    int demande = 0;

    int demandeur_travail;

    unsigned char numeroTravail;

    MPI_Request request;

    while (nbComputed < nbTiles && !noMoreWork)
    {
        if (mon_rang != 0)
        {
            // demande travail a 0
            MPI_Isend(&mon_rang, 1, MPI_INT32_T, 0, 0, MPI_COMM_WORLD, &request);
            //printf(" %d demande un tarvail\n",mon_rang);
            // recupere le numero du travail

            MPI_Irecv(&numeroTravail, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);

            //printf("%d a recu le travail %d\n",mon_rang,numeroTravail);

            if (status.MPI_TAG == 3)
            {   
                //printf("%d a recu qu'il n'y a plus de travail\n",mon_rang);
                noMoreWork = 1;
                continue;
            }

            // calcul

            //printf("%d calcul\n",mon_rang);

            for(int i=0; i<ligne_par_process*w;i++){

                localBuff[i+1]=mon_rang+1;
            }

            // envoit le resultat a 0
            localBuff[0] = numeroTravail;
            MPI_Isend(localBuff, ligne_par_process*w+ 1, MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD, &request);
            //printf("%d a envoye sont travail a 0\n",mon_rang);
        }

        if (mon_rang == 0)
        {
            // attend demande travail
            MPI_Irecv(&demandeur_travail, 1, MPI_INT32_T, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            //printf("0 a recue une demande de %d\n",demandeur_travail);

            // envoit travail numero nbComputed et nbComputed++ si il y en a
            MPI_Isend(&nbComputed, 1, MPI_INT32_T, demandeur_travail, 0, MPI_COMM_WORLD, &request);
            nbComputed++;
            //printf("0 a envoye le travail %d a %d\n",nbComputed-1,demandeur_travail);
            // recuperere le resultat

            MPI_Irecv(localBuffRecepetion, ligne_par_process*w + 1, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            printf("0 a recue le travail numero %d\n",localBuffRecepetion[0]);
            
            unsigned char numeroTravailRecue = localBuffRecepetion[0];
            // printf("numero travail recue : %d\n",numeroTravailRecue);

            for (int i = 0; i < w * ligne_par_process+1; i++)
            {
                
                printf("%d ",localBuffRecepetion[i]);

                
                

                // printf("copie en %d pour le travail numero %d\n", numeroTravailRecue * (w * h / nbTiles) + i, numeroTravailRecue);
                image[numeroTravailRecue * ligne_par_process*w + i] = localBuffRecepetion[i + 1];
                //printf("%d ",numeroTravailRecue * ligne_par_process*w + i);
            }
            printf("w*ligneparprocess : %d",w*ligne_par_process);
            printf("\n");
        }
    }

    if (mon_rang == 0)
    {
        noMoreWork = 1;
        for (int i = 1; i < p; i++)
        {
            MPI_Isend(&noMoreWork, 0, MPI_INT, i, 3, MPI_COMM_WORLD, &request);
        }
        for (int i = 0; i < h * w; i++)
        {
            //printf("%d ", image[i]);
        }
    }

    // faire un broadcast pour dire que c'est finit

    /* start timer */
    // double start = wallclock_time();

    /* Process the grid point by point */

    /* stop timer */
    // double end = wallclock_time();
    // fprintf(stderr, "Total computing time: %g sec\n", end - start);

    /* Save the image in the output file "mandel.ras" */
    if (my_rank == 0)
    {
        for(int i=0; i<h*w; i++){
            printf("%d ",image[i]);
        }
        save_rasterfile("mandel.ras", w, h, image);
    }

    MPI_Finalize();
    return 0;
}
