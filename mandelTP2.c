#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h> /* compile with -lm */
#include <sys/time.h>
#include <arpa/inet.h> /* htonl */

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

unsigned char cos_component(int i, double freq)
{
    double iD = i;
    iD = cos(iD / 255.0 * 2 * M_PI * freq);
    iD += 1;
    iD *= 128;
    return iD;
}

double wallclock_time()
{
    struct timeval tmp_time;
    gettimeofday(&tmp_time, NULL);
    return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6);
}

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

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    if (argc == 1)
    {
        fprintf(stderr, "%s\n", info);
        return 1;
    }

    int my_rank;
    int nbProcess;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbProcess);

    double xmin = -2; /* Domain of computation, in the complex plane */
    double ymin = -2;
    double xmax = 2;
    double ymax = 2;
    int w = 1000;
    int h = 1000;

    int depth = 100;

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

    // decoupage de l'image
    int nbTiles = h;

    int paddedH = h;

    // on gere le padding
    if (h % nbTiles != 0)
    {
        paddedH = h + (nbTiles - h % nbTiles);
    }

    int linesPerTiles = paddedH / nbTiles;
    int pixelPerTiles = linesPerTiles * w;

    int tilesSent = 0;
    int tilesReceive = 0;

    double xinc = (xmax - xmin) / (w - 1);
    double yinc = (ymax - ymin) / (h - 1);

    int stop = -1;

    // numero de ligne dans tag

    if (my_rank == 0)
    {
        MPI_Status status;
        MPI_Request request;

        /* start timer */
        double start = wallclock_time();
        unsigned char *image = malloc(paddedH * w * sizeof(unsigned char));
        unsigned char *tempImage = malloc(pixelPerTiles * sizeof(unsigned char));

        // on envoi a tout le monde au debut
        for (int process = 1; process < nbProcess; process++)
        {
            if (tilesSent < nbTiles)
            {
                MPI_Isend(&tilesSent, 1, MPI_INT, process, 0, MPI_COMM_WORLD, &request);
                tilesSent++;
            }
        }

        // on continue a envoyer le reste quand on a des reponses
        while (tilesReceive < nbTiles)
        {
            MPI_Irecv(tempImage, pixelPerTiles, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);

            MPI_Wait(&request, &status);
            tilesReceive++;

            memcpy(&image[status.MPI_TAG * pixelPerTiles], tempImage, pixelPerTiles);

            if (tilesSent < nbTiles)
            {
                MPI_Isend(&tilesSent, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &request);
                tilesSent++;
            }
            else
            {
                MPI_Isend(&stop, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &request);
                MPI_Wait(&request,&status);
            }
        }

        save_rasterfile("mandelTP2.ras", w, h, image);
        /* stop timer */
        double end = wallclock_time();
        fprintf(stderr, "Total computing time: %g sec\n", end - start);
    }
    else
    {
        int lineToWorkOn;
        unsigned char *localImage = malloc(pixelPerTiles * sizeof(unsigned char));
        MPI_Status statusSend;
        MPI_Request requestSend;
        MPI_Status statusRcv;
        MPI_Request requestRcv;
        while (1)
        {
            MPI_Irecv(&lineToWorkOn, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &requestRcv);

            MPI_Wait(&requestRcv, &statusRcv);

            if (lineToWorkOn == -1)
            {
                break;
            }

            double y = ymin + (linesPerTiles * lineToWorkOn) * yinc;
            for (int i = 0; i < linesPerTiles; i++)
            {
                double x = xmin;

                if (linesPerTiles * lineToWorkOn + i >= h)
                {
                    break;
                }
                for (int j = 0; j < w; j++)
                {
                    localImage[j + i * w] = xy2color(x, y, depth);
                    x += xinc;
                }
                y += yinc;
            }

            MPI_Isend(localImage, pixelPerTiles, MPI_UNSIGNED_CHAR, 0, lineToWorkOn, MPI_COMM_WORLD, &requestSend);
        }
    }

    MPI_Finalize();

    return 0;
}
