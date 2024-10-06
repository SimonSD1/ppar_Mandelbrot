#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int my_rank;
    int nbProcess;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbProcess);

    MPI_Status status;
    MPI_Request request;

    int w = 100;
    int h = 100;

    int nbTiles = 50;

    int LinesPerProcess = h / nbTiles;
    int PixelPerProcess = LinesPerProcess * w;

    int linesSent = 0;

    int stop = -1;

    // tag 0 = envoi de ligne a calculer
    // tag 1 = reponse de worker

    if (my_rank == 0)
    {
        unsigned char *image = malloc(h * w * sizeof(unsigned char));
        unsigned char *tempImage = malloc(LinesPerProcess*w*sizeof(unsigned char));
        // on envoi a tout le monde au debut
        for (int process = 1; process < nbProcess && linesSent < h; process++, linesSent++)
        {
            MPI_Isend(&linesSent, 1, MPI_INT, process, 0, MPI_COMM_WORLD, &request);
            printf("0 a envoye la ligne %d\n", linesSent);
        }

        // on continue a envoyer le reste quand on a des reponses
        while (linesSent < h)
        {
            MPI_Irecv(tempImage,LinesPerProcess*w,MPI_UNSIGNED_CHAR,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&request);
        }

        for (int i = 1; i < nbProcess; i++)
        {
            printf("envoi stop\n");
            MPI_Isend(&stop, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
        }
    }
    else
    {
        int lineToWorkOn;
        unsigned char *localImage = malloc(LinesPerProcess * w * sizeof(unsigned char));
        while (1)
        {
            MPI_Irecv(&lineToWorkOn, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);

            MPI_Wait(&request, &status);

            printf("le processus n° %d a recue la ligne %d\n", my_rank, lineToWorkOn);

            if (lineToWorkOn == -1)
            {
                break;
            }
        }
    }

    MPI_Finalize();

    return 0;
}
