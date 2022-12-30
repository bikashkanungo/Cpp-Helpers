#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "mpi.h"

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *buffer;
    long BUFFER_SIZE = (512 * 1024); /* # elements in buffer *for each MPI proc */
    if (argc > 1) {
        BUFFER_SIZE = atol(argv[1]);
    }
    int i;
    buffer = (int *)malloc(BUFFER_SIZE * sizeof(int));
    for (i = 0; i < BUFFER_SIZE; i++) {
        buffer[i] = rank * BUFFER_SIZE + i;
    }

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, "data.bin", MPI_MODE_CREATE | MPI_MODE_RDWR,
		    MPI_INFO_NULL, &fh);

    int count = 1;
    int blocklength = BUFFER_SIZE;
    int stride = BUFFER_SIZE;
    MPI_Datatype etype = MPI_INT;
    MPI_Datatype filetype;
    MPI_Type_vector(count, blocklength, stride, etype, &filetype);
    MPI_Type_commit(&filetype);

    MPI_Status status;
    /* Both methods end in the same error; Changing to non collective versions works
     * fine. */
    //#if 0
    MPI_Offset disp = (MPI_Offset)rank * BUFFER_SIZE * sizeof(int);
    MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, (void *)buffer, BUFFER_SIZE, etype, &status);
    //#else
    //MPI_File_write_at_all(fh, BUFFER_SIZE * rank, (void *)buffer, BUFFER_SIZE,
    //                      etype, &status);
    //#endif

    MPI_File_close(&fh);
    free(buffer);

    MPI_Finalize();

    return 0;
}
