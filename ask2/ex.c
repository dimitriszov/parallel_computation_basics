/************************************************************/
/* Εργαστήριο εισαγωγής στον παράλληλο υπολογισμό           */
/* Άσκηση 2 / 2021-22                                       */
/* Ονοματεπώνυμο: Ζωβοϊλης Δημήτριος-Μάριος                 */
/* ΑΜ: 19390064                                             */
/* Το πρόγραμμα αυτό τρέχει μόνο όταν το ‘N’ είναι ακέραιο  */
/* πολλαπλάσιο του ‘p’                                      */
/************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int* read_data(int*, int);
int is_strictly_diagonally_dominant(int*, int, int, int);
int find_max_Aii(int*, int, int, int);
int* create_array_B(int*, int, int, int, int);
int find_min(int*, int, int, int*);
void print_2d_mat(int*, int, int);


/**********************************/
/*          driver code           */
/**********************************/
int main(int argc, char** argv)
{
    int rank, root;
    int n, p, size, rows;
    int *A, *A_loc;
    int loc_result, global_result;
    int max, max_loc;
    int *B, *B_loc;
    int data_pair[2], reduction_result[2];
    int loc_min_location;
    int min, min_i, min_j;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    if (rank == 0)
    {
        /**********************************/
        /*        read input data         */
        /**********************************/
        A = read_data(&n, p);
        // print_2d_mat(A, n, n);
    }
    
    /**********************************/
    /* Αποστολή στις διεργασίες του   */
    /* πλήθους των στοιχείων          */
    /**********************************/
    root = 0; // The processor that handles IO has rank == 0
    MPI_Bcast(&n, 1, MPI_INT, root, MPI_COMM_WORLD);
    
    size = n * n / p; // amount of integers in A_loc
    rows = n / p; // rows of array A_loc
    
    // allocate memory for array A_loc
    A_loc = (int *) malloc (size * sizeof(int));
    
    /**********************************/
    /* scatter array A from root to   */
    /* all the other processors into  */
    /* array A_loc                    */
    /**********************************/
    MPI_Scatter(A, size, MPI_INT, A_loc, size, MPI_INT, root, MPI_COMM_WORLD);
    
    /**********************************/
    /* Find for every A_loc if it is  */
    /* strictly diagonally dominant.  */
    /*                                */
    /* Then MPI_Allreduce the result  */
    /* with MPI_MIN so that even if   */
    /* one A_loc wasn't strictly      */
    /* diagonally dominant            */
    /* (loc_result == 0) then all the */
    /* processors will have           */
    /* global_result == 0             */
    /**********************************/
    loc_result = is_strictly_diagonally_dominant(A_loc, rows, n, rank);
    MPI_Allreduce(&loc_result, &global_result, 1, MPI_INT, MPI_MIN,
              MPI_COMM_WORLD);
    
    /**********************************/
    /* if one A_loc isn't strictly    */
    /* diagonally_dominant (we have   */
    /* global_result == 0) print no   */
    /* and then free the space        */
    /* allocated and exit.            */
    /**********************************/
    if(global_result == 0)
    {
        if (rank == 0)
        {
            printf("no\n");
            free(A);
        }
        
        // free the space and exit
        free(A_loc);
        MPI_Finalize();
        exit(0);
    }
    
    /**********************************/
    /*  Calculate max_Aii for every   */
    /*  one A_loc and then calculate  */
    /*   global using MPI_Allreduce   */
    /*   with MPI_MAX so that every   */
    /*    processor has this value    */
    /**********************************/
    max_loc = find_max_Aii(A_loc, rows, n, rank);
    MPI_Allreduce(&max_loc, &max, 1, MPI_INT, MPI_MAX,
              MPI_COMM_WORLD);
    
    /**********************************/
    /* if A is strictly diagonally    */
    /* dominant (global_result == 1)  */
    /* print yes and the value of m   */
    /**********************************/
    if (rank == 0)
    {
        printf("yes\n");
        printf("m = %d\n", max);
        // allocate memory for array B
        B = (int *) malloc(n*n * sizeof(int));
    }
    
    /**********************************/
    /* Create local arrays B_loc with */
    /* A_loc and max.                 */
    /*                                */
    /* After, use MPI_Gather to       */
    /* create array B on processor 0  */
    /* and then print array B.        */
    /**********************************/
    B_loc = create_array_B(A_loc, rows, n, rank, max);
    MPI_Gather(B_loc, size, MPI_INT, B, size, MPI_INT, root, MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("\nArray B:\n");
        print_2d_mat(B, n, n);
    }
    
    /**********************************/
    /*  Use find_min to find the min  */
    /*    element of every B_loc (    */
    /*     data_pair[0]) and it's     */
    /*     location (data_pair[1])    */
    /*                                */
    /* After, use MPI_Reduce with     */
    /* MPI_MINLOC so that proc 0 has  */
    /* the global min value and it's  */
    /* location (reduction_result[0]  */
    /*    and reduction_result[1])    */
    /**********************************/
    data_pair[0] = find_min(B_loc, size, rank, &loc_min_location);
    data_pair[1] = loc_min_location;
    MPI_Reduce(data_pair, reduction_result, 1, MPI_2INT, MPI_MINLOC, root, MPI_COMM_WORLD);
    
    /**********************************/
    /* Print min and its location and */
    /* EXIT                           */
    /**********************************/
    if (rank == 0)
    {
        // calculate and print min - min location
        min = reduction_result[0];
        min_i = reduction_result[1]/n;
        min_j = reduction_result[1]%n;
        printf("min = %d and it is located", min);
        printf(" on the %d row and the %d col\n", min_i, min_j);
        // free space that was allocated only on root
        free(A);
        free(B);
    }
    
    /**********************************/
    /* free the space dynamically     */
    /* allocated and exit             */
    /**********************************/
    free(A_loc);
    free(B_loc);
    free(array_min_location);
    MPI_Finalize();
    return 0;
}



/**********************************/
/*        useful functions        */
/**********************************/

/***********************************************************/
/* This function is used to read the dimensions of array A */
/* and read the numbers and then return them.              */
/*    - n == cols and rows of array A                      */
/*    - p == number of processors                          */
/*    - data == array containing the numbers read          */
/***********************************************************/
int* read_data(int *n, int p)
{
    int i, *data;
    /*
     * in order this program to work 
     * we need n % p == 0 and n > 1. The do -
     * while loop ensures that these criteria are 
     * met in order to continue the programm
     */
    do
    {
        // read n
        printf("Dose to N:\n");
        scanf("%d", n);
    } while((*n % p != 0) || *n < 1);
    data = (int *) malloc(*n*(*n) * sizeof(int));
    if(data == NULL)
    {
        printf("Error allocating memory");
        exit(-1);
    }
    printf("Dose tous %d arithmous:\n", *n*(*n));
    for (i=0; i<*n*(*n); i++)
        scanf("%d", &data[i]);
    
    return data;
}


/***********************************************************/
/* This function is used to check if an array is strictly  */
/* diagonally dominant. If it is 1 is returned, else 0.    */
/*    - mat == input array                                 */
/*    - rows == rows of mat                                */
/*    - cols == cols of mat                                */
/*    - rank == rank of the proc that called the function  */
/***********************************************************/
int is_strictly_diagonally_dominant(int *mat, int rows, int cols, int rank)
{
    int i, j, sum, real_i;
    for (i = 0; i < rows; i++)
    {
        real_i = i+rank*rows;
        // calculate sum
        for (j = 0, sum = 0; j < cols; j++)
            if(real_i != j)
                sum += abs(mat[i*cols + j]);
        // check for this row if the mat is strictly diagonally dominant
        if(sum >= abs(mat[i*cols+rank*rows + i]))
            return 0;
    }
    
    return 1;
}


/***********************************************************/
/* This function is used to find max Aii into an array and */
/* then returns it.                                        */
/*    - mat == input array                                 */
/*    - rows == rows of mat                                */
/*    - cols == cols of mat                                */
/*    - rank == rank of the proc that called the function  */
/*    - max == max A_ii on the array                       */
/***********************************************************/
int find_max_Aii(int *mat, int rows, int cols, int rank)
{
    int i, j, max;
    for (i = 0, max = abs(mat[i*cols+rank*rows + i]); i < rows; i++)
    {
        // find max A_ii
        if(abs(mat[i*cols+rank*rows + i]) > max)
            max = abs(mat[i*cols+rank*rows + i]);
    }
    return max;
}


/***********************************************************/
/* This function is used to create array B_loc based on    */
/* array A_loc, max and then return it.                    */
/*    - A_loc == input array                               */
/*    - rows == rows of mat                                */
/*    - cols == cols of mat                                */
/*    - rank == rank of the proc that called the function  */
/*    - max == max_Aii                                     */
/*    - B_loc == output array                              */
/***********************************************************/
int* create_array_B(int *A_loc, int rows, int cols, int rank, int max)
{
    int i, j, *B_loc;
    B_loc = (int *) malloc (rows*cols * sizeof(int));
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            if(i+rank*rows != j)
                B_loc[i*cols + j] = max - abs(A_loc[i*cols + j]);
            else
                B_loc[i*cols + j] = max;
        }
    }
    return B_loc;
}


/***********************************************************/
/* This function is used to find the minimum element of a  */
/* 2d array that is stored into an 1d array                */
/*    - array == 1d array containing a 2d array            */
/*    - size == number of elements on the array            */
/*    - rank == rank of the processor that called          */
/*              the function                               */
/*    - location == location of the minimum element        */
/*                                                         */
/* The location that will be returned contains the actual  */
/* location of the element on the global mat since we can  */
/* add the offset from the begining of the array adding on */
/* i the size of each subarray multiplied by the rank      */
/* (location = size*rank + i). After that in order to      */
/* convert the location to 2d we just have to find:        */
/* i == location / N                                       */
/* j == location % N                                       */
/***********************************************************/
int find_min(int *array, int size, int rank, int *location)
{
    int i, min;
    for(i = 1, min = array[0], *location = size*rank; i < size; i++)
        if(array[i] < min)
        {
            min = array[i];
            /**********************************/
            /* find the actual location of    */
            /* the element on the global mat  */
            /**********************************/
            *location = size*rank + i;
        }
    return min;
}


/***********************************************************/
/* This function is used to print a 2d array that is       */
/* stored into an 1d array(mat)                            */
/*    - mat == 1d array containing a 2d array              */
/*    - rows == rows of the array                          */
/*    - cols == columns of the array                       */
/***********************************************************/
void print_2d_mat(int *mat, int rows, int cols)
{
    int i, j, k, l;
    
    /**********************************/
    /* these print are for decorating */
    /* purposes                       */
    /**********************************/
    printf("+");
    for (i = 0; i<cols; i++)
    {
        for (j = 0; j < 6; j++)
        {
            printf("-");
        }
        printf("+");
    }
    printf("\n");
    
    for (i = 0; i<rows; i++)
    {
        for (j = 0; j < cols; j++)
        { // actual printing of the numbers
            printf("| %04d ", mat[i*cols + j]);
        }
        printf("|\n+");
        
        /**********************************/
        /* these print are for decorating */
        /* purposes                       */
        /**********************************/
        for (k = 0; k<cols; k++)
        {
            for (l = 0; l < 6; l++)
            {
                printf("-");
            }
            printf("+");
        }
        printf("\n");
    }
}
