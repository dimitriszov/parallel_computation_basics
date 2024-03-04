/************************************************************/
/* Εργαστήριο εισαγωγής στον παράλληλο υπολογισμό           */
/* Άσκηση 2 / 2021-22                                       */
/* Ονοματεπώνυμο: Ζωβοϊλης Δημήτριος-Μάριος                 */
/* ΑΜ: 19390064                                             */
/* Το πρόγραμμα αυτό τρέχει για οποιονδήποτε συνδυασμό      */
/* τιμών ‘Ν’ και ‘p’ (με χρήση των συναρτήσεων scatterv/    */
/* gatherv).                                                */
/************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int* read_data(int*, int);
void calc_scounts_displs(int**, int**, int, int);
int* find_sums_Aii(int[], int, int, int, int **);
int is_strictly_diagonally_dominant(int[], int[], int);
int find_max_Aii(int[], int);
int* create_array_B(int[], int, int, int, int);
int find_min(int[], int, int, int*);
void print_2d_mat(int[], int, int);


/**********************************/
/*          driver code           */
/**********************************/
int main(int argc, char** argv)
{
    int rank, root;
    int n, p, size, rows;
    int *sendcounts;    // array describing how many elements to send to each process
    int *displs;        // array describing the displacements where each segment begins
    int *A, *A_loc;
    int *sums, *global_sums, *Aii, *global_Aii;
    int result;
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
        //printf("Array A:\n");
        //print_2d_mat(A, n, n);
    }
    
    /**********************************/
    /* Αποστολή στις διεργασίες του   */
    /* πλήθους των στοιχείων          */
    /**********************************/
    root = 0; // The processor that handles IO has rank == 0
    MPI_Bcast(&n, 1, MPI_INT, root, MPI_COMM_WORLD);
    
    // calculate sendcounts and displs
    calc_scounts_displs(&sendcounts, &displs, n, p);
    
    // allocate memory for array A_loc
    A_loc = (int *) malloc (sendcounts[rank] * sizeof(int));
    
    /**********************************/
    /* scatter array A from root to   */
    /* all the other processors into  */
    /* array A_loc                    */
    /**********************************/
    MPI_Scatterv(A, sendcounts, displs, MPI_INT, A_loc, sendcounts[rank], 
                 MPI_INT, root, MPI_COMM_WORLD);
    
    // find local_sums and Aii
    sums = find_sums_Aii(A_loc, sendcounts[rank], displs[rank], n, &Aii);
    
    // get the sums and Aii for each row of Α on
    // global_sums and global_Aii using MPI_Reduce
    if (root == rank)
    {
        global_sums = (int *) malloc(n *sizeof(int));
        global_Aii = (int *) malloc(n *sizeof(int));
    }
    MPI_Reduce(sums, global_sums, n, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
    MPI_Reduce(Aii, global_Aii, n, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
    
    // check if A is strictly diagonally dominant based on sums and Aii
    if (rank == root)
    {
        result = is_strictly_diagonally_dominant(global_sums, global_Aii, n);
    }
    
    /**********************************/
    /* Broadcast the result so that   */
    /* all processors know if A is    */
    /* strictly diagonally dominant   */
    /**********************************/
    MPI_Bcast(&result, 1, MPI_INT, root, MPI_COMM_WORLD);
    
    /**********************************/
    /* if A isn't strictly diagonally */
    /* dominant (we have result == 0) */
    /* print no, then free the space  */
    /* allocated and exit.            */
    /**********************************/
    if(result == 0)
    {
        if (rank == 0)
        {
            printf("no\n");
            free(A);
            free(global_sums);
            free(global_Aii);
        }
        
        // free the space and exit
        free(A_loc);
        free(sendcounts);
        free(displs);
        free(sums);
        free(Aii);
        MPI_Finalize();
        exit(0);
    }

    /**********************************/
    /*  Calculate max_Aii for every   */
    /*  one Aii and then calculate    */
    /*   global using MPI_Allreduce   */
    /*   with MPI_MAX so that every   */
    /*    processor has this value    */
    /**********************************/
    max_loc = find_max_Aii(Aii, n);
    MPI_Allreduce(&max_loc, &max, 1, MPI_INT, MPI_MAX,
              MPI_COMM_WORLD);
    
    /**********************************/
    /* if A is strictly diagonally    */
    /* dominant (result == 1) print   */
    /* yes and the value of m.        */
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
    /* After, use MPI_Gatherv to      */
    /* create array B on processor 0  */
    /* and then print array B.        */
    /**********************************/
    B_loc = create_array_B(A_loc, sendcounts[rank], displs[rank], n, max);
    MPI_Gatherv(B_loc, sendcounts[rank], MPI_INT, B, sendcounts, displs, MPI_INT, root, MPI_COMM_WORLD);
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
    data_pair[0] = find_min(B_loc, sendcounts[rank], displs[rank], &loc_min_location);
    data_pair[1] = loc_min_location;
    MPI_Reduce(data_pair, reduction_result, 1, MPI_2INT, MPI_MINLOC, root, MPI_COMM_WORLD);
    
    /**********************************/
    /* Display the min element of     */
    /* array B and its location. Then */
    /* exit the program.              */
    /**********************************/
    if (rank == 0)
    {
        // calculate and print min, min location
        min = reduction_result[0];
        min_i = reduction_result[1]/n;
        min_j = reduction_result[1]%n;
        printf("min = %d and it is located", min);
        printf(" on the %d row and the %d col\n", min_i, min_j);
        
        // free space that was allocated only on root
        free(A);
        free(global_sums);
        free(global_Aii);
        free(B);
    }
    
    /**********************************/
    /* free the space dynamically     */
    /* allocated and exit             */
    /**********************************/
    free(A_loc);
    free(sendcounts);
    free(displs);
    free(sums);
    free(Aii);
    free(B_loc);
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
    } while(*n*(*n) < p);
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
/* This function is used to calculate the sendcounts and   */
/* displs and then return them                             */
/*    - sendcounts == array describing how many elements   */
/*                    to send to each process              */
/*    - displs == array describing the displacements where */
/*                each segment begins                      */
/*    - array_size == number of elements on array A        */
/*    - proc_size == number of procs                       */
/***********************************************************/
void calc_scounts_displs(int *sendcounts[], int *displs[], int array_size, int proc_size)
{
    int i, rem, sum;
    *sendcounts = malloc(sizeof(int)*proc_size);
    *displs = malloc(sizeof(int)*proc_size);
    // rem == elements remaining after division among processes
    rem = (array_size*array_size)%proc_size;
    sum = 0;  // Sum of counts. Used to calculate displacements
    for (i = 0; i < proc_size; i++) {
        (*sendcounts)[i] = (array_size*array_size)/proc_size;
        if (rem > 0) {
            (*sendcounts)[i]++;
            rem--;
        }

        (*displs)[i] = sum;
        sum += (*sendcounts)[i];
    }
}



/***********************************************************/
/* This function is used to read the dimensions of array A */
/* and read the numbers and then return them.              */
/*    - n == cols and rows of array A                      */
/*    - p == number of processors                          */
/*    - data == array containing the numbers read          */
/***********************************************************/
int* find_sums_Aii(int array[], int size, int displ, int n, int *Aii[])
{
    int i, j, dist, sum, *sums;
    int real_i, real_j;
    sums = (int *) malloc(n*sizeof(int));
    for(i = 0; i < n; i++)
        sums[i] = 0;
    *Aii = (int *) malloc(n*sizeof(int));
    for(i = 0; i <= n; i++)
        (*Aii)[i] = 0;
    dist = 0;
    for (i = 0; i <= (displ+size)/n - displ/n; i++)
    {
        real_i = i+displ/n;
        real_j = (displ+dist)%n;
        // calculate sum
        for (j = 0; j+real_j < n && dist < size; j++)
        {
            if(real_i != j+real_j)
                sums[real_i] += abs(array[dist]);
            else
                (*Aii)[real_i] = abs(array[dist]);
            dist++;
        }
    }
    return sums;
}


/***********************************************************/
/* This function is used to check if an array is strictly  */
/* diagonally dominant. If it is 1 is returned, else 0.    */
/*    - sums == array with the sums for each row of A      */
/*    - Aii == array with the Aii for each row of A        */
/*    - n == rows of array A                               */
/***********************************************************/
int is_strictly_diagonally_dominant(int sums[], int Aii[], int n)
{
    int i;
    for (i = 0; i < n; i++)
        if(sums[i] >= Aii[i])
            return 0;
    
    return 1;
}


/***********************************************************/
/* This function is used to find max Aii into an array of  */
/* Aiis and then returns it.                               */
/*    - Aii == array with the local Aii for each row of A  */
/*    - n == rows of array Aii                             */
/***********************************************************/
int find_max_Aii(int Aii[], int n)
{
    int i, max;
    for (i = 1, max = Aii[0]; i < n; i++)
        if(Aii[i] >= max)
            max = Aii[i];
    
    return max;
}


/***********************************************************/
/* This function is used to create array B_loc based on    */
/* array A_loc, max and then return it.                    */
/*    - A_loc == input array                               */
/*    - size == size of input array                        */
/*    - displ == displs of the proc                        */
/*    - max == max_Aii                                     */
/*    - B_loc == output array                              */
/***********************************************************/
int* create_array_B(int A_loc[], int size, int displ, int n, int max)
{
    int i, j, *B_loc;
    int dist, real_i, real_j;
    B_loc = (int *) malloc (size * sizeof(int));
    dist = 0;
    for (i = 0; i <= (displ+size)/n - displ/n; i++)
    {
        real_i = i+displ/n;
        real_j = (displ+dist)%n;
        for (j = 0; j+real_j < n && dist < size; j++)
        {
            if(real_i != j+real_j)
                B_loc[dist] = max - abs(A_loc[dist]);
            else
                B_loc[dist] = max;
            dist++;
        }
    }
    return B_loc;
}


/***********************************************************/
/* This function is used to find the minimum element of a  */
/* 2d array that is stored into an 1d array                */
/*    - array == 1d array containing a 2d array            */
/*    - size == number of elements on the array            */
/*    - displ == displs of the processor that called       */
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
int find_min(int array[], int size, int displ, int *location)
{
    int i, min;
    min = array[0];
    *location = displ;
    for(i = 1; i < size; i++)
        if(array[i] < min)
        {
            min = array[i];
            /**********************************/
            /* find the actual location of    */
            /* the element on the global mat  */
            /**********************************/
            *location = displ + i;
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
void print_2d_mat(int mat[], int rows, int cols)
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
