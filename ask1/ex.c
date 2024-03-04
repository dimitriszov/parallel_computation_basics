/************************************************************/
/* Εργαστήριο εισαγωγής στον παράλληλο υπολογισμό           */
/* Ασκηση 1 / 2021-22                                       */
/* Ονοματεπώνυμο: Ζωβοϊλης Δημήτριος-Μάριος                 */
/* ΑΜ: 19390064                                             */
/************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


int* read_data(int*, int);
int* calculate_info(int*, int);
int menu();


/**********************************/
/*          driver code           */
/**********************************/
int main(int argc, char** argv)
{
    int rank, option, source, target;
    int n, p, k, min_final, max_final, num;
    double m_final, var_local, var_final;
    int tag1=50, tag2=60, tag3=70, tag4=80, tag5=90, tag6=100, tag7 = 110;
    int *info, *info_loc;
    int *data, *data_loc;
    double *d, *d_loc;
    MPI_Status status;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    do {
        
        if (rank == 0)
        {
            /**********************************/
            /*        read input data         */
            /**********************************/
            data = read_data(&n, p);
            
            /**********************************/
            /*    send n to all processors    */
            /**********************************/
            for (target = 1; target < p; target++)
                MPI_Send(&n, 1, MPI_INT, target, tag1, MPI_COMM_WORLD);
            
            /**********************************/
            /* scatter initial matrix to the  */
            /* other processors               */
            /**********************************/
            num = n/p;
            k=num;
            for (target = 1; target < p; target++)
            {
                MPI_Send(&data[k], num, MPI_INT, target, tag2, MPI_COMM_WORLD);
                k+=num;
            }
            
            /**********************************/
            /* copy the initial part of       */
            /* matrix to data_loc             */
            /**********************************/
            data_loc = (int *) malloc (num * sizeof(int));
            for (k=0; k<num; k++)
                data_loc[k]=data[k];
        }
        else
        {
            /**********************************/
            /*   receive n from processor 0   */
            /**********************************/
            MPI_Recv(&n, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD, &status);
            
            /**********************************/
            /* receive the scattered part of  */
            /* the initial data and place     */
            /* them to  data_loc              */
            /**********************************/
            num = n/p;
            data_loc = (int *) malloc (num * sizeof(int));
            MPI_Recv(&data_loc[0], num, MPI_INT, 0, tag2, MPI_COMM_WORLD, &status);
        }
        
        /**********************************/
        /* calculate local m, min and max */
        /* and return it to matrix        */
        /* info_loc For the elements of   */
        /* info_mat we have that:         */
        /*    - info_mat[0] == m_local    */
        /*    - info_mat[1] == min_local  */
        /*    - info_mat[2] == max_local  */
        /**********************************/
        info_loc = calculate_info(data_loc, num);
        
        /**********************************/
        /* broadcast local m, min and max */
        /* and calculate final m, min and */
        /* max                            */
        /**********************************/
        if (rank != 0)
        {
            MPI_Send(&info_loc[0], 3, MPI_INT, 0, tag3, MPI_COMM_WORLD);
        }
        else
        {
            m_final = info_loc[0];
            min_final = info_loc[1];
            max_final = info_loc[2];
            // printf("\n m_local of process %d: %d\n", rank, info_loc[0]);
            for (source = 1; source < p; source++)
            {
                MPI_Recv(&info_loc[0], 3, MPI_INT, source, tag3, MPI_COMM_WORLD, &status);
                
                // calculate m_final
                m_final = m_final + info_loc[0];
                // printf("\n m_local of process %d: %d\n", source, info_loc[0]);
                
                // calculate min_final
                if(min_final > info_loc[1])
                    min_final = info_loc[1];
                
                // calculate max_final
                if(max_final < info_loc[2])
                    max_final = info_loc[2];
            }
            m_final = m_final/n;
            
            // print the m_final, min_final, max_final
            printf("\n\n\n m_final: %f", m_final);
            printf("\n min_final: %d", min_final);
            printf("\n max_final: %d", max_final);
        }
        
        info = (int*) malloc(2 * sizeof(int));
        if(info == NULL)
        {
            printf("Error allocating memory");
            exit(-1);
        }
        
        /**********************************/
        /* send back to all processors    */
        /* m_final, min and max           */
        /**********************************/
        if (rank == 0)
        {
            info[0] = min_final;
            info[1] = max_final;
            
            for (target = 1; target < p; target++)
            {
                MPI_Send(&m_final, 1, MPI_DOUBLE, target, tag4, MPI_COMM_WORLD);
                MPI_Send(&info[0], 2, MPI_INT, target, tag5, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Recv(&m_final, 1, MPI_DOUBLE, 0, tag4, MPI_COMM_WORLD, &status);
            MPI_Recv(&info[0], 2, MPI_INT, 0, tag5, MPI_COMM_WORLD, &status);
        }
        
        /**********************************/
        /* allocate memory for mat d_loc  */
        /**********************************/
        d_loc = (double *) malloc(num * sizeof(double));
        
        /**********************************/
        /* calculate var_local and mat d  */
        /**********************************/
        var_local = 0;
        for (k=0; k<num; k++)
        {
            // calculate var_local
            var_local = var_local + (data_loc[k]-m_final)*(data_loc[k]-m_final);
            
            // caclulate table d
            d_loc[k] = ((double)data_loc[k]-info[0])/(info[1]-info[0]) * 100;
        }
        
        /**********************************/
        /*    send var_local from each    */
        /*    processor to the main one   */
        /**********************************/
        if (rank != 0)
        {
            MPI_Send(&var_local, 1, MPI_DOUBLE, 0, tag6, MPI_COMM_WORLD);
            MPI_Send(&d_loc[0], num, MPI_DOUBLE, 0, tag7, MPI_COMM_WORLD);
        }
        else
        {
            // allocate memory for table d
            d = (double *) malloc(n * sizeof(double));
            
            for(k=0; k<num; k++)
                d[k] = d_loc[k];
            var_final = var_local;
            // printf("\n var_local of process %d: %f\n", rank, var_local);
            for (source = 1; source < p; source++)
            {
                MPI_Recv(&var_local, 1, MPI_DOUBLE, source, tag6, MPI_COMM_WORLD, &status);
                // calculate var_final
                var_final = var_final + var_local;
                // printf("\n var_local of process %d: %f\n", source, var_local);
                MPI_Recv(&d[source*num], num, MPI_DOUBLE, source, tag7, MPI_COMM_WORLD, &status);
            }
            var_final = (double)var_final/n;
            printf("\n\n\n var_final: %f\n", var_final);
            for(k=0; k<n; k++)
                printf("d[%d]=%f\n", k, d[k]);
        }
        
        
        /**********************************/
        /*      free resources used       */
        /**********************************/
        if (rank == 0)
        {
            free(data);
            free(d);
        }
        free(data_loc);
        free(d_loc);
        free(info_loc);
        
        /**********************************/
        /* To exit for loop we need to    */
        /* send to all processors the     */
        /* value of option                */
        /**********************************/
        if (rank == 0)
        {
            // display the menu and get option
            option = menu();
            
            // send the selection to other processors
            for (target = 1; target < p; target++)
                MPI_Send(&option, 1, MPI_INT, target, tag1, MPI_COMM_WORLD);
        } 
        else
        {
            MPI_Recv(&option, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD, &status);
        }
    } while(option == 1);
    MPI_Finalize();
}


/**********************************/
/*        useful functions        */
/**********************************/


/***********************************************************/
/* This function is used to read how many numbers the user */
/* will input, read the numbers and then return them.      */
/*    - n == length of matrix                              */
/*    - p == number of processors                          */
/*    - data == matrix with the numbers read               */
/***********************************************************/
int* read_data(int *n, int p)
{
    int i, *data;
    /*
     * in order to work with multiple processors 
     * we need n % p == 0 and n > 1. The do -
     * while loop ensures that these criteria are 
     * met in order to continue the programm
     */
    do
    {
        // read n
        printf("Dose plithos arithmon:\n");
        scanf("%d", n);
    } while((*n % p != 0) || *n < 1);
    data = (int *) malloc(*n * sizeof(int));
    if(data == NULL)
    {
        printf("Error allocating memory");
        exit(-1);
    }
    printf("Dose tous %d arithmous:\n", *n);
    for (i=0; i<*n; i++)
        scanf("%d", &data[i]);
    
    return data;
}


/***********************************************************/
/* This function is used by every processor and calculates */
/* m, min and max on the section of the initial mat that   */
/* they have. Then it returns a matrix containing those    */
/* elements.                                               */
/*    - data_loc is the matrix containing a section        */
/* of the initial mat,                                     */
/*    - size is the length of data_loc                     */
/* For the elements of info_mat we have that:              */
/*    - info_mat[0] == m_local                             */
/*    - info_mat[1] == min_local                           */
/*    - info_mat[2] == max_local                           */
/***********************************************************/
int* calculate_info(int *data_loc, int size)
{
    int i, *info_mat;
    
    /**********************************/
    /*  allocate memory for info_mat  */
    /**********************************/
    info_mat = (int*) malloc(3 * sizeof(int));
    if(info_mat == NULL)
    {
        printf("Error allocating memory");
        exit(-1);
    }
    
    /**********************************/
    /* calculate local m, min and max */
    /* For the elements of info_mat   */
    /* we have that:                  */
    /*    - info_mat[0] == m_local    */
    /*    - info_mat[1] == min_local  */
    /*    - info_mat[2] == max_local  */
    /**********************************/
    info_mat[0] = data_loc[0];
    info_mat[1] = info_mat[2] = data_loc[0];
    for (i = 1; i < size; i++) 
    {
        // calculate m_local
        info_mat[0] = info_mat[0] + data_loc[i];
        
        // calculate min_local
        if(info_mat[1] > data_loc[i])
            info_mat[1] = data_loc[i];
        
        // calculate max_local
        if(info_mat[2] < data_loc[i])
            info_mat[2] = data_loc[i];
    }
    
    return info_mat;
}


/***********************************************************/
/* This function displays the menu and returns the users   */
/* option.                                                 */
/***********************************************************/
int menu()
{
    int option;
    printf("\n\n----Menu----\n");
    printf("-To exit press 0\n-To continue press 1\nOption:");
    scanf("%d",&option);
    return option;
}
