#include <mpi.h>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) { 

    int rc, send,recv;


    int parent_size, parent_id, my_id, numprocs; 
    // parent intercomm
    MPI_Comm parent; 
    MPI_Init(&argc, &argv); 

    MPI_Comm_get_parent(&parent); 
    if (parent == MPI_COMM_NULL) {
        cout << "No parent!" << endl;
    }
    MPI_Comm_remote_size(parent, &parent_size); 
    MPI_Comm_rank(parent, &parent_id) ; 
    //cout << "Parent is of size: " << size << endl;
    if (parent_size != 1) {
        cout << "Something's wrong with the parent" << endl;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &my_id) ;     
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs) ;  

    cout << "I'm child process rank "<< my_id << " and we are " << numprocs << endl;
    cout << "The parent process rank "<< parent_id << " and we are " << parent_size << endl;

    /*
    // get value of send
    send = 7; // just an example
    recv = 0;

    rc = MPI_Reduce(&send, &recv, 1, MPI_INT, MPI_SUM, parent_id, parent);
    if (rc != MPI_SUCCESS)
        cout << my_id << " failure on mpi_reduce in WORKER" << endl;

    */
    MPI_Finalize(); 
    return 0; 
} 
