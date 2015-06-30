#include "mpi.h" 
int main(int argc, char *argv[]) 
{ 
   int size, rankworker, sizeworker; 
   MPI_Comm parent; 
   MPI_Init(&argc, &argv); 
   MPI_Comm_get_parent(&parent); 
   if (parent == MPI_COMM_NULL) error("No parent!"); 
   MPI_Comm_remote_size(parent, &size);
   MPI_Comm_size(MPI_COMM_WORLD, &sizeworker);  
   if (size != 1) error("Something's wrong with the parent"); 
 
   /* 
    * Parallel code here.  
    * The manager is represented as the process with rank 0 in (the remote 
    * group of) MPI_COMM_PARENT.  If the workers need to communicate among 
    * themselves, they can use MPI_COMM_WORLD. 
    */ 
    //printf("I am %d of %d\n",
    //printf("size: %d  - sizeworker: %d\n",size,sizeworker);
    MPI_Comm_rank(MPI_COMM_WORLD,&rankworker);
    printf("I am %d of %d\n",rankworker,sizeworker);
     
     
    //receive()
    
    //calc()
    
    //send()
   MPI_Finalize(); 
   return 0; 
} 
 
