#include <omp.h> 
#include <stdio.h> 
#include <stdlib.h> 
  
int main(int argc, char* argv[]) 
{
    #pragma omp parallel 
    { 
       int id = omp_get_thread_num(); 
       printf("Hello World from thread id %d \n",id);
    } 
} 