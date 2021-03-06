#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <vector>
#include <fstream>
using namespace std;
void show_vector( vector<int>&a)
{
    for (vector<int>::iterator it = a.begin() ; it!=a.end() ; ++it)
        cout<<*it<<"   ";
}
int main() {

    std::vector<int> vecx;
    std::vector<int> vecy;
   // omp_set_num_threads(2);
    const long samples=1000000l;
    long in = 0l;
#pragma omp parallel num_threads(4) reduction(+:in)
    {
        ofstream fout;
        int i, n, k=0;
        string name;
        i= omp_get_thread_num();
        std::string ii = std::to_string(i+1);
        n=omp_get_num_threads();
        std::string nn = std::to_string(n);
        name="/home/egor/quest/trng/bin/all_"+nn+"_curr_"+ii+".txt";
        //cout<< name<<endl;
        fout.open(name);

        trng::yarn2 r;
        int size=omp_get_num_threads();
        int rank=omp_get_thread_num();
        //std::cout<<"size:  "<<omp_get_num_threads()<<"  "<<"rank:  "<<omp_get_thread_num()<<std::endl;
        trng::uniform01_dist<> u;
        //r.jump(2*(rank*samples/size));

        if (rank==0)
            r.jump(0);
        else if(rank==1)
            r.jump(2);
        else if (rank==2)
            r.jump(4);
        else if(rank==3)
            r.jump(6);
        for (long i=rank*samples/size; i<(rank+1)*samples/size; ++i) {

            //cout<<rank<<endl;
//            r.jump(2*(size*k+rank));



            k++;
            double x=u(r), y=u(r);
            // choose random x − and y − coordinates
            fout <<x<<"    "<<y<<"\n" ;
            if (x*x+y*y<=1.0){
                ++in;
            }
            if (size==4)
            r.jump(size+2);
            else if (size==2)
                r.jump(size);
            else if (size == 3)
                r.jump(size+1);
            else if(size==1)
                r.jump(0);

        }
        fout.close();
    }
    std::cout << " pi = " << 4.0*in/samples << std::endl;
    return EXIT_SUCCESS;
}







/*
#include <cstdlib>
#include <iostream>
#include "/home/egor/LIBS/openmpi/include/mpi.h"
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
int main(int argc, char *argv[]) {
    const long samples=1000000l;
    trng::yarn2 r;
    MPI::Init(argc, argv);
    int size=MPI::COMM_WORLD.Get_size();
    int rank=MPI::COMM_WORLD.Get_rank();
    long in=0l;
    trng::uniform01_dist<> u;
    r.jump(2*(rank*samples/size));
    // total number of points in square
    // random number engine
    // i n i t i a l i s e MPI environment
    // get total number of processes
    // get rank of current process
    // number of points in c i r c l e
    // random number distribution
    // jump ahead
    // throw random points into square and distribute workload over a l l processes
    for (long i=rank*samples/size; i<(rank+1)*samples/size; ++i) {
        double x=u(r), y=u(r);
        // choose random x − and y − coordinates
        if (x*x+y*y<=1.0)
            // i s point in c i r c l e ?
            ++in;
        // increase counter
    }
    // calculate sum of a l l local variables ’ in ’ and storre result in ’ in_all ’ on process 0
    long in_all;
    MPI::COMM_WORLD.Reduce(&in, &in_all, 1, MPI::LONG, MPI::SUM, 0);
    if (rank==0)
        // print result
        std::cout << " pi = " << 4.0*in_all/samples << std::endl;
    MPI::Finalize();
    // quit MPI
    return EXIT_SUCCESS;
}
*/


/*
#include <iostream>
#include <cstdlib>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

int main() {
const long samples=100000000l;
long in=0l;
double PI, s;
trng::yarn2 r;
trng::uniform01_dist<> u;
//total number of points in square
//no points in c i r c l e
//random number engine
//random number distribution
// throw random points into square
for (long i=0; i<samples; ++i) {
double x=u(r), y=u(r);
// choose random x − and y − coordinates
if (x*x+y*y<=1.0)
// i s point in c i r c l e ?
++in;
// increase counter
}
s=samples;
PI=4.0*(double(in))/(double(s));
std::cout<<in<<"  "<<samples<<std::endl;
std::cout << "pi = " << PI << std::endl;
printf("%.10f \n", PI);
return EXIT_SUCCESS;
}*/
