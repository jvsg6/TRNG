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
    omp_set_num_threads(1);
    const long samples=1000000l;
    long in=0l;
    std::cout<<"size:  "<<omp_get_num_threads()<<"  "<<"rank:  "<<omp_get_thread_num()<<std::endl;
    // total number of points in square
    // number of points in c i r c l e
    // distribute workload over a l l processes and make a global reduction
#pragma omp parallel /*num_threads(2)*/ reduction(+:in)
    {
        ofstream fout;
        int i, n ;
        string name;
        i= omp_get_thread_num();
        std::string ii = std::to_string(i);
        n=omp_get_num_threads();
        std::string nn = std::to_string(n);
        name="thread_"+ii+"from_"+nn+".txt";
        fout.open(name);
        fout << "Работа с файлами в С++";
        fout.close();
        trng::yarn2 r;
        int size=omp_get_num_threads();
        int rank=omp_get_thread_num();
        std::cout<<"size:  "<<omp_get_num_threads()<<"  "<<"rank:  "<<omp_get_thread_num()<<std::endl;
        trng::uniform01_dist<> u;
        r.jump(2*(rank*samples/size));
        for (long i=rank*samples/size; i<(rank+1)*samples/size; ++i) {
            double x=u(r), y=u(r);
            // choose random x − and y − coordinates
            vecx.push_back (x);
            vecy.push_back (y);
            if (x*x+y*y<=1.0){
                // i s point in c i r c l e ?
                ++in;
                //std::cout<<"in:  "<<in<<std::endl;
            }
            // increase thread − local counter
        }
    }
    // print result

    //for (vector<int>::iterator it = vecx.begin() ; it!=vecx.end() ; ++it)
        //cout<<*it<<"   "<<endl;
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
