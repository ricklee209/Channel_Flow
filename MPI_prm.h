



#define np 6

int nproc = np; // how many CPUs //

int istart, iend, icount, r_nbr, l_nbr, lastp, iroot;
int itag, isrc, idest, istart1, icount1,icount2, istart2, iend1, iend2;

int gstart0[np], gend0[np],gstart[np], gend[np], gcount[np];

double tstart,tend;

MPI_Request l_reqs1,l_reqs2,l_reqs3,l_reqs4,l_reqs5,
			r_reqs1,r_reqs2,r_reqs3,r_reqs4,r_reqs5; 

MPI_Request r1,r2,r3,r4,r5,r6,r7,r8,r9,r10;



