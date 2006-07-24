#include <malloc.h>
#include <pthread.h>
#include <mpi.h>
#include <stdio.h>

/* The new variant: uses allgatherv instead of alltoallv
   wait does the free too, calls can be chained */
/* When chaining, note you have to WAIT for only the last one.
   The program behaviour will be undefined, if the intermediate values are
   WAIT-ed as well */
/* For an example of usage, see TransmitCoarse */

/* NB! all the structures used on sending need to be in the memory during
       the whole sending process, not just its init */

struct mpi_allgatherv_args
{
	void *sendbuf;
	void *sendcnt;
	void *sendtype;
	void *recvbuf;
	void *recvcnts;
	void *rdispls;
	void *recvtype;
	void *comm;
	volatile int finished;
        struct mpi_allgatherv_args *prev;
	pthread_t thread;
};

void *mpi_allgatherv_nb_threadproc(void *ptr)
{
	long long ierr = 0;
	struct mpi_allgatherv_args *args = (struct mpi_allgatherv_args *) ptr;
        
        // Wait for the previous to finish, so we could sequence threads
        if (args->prev!=NULL)
        {
            pthread_join(args->prev->thread, NULL);
	    free(args->prev);
        }
        
	mpi_allgatherv_(args->sendbuf, args->sendcnt, args->sendtype, args->recvbuf, args->recvcnts, args->rdispls, args->recvtype, args->comm, &ierr);

	args->finished = 1;
        
        return 0;
}

void getnull_(void **res)
{   res=NULL; }

void mpi_allgatherv_nb_chain_(
	void *sendbuf, void *sendcnt, void *sendtype,
	void *recvbuf, void *recvcnts, void *rdispls, void *recvtype,
	void *comm, void ** prev, void **res)
{
	long long ierr = 0;
	long long exitcode = -1;
	struct mpi_allgatherv_args *args;
	args = (struct mpi_allgatherv_args *) malloc(sizeof(struct mpi_allgatherv_args));

	if (args == NULL)
	{
		mpi_abort_(comm, &exitcode, &ierr);
		return;
	}

	args->sendbuf = sendbuf;
	args->sendcnt = sendcnt;
	args->sendtype = sendtype;
	args->recvbuf = recvbuf;
	args->recvcnts = recvcnts;
	args->rdispls = rdispls;
	args->recvtype = recvtype;
	args->comm = comm;
	args->finished = 0;

        if (prev==NULL) args->prev=NULL;
        else args->prev = *((struct mpi_allgatherv_args **) prev);

	if (pthread_create(&args->thread, NULL, mpi_allgatherv_nb_threadproc, (void *) args) != 0)
	{
		mpi_abort_(comm, &exitcode, &ierr);
		return;
	}
	*res = args;
}

void mpi_allgatherv_nb_init_(
	void *sendbuf, void *sendcnt, void *sendtype,
	void *recvbuf, void *recvcnts, void *rdispls, void *recvtype,
	void *comm, void **res)
{
    mpi_allgatherv_nb_chain_(sendbuf,sendcnt,sendtype,
                             recvbuf,recvcnts,rdispls,recvtype,
                             comm, NULL, res);
}

void mpi_allgatherv_nb_wait_(void **ptr)
{
	struct mpi_allgatherv_args *args;
	args = *((struct mpi_allgatherv_args **) ptr);
//	if (!args->finished) 
                pthread_join(args->thread, NULL);

	free(args);
}

void mpi_allgatherv_nb_test_(void **ptr, int *res)
{
	struct mpi_allgatherv_args *args;
	args = *((struct mpi_allgatherv_args **) ptr);
	*res = !args->finished;
}

/* The old variant. ALLTOALL seems to be a bit too much,
   as it is meant for cases where everyone has different data for
   all the others instead of sending the same thing to everyone */

struct mpi_alltoall_args
{
	void *sendbuf;
	void *sendcnts;
	void *sdispls;
	void *sendtype;
	void *recvbuf;
	void *recvcnts;
	void *rdispls;
	void *recvtype;
	void *comm;
	volatile int finished;
	pthread_t thread;
};

void *mpi_alltoall_nb_threadproc(void *ptr)
{
	long long ierr = 0;
	struct mpi_alltoall_args *args = (struct mpi_alltoall_args *) ptr;
	mpi_alltoallv_(args->sendbuf, args->sendcnts, args->sdispls, args->sendtype, args->recvbuf, args->recvcnts, args->rdispls, args->recvtype, args->comm, &ierr);
	args->finished = 1;
}

void mpi_alltoall_nb_init_(
	void *sendbuf, void *sendcnts, void *sdispls, void *sendtype,
	void *recvbuf, void *recvcnts, void *rdispls, void *recvtype,
	void *comm, void **res)
{
	long long ierr = 0;
	long long exitcode = -1;
	struct mpi_alltoall_args *args;
	args = (struct mpi_alltoall_args *) malloc(sizeof(struct mpi_alltoall_args));
	if (args == NULL)
	{
		mpi_abort_(comm, &exitcode, &ierr);
		return;
	}
	args->sendbuf = sendbuf;
	args->sendcnts = sendcnts;
	args->sdispls = sdispls;
	args->sendtype = sendtype;
	args->recvbuf = recvbuf;
	args->recvcnts = recvcnts;
	args->rdispls = rdispls;
	args->recvtype = recvtype;
	args->comm = comm;
	args->finished = 0;
	if (pthread_create(&args->thread, NULL, mpi_alltoall_nb_threadproc, (void *) args) != 0)
	{
		mpi_abort_(comm, &exitcode, &ierr);
		return;
	}
	*res = args;
}

void mpi_alltoall_nb_wait_(void **ptr)
{
	struct mpi_alltoall_args *args;
	args = *((struct mpi_alltoall_args **) ptr);
	if (!args->finished) pthread_join(args->thread, NULL);
}

void mpi_alltoall_nb_test_(void **ptr, int *res)
{
	struct mpi_alltoall_args *args;
	args = *((struct mpi_alltoall_args **) ptr);
	*res = !args->finished;
}

void mpi_alltoall_nb_free_(void **ptr)
{
	mpi_alltoall_nb_wait_(ptr);
	struct mpi_alltoall_args *args;
	args = *((struct mpi_alltoall_args **) ptr);
	free(args);
}