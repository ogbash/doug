#include <malloc.h>
#include <pthread.h>
#include <mpi.h>

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
