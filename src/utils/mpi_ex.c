// DOUG - Domain decomposition On Unstructured Grids
// Copyright (C) 1998-2006 Faculty of Computer Science, University of Tartu and
// Department of Mathematics, University of Bath
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// or contact the author (University of Tartu, Faculty of Computer Science, Chair
// of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
// mailto:info(at)dougdevel.org)

#include <doug_config.h>

#include <malloc.h>
#include <pthread.h>
#include <mpi.h>
#include <stdio.h>

#if FORTRANDOUBLEUNDERSCORE
#define getnull_ getnull__
#define mpi_allgatherv_nb_chain_ mpi_allgatherv_nb_chain__
#define mpi_allgatherv_nb_init_ mpi_allgatherv_nb_init__
#define mpi_allgatherv_nb_wait_ mpi_allgatherv_nb_wait__
#define mpi_allgatherv_nb_test_ mpi_allgatherv_nb_test__
#define mpi_alltoall_nb_init_ mpi_alltoall_nb_init__
#define mpi_alltoall_nb_wait_ mpi_alltoall_nb_wait__
#define mpi_alltoall_nb_test_ mpi_alltoall_nb_test__
#define mpi_alltoall_nb_free_ mpi_alltoall_nb_free__

#define mpi_abort_ mpi_abort__
#define mpi_allgatherv_ mpi_allgatherv__
#define mpi_alltoallv_ mpi_alltoallv__
#endif

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
