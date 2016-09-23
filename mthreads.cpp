/* 
TODO:
 - Implement using Sys V instead and use type identifier to select among threads
 - Make sure threads are not restarted
 - If POSIX queue, use name identifier inside messages
 - Use msgctl to get PID of last mq_send

 */
#define _NTHREADS 6

#include "common.h"

#include <csignal>
#include <pthread.h>
#include <mqueue.h>
#include <unistd.h>
#include <sched.h>
#include <linux/unistd.h>
#include <sys/syscall.h>

//#include <sys/msg.h>

#include <cstdio>
#include <cstdlib>

// Error handling
#include <cerrno>
#include <cstring> // strerror

#include "mexarg.h"


// Definitions from /usr/include/bits/local_lim.h

#define 	MQ_OPEN_MAX   8
//#define 	MQ_PRIO_MAX   32
#define 	MQ_BLOCK   0
#define 	MQ_NAME_MAX   80
#define 	MQ_MIN_MSG_PRIORITY   0
#define 	MQ_MAX_MSG_PRIORITY   MQ_PRIO_MAX
#define 	MAX_PQUEUES   4
#define 	MAX_MSGSIZE   50
#define 	MAX_MSGS   10
#define 	INVALID_PQUEUE   0
#define 	MQIDX   0

#define MSG_SIZE       8192

#define FailErrExit(msg,ret) {                                  \
    (void)fprintf(stderr, "FAILED: %s(errno=%d strerror=%s)\n", \
                  msg, errno, strerror(errno));                 \
    (void)fflush(stderr);                                       \
    return ret; }

#define FailErr(msg) {                                          \
    (void)fprintf(stderr, "FAILED: %s(errno=%d strerror=%s)\n", \
                  msg, errno, strerror(errno));                 \
    (void)fflush(stderr);}

#define Usage   "Usage error. see above"

static mqd_t mqd[_NTHREADS];
pthread_t pids[_NTHREADS];

static mqd_t mqd_master;

static bool threads_initialized = false;

typedef struct thread_arg {
  int cpu_id;
  int thread_id;
  mqd_t mqd_master;
  mqd_t mqd;
} thread_arg_t;

static void *thread_function(void *arg);
static pid_t gettid( void );

static int32_t nthreads = 0;

bool mthreads_work(mxArray* plhs[]) {

  bool retval = true;
  struct mq_attr attr;
  unsigned int prio = 0;
  char* buf = NULL;

  const int n_odim = 2;
  bool* p_blhs = NULL;
  mwSize o_dims[n_odim];

  o_dims[0] = mwSize(1);
  o_dims[1] = mwSize(1);

  Call(plhs[0] = mxCreateLogicalArray,
       (n_odim, (const mwSize*)o_dims));

  p_blhs = ((bool*) mxGetData(plhs[0]));

  for (int32_t i=0;i<nthreads;i++) {
    char name[9];
    sprintf(name,"/jmh%04d",i);
    CallErrExit(mqd[i] = mq_open,
		(name, O_WRONLY | O_CREAT ,
		 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
		 NULL), false);
    if (mq_send(mqd[i],"RUN",3,prio)==-1)
      perror ("mq_send()");
  }

  // Wait for threads to finish, listen to own msg queue
  CallErrExit(mqd_master = mq_open,
	      ("/jmh-master", O_RDONLY | O_CREAT ,
	       S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
	       NULL), false);

  mq_getattr(mqd_master, &attr);

  buf = (char*) malloc(attr.mq_msgsize*sizeof(char));

  size_t n_to_go = nthreads;


  // TODO: Verify thread id of sending process
  while (true) {
    mq_receive(mqd_master, &buf[0], attr.mq_msgsize, &prio);
    if (!strncmp(buf,"DONE",4)) {
      printf("Master: Thread done\n");
      n_to_go--;
    }
    else if (!strncmp(buf,"READY",5)) {
      printf("Master: Thread ready\n");
    }
    else {
      printf("Unexpected message: %s\n",buf);
      retval = false;
      break;
    }
    if (n_to_go==0)
      break;
  }
  
  printf("Master: All threads done\n");

  free(buf);

  *p_blhs = retval;

  return retval;
}

bool mthreads_end(mxArray* plhs[]) {
  bool retval = false;
  unsigned int prio = 0;

  const int n_odim = 2;
  bool* p_blhs = NULL;
  mwSize o_dims[n_odim];

  o_dims[0] = mwSize(1);
  o_dims[1] = mwSize(1);

  Call(plhs[0] = mxCreateLogicalArray,
       (n_odim, (const mwSize*)o_dims));

  p_blhs = ((bool*) mxGetData(plhs[0]));

  if (threads_initialized) {
    char name[9];
    
    printf("nthreads: %ld\n",nthreads);

    for (int32_t i=0;i<nthreads;i++) {
      sprintf(name,"/jmh%04d",i);
      CallErrExit(mqd[i] = mq_open,
		  (name, O_WRONLY | O_CREAT ,
		   S_IRUSR | S_IWUSR | S_IRGRP  | S_IROTH,
		   NULL), false);
      if (mq_send(mqd[i],"EXIT",4,prio)==-1)
	perror ("mq_send()");
      mq_close(mqd[i]);
    }
    threads_initialized = false;
    retval = true;
  }

  *p_blhs = retval;

  mexUnlock();

  return retval;
}

thread_arg_t thread_args[_NTHREADS];

bool mthreads_init(mxArray* plhs[],
		   const mxArray* mx_nthreads) {
  struct mq_attr attr, old_attr;
  pthread_attr_t p_attr;
  cpu_set_t cpuset;
  unsigned int prio;

  char* buf = NULL;
  bool* p_blhs = NULL;
  bool retval = false;

  int nproc = 1;

  const int n_odim = 2;
  mwSize o_dims[n_odim];

  mexLock();

  Call(mxIsScalarInt32, (mx_nthreads));
  nthreads = mxGetInt(mx_nthreads);

  o_dims[0] = mwSize(1);
  o_dims[1] = mwSize(1);
  
  Call(plhs[0] = mxCreateLogicalArray,
       (n_odim, (const mwSize*)o_dims));
  
  if ((!threads_initialized) && !(nthreads>_NTHREADS)) {
    printf("nthreads: %ld\n",nthreads);
    
    CPU_ZERO(&cpuset);
    CallErr(sched_getaffinity,
	    (gettid(), sizeof( cpu_set_t ), &cpuset));
  
    nproc = CPU_COUNT(&cpuset);

    printf("nproc: %lu\n",nproc);

    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    // Open message queue for main thread
    CallErrExit(mqd_master = mq_open,
		("/jmh-master",O_RDONLY | O_CREAT ,
		 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
		 NULL), EXIT_FAILURE);
    
    // Remove previous messages on the queue
    mq_getattr(mqd_master, &attr);
    
    //  printf("Allocating buffer: %ld\n",attr.mq_msgsize);
    buf = (char*) malloc((attr.mq_msgsize+1)*sizeof(char));
  
    if (attr.mq_curmsgs !=0) {
      printf("Removing old master messages");
      attr.mq_flags = O_NONBLOCK;
      mq_setattr (mqd_master, &attr, &old_attr);    
      while (mq_receive (mqd_master, &buf[0], attr.mq_msgsize, &prio) != -1) 
	printf(".");
      printf("\n");
      if (errno != EAGAIN) { 
	perror ("mq_receive()");
	exit (EXIT_FAILURE);
      }
      // Restore the attributes
      mq_setattr (mqd_master, &old_attr, 0);            
    }
    
    mq_close(mqd_master);
    
    pthread_attr_init(&p_attr);
    pthread_attr_setdetachstate(&p_attr, PTHREAD_CREATE_DETACHED);
    
    printf("Creating message queues:\n");
    for (int32_t i=0;i<nthreads;i++) {
      char name[9];
      sprintf(name,"/jmh%04d",i);
      CallErrExit(mqd[i] = mq_open,
		  (name, O_RDONLY | O_CREAT ,
		   S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
		   NULL), EXIT_FAILURE);
      printf("\t%s\n",name);
      mq_getattr(mqd[i], &attr);
      
      if (attr.mq_curmsgs !=0) {
	
	attr.mq_flags = O_NONBLOCK;
	mq_setattr (mqd[i], &attr, &old_attr);    
	// Remove previously unfinished jobs
	while (mq_receive (mqd[i], &buf[0], attr.mq_msgsize, &prio) != -1) 
	  printf(".");
	printf("\n");
	if (errno != EAGAIN) { 
	  perror ("mq_receive()");
	  exit (EXIT_FAILURE);
	}
	// Restore the attributes
	mq_setattr (mqd[i], &old_attr, 0);            
      }
      // Must be global if used by threads (after this function has exited)
      thread_args[i].cpu_id = ((int)i) % nproc;
      thread_args[i].thread_id = i;
      thread_args[i].mqd_master = mqd_master;
      thread_args[i].mqd        = mqd[i];
      mq_close(mqd[i]);
    }
    for (int32_t i=0;i<nthreads;i++) {
      CallErr(pthread_create,
	      (&pids[i], &p_attr, thread_function,
	       (void *)&thread_args[i]));
      printf("pid[%d] = %u\n",i,pids[i]);
    }
    
    threads_initialized = true;
    retval = true;
    // Wait for threads ready
    free(buf);
  }
  else {
    mexPrintf("Error initializing threads\n");
  }

  p_blhs = ((bool*) mxGetData(plhs[0]));
  *p_blhs = retval;
  
  return retval;
}

bool mthreads_mex(const int nlhs, mxArray *plhs[],
		  const int nrhs, const mxArray *prhs[]) {

  char type[256];
  //  char* attribute = NULL;

  /* Check number of arguments */    
  if (nrhs < 1 || nlhs < 0 || !mxIsChar(prhs[0])) {
    Fail(Usage);
  }
  Call(mxu_fixed_string, (type, 256, prhs[0], "1st argument"));

  if (!strcmp(type, "mthreads,init")) {
    if ((nrhs != 2) || (nlhs != 1)) {
      Fail(Usage);
    }
    Call(mthreads_init, (plhs, prhs[1]));
  }
  else if (!strcmp(type, "mthreads,work")) {
    if ((nrhs != 1) || (nlhs != 1)) {
      Fail(Usage);
    }
    Call(mthreads_work, (plhs));
  }
  else if (!strcmp(type, "mthreads,end")) {
    if ((nrhs != 1) || (nlhs != 1)) {
      Fail(Usage);
    }
    Call(mthreads_end, (plhs));
  }
  else {
    Fail(Usage);
  }
  return true;
}

void* thread_function(void *arg) {

  char buf[MSG_SIZE];

  thread_arg_t* threadarg = NULL;
  cpu_set_t set;
  struct mq_attr attr;
  unsigned int prio;

  threadarg = (thread_arg_t*) arg;

  int thread_id = threadarg->thread_id;

  // Lowest priority
  prio = 0;

  CPU_ZERO( &set );
  CPU_SET( threadarg->cpu_id, &set );

  CallErrExit(sched_setaffinity,
	      (gettid(), sizeof(cpu_set_t), &set), NULL);

  mqd_t mqd_master2;
  mqd_t mqd_mine;

  CallErrExit(mqd_master2 = mq_open,
	      ("/jmh-master",O_WRONLY,
	       S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
	       NULL), NULL);

  // Post Message
  if (mq_send(mqd_master2,"READY",5,prio)==-1)
    perror ("mq_send(): READY");

  printf("Thread %d: READY\n",thread_id);

  int state = 0; // Ready

  char name[9];
  sprintf(name,"/jmh%04d",thread_id);
  

  CallErrExit(mqd_mine = mq_open,
	      (name,O_RDONLY,
	       S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
	       NULL), NULL);

  mq_getattr (mqd_mine, &attr);

  attr.mq_flags &= ~O_NONBLOCK;
  mq_setattr(mqd_mine,&attr,0);
  int rett = 0;
  while(mq_receive (mqd_mine, &buf[0], attr.mq_msgsize, &prio) != -1) {
    printf ("Thread %d: Received a message with priority %d.\n", thread_id,prio);
  
    if (errno == EAGAIN) { 
      printf("Apparently non-blocking\n");
      pthread_exit(NULL);
    }

    if (!strcmp(buf,"RESET"))
      state = 0;
    else if (!strncmp(buf,"RUN",3))
      state = 1;
    else if (!strncmp(buf,"EXIT",4)) {
      state = 3;
      printf("Thread %d: EXIT\n",thread_id);
      break;
    }
    else {
      printf("Unknown message: %s", buf);
    }

    switch (state) {
    case 0:
    case 1:
      printf("Thread %d: Computing\n",thread_id);
      state = 2;
    case 2:
      printf("Thread %d: Done\n",thread_id);
      rett = mq_send(mqd_master2,"DONE",4,prio);
      if (rett==-1) {
	printf("Error sending DONE\n");
	perror ("mq_send()");
      }
      state = 0;
      break;
    }
  };

  mq_close(mqd_mine);
  mq_close(mqd_master2);

  pthread_exit(NULL);
}


void mexFunction(int nlhs,mxArray *plhs[],
                 int nrhs,const mxArray *prhs[]) {
  if (!nlhs && !nrhs) {
    return;
  }
  if (!mthreads_mex(nlhs, plhs, nrhs, prhs)) {
    mexErrMsgTxt("LORT");
  }
}

static pid_t gettid( void ) {
  pid_t pid;
  CallErr(pid = syscall, (__NR_gettid));
  return pid;
}
