#define _MULTITHREADED
#define _NTHREADS 8

#include <csignal>
#include <pthread.h>
#include <unistd.h>
#include <sched.h>
#include <linux/unistd.h>

#include <cstdio>
#include <cstdlib>

// Error handling
#include <cerrno>
#include <cstring> // strerror

#include "mex.h"

#define CallErr(fun, arg)  { if ((fun arg)<0)   \
      FailErr(#fun) }

#define CallErrExit(fun, arg, ret)  { if ((fun arg)<0)  \
      FailErrExit(#fun,ret) }

#define FailErrExit(msg,ret) {                                  \
    (void)fprintf(stderr, "FAILED: %s(errno=%d strerror=%s)\n", \
                  msg, errno, strerror(errno));                 \
    (void)fflush(stderr);                                       \
    return ret; }

#define FailErr(msg) {                                          \
    (void)fprintf(stderr, "FAILED: %s(errno=%d strerror=%s)\n", \
                  msg, errno, strerror(errno));                 \
    (void)fflush(stderr);}

typedef struct thread_arg {
  int cpu_id;
  int thread_id;
} thread_arg_t;

static volatile int quitflag;

// Experiment waiting for other threads
pthread_mutex_t lock   = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t waitloc = PTHREAD_COND_INITIALIZER;

pthread_attr_t attr;

pthread_t pids[_NTHREADS];

sigset_t mask; // No effect to use multiple
sigset_t oldmask;

static pid_t gettid( void );
static void *thread_function(void *arg);
static void *thread_irq_function(void *arg);

void mexFunction(int nlhs,mxArray *plhs[],
                 int nrhs,const mxArray *prhs[]) {

  cpu_set_t cpuset;
  int nproc = 1;
  int i;
  thread_arg_t thread_args[_NTHREADS];

  pthread_t irqpid;

  quitflag = 0;

  CPU_ZERO(&cpuset);
  CallErr(sched_getaffinity,
          (gettid(), sizeof( cpu_set_t ), &cpuset));
  
  nproc = CPU_COUNT(&cpuset);

  for (i=0 ; i < _NTHREADS ; i++) {
    thread_args[i].cpu_id = i % nproc;
    thread_args[i].thread_id = i;
  }

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  // We pray for no locks on buffers and setbuf will work, if not we
  // need to use filelock() on on FILE* access, tricky
  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  sigemptyset(&mask);
  sigaddset(&mask, SIGINT);
  sigaddset(&mask, SIGQUIT);
  sigaddset(&mask, SIGUSR2);
  CallErr(pthread_sigmask, (SIG_BLOCK, &mask, &oldmask));

  srand(time(NULL));

  thread_arg_t irq_thread_arg;
  irq_thread_arg.cpu_id = 0;
  
  CallErr(pthread_create,
          (&irqpid, &attr, thread_irq_function,
           (void*)&irq_thread_arg));

  for (i = 0; i < _NTHREADS; i++)
    CallErr(pthread_create,
            (&pids[i], &attr, thread_function,
             (void *)&thread_args[i]));

  for (i = 0; i < _NTHREADS; i++)
    CallErr(pthread_join, (pids[i], NULL));

  
  //    kill(getpid(), sig);
  raise(SIGUSR2);

  // Restore signal handler for main thread
  CallErr(sigprocmask, (SIG_SETMASK, &oldmask, NULL));
  
  // Clean up and exit 
  CallErr(pthread_attr_destroy, (&attr));

}

// Thread-safe, since pid is on stack and syscall is thread-safe
static pid_t gettid( void ) {
  pid_t pid;
  CallErr(pid = syscall, (__NR_gettid));
  return pid;
}

void *thread_irq_function(void *arg) {

  int signo;
  
  cpu_set_t set;

  thread_arg_t* thread_arg = (thread_arg_t*) arg;
  
  CPU_ZERO( &set );
  CPU_SET(thread_arg->cpu_id, &set);

  CallErr(sched_setaffinity,
          (gettid(), sizeof(cpu_set_t), &set));

  while(1) {
    CallErrExit(sigwait,(&mask, &signo),NULL);
    switch (signo) {
    case SIGINT:
    case SIGQUIT:
      pthread_mutex_lock(&lock);
      quitflag = 1;
      pthread_mutex_unlock(&lock);
      pthread_cond_signal(&waitloc);
      fprintf(stderr,"\n");
    case SIGUSR2:
      pthread_exit(NULL);
    default:
      printf("unexpected signal %d\n", signo);
      exit(1);
    }
  }
}

void *thread_function(void *arg) {
  thread_arg_t* threadarg;
  int thread_id;
  cpu_set_t set;

  threadarg = (thread_arg_t*) arg;

  thread_id = threadarg->thread_id+1;

  CPU_ZERO( &set );
  CPU_SET( threadarg->cpu_id, &set );
  CallErrExit(sched_setaffinity, (gettid(), sizeof(cpu_set_t), &set ),
              NULL);

  sigemptyset(&mask);
  sigaddset(&mask, SIGINT);
  sigaddset(&mask, SIGQUIT);
  sigaddset(&mask, SIGUSR2);

  CallErrExit(pthread_sigmask, (SIG_BLOCK, &mask, NULL), NULL);

  // While loop waiting for exit condition
  size_t i;

  for (i=0;i<10;i++) {
    sleep(rand() % 2);
    pthread_mutex_lock(&lock);
    if (quitflag != 0) {
      pthread_mutex_unlock(&lock);
      pthread_exit(NULL);
    }
    pthread_mutex_unlock(&lock);
  }
  printf("thread %d: finished\n",thread_id);
  pthread_exit(NULL);
}

