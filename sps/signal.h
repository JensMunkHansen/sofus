/**
 * @file   signal.h
 * @author Jens Munk Hansen <jmh@jmhlaptop.parknet.dk>
 * @date   Sat Jul 18 15:50:10 2015
 *
 * @brief  Wrapper signal.h (C header)
 *
 * TODO: Export functions:
 *       sigaction, sigemptyset, sigaddset, sigprocmask, pthread_sigmask
 *
 * Needed are sigemptyset, sigaddset, pthread_sigmask
 */
/*
 *  This file is part of SOFUS.
 *
 *  SOFUS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SOFUS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SOFUS.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <sps/cenv.h>

#ifdef HAVE_CONFIG
# include <sps/config.h>
#endif

#ifdef __GNUC__
# include <signal.h>
#else

# include <limits.h>
# include <stdint.h>
# include <time.h>
# include <signal.h>
#if __WORDSIZE == 64
typedef uint64_t sigset_t;
#else
typedef uint32_t sigset_t;
#endif

// Alternative definition
//#define sigemptyset(what)   (*(what) = 0, 0)
//#define sigaddset(what,sig) (*(what) |= (1<<(sig)))

#define SA_NOCLDSTOP 1          /* Do not generate SIGCHLD when children
                       stop */
#define SA_SIGINFO   2          /* Invoke the signal catching function
                               with three arguments instead of one
                            */
#define SA_RESTART   0x10000000     /* Restart syscall on signal return */
#define SA_ONSTACK   0x20000000     /* Call signal handler on alternate
                                     signal stack provided by
                                     sigaltstack(2). */
#define SA_NODEFER   0x40000000     /* Don't automatically block the signal
                       when its handler is being executed  */
#define SA_RESETHAND 0x80000000     /* Reset to SIG_DFL on entry to handler */
#define SA_ONESHOT   SA_RESETHAND   /* Historical linux name */
#define SA_NOMASK    SA_NODEFER     /* Historical linux name */

/* Used internally by cygwin.  Included here to group everything in one place.
   Do not use.  */
#define _SA_INTERNAL_MASK 0xf000    /* bits in this range are internal */

#undef  MINSIGSTKSZ
#define MINSIGSTKSZ  8192
#undef  SIGSTKSZ
#define SIGSTKSZ    32768

#define SIGHUP  1   /* hangup */
#define SIGINT  2   /* interrupt */
#define SIGQUIT 3   /* quit */
#define SIGILL  4   /* illegal instruction (not reset when caught) */
#define SIGTRAP 5   /* trace trap (not reset when caught) */
#define SIGEMT  7   /* EMT instruction */
#define SIGFPE  8   /* floating point exception */
#define SIGKILL 9   /* kill (cannot be caught or ignored) */
#define SIGBUS  10  /* bus error */
#define SIGSEGV 11  /* segmentation violation */
#define SIGSYS  12  /* bad argument to system call */
#define SIGPIPE 13  /* write on a pipe with no one to read it */
#define SIGALRM 14  /* alarm clock */
#define SIGTERM 15  /* software termination signal from kill */
#define SIGURG  16  /* urgent condition on IO channel */
#define SIGSTOP 17  /* sendable stop signal not from tty */
#define SIGTSTP 18  /* stop signal from tty */
#define SIGCONT 19  /* continue a stopped process */
#define SIGCHLD 20  /* to parent on child stop or exit */
#define SIGCLD  20  /* System V name for SIGCHLD */
#define SIGTTIN 21  /* to readers pgrp upon background tty read */
#define SIGTTOU 22  /* like TTIN for output if (tp->t_local&LTOSTOP) */
#define SIGIO   23  /* input/output possible signal */
#define SIGPOLL SIGIO   /* System V name for SIGIO */
#define SIGXCPU 24  /* exceeded CPU time limit */
#define SIGXFSZ 25  /* exceeded file size limit */
#define SIGVTALRM 26    /* virtual time alarm */
#define SIGPROF 27  /* profiling time alarm */
#define SIGWINCH 28 /* window changed */
#define SIGLOST 29  /* resource lost (eg, record-lock lost) */
#define SIGPWR  SIGLOST /* power failure */
#define SIGUSR1 30  /* user defined signal 1 */
#define SIGUSR2 31  /* user defined signal 2 */

typedef unsigned long timer_t;
typedef int pid_t;

typedef uint32_t uid_t;

typedef void (*_sig_func_ptr)(int);

#define __SI_PAD_SIZE 32

#define __uint32_size(__x) (max(sizeof (__x) / sizeof (uint32_t), 1))

#define __SI_CYG_PAD (__SI_PAD_SIZE - __uint32_size (void *))

typedef union sigval {
  int sival_int;            /* integer signal value */
  void  *sival_ptr;         /* pointer signal value */
} sigval_t;

#pragma pack(push,4)
struct _sigcommune {
  uint32_t _si_code;
  void *_si_read_handle;
  void *_si_write_handle;
  void *_si_process_handle;
  union {
    int _si_fd;
    void *_si_pipe_fhandler;
    char *_si_str;
  };
};

#define __SI_PAD_SIZE 32
#ifdef __INSIDE_CYGWIN__
# ifndef max
#   define max(a,b) (((a) > (b)) ? (a) : (b))
# endif /*max*/
# define __uint32_size(__x) (max(sizeof (__x) / sizeof (uint32_t), 1))

/* This padding represents the elements of the last struct in siginfo_t,
   aligning the elements to the end to avoid conflicts with other struct
   members. */
# define __SI_CYG_PAD (__SI_PAD_SIZE - __uint32_size (void *))
#endif /*__INSIDE_CYGWIN__*/

typedef struct {
  int si_signo;                         /* signal number */
  int si_code;                          /* signal code */
  pid_t si_pid;                         /* sender's pid */
  uid_t si_uid;                         /* sender's uid */
  int si_errno;                         /* errno associated with signal */

  union {
    uint32_t __pad[__SI_PAD_SIZE];    /* plan for future growth */
    struct _sigcommune _si_commune;     /* cygwin ipc */
    struct {
      union {
        sigval_t si_sigval;             /* signal value */
        sigval_t si_value;              /* signal value */
      };
      struct {
        timer_t si_tid;                 /* timer id */
        unsigned int si_overrun;        /* overrun count */
      };
    };
    /* SIGCHLD */
    struct {
      int si_status;                    /* exit code */
      clock_t si_utime;                 /* user time */
      clock_t si_stime;                 /* system time */
    };

    void *si_addr;                      /* faulting address for core dumping
                                           signals */
    /* Cygwin internal fields */
#ifdef __INSIDE_CYGWIN__
    struct {
      __uint32_t __pad2[__SI_CYG_PAD];  /* Locate at end of struct */
      void *si_cyg;                     /* pointer to block containing
                                           cygwin-special info */
    };
#endif /*__INSIDE_CYGWIN__*/
  };
} siginfo_t;
#pragma pack(pop)

#define SIG_SETMASK 0   /* set mask with sigprocmask() */
#define SIG_BLOCK 1 /* set of signals to block */
#define SIG_UNBLOCK 2   /* set of signals to, well, unblock */



#define _EXFUN(name, proto)     __cdecl name proto

int _EXFUN(sigaction, (int, const struct sigaction *, struct sigaction *));

int _EXFUN(sigemptyset, (sigset_t *));

int _EXFUN(sigaddset, (sigset_t *, const int));

struct sigaction {
  union {
    _sig_func_ptr sa_handler; /* SIG_DFL, SIG_IGN, or pointer to a function */
    void  (*sa_sigaction) ( int, siginfo_t *, void * );
  };
  sigset_t sa_mask;
  int sa_flags;
};


extern "C" int
sigprocmask (int how, const sigset_t *set, sigset_t *oldset);

extern "C" int
pthread_sigmask (int operation, const sigset_t *set, sigset_t *old_set);

#endif


#ifndef SA_NOCLDWAIT
# define SA_NOCLDWAIT 0
#endif

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
