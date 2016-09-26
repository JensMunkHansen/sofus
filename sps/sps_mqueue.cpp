#include <sps/config.h>
#include <sps/sps_mqueue.hpp>
#include <sps/stdlib.h>

#ifndef MSGMAX
# define MSGMAX 8192
#endif

#ifdef HAVE_MQUEUE_H
int mq_clear(const char* qname)
{
  mqd_t mqd;
  struct mq_attr qattr, old_qattr;
  unsigned int prio = 0;
  char buf[MSGMAX];

  // TODO: Consider changing O_RDWR | O_CREAT to O_RDONLY. At the moment, we create if non-existing
  CallErrReturn(mqd = mq_open,
                (qname, O_RDWR | O_CREAT,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
                 NULL), EXIT_SUCCESS);

  // Remove messages on the queue
  mq_getattr(mqd, &qattr);

  if (qattr.mq_curmsgs !=0) {
    qattr.mq_flags = O_NONBLOCK;
    mq_setattr (mqd, &qattr, &old_qattr);
    while (mq_receive(mqd, &buf[0], qattr.mq_msgsize, &prio) != -1) {
    }
    if (errno != EAGAIN) {
      perror("mq_receive()");
      return EXIT_FAILURE;
    }
    // Restore attributes
    mq_setattr(mqd, &old_qattr, 0);
  }
  CallErrReturn(mq_close,(mqd),EXIT_FAILURE);
  return EXIT_SUCCESS;
}
#else
int mq_clear(const char* qname)
{
  return EXIT_FAILURE;
}
#endif

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
