#ifndef BDJR_H_
#define BDJR_H_

#include "pcmax.h"

namespace PCmax { namespace BDJR {

  /*
   * Solves instance I to schedule S and returns the makespan of S.
   */
  nat ComputeSchedule(const Instance& I, Schedule& S);

}}

#endif
