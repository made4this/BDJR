#ifndef MULTIFIT_H_
#define MULTIFIT_H_

#include "pcmax.h"

namespace PCmax {

  namespace LPT {
    /*
     * Compute the makespan without computing the schedule.
     */
    nat ComputeMakespan(const Instance& I);

    /*
     * Schedules instance I on top of schedule S and returns the makespan.
     */
    nat ComputeSchedule(const Instance& I, Schedule& S);

  }

  namespace MF {

    /*
    * Compute the makespan without computing the schedule.
    */
    nat ComputeMakespan(const Instance& I);

    /*
    * Compute the schedule and return the makespan.
    */
    nat ComputeSchedule(const Instance& I, Schedule& S);

  }

  namespace DJMS {
    /*
    * Compute the makespan without computing the schedule.
    */
    nat ComputeMakespan(const Instance& I);

    /*
    * Compute the schedule and return the makespan.
    */
    nat ComputeSchedule(const Instance& I, Schedule& S);

  }

}

#endif
