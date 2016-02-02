/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// general Timer class; allows nested timers
// written by Peter Strazdins, Jun 13

#ifndef TIMER_INCLUDED
#define TIMER_INCLUDED

#include <mpi.h>
#include <string> //std::string
#include <assert.h>

class Timer {
 private:

  typedef std::pair<int, double> TimerData; // number of calls & acc time
  std::map<std::string, TimerData> timers;  
  std::map<std::string, double> startTime;
  std::map<std::string, int> paramSize;  
  // level of timer (0 = highest) used to control output
  std::map<std::string, int> timerLevel; 
   
  public:
  Timer() {}

  void start(std::string section, int size, int level) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // printf("%d: start timer: %s\n", rank, section.c_str());

    if (timers.count(section) == 0) { // new section
      paramSize[section] = size;
      timerLevel[section] = level;
      timers[section] = TimerData(1, 0.0);
    } else {
      if (startTime[section] != 0.0)
        printf("%d start timer: started timer %s %e\n", rank, 
	       section.c_str(), startTime[section]);
      fflush(stdout);
      assert(startTime[section] == 0.0); // previous use has ended
    }
    startTime[section] = MPI_Wtime();
  } //start()

  void stop(std::string section) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // printf("%d: stop timer: %s\n", rank,section.c_str());
   // check timer is in use
    assert(startTime.count(section) > 0 && startTime[section] > 0.0); 
    TimerData thisTimer = timers[section];
    timers[section] = TimerData(thisTimer.first + 1, 
				thisTimer.second + MPI_Wtime() -
				startTime[section]);
    startTime[section] = 0.0;
  } //stop()

  // dump timer summaries by process 0 of comm
  // note: undesirable to assume that timers is identical on each process
  void dump(MPI_Comm comm, int maxLevel) {
    const unsigned int KeySzEnd = 0; // signifies no more keys to send
    int rank, grank, gsize;
    MPI_Comm_rank(comm, &grank);
    MPI_Comm_size(comm, &gsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (grank == 0) { // get all timers on this rank
      for (std::map<std::string, TimerData>::iterator it = timers.begin();
	   it != timers.end(); it++) {
	TimerData thisTimer = it->second; 
	std::string key = it->first;
	assert(timerLevel.count(key) > 0);
	if (timerLevel[key] <= maxLevel) {
	  int keySz = key.length() + 1, nVals, nValsI = 1;        
	  double myTime = thisTimer.second / std::max(thisTimer.first, 1);
	  double timeSum, timeMax;
	  MPI_Bcast(&keySz, 1, MPI_INT, 0, comm);
	  MPI_Bcast((void *) key.c_str(), keySz, MPI_CHAR, 0, comm);
	  MPI_Reduce(&nValsI, &nVals, 1, MPI_INT, MPI_SUM, 0, comm);
	  MPI_Reduce(&myTime, &timeSum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	  MPI_Reduce(&myTime, &timeMax, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	  assert (paramSize.count(key) > 0);
	  printf("%d: timer %s: sz=%d, avg=%.2e, max=%.2e (%d/%d)\n",
		 rank, it->first.c_str(), paramSize[key], 
		 timeSum / nVals, timeMax, nVals, gsize);
	}
      } //for (it...) 
      MPI_Bcast((void *) &KeySzEnd, 1, MPI_INT, 0, comm);
    } else {
      unsigned int keySz;
      MPI_Bcast(&keySz, 1, MPI_INT, 0, comm);
      while (keySz != KeySzEnd) {
	std::string key; char keyBuf[256];  // section name
	assert(keySz < sizeof(keyBuf));
	MPI_Bcast(keyBuf, keySz, MPI_CHAR, 0, comm);
	key.assign(keyBuf);
	int hasKey = timers.count(key) > 0, tmpi;
	double tmp;
	double myTime = !hasKey? 0.0: timers[key].second /
	                              std::max(timers[key].first, 1);
	MPI_Reduce(&hasKey, &tmpi, 1, MPI_INT, MPI_SUM, 0, comm);
	MPI_Reduce(&myTime, &tmp, 1, MPI_DOUBLE,  MPI_SUM, 0, comm);
	MPI_Reduce(&myTime, &tmp, 1, MPI_DOUBLE,  MPI_MAX, 0, comm);	
	MPI_Bcast(&keySz, 1, MPI_INT, 0, comm);
      }      
    }   
  } //dump()

};

#endif /*TIMER_INCLUDED*/
