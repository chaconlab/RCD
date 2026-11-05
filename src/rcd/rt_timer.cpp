#include "rt_timer.h"
//#include "stdio.h"

#if defined  _WIN32 || defined _WIN64
// Windows version of the timer
void timerReal::startTimer()
{
    QueryPerformanceCounter(&start);
}

void timerReal::stopTimer()
{
    QueryPerformanceCounter(&stop);
}

double timerReal::getElapsedTime()
{
	LARGE_INTEGER time;
	time.QuadPart = stop.QuadPart - start.QuadPart;
    return ((double)time.QuadPart /(double)frequency.QuadPart);
}



#else
	#ifdef __MACH__
		// OS X does not have clock_gettime, use clock_get_time
		// MAC version of the timer
		void timerReal::startTimer()
		{
			clock_serv_t cclock;
			mach_timespec_t mts;
			host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
			clock_get_time(cclock, &mts);
			mach_port_deallocate(mach_task_self(), cclock);
			start.tv_sec = mts.tv_sec;
			start.tv_nsec = mts.tv_nsec;
		}

		void timerReal::stopTimer()
		{
			clock_serv_t cclock;
			mach_timespec_t mts;
			host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
			clock_get_time(cclock, &mts);
			mach_port_deallocate(mach_task_self(), cclock);
			stop.tv_sec = mts.tv_sec;
			stop.tv_nsec = mts.tv_nsec;
		}

		double timerReal::getElapsedTime()
		{
			double elapsed;
			elapsed = (double)(stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
			return elapsed;
		}
	#else
		// Linux version of the timer
		void timerReal::startTimer()
		{
			clock_gettime(CLOCK_REALTIME,&start);
		}

		void timerReal::stopTimer()
		{
			clock_gettime(CLOCK_REALTIME,&stop);
		}

		double timerReal::getElapsedTime()
		{
			double elapsed;
//			fprintf(stderr,"stop.tv_sec= %ld  start.tv_sec=%ld  stop.tv_nsec=%ld  start.tv_nsec=%ld \n",stop.tv_sec,start.tv_sec,stop.tv_nsec,start.tv_nsec);
			elapsed = (double)(stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
			return elapsed;
		}
	#endif
#endif
