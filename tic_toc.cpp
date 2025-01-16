// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "tic_toc.hpp"

#include <ctime>

namespace mfem
{

namespace internal
{

class StopWatch
{
private:
   std::clock_t user_time, start_utime;
   short Running;

public:
   StopWatch();
   inline void Clear();
   inline void Start();
   inline void Stop();
   inline double Resolution();
   inline double RealTime();
   inline double UserTime();
   inline double SystTime();
};

StopWatch::StopWatch()
{
   user_time = 0;
   Running = 0;
}


inline void StopWatch::Clear()
{
   user_time = 0;
   if (Running)
   {
      start_utime = std::clock();
   }
}

inline void StopWatch::Start()
{
   if (Running) { return; }
   start_utime = std::clock();
   Running = 1;
}

inline void StopWatch::Stop()
{
   if (!Running) { return; }
   user_time += ( std::clock() - start_utime );
   Running = 0;
}

inline double StopWatch::Resolution()
{
   return 1.0 / CLOCKS_PER_SEC; // potential resolution
}

inline double StopWatch::RealTime()
{
   return UserTime();
}

inline double StopWatch::UserTime()
{
   std::clock_t utime = user_time;
   if (Running)
   {
      utime += (std::clock() - start_utime);
   }
   return (double)(utime) / CLOCKS_PER_SEC;
}

inline double StopWatch::SystTime()
{
   return 0.0;
}

} // namespace internal


StopWatch::StopWatch() : M(new internal::StopWatch) { }

StopWatch::StopWatch(const StopWatch &sw)
   : M(new internal::StopWatch(*(sw.M))) { }

void StopWatch::Clear()
{
   M->Clear();
}

void StopWatch::Start()
{
   M->Start();
}

void StopWatch::Restart()
{
   Clear();
   Start();
}

void StopWatch::Stop()
{
   M->Stop();
}

double StopWatch::Resolution()
{
   return M->Resolution();
}

double StopWatch::RealTime()
{
   return M->RealTime();
}

double StopWatch::UserTime()
{
   return M->UserTime();
}

double StopWatch::SystTime()
{
   return M->SystTime();
}

StopWatch::~StopWatch() = default;


StopWatch tic_toc;

void tic()
{
   tic_toc.Clear();
   tic_toc.Start();
}

double toc()
{
   return tic_toc.UserTime();
}

}
