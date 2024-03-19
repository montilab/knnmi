#ifndef EIGEN_WARNINGS_DISABLED
#define EIGEN_WARNINGS_DISABLED

// This file has been gutted to avoid notes in the devtools::check() function.
// This has no impact on the library as this just disables compiler warnings,
// and none of these are raised by the knnmi compilation.


#endif
  
#else
  // warnings already disabled:
# ifndef EIGEN_WARNINGS_DISABLED_2
#  define EIGEN_WARNINGS_DISABLED_2
# elif defined(EIGEN_INTERNAL_DEBUGGING)
#  error "Do not include \"DisableStupidWarnings.h\" recursively more than twice!"
# endif
  


#endif // not EIGEN_WARNINGS_DISABLED
