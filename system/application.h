#ifndef APPLICATION_H
#define APPLICATION_H

#include <aocommon/logger.h>

#include <sys/types.h>
#include <sys/wait.h>

#include <unistd.h>

namespace wsclean {

class Application {
 public:
  static void Run(const std::string& commandLine) {
    aocommon::Logger::Info << "Running: " << commandLine << '\n';
    const char* commandLineCStr = commandLine.c_str();
    int pid = vfork();
    switch (pid) {
      case -1:  // Error
        throw std::runtime_error(
            "Could not vfork() new process for executing command line "
            "application");
      case 0:  // Child
        execl("/bin/sh", "sh", "-c", commandLineCStr, NULL);
        _exit(127);
    }
    // Wait for process to terminate
    int pStatus;
    do {
      int pidReturn;
      do {
        pidReturn = waitpid(pid, &pStatus, 0);
      } while (pidReturn == -1 && errno == EINTR);
    } while (!WIFEXITED(pStatus) && !WIFSIGNALED(pStatus));
    if (WIFEXITED(pStatus)) {
      // all good
      // const int exitStatus = WEXITSTATUS(pStatus);
    } else {
      throw std::runtime_error(
          "Running command line application returned an error");
    }
  }
};

}  // namespace wsclean

#endif
