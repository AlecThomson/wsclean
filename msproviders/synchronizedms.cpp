#include "synchronizedms.h"

namespace wsclean {

std::set<std::string> SynchronizedMS::MSLock::_openFiles;

std::condition_variable SynchronizedMS::MSLock::_condition;

std::mutex SynchronizedMS::MSLock::_mutex;

}  // namespace wsclean
