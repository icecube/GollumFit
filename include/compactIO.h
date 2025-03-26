#ifndef COMPACTIO_H_
#define COMPACTIO_H_

#include <deque>
#include "Event.h"
#include "GollumMCSet.h"

/**
* @file CompactIO.h
* @brief File checks
*/

namespace gollumfit {

namespace dump {
  void splatData(const std::string& filename, const uint32_t progChecksum, const std::deque<Event>& exp, const std::deque<Event>& sim);
  void unsplatData(const std::string& filename, const uint32_t progChecksum, std::deque<Event>& exp, std::deque<Event>& sim);
}

}

#endif /* COMPACTIO_H_ */
