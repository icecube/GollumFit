#ifndef COMPACTIO_H_
#define COMPACTIO_H_

#include <vector>
#include "Event.h"
#include "GollumMCSet.h"

/**
* @file CompactIO.h
* @brief File checks
*/

namespace gollumfit {

namespace dump {
  void splatData(const std::string& filename, const uint32_t progChecksum, const std::vector<Event>& exp, const std::vector<Event>& sim);
  void unsplatData(const std::string& filename, const uint32_t progChecksum, std::vector<Event>& exp, std::vector<Event>& sim);
}

}

#endif /* COMPACTIO_H_ */
