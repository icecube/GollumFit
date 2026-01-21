#ifndef COMPACTIO_H_
#define COMPACTIO_H_

#include <deque>
#include <fstream>
#include <boost/crc.hpp>
#include <PhysTools/tableio.h>

#include "analysisWeighting.h"
#include "build_info.h"
#include "Event.h"
#include "GollumMCSet.h"

/**
* @file CompactIO.h
* @brief File checks
*/

namespace gollumfit {

namespace dump {
  static void crc_write(std::ofstream& out, boost::crc_32_type& crc, const void* ptr, size_t n);
  static void crc_read(std::ifstream& in, boost::crc_32_type& crc, void* ptr, size_t n);
  void splatData(const std::string& filename, const uint32_t progChecksum, const std::deque<Event>& exp, const std::deque<Event>& sim);
  void unsplatData(const std::string& filename, const uint32_t progChecksum, std::deque<Event>& exp, std::deque<Event>& sim);
}

}

#endif /* COMPACTIO_H_ */
