#include "compactIO.h"

/**
* @file CompactIO.cpp
* @brief File checks
*/

namespace gollumfit {

namespace dump {

/**
* @brief Computes a CRC32 checksum for the contents of a given file.
*
* This function opens the file specified by 'filename' in binary mode and computes
* a CRC32 checksum of the contents. It uses a buffer to process the file in chunks,
* which helps manage memory usage for large files.
*
* @param filename The path to the file for which the checksum is to be computed.
* @return uint32_t The computed CRC32 checksum of the file's contents.
* @throw std::runtime_error If the file cannot be opened for reading.
*/
uint32_t getFileChecksum(const std::string& filename){
	const std::streamsize bufferSize=1u<<16; //1MB
	
	std::ifstream file(filename, std::ios_base::binary);
	if(!file)
		throw std::runtime_error("Unable to open "+filename+" for reading");
	
	boost::crc_32_type fileCRC;
	char buffer[bufferSize];
	do{
		file.read(buffer, bufferSize);
		fileCRC.process_bytes(buffer, file.gcount());
	} while(file);
	
	return(fileCRC.checksum());
}


// Functions to simplify reading/writing the compact data files
static void crc_write(std::ofstream& out, 
                      boost::crc_32_type& crc,
                      const void* ptr, 
                      size_t n) {
    out.write(reinterpret_cast<const char*>(ptr), n);
    crc.process_bytes(ptr, n);
}

static void crc_read(std::ifstream& in, 
                     boost::crc_32_type& crc,
                     void* ptr, 
                     size_t n)
{
    in.read(reinterpret_cast<char*>(ptr), n);
    if (!in) throw std::runtime_error("Unexpected EOF / read failure");
    crc.process_bytes(ptr, n);
}

/**
* @brief Writes experimental and simulation event data to a file, along with checksums.
*
* The 'splatData' function creates a file specified by 'filename' and writes the
* CRC32 checksum of the program, followed by the experimental and simulation event
* data. Each data block is preceded by its size and followed by a checksum.
*
* @param filename The path to the file where the data will be written.
* @param progChecksum A CRC32 checksum representing the state of the program.
* @param exp A deque of Event objects representing experimental data.
* @param sim A deque of Event objects representing simulated data.
* @throw std::runtime_error If the file cannot be opened for writing.
*/
void splatData(const std::string& filename, const uint32_t progChecksum, const std::deque<Event>& exp, const std::deque<Event>& sim){
	std::ofstream datafile(filename);
	if(!datafile)
		throw std::runtime_error("Unable to open "+filename+" for writing");
	boost::crc_32_type fileCRC;

    const char format_header[8] = "FASTMC";
    const uint32_t format_version = 1;

    const std::string metadata =
        std::string("{")
        + "\"git_hash\":\"" + buildinfo::GIT_HASH + "\","
        + "\"git_dirty\":" + (buildinfo::GIT_DIRTY ? "true" : "false") + ","
        + "\"build_time_utc\":\"" + buildinfo::BUILD_TIME_UTC + "\","
        + "\"build_user\":\"" + std::string(buildinfo::BUILD_USER) + "\""
        + "}";
	const uint32_t meta_length = static_cast<uint32_t>(metadata.size());

    // write the different types of meta data
    crc_write(datafile, fileCRC, format_header, sizeof(format_header));
    crc_write(datafile, fileCRC, &format_version, sizeof(format_version));
    crc_write(datafile, fileCRC, &meta_length, sizeof(meta_length));
    crc_write(datafile, fileCRC, metadata.data(), metadata.size());

    crc_write(datafile, fileCRC, &progChecksum, sizeof(progChecksum));

    // write exp data
    uint64_t n = static_cast<uint64_t>(exp.size());
    crc_write(datafile, fileCRC, &n, sizeof(n));
    for (const Event& e : exp) {
        crc_write(datafile, fileCRC, &e, sizeof(e));
    }

    // write sim
    n = static_cast<uint64_t>(sim.size());
    crc_write(datafile, fileCRC, &n, sizeof(n));
    for (const Event& e : sim) {
        crc_write(datafile, fileCRC, &e, sizeof(e));
    }

    const uint32_t checksum = fileCRC.checksum();
    datafile.write(reinterpret_cast<const char*>(&checksum), sizeof(checksum));
}


/**
* @brief Reads experimental and simulation event data from a file, verifying checksums.
*
* The 'unsplatData' function opens a file specified by 'filename' and reads the
* stored program checksum, verifying it against an expected value. It then reads
* the experimental and simulation event data, checking the integrity of the data
* with CRC32 checksums. If the data exceeds a safety limit or checksums do not
* match, an exception is thrown.
*
* @param filename The path to the file from which to read the data.
* @param expectedChecksum The expected CRC32 checksum for program verification.
* @param exp A deque to be populated with Event objects representing experimental data.
* @param sim A deque to be populated with Event objects representing simulated data.
* @throw std::runtime_error If the file cannot be opened for reading, if the safety limit
*       for the number of events is exceeded, or if the checksums do not match the data.
*/
void unsplatData(const std::string& filename, const uint32_t expectedChecksum, std::deque<Event>& exp, std::deque<Event>& sim){
	std::ifstream datafile(filename);
	if(!datafile)
		throw std::runtime_error("Unable to open "+filename+" for reading");
	boost::crc_32_type fileCRC;

    // look for the "FASTMC" header
    char header[8] = {};
    datafile.read(header, sizeof(header));
    const char expected_header[8] = "FASTMC"; // expands to FASTMC\0\0
    const bool found_header = (std::memcmp(header, expected_header, sizeof(expected_header)) == 0);
	
    if (!found_header){ // Old FastMC format
        datafile.clear();
        datafile.seekg(0, std::ios::beg); // move back
        fileCRC.reset();

        uint32_t storedChecksum = 0;
        crc_read(datafile, fileCRC, &storedChecksum, sizeof(storedChecksum));
        if(storedChecksum!=expectedChecksum){
            std::ostringstream ss;
            ss << "Program checksum stored in " << filename << ", " << std::hex << storedChecksum << ", does not match expected (current) checksum, " << expectedChecksum;
            throw std::runtime_error(ss.str());
        }
        
        size_t size = 0;
        const size_t maxEvents=static_cast<size_t>(5e7); //as a vague sort-of safety check assume that there will never be more than 50 million events
        crc_read(datafile, fileCRC, &size, sizeof(size));
        if(size>maxEvents)
            throw std::runtime_error(filename+" claims to contain "+boost::lexical_cast<std::string>(size)
                                    +" experimental events, which is larger than the safety limit of "
                                    +boost::lexical_cast<std::string>(maxEvents));
        exp.resize(size);
        for(Event& e : exp){
            crc_read(datafile, fileCRC, &e, sizeof(e));
        }
        
        datafile.read((char*)&size,sizeof(size));
        fileCRC.process_bytes((char*)&size,sizeof(size));
        if(size>maxEvents)
            throw std::runtime_error(filename+" claims to contain "+boost::lexical_cast<std::string>(size)
                                    +" simulated events, which is larger than the safety limit of "
                                    +boost::lexical_cast<std::string>(maxEvents));
        sim.resize(size);
        for(Event& e : sim){
            crc_read(datafile, fileCRC, &e, sizeof(e));
        }
        
        uint32_t storedFileChecksum = 0;
        datafile.read(reinterpret_cast<char*>(&storedFileChecksum), sizeof(storedFileChecksum));
        if(storedFileChecksum!=fileCRC.checksum()){
            std::ostringstream ss;
            ss << filename << " appears to be corrupted: stored checksum, " << std::hex << storedFileChecksum << ", does not match recomputed checksum, " << fileCRC.checksum();
            throw std::runtime_error(ss.str());
        }
        return;
    } else { // New FastMC format
        fileCRC.process_bytes(header, sizeof(header));

        uint32_t format_version = 0;
        crc_read(datafile, fileCRC, &format_version, sizeof(format_version));

        uint32_t meta_length = 0;
        crc_read(datafile, fileCRC, &meta_length, sizeof(meta_length));

        // Read and ignore metadata for now (you can store it if you want)
        std::string metadata(meta_length, '\0');
        if (meta_length > 0)
            crc_read(datafile, fileCRC, metadata.data(), metadata.size());

        // Now the same as old format, but with version-dependent count type
        uint32_t storedProgChecksum = 0;
        crc_read(datafile, fileCRC, &storedProgChecksum, sizeof(storedProgChecksum));
        if (storedProgChecksum != expectedChecksum) {
            std::ostringstream ss;
            ss << "Program checksum stored in " << filename << ", " << std::hex << storedProgChecksum
               << ", does not match expected (current) checksum, " << expectedChecksum;
            throw std::runtime_error(ss.str());
        }

        const uint64_t maxEvents = 999'000'000ULL;

        uint64_t nExp = 0;
        crc_read(datafile, fileCRC, &nExp, sizeof(nExp));
        if (nExp > maxEvents)
            throw std::runtime_error(filename + " claims to contain " +
                                     boost::lexical_cast<std::string>(nExp) +
                                     " experimental events, which exceeds safety limit " +
                                     boost::lexical_cast<std::string>(maxEvents));

        exp.resize(static_cast<size_t>(nExp));
        for (Event& e : exp) {
            crc_read(datafile, fileCRC, &e, sizeof(e));
        }

        uint64_t nSim = 0;
        crc_read(datafile, fileCRC, &nSim, sizeof(nSim));
        if (nSim > maxEvents)
            throw std::runtime_error(filename + " claims to contain " +
                                     boost::lexical_cast<std::string>(nSim) +
                                     " simulated events, which exceeds safety limit " +
                                     boost::lexical_cast<std::string>(maxEvents));

        sim.resize(static_cast<size_t>(nSim));
        for (Event& e : sim) {
            crc_read(datafile, fileCRC, &e, sizeof(e));
        }

        // Final CRC (not included in CRC calculation)
        uint32_t storedFileChecksum = 0;
        datafile.read(reinterpret_cast<char*>(&storedFileChecksum), sizeof(storedFileChecksum));
        if (!datafile) throw std::runtime_error("Unable to read final checksum");
        if (storedFileChecksum != fileCRC.checksum()) {
            std::ostringstream ss;
            ss << filename << " appears corrupted: stored checksum " << std::hex << storedFileChecksum
               << " does not match recomputed checksum " << fileCRC.checksum();
            throw std::runtime_error(ss.str());
        }

        return;
    }

}
} // closing dump namespace

} // closing gollumfit namespace
