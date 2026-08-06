#include "utils.h"

/**
* @brief Calculate the center points of the intervals defined by a vector of edges.
* @param edges Rvalue reference to a vector of doubles representing interval edges.
* @return A vector of doubles representing the center points of the intervals.
*
* This function takes a vector of edges where each pair of adjacent values define the 
* boundaries of intervals. It calculates the center of each interval and returns a new vector
* containing these center values. The returned vector will have a size of one less than the input
* vector, as it calculates the centers between adjacent edges.
*
* @note This version of the function accepts an rvalue reference, which allows the use of temporaries.
*/
std::vector<double> edgesToCenters(const std::vector<double> && edges) {
    std::vector<double> res(edges.size()-1);
    for(unsigned int i=0; i<edges.size()-1; ++i) {
        res[i] = (edges[i] + edges[i+1])/2.0;
    }
    return res;
}

/**
* @brief Calculate the center points of the intervals defined by a vector of edges.
* @param edges Lvalue reference to a constant vector of doubles representing interval edges.
* @return A vector of doubles representing the center points of the intervals.
*
* This function operates like the rvalue reference overload, but is designed to work with
* lvalue references. It takes a vector of edges where each pair of adjacent values define
* the boundaries of intervals. It calculates the center of each interval and returns a new vector
* containing these center values. The returned vector will have a size of one less than the input
* vector, as it calculates the centers between adjacent edges.
*
* @note This version of the function accepts an lvalue reference, making it suitable for named variables.
*/
std::vector<double> edgesToCenters(const std::vector<double> & edges) {
    std::vector<double> res(edges.size()-1);
    for(unsigned int i=0; i<edges.size()-1; ++i) {
        res[i] = (edges[i] + edges[i+1])/2.0;
    }
    return res;
}

std::string get_current_datetime_string() {
    using namespace std::chrono;

    const auto now = system_clock::now();
    const std::time_t t = system_clock::to_time_t(now);

    std::tm tm{};
    gmtime_r(&t, &tm);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return oss.str();
}
