#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <boost/crc.hpp>
#include <boost/shared_ptr.hpp>
#include <hdf5.h>
#include <nuSQuIDS/tools.h>
#include <PhysTools/autodiff.h>

//Stuff to enable converting between types of shared_ptrs
/**
* @brief Anonymous namespace to prevent external linkage.
*/
namespace {

	/**
    * @brief A holder struct template that encapsulates a shared pointer.
    * 
    * @tparam SharedPointer The type of the shared pointer.
    */
	template<class SharedPointer> struct Holder {
		SharedPointer p;

		/**
        * @brief Construct a new Holder object from a shared pointer.
        * 
    	* @param p A constant reference to a shared pointer of type SharedPointer.
        */
		Holder(const SharedPointer &p) : p(p) {}

		/**
        * @brief Copy constructor for the Holder struct.
        * 
        * @param other A constant reference to another Holder object.
        */
		Holder(const Holder &other) : p(other.p) {}
		//Holder(Holder &&other) : p(std::move<SharedPointer>(other.p)) {}

		/**
        * @brief Function call operator overload, can be customized for specific deleter functionality.
        */
		void operator () (...) const {}
    };
	
	/**
    * @brief Converts a boost::shared_ptr to a std::shared_ptr.
    * 
    * @tparam T The object type managed by the shared pointer.
    * @param p A constant reference to a boost::shared_ptr of type T.
    * @return std::shared_ptr<T> A std::shared_ptr of type T that shares ownership with the input.
    */
	template<class T> std::shared_ptr<T> to_std_ptr(const boost::shared_ptr<T> &p) {
		typedef Holder<std::shared_ptr<T>> H;
		if(H* h = boost::get_deleter<H, T>(p))
			return h->p;
		return std::shared_ptr<T>(p.get(), Holder<boost::shared_ptr<T>>(p));
	}

	/**
    * @brief Converts a std::shared_ptr to a boost::shared_ptr.
    * 
    * @tparam T The object type managed by the shared pointer.
    * @param p A constant reference to a std::shared_ptr of type T.
    * @return boost::shared_ptr<T> A boost::shared_ptr of type T that shares ownership with the input.
    */
	template<class T> boost::shared_ptr<T> to_boost_ptr(const std::shared_ptr<T> &p){
		typedef Holder<boost::shared_ptr<T>> H;
		if(H* h = std::get_deleter<H, T>(p))
			return h->p;
		return boost::shared_ptr<T>(p.get(), Holder<std::shared_ptr<T>>(p));
	}
}

std::vector<double> edgesToCenters(const std::vector<double> && edges);
std::vector<double> edgesToCenters(const std::vector<double> & edges);

#endif //UTILS_H
