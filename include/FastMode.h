#ifndef FASTMODE_H_
#define FASTMODE_H_

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>
#include <PhysTools/histogram.h>
#include <PhysTools/likelihood/likelihood.h>


/**
* @file FastMode.h
* @brief Header file for fast mode functions.
* 
*/

namespace gollumfit {
namespace fastmode {
namespace detail {

    
//some useful type traits
/**
* @brief Type trait to remove std::reference_wrapper from a type if present.
* 
* @tparam T The type to be processed.
*/
template<typename T>
struct remove_reference_wrapper{ using type=T; };

/**
* @brief Specialization of remove_reference_wrapper for std::reference_wrapper types.
* 
* @tparam T The type wrapped by std::reference_wrapper.
*/
template<typename T>
struct remove_reference_wrapper<std::reference_wrapper<T>>{ using type=T; };

/**
* @brief Type trait to ensure a type is wrapped in std::reference_wrapper.
* 
* @tparam T The type to be wrapped.
*/
template<typename T>
struct ensure_reference_wrapper{ using type=typename std::reference_wrapper<T>; };                                                                                                                                                

/**
* @brief Specialization of ensure_reference_wrapper that maintains std::reference_wrapper type.
* 
* @tparam T The type already wrapped by std::reference_wrapper.
*/
template<typename T>
struct ensure_reference_wrapper<std::reference_wrapper<T>>{
    using type=typename std::reference_wrapper<T>;
};

/**
* @brief Functor to dereference an object if it is wrapped in std::reference_wrapper or return it as-is.
* 
* @tparam T The type of the object to be dereferenced.
*/
template<typename T>
struct dereference {
    T& operator()(T& t) {
        return t;
    }
};

/**
* @brief Specialization of dereference for std::reference_wrapper types.
* 
* @tparam T The type wrapped by std::reference_wrapper.
*/
template<typename T>
struct dereference<std::reference_wrapper<T>> {
    T & operator()(std::reference_wrapper<T> t) {
        return t.get();
    }
};


/**
* @brief Base case for recursive template function apply_to_tuple.
* 
* @tparam I The index of the tuple to apply the function to.
* @tparam f The function to apply.
* @tparam Tp Variadic template parameter pack for the tuple elements.
*/
template<std::size_t I = 0, typename f, typename... Tp> 
inline typename std::enable_if<I == sizeof...(Tp), void>::type
apply_to_tuple(f& do_thing, std::tuple<Tp...>& t)
{ }

/**
* @brief Recursive template function to apply a function to each element of a tuple.
* 
* @tparam I The index of the tuple to apply the function to.
* @tparam f The function to apply.
* @tparam Tp Variadic template parameter pack for the tuple elements.
*/
template<std::size_t I = 0, typename f, typename... Tp> 
inline typename std::enable_if<I < sizeof...(Tp), void>::type
apply_to_tuple(f& do_thing, std::tuple<Tp...>& t)                                                                           
{
    do_thing(std::get<I>(t));
    apply_to_tuple<I + 1, f, Tp...>(do_thing, t); 
}

// Forward declaration for class P
template<typename X>
class P;

/**
* @brief Function object to apply an action to all bins in a histogram.
* 
* @tparam Action The type of the action to apply.
*/
template<typename Action>
struct apply_to_bins {
    Action& action;
    /**
    * @brief Constructs an apply_to_bins object with the given action.
    * 
    * @param action The action to apply to each bin.
    */
    apply_to_bins(Action& action):action(action){}
    /**
    * @brief Applies the action to each bin in the histogram.
    * 
    * @tparam HType The type of the histogram.
    * @param hist The histogram to which the action will be applied.
    */
    template<typename HType>
    void operator()(HType& hist) {
        for(auto it=hist.begin(), end=hist.end(); it!=end; ++it) {
            //const auto& itc = static_cast<const phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(*it);
            //action(itc.entries());
            //P<decltype(*it)> itType;
            action(*it);
        }
    }
};

/**
* @brief Function object to apply an action to a histogram or tuple of histograms.
* 
* @tparam Action The type of the action to apply.
*/
template<typename Action>
struct apply_to_histogram {
    Action& action;
    /**
    * @brief Constructs an apply_to_histogram object with the given action.
    * 
    * @param action The action to apply to the histogram(s).
    */
    apply_to_histogram(Action& action):action(action){}
    template<typename T>
    struct apply_general {
        Action& action;
        apply_general(Action& action):action(action){}
        void operator()(T& t) {
            action(t);
        }
    };

    /**
    * @brief Specialization of the apply_general struct for applying an action to a std::tuple.
    *
    * This specialization of apply_general is designed to work with std::tuple objects. 
    * It takes an action and applies it to each element of the tuple using the apply_to_tuple function.
    *
    * @tparam Tp... The types of the elements contained within the std::tuple.
    */
    template<typename... Tp>
    struct apply_general<std::tuple<Tp...>> {
        Action& action; ///< Reference to the action to be applied to each element of the tuple.

        /**
        * @brief Constructs an apply_general object for a tuple with the given action.
        *
        * @param action The action to be applied to each element of the tuple.
        */
        apply_general(Action& action):action(action){}

        /**
        * @brief Applies the stored action to each element of the tuple.
        *
        * Iterates over each element of the tuple and applies the action using apply_to_tuple, 
        * which recursively processes each tuple element with the given action.
        *
        * @param t The tuple to which the action will be applied.
        */
        void operator()(std::tuple<Tp...>& t) {
            apply_to_tuple(action, t);
        }
    };

    /**
    * @brief Applies the action to a histogram or recursively to each histogram in a tuple.
    * 
    * @tparam T The type of the histogram or tuple.
    * @param t The histogram or tuple to which the action will be applied.
    */
    template<typename T>
    void operator()(T& t) {
        apply_general<T> ag(action);
        ag(t);
    }
    
};

/**
* @brief Function object to apply an action to all bins in a histogram.
*
* @tparam Action The type of the action to apply.
*/
template<typename Action>
struct apply_to_histogram_bins {
    Action& action;

    /**
    * @brief Constructs an apply_to_histogram_bins object with the given action.
    *
    * @param action The action to apply to each bin in the histogram.
    */
    apply_to_histogram_bins(Action& action):action(action){}

    /**
    * @brief Applies the action to each bin in the histogram.
    *
    * @tparam T The type of the histogram.
    * @param t The histogram to which the action will be applied.
    */
    template<typename T>
    void operator()(T& t) {
        apply_to_bins<Action> atb(action);
        apply_to_histogram<apply_to_bins<Action>> ath(atb);
        ath(t);
    }
};

/**
* @brief Function object to clear a bin by setting its value to zero.
*
* @tparam binType The data type of the bin's content.
*/
template<typename binType>
struct clear_bin {
    
    /**
    * @brief Sets the value of a bin to zero.
    *
    * @tparam T The type of the bin.
    * @param t The bin to be cleared.
    */
    template<typename T>
    void operator()(T& t) {
        t = binType(0.);
    }
};

/**
* @brief Type trait to extract the data type from a histogram or a tuple of histograms.
*
* @tparam T The type of the histogram or tuple.
*/
template<typename T> struct extract_data_type { using type=typename T::dataType;};

/**
* @brief Specialization of extract_data_type for tuples that extracts the data type 
* from the first element of the tuple.
*
* @tparam Tp... The types of the elements in the tuple.
*/
template<typename... Tp> struct extract_data_type<std::tuple<Tp...>> { using type=typename std::tuple_element<0, std::tuple<Tp...>>::type::dataType;};


/**
* @brief Function object to clear all bins in a histogram.
*
* @tparam T The type of the histogram.
*/
template<typename T>
struct clear_histogram {

    /**
    * @brief Clears all bins in the histogram by setting their values to zero.
    *
    * @param t The histogram to be cleared.
    */
    void operator()(T& t) {
        clear_bin<typename extract_data_type<T>::type> cb;
        apply_to_histogram_bins<decltype(cb)> clear(cb);
        clear(t);
    }
};

/**
* @brief Accumulates events into a meta histogram using a binner and combiner.
*
* @tparam MetaHistType The type of the meta histogram.
* @tparam MetaBinner The type of the binner function object.
* @tparam Combiner The type of the combiner function object.
* @tparam Event The type of the event data.
*/
template<typename MetaHistType, typename MetaBinner, typename Combiner, typename Event>
struct meta_event_accumulator {
    MetaHistType& meta_hist; ///< Reference to the meta histogram.
    MetaBinner binner;  ///< Binner function object instance.
    Combiner combiner;  ///< Combiner function object instance.

    std::deque<Event> meta_events;  ///< Queue of accumulated meta events.
    
    /**
    * @brief Constructs a meta_event_accumulator object.
    *
    * @param mh Reference to the meta histogram.
    * @param b The binner function object instance.
    * @param c The combiner function object instance.
    */
    meta_event_accumulator(MetaHistType& mh, MetaBinner b, Combiner c):meta_hist(mh), binner(b), combiner(c){}
    
    /**
    * @brief Processes a collection of events and populates the meta event queue.
    *
    * @tparam dataType The type of the data in the events vector.
    * @param events The vector of events to process.
    */
    template<typename dataType>
    void pops(std::vector<dataType>& events) {
        //Combiner& c = combiner;
        //std::deque<Event>& me = meta_events;
        
        auto combine_and_store = [&](const phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event> > m_events){
            if(m_events.entries().size() > 0) {
                Event e = this->combiner(m_events.entries());
                this->meta_events.push_back(e);
            }
        };
        for(const Event& e : events) {
            binner(meta_hist, e);
        }
        apply_to_histogram_bins<decltype(combine_and_store)> athb(combine_and_store);
        athb(meta_hist);
        clear_histogram<MetaHistType> clear;
        clear(meta_hist);
    }

    /**
    * @brief Function object to process elements of the histogram.
    * 
    * This struct provides the
    * functionality to process elements of a histogram by applying the pops
    * method of the meta_event_accumulator.
    * 
    * @tparam T The type of the element in the histogram to process.
    */
    template<typename T>
    struct ops {
        meta_event_accumulator<MetaHistType, MetaBinner, Combiner,Event>* this_;
        /**
        * @brief Constructs an ops object associated with a meta_event_accumulator.
        * 
        * @param t A pointer to the meta_event_accumulator that this ops will use to process elements.
        */
        ops(meta_event_accumulator<MetaHistType, MetaBinner, Combiner,Event>* t):this_(t){}
        /**
        * @brief Processes an element of the histogram using the associated meta_event_accumulator.
        * 
        * Casts the element to a phys_tools::likelihood::entryStoringBin containing reference_wrapped
        * const Events, extracts its entries, and passes them to meta_event_accumulator's pops method.
        * 
        * @param t The histogram element to process.
        */
        void operator()(T& t) {
            auto& itc = static_cast<const phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(t);
            auto events = itc.entries();
            //typedef decltype(events) eventsType;
            this_->pops(events);
        }
    };
    
    /**
    * @brief Specialization of the ops function object for vectors of data.
    * 
    * This specialization handles the case where the element to be processed is
    * a vector of data, allowing direct processing without casting.
    * 
    * @tparam dataType The type of data contained in the vector.
    */
    template<typename dataType>
    struct ops<std::vector<dataType>> {
        meta_event_accumulator<MetaHistType, MetaBinner, Combiner,Event>* this_;
        /**
        * @brief Constructs an ops object associated with a meta_event_accumulator.
        * 
        * @param t A pointer to the meta_event_accumulator that this ops will use to process elements.
        */
        ops(meta_event_accumulator<MetaHistType, MetaBinner, Combiner,Event>* t):this_(t){}
        /**
        * @brief Processes a vector of data using the associated meta_event_accumulator.
        * 
        * Passes the vector directly to the meta_event_accumulator's pops method.
        * 
        * @param events The vector of data to process.
        */
        void operator()(std::vector<dataType>& events) {
            this_->pops(events);
        }
    };

    /**
    * @brief Invokes the appropriate ops function object based on the type T.
    *
    * @tparam T The type of the histogram element.
    * @param t The element to process.
    */
    template<typename T>
    void operator()(T& t) {
        ops<T>(this)(t);
    }
};


} // namespace detail

/**
* @brief Processes the given histogram to accumulate events into a meta histogram.
* 
* This function takes a histogram and processes it using the provided meta histogram,
* binner, and combiner to accumulate the results into a set of meta events. These
* meta events are then returned as a deque.
* 
* @tparam Event The event type to be accumulated into the meta histogram.
* @tparam MetaHistType The type of the meta histogram that accumulates the meta events.
* @tparam MetaBinner The type of the binner function object used to bin events into the meta histogram.
* @tparam Combiner The type of the combiner function object used to combine events within each bin of the meta histogram.
* @tparam Histogram The type of the histogram that contains the events to be processed.
* 
* @param meta_hist A reference to the meta histogram where the accumulated events will be stored.
* @param meta_binner A reference to the binner function object which dictates how events are binned in the meta histogram.
* @param combiner An r-value reference to the combiner function object which dictates how events are combined within the bins of the meta histogram.
* @param hist A reference to the histogram containing the events to be processed.
* 
* @return A deque containing the accumulated meta events from processing the histogram.
*/
template<typename Event, typename MetaHistType, typename MetaBinner, typename Combiner, typename Histogram>
std::deque<Event> get_fastmode_events(MetaHistType& meta_hist, MetaBinner& meta_binner, Combiner&& combiner, Histogram& hist){
    detail::meta_event_accumulator<MetaHistType, MetaBinner, Combiner, Event> fm(meta_hist, meta_binner, combiner);
    detail::apply_to_histogram_bins<decltype(fm)> ath(fm);
    ath(hist);
    return fm.meta_events;
}

} // namespace fastmode
} // namespace gollumfit

#endif /* FASTMODE_H_ */
