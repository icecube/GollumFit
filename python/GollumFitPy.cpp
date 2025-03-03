#include <boost/python.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/overloads.hpp>
#include <functional>
#include "container_conversions.h"
#include "GollumFit.h"

#include <numpy/ndarrayobject.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

#include <nuSQuIDS/marray.h>

using namespace boost::python;
namespace bp = boost::python;
namespace GF = gollumfit;
namespace nsq = nusquids;

template<class T>
struct VecToList
{
  static PyObject* convert(const std::vector<T>& vec){
    boost::python::list* l = new boost::python::list();
    for(size_t i =0; i < vec.size(); i++)
      (*l).append(vec[i]);

    return l->ptr();
  }
};

template<class T>
struct DeqToList
{
  static PyObject* convert(const std::deque<T>& vec){
    boost::python::list* l = new boost::python::list();
    for(size_t i =0; i < vec.size(); i++)
      (*l).append(vec[i]);

    return l->ptr();
  }
};

void ListToVec(std::vector<unsigned int> &ret, bp::list l){
  for(int i=0;i<bp::len(l);i++)
    ret.push_back(bp::extract<unsigned int>(l[i]));
}

// converting marray to numpy array and back
template<unsigned int DIM>
struct marray_to_numpyarray {
  static PyObject* convert( nsq::marray<double,DIM> const & iarray){
    // get the data from the marray
    double * data = iarray.size() ? const_cast<double*>(iarray.get_data()) : static_cast<double*>(NULL);
    // construct numpy object
    npy_intp size[DIM];
    for(unsigned int i = 0; i < DIM; i++)
      size[i] = iarray.extent(i);
    PyArrayObject * pyObj = (PyArrayObject*) PyArray_SimpleNew(DIM,size,PyArray_DOUBLE);
    memcpy(pyObj->data, data, sizeof(double) * iarray.size());

    return PyArray_Return(pyObj);
  }
};

template<typename T,unsigned int DIM>
static nsq::marray<T,DIM> numpyarray_to_marray(PyObject * iarray, NPY_TYPES type_num){
  // es un array de numpy
  if (! PyArray_Check(iarray) )
  {
    PyErr_SetString(PyExc_TypeError, "numpyarray_to_marray: Input is not a numpy array.");
    boost::python::throw_error_already_set();
  }
  // si es que fuera un array de numpy castearlo
  //PyArrayObject* numpy_array = (PyArrayObject*) iarray;
  // lets get the contiguos C-style array
  PyArrayObject* numpy_array = PyArray_GETCONTIGUOUS((PyArrayObject*)iarray);

  // revisemos que los tipos del array sean dobles o que
  if ( PyArray_DESCR(numpy_array)->type_num != type_num )
  {
    if ( PyArray_DESCR(numpy_array)->type_num == NPY_LONG &&
         PyArray_ITEMSIZE(numpy_array) == 4 && type_num == NPY_INT)
    {
      // numpy on 32 bits sets numpy.int32 to NPY_LONG. So its all ok.
    }
    else
    {
      PyErr_SetString(PyExc_TypeError, "numpyarray_to_marray: numpy type is not the same as the input array type.");
      boost::python::throw_error_already_set();
    }
  }

  // arrays vacios
  if (PyArray_SIZE(numpy_array) == 0){
      PyErr_SetString(PyExc_TypeError,"numpyarray_to_marray: empty numpy array.");
      boost::python::throw_error_already_set();
  }

  // create numpy iterator
  NpyIter* iter = NpyIter_New(numpy_array, NPY_ITER_READONLY|
                             NPY_ITER_EXTERNAL_LOOP|
                             NPY_ITER_REFS_OK,
                             NPY_KEEPORDER, NPY_NO_CASTING,
                             NULL);

  unsigned int array_dim = PyArray_NDIM(numpy_array);
  assert(DIM == array_dim && "No matching dimensions.");

  // get numpy array shape and create marray object
#ifdef NPY_1_7_API_VERSION
  npy_intp* array_shape = PyArray_SHAPE(numpy_array);
#else
  npy_intp* array_shape = PyArray_DIMS(numpy_array);
#endif
  std::vector<size_t> dimensions;
  for(unsigned int i = 0; i < array_dim; i++)
    dimensions.push_back(array_shape[i]);

  // construct output object
  nsq::marray<T,DIM> oarray;
  oarray.resize(dimensions);
  auto it = oarray.begin();

  NpyIter_IterNextFunc *iternext = NpyIter_GetIterNext(iter, NULL);
  char** dataptr = NpyIter_GetDataPtrArray(iter);
  npy_intp* strideptr = NpyIter_GetInnerStrideArray(iter);
  npy_intp* sizeptr = NpyIter_GetInnerLoopSizePtr(iter);
  npy_intp iop, nop = NpyIter_GetNOp(iter);

  // magic to make the int work
  bool magic = false;
  if ( type_num == NPY_INT or type_num == NPY_LONG )
    magic = true;

  do{
    char* data = *dataptr;
    npy_intp count = *sizeptr;
    npy_intp stride = *strideptr;

    while (count--)
    {
      for (iop = 0; iop < nop; ++iop, data+=stride){
        if (magic)
          *it++ = *(T*)(reinterpret_cast<int*>(data));
        else
          *it++ = *(T*)(data);
      }
    }
  } while(iternext(iter));

  NpyIter_Deallocate(iter);

  return oarray;
}

template<typename key_t, typename val_t>
struct Dict2Map {

    /// The type of the map we convert the Python dict into
    typedef std::map<key_t, val_t> map_t;

    /// constructor
    /// registers the converter with the Boost.Python runtime
    Dict2Map() {
        bp::converter::registry::push_back(
            &convertible,
            &construct,
            bp::type_id<map_t>()
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
            , &bp::converter::wrap_pytype<&PyDict_Type>::get_pytype
#endif
        );
    }

    /// Check if conversion is possible
    static void* convertible(PyObject* objptr) {
        return PyDict_Check(objptr)? objptr: nullptr;
    }

    /// Perform the conversion
    static void construct(
        PyObject* objptr,
        bp::converter::rvalue_from_python_stage1_data* data
    ) {
        // convert the PyObject pointed to by `objptr` to a bp::dict
        bp::handle<> objhandle{ bp::borrowed(objptr) };   // "smart ptr"
        bp::dict d{ objhandle };

        // get a pointer to memory into which we construct the map
        // this is provided by the Python runtime
        void* storage =
            reinterpret_cast<
                bp::converter::rvalue_from_python_storage<map_t>*
            >(data)->storage.bytes;

        // placement-new allocate the result
        new(storage) map_t{};

        // iterate over the dictionary `d`, fill up the map `m`
        //map_t& m;
        //m = *(static_cast<map_t *>(storage));
        map_t& m = *(static_cast<map_t *>(storage));
        bp::list keys{ d.keys() };
        int keycount{ static_cast<int>(bp::len(keys)) };
        for (int i = 0; i < keycount; ++i) {
            // get the key
            bp::object keyobj{ keys[i] };
            bp::extract<key_t> keyproxy{ keyobj };
            if (! keyproxy.check()) {
                PyErr_SetString(PyExc_KeyError, "Bad key type");
                bp::throw_error_already_set();
            }
            key_t key = keyproxy();

            // get the corresponding value
            bp::object valobj{ d[keyobj] };
            bp::extract<val_t> valproxy{ valobj };
            if (! valproxy.check()) {
                PyErr_SetString(PyExc_ValueError, "Bad value type");
                bp::throw_error_already_set();
            }
            val_t val = valproxy();
            m[key] = val;
        }

        // remember the location for later
        data->convertible = storage;
    }
};

// auxiliary wrapper functions // evil // Carlos

static double wrap_SetData(GF::GollumFit* st, PyObject * array){
  if (! PyArray_Check(array) )
  {
    throw std::runtime_error("GollumFit::Error:Input array is not a numpy array.");
  }

  PyArrayObject* numpy_array = (PyArrayObject*)array;
  unsigned int array_dim = PyArray_NDIM(numpy_array);
  NPY_TYPES type = (NPY_TYPES) PyArray_DESCR(numpy_array)->type_num;

  // things i think can cast ok to doubles
  if (!( type == NPY_LONG or type == NPY_INT or type == NPY_SHORT or type == NPY_FLOAT or
      type == NPY_DOUBLE or type == NPY_LONGDOUBLE or type == NPY_CFLOAT or type == NPY_CDOUBLE))
    throw std::runtime_error("GollumFit::Error:Input numpy array cannot be meaninfully casted into double.");

  if ( array_dim == 2 ) {
    nsq::marray<double,2> state = numpyarray_to_marray<double,2>(array, type);
    double value = st->SetData(state);
    return value;
  } else
    throw std::runtime_error("GollumFit::Input array has wrong dimenions.");
}


static Dict2Map<unsigned int, double> reg{};

BOOST_PYTHON_MODULE(GollumFitPy)
{
  // import numpy array definitions
  import_array();
  import_ufunc();

  enum_<GF::FluxComponent>("FluxComponent")
    .value("atmConv",GF::FluxComponent::atmConv)
    .value("atmPrompt",GF::FluxComponent::atmPrompt)
    .value("diffuseAstro",GF::FluxComponent::diffuseAstro)
  ;

  class_<GF::FitResult, std::shared_ptr<GF::FitResult> >("FitResult",init<>())
    .def_readwrite("params",&GF::FitResult::params)
    .def_readwrite("likelihood",&GF::FitResult::likelihood)
    .def_readwrite("aux_likelihood",&GF::FitResult::aux_likelihood)
    .def_readwrite("nEval",&GF::FitResult::nEval)
    .def_readwrite("nGrad",&GF::FitResult::nGrad)
    .def_readwrite("succeeded",&GF::FitResult::succeeded)
  ;

  class_<GF::hist_marray>("hist_marray");

  class_<GF::DataPaths, boost::noncopyable,std::shared_ptr<GF::DataPaths> >("DataPaths",init<>())
    .def_readwrite("compact_file_path",&GF::DataPaths::compact_file_path)
    .def_readwrite("neutrino_cc_xs_spline_path",&GF::DataPaths::neutrino_cc_xs_spline_path)
    .def_readwrite("antineutrino_cc_xs_spline_path",&GF::DataPaths::antineutrino_cc_xs_spline_path)
    .def_readwrite("neutrino_nc_xs_spline_path",&GF::DataPaths::neutrino_nc_xs_spline_path)
    .def_readwrite("antineutrino_nc_xs_spline_path",&GF::DataPaths::antineutrino_nc_xs_spline_path)
    .def_readwrite("diff_neutrino_cc_xs_spline_path",&GF::DataPaths::diff_neutrino_cc_xs_spline_path)
    .def_readwrite("diff_antineutrino_cc_xs_spline_path",&GF::DataPaths::diff_antineutrino_cc_xs_spline_path)
    .def_readwrite("diff_neutrino_nc_xs_spline_path",&GF::DataPaths::diff_neutrino_nc_xs_spline_path)
    .def_readwrite("diff_antineutrino_nc_xs_spline_path",&GF::DataPaths::diff_antineutrino_nc_xs_spline_path)
    .def_readwrite("data_path",&GF::DataPaths::data_path)
    .def_readwrite("mc_path",&GF::DataPaths::mc_path)
    .def_readwrite("conventional_nusquids_atmospheric_file",&GF::DataPaths::conventional_nusquids_atmospheric_file)
    .def_readwrite("prompt_nusquids_atmospheric_file",&GF::DataPaths::prompt_nusquids_atmospheric_file)
    .def_readwrite("astro_nusquids_file",&GF::DataPaths::astro_nusquids_file)
    .def_readwrite("domeff_spline_path",&GF::DataPaths::domeff_spline_path)
    .def_readwrite("holeice_spline_path",&GF::DataPaths::holeice_spline_path)
    .def_readwrite("attenuation_spline_path",&GF::DataPaths::attenuation_spline_path)
    .def_readwrite("ice_gradient_spline_path",&GF::DataPaths::ice_gradient_spline_path)
    .def_readwrite("atmospheric_density_spline_path",&GF::DataPaths::atmospheric_density_spline_path)
    .def_readwrite("atmospheric_kaonlosses_spline_path",&GF::DataPaths::atmospheric_kaonlosses_spline_path)
    .def_readwrite("hadronic_spline_path",&GF::DataPaths::hadronic_spline_path)
    .def_readwrite("cosmic_ray_spline_path",&GF::DataPaths::cosmic_ray_spline_path)
  ;

  class_<GF::SteeringParams, boost::noncopyable,std::shared_ptr<GF::SteeringParams> >("SteeringParams",init<>())
    .def_readwrite("minFitEnergy",&GF::SteeringParams::minFitEnergy, "Minimum energy in the fit")
    .def_readwrite("maxFitEnergy",&GF::SteeringParams::maxFitEnergy, "Maximum energy in the fit")
    .def_readwrite("minCosth",&GF::SteeringParams::minCosth)
    .def_readwrite("maxCosth",&GF::SteeringParams::maxCosth)
    .def_readwrite("logEbinEdge",&GF::SteeringParams::logEbinEdge)
    .def_readwrite("logEbinWidth",&GF::SteeringParams::logEbinWidth)
    .def_readwrite("cosThbinEdge",&GF::SteeringParams::cosThbinEdge)
    .def_readwrite("cosThbinWidth",&GF::SteeringParams::cosThbinWidth)
    .def_readwrite("ice_gradient_filename",&GF::SteeringParams::ice_gradient_filename)
    .def_readwrite("active_hadronic_parameters",&GF::SteeringParams::active_hadronic_parameters)
    .def_readwrite("active_cosmicray_parameters",&GF::SteeringParams::active_cosmicray_parameters)
    .def_readwrite("uncertaintyModSigmaOverMu",&GF::SteeringParams::uncertaintyModSigmaOverMu, "Adds a relative error to the effective likelihood, where the input parameter is the relative error. The relative error is the same for all the bins.")
    .def_readwrite("model_label",&GF::SteeringParams::model_label)
    .def_readwrite("simToLoad",&GF::SteeringParams::simToLoad)
    .def_readwrite("evalThreads",&GF::SteeringParams::evalThreads)
    .def_readwrite("readCompact",&GF::SteeringParams::readCompact)
    .def_readwrite("fullLivetime",&GF::SteeringParams::fullLivetime)
    .def_readwrite("change_tol",&GF::SteeringParams::change_tol)
    .def_readwrite("grad_tol",&GF::SteeringParams::grad_tol)
    .def_readwrite("energyName",&GF::SteeringParams::energyName)
    .def_readwrite("selectionStart",&GF::SteeringParams::selectionStart)
  ;


  class_<GF::FitParameters, boost::noncopyable,std::shared_ptr<GF::FitParameters> >("FitParameters",init<>())
    .def_readwrite("convNorm",&GF::FitParameters::convNorm)
    .def_readwrite("promptNorm",&GF::FitParameters::promptNorm)
    .def_readwrite("zenithCorrection",&GF::FitParameters::zenithCorrection)
    .def_readwrite("kaonLosses",&GF::FitParameters::kaonLosses)
    .def_readwrite("hadronicHEkp",&GF::FitParameters::hadronicHEkp)
    .def_readwrite("hadronicHEkm",&GF::FitParameters::hadronicHEkm)
    .def_readwrite("hadronicVHE1pip",&GF::FitParameters::hadronicVHE1pip)
    .def_readwrite("hadronicVHE1pim",&GF::FitParameters::hadronicVHE1pim)
    .def_readwrite("hadronicVHE3kp",&GF::FitParameters::hadronicVHE3kp)
    .def_readwrite("hadronicVHE3km",&GF::FitParameters::hadronicVHE3km)
    .def_readwrite("hadronicVHE3pip",&GF::FitParameters::hadronicVHE3pip)
    .def_readwrite("hadronicVHE3pim",&GF::FitParameters::hadronicVHE3pim)
    .def_readwrite("hadronicVHE3p",&GF::FitParameters::hadronicVHE3p)
    .def_readwrite("hadronicVHE3n",&GF::FitParameters::hadronicVHE3n)
    .def_readwrite("cosmicRay1",&GF::FitParameters::cosmicRay1)
    .def_readwrite("cosmicRay2",&GF::FitParameters::cosmicRay2)
    .def_readwrite("cosmicRay3",&GF::FitParameters::cosmicRay3)
    .def_readwrite("cosmicRay4",&GF::FitParameters::cosmicRay4)
    .def_readwrite("cosmicRay5",&GF::FitParameters::cosmicRay5)
    .def_readwrite("cosmicRay6",&GF::FitParameters::cosmicRay6)
    .def_readwrite("icegrad0",&GF::FitParameters::icegrad0)
    .def_readwrite("icegrad1",&GF::FitParameters::icegrad1)
    .def_readwrite("icegrad2",&GF::FitParameters::icegrad2)
    .def_readwrite("icegrad3",&GF::FitParameters::icegrad3)
    .def_readwrite("icegrad4",&GF::FitParameters::icegrad4)
    .def_readwrite("icegrad5",&GF::FitParameters::icegrad5)
    .def_readwrite("icegrad6",&GF::FitParameters::icegrad6)
    .def_readwrite("icegrad7",&GF::FitParameters::icegrad7)
    .def_readwrite("icegrad8",&GF::FitParameters::icegrad8)
    .def_readwrite("domEfficiency",&GF::FitParameters::domEfficiency)
    .def_readwrite("holeiceForward",&GF::FitParameters::holeiceForward)
    .def_readwrite("astroNorm",&GF::FitParameters::astroNorm)
    .def_readwrite("astroDeltaGamma",&GF::FitParameters::astroDeltaGamma)
    .def_readwrite("astroDeltaGammaSec",&GF::FitParameters::astroDeltaGammaSec)
    .def_readwrite("astroPivot",&GF::FitParameters::astroPivot)
    .def_readwrite("NeutrinoAntineutrinoRatio",&GF::FitParameters::NeutrinoAntineutrinoRatio)
    .def_readwrite("nuxs",&GF::FitParameters::nuxs)
    .def_readwrite("nubarxs",&GF::FitParameters::nubarxs)
  ;

  class_<GF::FitParametersBound, boost::noncopyable,std::shared_ptr<GF::FitParametersBound> >("FitParametersBound",init<>())
    .def_readwrite("convNormMin",&GF::FitParametersBound::convNormMin)
    .def_readwrite("convNormMax",&GF::FitParametersBound::convNormMax)    
    .def_readwrite("promptNormMin",&GF::FitParametersBound::promptNormMin)
    .def_readwrite("promptNormMax",&GF::FitParametersBound::promptNormMax)    
    .def_readwrite("zenithCorrectionMin",&GF::FitParametersBound::zenithCorrectionMin)
    .def_readwrite("zenithCorrectionMax",&GF::FitParametersBound::zenithCorrectionMax)    
    .def_readwrite("kaonLossesMin",&GF::FitParametersBound::kaonLossesMin)
    .def_readwrite("kaonLossesMax",&GF::FitParametersBound::kaonLossesMax)    
    .def_readwrite("hadronicHEkpMin",&GF::FitParametersBound::hadronicHEkpMin)
    .def_readwrite("hadronicHEkpMax",&GF::FitParametersBound::hadronicHEkpMax)    
    .def_readwrite("hadronicHEkmMin",&GF::FitParametersBound::hadronicHEkmMin)
    .def_readwrite("hadronicHEkmMax",&GF::FitParametersBound::hadronicHEkmMax)    
    .def_readwrite("hadronicVHE1pipMin",&GF::FitParametersBound::hadronicVHE1pipMin)
    .def_readwrite("hadronicVHE1pipMax",&GF::FitParametersBound::hadronicVHE1pipMax)    
    .def_readwrite("hadronicVHE1pimMin",&GF::FitParametersBound::hadronicVHE1pimMin)
    .def_readwrite("hadronicVHE1pimMax",&GF::FitParametersBound::hadronicVHE1pimMax)    
    .def_readwrite("hadronicVHE3kpMin",&GF::FitParametersBound::hadronicVHE3kpMin)
    .def_readwrite("hadronicVHE3kpMax",&GF::FitParametersBound::hadronicVHE3kpMax)    
    .def_readwrite("hadronicVHE3kmMin",&GF::FitParametersBound::hadronicVHE3kmMin)
    .def_readwrite("hadronicVHE3kmMax",&GF::FitParametersBound::hadronicVHE3kmMax)    
    .def_readwrite("hadronicVHE3pipMin",&GF::FitParametersBound::hadronicVHE3pipMin)
    .def_readwrite("hadronicVHE3pipMax",&GF::FitParametersBound::hadronicVHE3pipMax)    
    .def_readwrite("hadronicVHE3pimMin",&GF::FitParametersBound::hadronicVHE3pimMin)
    .def_readwrite("hadronicVHE3pimMax",&GF::FitParametersBound::hadronicVHE3pimMax)    
    .def_readwrite("hadronicVHE3pMin",&GF::FitParametersBound::hadronicVHE3pMin)
    .def_readwrite("hadronicVHE3pMax",&GF::FitParametersBound::hadronicVHE3pMax)    
    .def_readwrite("hadronicVHE3nMin",&GF::FitParametersBound::hadronicVHE3nMin)
    .def_readwrite("hadronicVHE3nMax",&GF::FitParametersBound::hadronicVHE3nMax)    
    .def_readwrite("cosmicRay1Min",&GF::FitParametersBound::cosmicRay1Min)
    .def_readwrite("cosmicRay1Max",&GF::FitParametersBound::cosmicRay1Max)    
    .def_readwrite("cosmicRay2Min",&GF::FitParametersBound::cosmicRay2Min)
    .def_readwrite("cosmicRay2Max",&GF::FitParametersBound::cosmicRay2Max)    
    .def_readwrite("cosmicRay3Min",&GF::FitParametersBound::cosmicRay3Min)
    .def_readwrite("cosmicRay3Max",&GF::FitParametersBound::cosmicRay3Max)    
    .def_readwrite("cosmicRay4Min",&GF::FitParametersBound::cosmicRay4Min)
    .def_readwrite("cosmicRay4Max",&GF::FitParametersBound::cosmicRay4Max)    
    .def_readwrite("cosmicRay5Min",&GF::FitParametersBound::cosmicRay5Min)
    .def_readwrite("cosmicRay5Max",&GF::FitParametersBound::cosmicRay5Max)    
    .def_readwrite("cosmicRay6Min",&GF::FitParametersBound::cosmicRay6Min)
    .def_readwrite("cosmicRay6Max",&GF::FitParametersBound::cosmicRay6Max)    
    .def_readwrite("icegrad0Min",&GF::FitParametersBound::icegrad0Min)
    .def_readwrite("icegrad0Max",&GF::FitParametersBound::icegrad0Max)    
    .def_readwrite("icegrad1Min",&GF::FitParametersBound::icegrad1Min)
    .def_readwrite("icegrad1Max",&GF::FitParametersBound::icegrad1Max)    
    .def_readwrite("icegrad2Min",&GF::FitParametersBound::icegrad2Min)
    .def_readwrite("icegrad2Max",&GF::FitParametersBound::icegrad2Max)    
    .def_readwrite("icegrad3Min",&GF::FitParametersBound::icegrad3Min)
    .def_readwrite("icegrad3Max",&GF::FitParametersBound::icegrad3Max)    
    .def_readwrite("icegrad4Min",&GF::FitParametersBound::icegrad4Min)
    .def_readwrite("icegrad4Max",&GF::FitParametersBound::icegrad4Max)    
    .def_readwrite("icegrad5Min",&GF::FitParametersBound::icegrad5Min)
    .def_readwrite("icegrad5Max",&GF::FitParametersBound::icegrad5Max)    
    .def_readwrite("icegrad6Min",&GF::FitParametersBound::icegrad6Min)
    .def_readwrite("icegrad6Max",&GF::FitParametersBound::icegrad6Max)    
    .def_readwrite("icegrad7Min",&GF::FitParametersBound::icegrad7Min)
    .def_readwrite("icegrad7Max",&GF::FitParametersBound::icegrad7Max)    
    .def_readwrite("icegrad8Min",&GF::FitParametersBound::icegrad8Min)
    .def_readwrite("icegrad8Max",&GF::FitParametersBound::icegrad8Max)    
    .def_readwrite("domEfficiencyMin",&GF::FitParametersBound::domEfficiencyMin)
    .def_readwrite("domEfficiencyMax",&GF::FitParametersBound::domEfficiencyMax)    
    .def_readwrite("holeiceForwardMin",&GF::FitParametersBound::holeiceForwardMin)
    .def_readwrite("holeiceForwardMax",&GF::FitParametersBound::holeiceForwardMax)    
    .def_readwrite("astroNormMin",&GF::FitParametersBound::astroNormMin)
    .def_readwrite("astroNormMax",&GF::FitParametersBound::astroNormMax)    
    .def_readwrite("astroDeltaGammaMin",&GF::FitParametersBound::astroDeltaGammaMin)
    .def_readwrite("astroDeltaGammaMax",&GF::FitParametersBound::astroDeltaGammaMax)    
    .def_readwrite("astroDeltaGammaSecMin",&GF::FitParametersBound::astroDeltaGammaSecMin)
    .def_readwrite("astroDeltaGammaSecMax",&GF::FitParametersBound::astroDeltaGammaSecMax)    
    .def_readwrite("astroPivotMin",&GF::FitParametersBound::astroPivotMin)
    .def_readwrite("astroPivotMax",&GF::FitParametersBound::astroPivotMax)    
    .def_readwrite("NeutrinoAntineutrinoRatioMin",&GF::FitParametersBound::NeutrinoAntineutrinoRatioMin)
    .def_readwrite("NeutrinoAntineutrinoRatioMax",&GF::FitParametersBound::NeutrinoAntineutrinoRatioMax)    
    .def_readwrite("nuxsMin",&GF::FitParametersBound::nuxsMin)
    .def_readwrite("nuxsMax",&GF::FitParametersBound::nuxsMax)    
    .def_readwrite("nubarxsMin",&GF::FitParametersBound::nubarxsMin)
    .def_readwrite("nubarxsMax",&GF::FitParametersBound::nubarxsMax)    
  ;

  class_<GF::FitParametersFlag, boost::noncopyable,std::shared_ptr<GF::FitParametersFlag> >("FitParametersFlag",init<>())
    .def(init<bool>())
    .def_readwrite("convNorm",&GF::FitParametersFlag::convNorm)
    .def_readwrite("promptNorm",&GF::FitParametersFlag::promptNorm)
    .def_readwrite("zenithCorrection",&GF::FitParametersFlag::zenithCorrection)
    .def_readwrite("kaonLosses",&GF::FitParametersFlag::kaonLosses)
    .def_readwrite("hadronicHEkp",&GF::FitParametersFlag::hadronicHEkp)
    .def_readwrite("hadronicHEkm",&GF::FitParametersFlag::hadronicHEkm)
    .def_readwrite("hadronicVHE1pip",&GF::FitParametersFlag::hadronicVHE1pip)
    .def_readwrite("hadronicVHE1pim",&GF::FitParametersFlag::hadronicVHE1pim)
    .def_readwrite("hadronicVHE3kp",&GF::FitParametersFlag::hadronicVHE3kp)
    .def_readwrite("hadronicVHE3km",&GF::FitParametersFlag::hadronicVHE3km)
    .def_readwrite("hadronicVHE3pip",&GF::FitParametersFlag::hadronicVHE3pip)
    .def_readwrite("hadronicVHE3pim",&GF::FitParametersFlag::hadronicVHE3pim)
    .def_readwrite("hadronicVHE3p",&GF::FitParametersFlag::hadronicVHE3p)
    .def_readwrite("hadronicVHE3n",&GF::FitParametersFlag::hadronicVHE3n)
    .def_readwrite("cosmicRay1",&GF::FitParametersFlag::cosmicRay1)
    .def_readwrite("cosmicRay2",&GF::FitParametersFlag::cosmicRay2)
    .def_readwrite("cosmicRay3",&GF::FitParametersFlag::cosmicRay3)
    .def_readwrite("cosmicRay4",&GF::FitParametersFlag::cosmicRay4)
    .def_readwrite("cosmicRay5",&GF::FitParametersFlag::cosmicRay5)
    .def_readwrite("cosmicRay6",&GF::FitParametersFlag::cosmicRay6)
    .def_readwrite("icegrad0",&GF::FitParametersFlag::icegrad0)
    .def_readwrite("icegrad1",&GF::FitParametersFlag::icegrad1)
    .def_readwrite("icegrad2",&GF::FitParametersFlag::icegrad2)
    .def_readwrite("icegrad3",&GF::FitParametersFlag::icegrad3)
    .def_readwrite("icegrad4",&GF::FitParametersFlag::icegrad4)
    .def_readwrite("icegrad5",&GF::FitParametersFlag::icegrad5)
    .def_readwrite("icegrad6",&GF::FitParametersFlag::icegrad6)
    .def_readwrite("icegrad7",&GF::FitParametersFlag::icegrad7)
    .def_readwrite("icegrad8",&GF::FitParametersFlag::icegrad8)
    .def_readwrite("domEfficiency",&GF::FitParametersFlag::domEfficiency)
    .def_readwrite("holeiceForward",&GF::FitParametersFlag::holeiceForward)
    .def_readwrite("astroNorm",&GF::FitParametersFlag::astroNorm)
    .def_readwrite("astroDeltaGamma",&GF::FitParametersFlag::astroDeltaGamma)
    .def_readwrite("astroDeltaGammaSec",&GF::FitParametersFlag::astroDeltaGammaSec)
    .def_readwrite("astroPivot",&GF::FitParametersFlag::astroPivot)
    .def_readwrite("NeutrinoAntineutrinoRatio",&GF::FitParametersFlag::NeutrinoAntineutrinoRatio)
    .def_readwrite("nuxs",&GF::FitParametersFlag::nuxs)
    .def_readwrite("nubarxs",&GF::FitParametersFlag::nubarxs)
  ;

  class_<GF::Priors, boost::noncopyable,std::shared_ptr<GF::Priors> >("Priors",init<>())
    .def_readwrite("convNormCenter",&GF::Priors::convNormCenter)
    .def_readwrite("convNormWidth",&GF::Priors::convNormWidth) 
    .def_readwrite("convNormCenter",&GF::Priors::convNormCenter)
    .def_readwrite("promptNormWidth",&GF::Priors::promptNormWidth) 
    .def_readwrite("promptNormCenter",&GF::Priors::promptNormCenter)
    .def_readwrite("zenithCorrectionWidth",&GF::Priors::zenithCorrectionWidth) 
    .def_readwrite("zenithCorrectionCenter",&GF::Priors::zenithCorrectionCenter)
    .def_readwrite("kaonLossesWidth",&GF::Priors::kaonLossesWidth) 
    .def_readwrite("kaonLossesCenter",&GF::Priors::kaonLossesCenter)
    .def_readwrite("hadronicHEkpWidth",&GF::Priors::hadronicHEkpWidth) 
    .def_readwrite("hadronicHEkpCenter",&GF::Priors::hadronicHEkpCenter)
    .def_readwrite("hadronicHEkmWidth",&GF::Priors::hadronicHEkmWidth) 
    .def_readwrite("hadronicHEkmCenter",&GF::Priors::hadronicHEkmCenter)
    .def_readwrite("hadronicVHE1pipWidth",&GF::Priors::hadronicVHE1pipWidth) 
    .def_readwrite("hadronicVHE1pipCenter",&GF::Priors::hadronicVHE1pipCenter)
    .def_readwrite("hadronicVHE1pimWidth",&GF::Priors::hadronicVHE1pimWidth) 
    .def_readwrite("hadronicVHE1pimCenter",&GF::Priors::hadronicVHE1pimCenter)
    .def_readwrite("hadronicVHE3kpWidth",&GF::Priors::hadronicVHE3kpWidth) 
    .def_readwrite("hadronicVHE3kpCenter",&GF::Priors::hadronicVHE3kpCenter)
    .def_readwrite("hadronicVHE3kmWidth",&GF::Priors::hadronicVHE3kmWidth) 
    .def_readwrite("hadronicVHE3kmCenter",&GF::Priors::hadronicVHE3kmCenter)
    .def_readwrite("hadronicVHE3pipWidth",&GF::Priors::hadronicVHE3pipWidth) 
    .def_readwrite("hadronicVHE3pipCenter",&GF::Priors::hadronicVHE3pipCenter)
    .def_readwrite("hadronicVHE3pimWidth",&GF::Priors::hadronicVHE3pimWidth) 
    .def_readwrite("hadronicVHE3pimCenter",&GF::Priors::hadronicVHE3pimCenter)
    .def_readwrite("hadronicVHE3pWidth",&GF::Priors::hadronicVHE3pWidth) 
    .def_readwrite("hadronicVHE3pCenter",&GF::Priors::hadronicVHE3pCenter)
    .def_readwrite("hadronicVHE3nWidth",&GF::Priors::hadronicVHE3nWidth) 
    .def_readwrite("hadronicVHE3nCenter",&GF::Priors::hadronicVHE3nCenter)
    .def_readwrite("cosmicRay1Width",&GF::Priors::cosmicRay1Width) 
    .def_readwrite("cosmicRay1Center",&GF::Priors::cosmicRay1Center)
    .def_readwrite("cosmicRay2Width",&GF::Priors::cosmicRay2Width) 
    .def_readwrite("cosmicRay2Center",&GF::Priors::cosmicRay2Center)
    .def_readwrite("cosmicRay3Width",&GF::Priors::cosmicRay3Width) 
    .def_readwrite("cosmicRay3Center",&GF::Priors::cosmicRay3Center)
    .def_readwrite("cosmicRay4Width",&GF::Priors::cosmicRay4Width) 
    .def_readwrite("cosmicRay4Center",&GF::Priors::cosmicRay4Center)
    .def_readwrite("cosmicRay5Width",&GF::Priors::cosmicRay5Width) 
    .def_readwrite("cosmicRay5Center",&GF::Priors::cosmicRay5Center)
    .def_readwrite("cosmicRay6Width",&GF::Priors::cosmicRay6Width) 
    .def_readwrite("cosmicRay6Center",&GF::Priors::cosmicRay6Center)
    .def_readwrite("icegrad0Width",&GF::Priors::icegrad0Width) 
    .def_readwrite("icegrad0Center",&GF::Priors::icegrad0Center)
    .def_readwrite("icegrad1Width",&GF::Priors::icegrad1Width) 
    .def_readwrite("icegrad1Center",&GF::Priors::icegrad1Center)
    .def_readwrite("icegrad2Width",&GF::Priors::icegrad2Width) 
    .def_readwrite("icegrad2Center",&GF::Priors::icegrad2Center)
    .def_readwrite("icegrad3Width",&GF::Priors::icegrad3Width) 
    .def_readwrite("icegrad3Center",&GF::Priors::icegrad3Center)
    .def_readwrite("icegrad4Width",&GF::Priors::icegrad4Width) 
    .def_readwrite("icegrad4Center",&GF::Priors::icegrad4Center)
    .def_readwrite("icegrad5Width",&GF::Priors::icegrad5Width) 
    .def_readwrite("icegrad5Center",&GF::Priors::icegrad5Center)
    .def_readwrite("icegrad6Width",&GF::Priors::icegrad6Width) 
    .def_readwrite("icegrad6Center",&GF::Priors::icegrad6Center)
    .def_readwrite("icegrad7Width",&GF::Priors::icegrad7Width) 
    .def_readwrite("icegrad7Center",&GF::Priors::icegrad7Center)
    .def_readwrite("icegrad8Width",&GF::Priors::icegrad8Width) 
    .def_readwrite("icegrad8Center",&GF::Priors::icegrad8Center)
    .def_readwrite("domEfficiencyWidth",&GF::Priors::domEfficiencyWidth) 
    .def_readwrite("domEfficiencyCenter",&GF::Priors::domEfficiencyCenter)
    .def_readwrite("holeiceForwardWidth",&GF::Priors::holeiceForwardWidth) 
    .def_readwrite("holeiceForwardCenter",&GF::Priors::holeiceForwardCenter)
    .def_readwrite("astroNormWidth",&GF::Priors::astroNormWidth) 
    .def_readwrite("astroNormCenter",&GF::Priors::astroNormCenter)
    .def_readwrite("astroDeltaGammaWidth",&GF::Priors::astroDeltaGammaWidth) 
    .def_readwrite("astroDeltaGammaCenter",&GF::Priors::astroDeltaGammaCenter)
    .def_readwrite("astroDeltaGammaSecWidth",&GF::Priors::astroDeltaGammaSecWidth) 
    .def_readwrite("astroDeltaGammaSecCenter",&GF::Priors::astroDeltaGammaSecCenter)
    .def_readwrite("astroPivotMin",&GF::Priors::astroPivotMin) 
    .def_readwrite("astroPivotMax",&GF::Priors::astroPivotMax) 
    .def_readwrite("NeutrinoAntineutrinoRatioWidth",&GF::Priors::NeutrinoAntineutrinoRatioWidth) 
    .def_readwrite("NeutrinoAntineutrinoRatioCenter",&GF::Priors::NeutrinoAntineutrinoRatioCenter)
    .def_readwrite("nuxsWidth",&GF::Priors::nuxsWidth) 
    .def_readwrite("nuxsCenter",&GF::Priors::nuxsCenter)
    .def_readwrite("nubarxsWidth",&GF::Priors::nubarxsWidth) 
    .def_readwrite("nubarxsCenter",&GF::Priors::nubarxsCenter)
    .def("SetFluxCorr",&GF::Priors::SetFluxCorr)
    .def("SetIceGradientsCorr",&GF::Priors::SetIceGradientsCorr)
 ;

  class_<GF::GollumFit, boost::noncopyable, std::shared_ptr<GF::GollumFit> >("GollumFit", init<GF::DataPaths,GF::SteeringParams>())
    .def("ReConfig",&GF::GollumFit::ReConfig<false>)
    .def("CheckDataLoaded",&GF::GollumFit::CheckDataLoaded)
    .def("CheckSimulationLoaded",&GF::GollumFit::CheckSimulationLoaded)
    .def("CheckCrossSectionWeighterConstructed",&GF::GollumFit::CheckCrossSectionWeighterConstructed)
    .def("CheckFluxWeighterConstructed",&GF::GollumFit::CheckFluxWeighterConstructed)
    .def("CheckLeptonWeighterConstructed",&GF::GollumFit::CheckLeptonWeighterConstructed)
    .def("CheckDataHistogramConstructed",&GF::GollumFit::CheckDataHistogramConstructed)
    .def("CheckLeptonWeighterConstructed",&GF::GollumFit::CheckLeptonWeighterConstructed)
    .def("CheckDataHistogramConstructed",&GF::GollumFit::CheckDataHistogramConstructed)
    .def("CheckSimulationHistogramConstructed",&GF::GollumFit::CheckSimulationHistogramConstructed)
    .def("CheckLikelihoodProblemConstruction",&GF::GollumFit::CheckLikelihoodProblemConstruction)
    .def("ReportStatus",&GF::GollumFit::ReportStatus)
    .def("GetDataDistribution",&GF::GollumFit::GetDataDistribution)
    .def("ConstructFastMode",&GF::GollumFit::ConstructFastMode)
    .def("ConstructLikelihoodProblem",&GF::GollumFit::ConstructLikelihoodProblem)
    .def("WriteCompact",&GF::GollumFit::WriteCompact)
    .def("GetEnergyBinsMC",&GF::GollumFit::GetEnergyBinsMC)
    .def("GetZenithBinsMC",&GF::GollumFit::GetZenithBinsMC)
    .def("GetTopologyBinsMC",&GF::GollumFit::GetTopologyBinsMC)
    .def("GetExpectation",(nsq::marray<double,3>(GF::GollumFit::*)(GF::FitParameters)const)&GF::GollumFit::GetExpectation)
    .def("GetSquareExpectation",(nsq::marray<double,3>(GF::GollumFit::*)(GF::FitParameters)const)&GF::GollumFit::GetSquareExpectation)
    .def("GetRealization",(nsq::marray<double,3>(GF::GollumFit::*)(GF::FitParameters,int)const)&GF::GollumFit::GetRealization)
    .def("GetDataEvents",&GF::GollumFit::GetDataEvents)
    .def("GetRealizationEvents",(nsq::marray<double,2>(GF::GollumFit::*)(GF::FitParameters,int)const)&GF::GollumFit::GetRealizationEvents)
    .def("GetExpectationEvents",(nsq::marray<double,2>(GF::GollumFit::*)(GF::FitParameters)const)&GF::GollumFit::GetExpectationEvents)
    .def("CheckExpectation",(int(GF::GollumFit::*)(GF::FitParameters)const)&GF::GollumFit::CheckExpectation)
    .def("EvalLLH",(double(GF::GollumFit::*)(GF::FitParameters,bool)const)&GF::GollumFit::EvalLLH)
    .def("SetData",wrap_SetData)
    .def("MinLLH",(GF::FitResult(GF::GollumFit::*)() const)&GF::GollumFit::MinLLH)
    .def("SetFitParametersSeed",(void(GF::GollumFit::*)(std::vector<GF::FitParameters>))&GF::GollumFit::SetFitParametersSeed)
    .def("SetFitParametersFlag",&GF::GollumFit::SetFitParametersFlag)
    .def("SetFitParametersBound",&GF::GollumFit::SetFitParametersBound)
    .def("SetFitParametersPriors",&GF::GollumFit::SetFitParametersPriors)
  ;

  enum_<LW::ParticleType>("ParticleType")
    .value("NuE",LW::ParticleType::NuE)
    .value("NuMu",LW::ParticleType::NuMu)
    .value("NuTau",LW::ParticleType::NuTau)
    .value("NuEBar",LW::ParticleType::NuEBar)
    .value("NuMuBar",LW::ParticleType::NuMuBar)
    .value("NuTauBar",LW::ParticleType::NuTauBar)

    .value("EMinus",LW::ParticleType::EMinus)
    .value("EPlus",LW::ParticleType::EPlus)
    .value("MuMinus",LW::ParticleType::MuMinus)
    .value("MuPlus",LW::ParticleType::MuPlus)
    .value("TauMinus",LW::ParticleType::TauMinus)
    .value("TauPlus",LW::ParticleType::TauPlus)

    .value("unknown",LW::ParticleType::unknown)
    .value("Hadrons",LW::ParticleType::Hadrons)
  ;

  class_<Event, std::shared_ptr<Event> >("Event", init<>())
    // MC Truth properties
    .def_readonly("final_state_particle_0",&Event::final_state_particle_0)
    .def_readonly("final_state_particle_1",&Event::final_state_particle_1)
    .def_readonly("primaryType",&Event::primaryType)
    .def_readonly("primaryEnergy",&Event::primaryEnergy)
    .def_readonly("primaryAzimuth",&Event::primaryAzimuth)
    .def_readonly("primaryZenith",&Event::primaryZenith)
    .def_readonly("totalColumnDepth",&Event::totalColumnDepth)
    .def_readonly("intX",&Event::intX)
    .def_readonly("intY",&Event::intY)
    // Reco properties
    .def_readonly("energy",&Event::energy)
    .def_readonly("zenith",&Event::zenith)
  ;

  // python container to vector<double> convertion
  using namespace scitbx::boost_python::container_conversions;
  from_python_sequence< std::vector<double>, variable_capacity_policy >();
  to_python_converter< std::vector<double, class std::allocator<double> >, VecToList<double> > ();
  from_python_sequence< std::vector<int>, variable_capacity_policy >();
  to_python_converter< std::vector<int, class std::allocator<int> >, VecToList<int> > ();
  from_python_sequence< std::vector<unsigned int>, variable_capacity_policy >();
  to_python_converter< std::vector<unsigned int, class std::allocator<unsigned int> >, VecToList<unsigned int> > ();

  from_python_sequence< std::deque<Event>, variable_capacity_policy >();
  to_python_converter< std::deque<Event, class std::allocator<Event> >, DeqToList<Event> > ();

  from_python_sequence< std::vector<GF::FitParameters>, variable_capacity_policy >();
  to_python_converter< std::vector<GF::FitParameters, class std::allocator<GF::FitParameters> >, VecToList<GF::FitParameters> > ();

  from_python_sequence< std::vector<string>, variable_capacity_policy >();
  to_python_converter< std::vector<string, class std::allocator<string> >, VecToList<string> > ();

  to_python_converter< nsq::marray<double,1> , marray_to_numpyarray<1> >();
  to_python_converter< nsq::marray<double,2> , marray_to_numpyarray<2> >();
  to_python_converter< nsq::marray<double,3> , marray_to_numpyarray<3> >();
  to_python_converter< GF::hist_marray , marray_to_numpyarray<3> >();

}
