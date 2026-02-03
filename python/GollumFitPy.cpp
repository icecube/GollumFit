#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>

#include <functional>
#include "GollumFit.h"

#include <numpy/ndarrayobject.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

#include <nuSQuIDS/marray.h>

using namespace pybind11;
namespace py = pybind11;
namespace GF = gollumfit;
namespace nsq = nusquids;

namespace pybind11 { namespace detail {
    template <unsigned int Dim>
    struct type_caster<nsq::marray<double, Dim>> {
    private:
        using T = nsq::marray<double, Dim>;
    public:
        PYBIND11_TYPE_CASTER(T, _("marray<double,Dim>"));

        // Python -> C++
        bool load(handle src, bool) {
            if (!py::isinstance<py::array>(src))
                return false;

            auto array = py::array::ensure(src);
            if (!array || array.ndim() != Dim || !py::isinstance<py::array_t<double>>(array))
                return false;

            std::array<size_t, Dim> shape;
            for (size_t i = 0; i < Dim; ++i)
                shape[i] = static_cast<size_t>(array.shape(i));

            double* ptr = static_cast<double*>(array.mutable_data());
            if (!ptr)
                return false;

            value = T();
            value.resize(shape);
            std::memcpy(value.get_data(), array.data(), sizeof(double) * array.size());
            return true;
        }

        // C++ -> Python
        static handle cast(const nsq::marray<double, Dim>& arr, return_value_policy, handle parent) {
            std::vector<ssize_t> shape(Dim);
            std::vector<ssize_t> strides(Dim);

            ssize_t stride = sizeof(double);
            for (ssize_t i = Dim - 1; i >= 0; --i) {
                shape[i] = arr.extent(i);
                strides[i] = stride;
                stride *= shape[i];
            }

            return py::array(py::buffer_info(
                const_cast<double*>(arr.get_data()), // assume mutable data
                sizeof(double),
                py::format_descriptor<double>::format(),
                Dim,
                shape,
                strides
            )).release();
        }
    };
}} // namespace pybind11::detail


// auxiliary wrapper functions // evil // Carlos

// static double wrap_SetData(GF::GollumFit* st, PyObject * array){
//   if (! PyArray_Check(array) )
//   {
//     throw std::runtime_error("GollumFit::Error:Input array is not a numpy array.");
//   }

//   PyArrayObject* numpy_array = (PyArrayObject*)array;
//   unsigned int array_dim = PyArray_NDIM(numpy_array);
//   NPY_TYPES type = (NPY_TYPES) PyArray_DESCR(numpy_array)->type_num;

//   // things i think can cast ok to doubles
//   if (!( type == NPY_LONG or type == NPY_INT or type == NPY_SHORT or type == NPY_FLOAT or
//       type == NPY_DOUBLE or type == NPY_LONGDOUBLE or type == NPY_CFLOAT or type == NPY_CDOUBLE))
//     throw std::runtime_error("GollumFit::Error:Input numpy array cannot be meaninfully casted into double.");

//   if ( array_dim == 2 ) {
//     nsq::marray<double,2> state = numpyarray_to_marray<double,2>(array, type);
//     double value = st->SetData(state);
//     return value;
//   } else
//     throw std::runtime_error("GollumFit::Input array has wrong dimenions.");
// }


//static Dict2Map<unsigned int, double> reg{};

PYBIND11_MODULE(GollumFitPy, m)
{
  m.doc() = "Python bindings for GollumFit via pybind11";

  // m.def("marray_to_numpy1", &marray_to_numpy<1>);
  // m.def("marray_to_numpy2", &marray_to_numpy<2>);
  // m.def("marray_to_numpy3", &marray_to_numpy<3>);
  // m.def("numpy_to_marray1", &numpy_to_marray<1>);

  py::enum_<GF::FluxComponent>(m, "FluxComponent")
    .value("atmConv",GF::FluxComponent::atmConv)
    .value("atmPrompt",GF::FluxComponent::atmPrompt)
    .value("diffuseAstro",GF::FluxComponent::diffuseAstro)
    .export_values()
  ;

  py::class_<GF::FitResult, std::shared_ptr<GF::FitResult> >(m, "FitResult")
    .def(py::init<>())
    .def_readwrite("params",&GF::FitResult::params)
    .def_readwrite("likelihood",&GF::FitResult::likelihood)
    .def_readwrite("aux_likelihood",&GF::FitResult::aux_likelihood)
    .def_readwrite("nEval",&GF::FitResult::nEval)
    .def_readwrite("nGrad",&GF::FitResult::nGrad)
    .def_readwrite("succeeded",&GF::FitResult::succeeded)
  ;

  py::class_<GF::hist_marray>(m, "hist_marray")
    .def(py::init<>())
  ;

  py::class_<GF::DataPaths,std::shared_ptr<GF::DataPaths> >(m, "DataPaths")
    .def(py::init<>())
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

  py::class_<GF::SteeringParams ,std::shared_ptr<GF::SteeringParams> >(m, "SteeringParams")
    .def(py::init<>())
    .def_readwrite("minFitEnergy",&GF::SteeringParams::minFitEnergy, "Minimum energy in the fit")
    .def_readwrite("maxFitEnergy",&GF::SteeringParams::maxFitEnergy, "Maximum energy in the fit")
    .def_readwrite("minCosth",&GF::SteeringParams::minCosth)
    .def_readwrite("maxCosth",&GF::SteeringParams::maxCosth)
    .def_readwrite("logEbinEdge",&GF::SteeringParams::logEbinEdge)
    .def_readwrite("logEbinWidth",&GF::SteeringParams::logEbinWidth)
    .def_readwrite("cosThbinEdge",&GF::SteeringParams::cosThbinEdge)
    .def_readwrite("cosThbinWidth",&GF::SteeringParams::cosThbinWidth)
    .def_readwrite("astroOscAvgScale",&GF::SteeringParams::astroOscAvgScale)
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
    .def_readwrite("enableTotalNorm",&GF::SteeringParams::enableTotalNorm)
  ;


  py::class_<GF::FitParameters, std::shared_ptr<GF::FitParameters> >(m, "FitParameters")
    .def(py::init<>())
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

  py::class_<GF::FitParametersBound, std::shared_ptr<GF::FitParametersBound> >(m, "FitParametersBound")
    .def(py::init<>())
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

  py::class_<GF::FitParametersFlag, std::shared_ptr<GF::FitParametersFlag> >(m, "FitParametersFlag")
    .def(py::init<>())
    .def(py::init<bool>())
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

  py::class_<GF::Priors, std::shared_ptr<GF::Priors> >(m, "Priors")
    .def(py::init<>())
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

  py::class_<GF::GollumFit, std::shared_ptr<GF::GollumFit> >(m, "GollumFit")
    .def(py::init<GF::DataPaths, GF::SteeringParams>())
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
    .def("SetData", &GF::GollumFit::SetData)
    .def("MinLLH",(GF::FitResult(GF::GollumFit::*)() const)&GF::GollumFit::MinLLH)
    .def("SetFitParametersSeed",(void(GF::GollumFit::*)(std::vector<GF::FitParameters>))&GF::GollumFit::SetFitParametersSeed)
    .def("SetFitParametersFlag",&GF::GollumFit::SetFitParametersFlag)
    .def("SetFitParametersBound",&GF::GollumFit::SetFitParametersBound)
    .def("SetFitParametersPriors",&GF::GollumFit::SetFitParametersPriors)
  ;

  py::enum_<LW::ParticleType>(m, "ParticleType")
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

    .export_values()
  ;

  py::class_<Event, std::shared_ptr<Event> >(m, "Event")
    // MC Truth properties
    .def(py::init<>())
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

  // to_python_converter< nsq::marray<double,1> , marray_to_numpy<1> >();
  // to_python_converter< nsq::marray<double,2> , marray_to_numpy<2> >();
  // to_python_converter< nsq::marray<double,3> , marray_to_numpy<3> >();
  // to_python_converter< GF::hist_marray , marray_to_numpy<3> >();
}
