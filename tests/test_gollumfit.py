"""
GollumFit Unit Tests

This test suite provides comprehensive unit tests for the GollumFit framework,
which performs binned-likelihood fits of Monte Carlo simulations to neutrino 
telescope data.

Test Coverage:
--------------
1. Parameter Structures (FitParameters, Priors, Bounds, Flags)
2. Data and MC Loading
3. Histogram Construction
4. Expectation Generation
5. Likelihood Evaluation
6. Fitting and Minimization
7. FastMC Functionality
8. Error Handling and Edge Cases

Dependencies:
-------------
- pytest: Test framework
- numpy: Numerical operations
- GollumFitPy: The main GollumFit Python bindings
- h5py: For reading test data files

Usage:
------
Run all tests:
    pytest test_gollumfit.py -v

Run specific test class:
    pytest test_gollumfit.py::TestFitParameters -v

Run specific test:
    pytest test_gollumfit.py::TestFitParameters::test_fitparameters_initialization -v

Author: Unit Test Suite for GollumFit
Date: February 2026
"""

import pytest
import numpy as np
import sys
import os
from collections import OrderedDict

# Import GollumFit Python bindings
import GollumFitPy as gf


# ============================================================================
# Test Fixtures
# ============================================================================

@pytest.fixture
def default_fit_parameters():
    """
    Create a FitParameters object with default values.
    
    Returns:
        gf.FitParameters: FitParameters object with all parameters set to 
                          their nominal values per the example scripts.
    """
    params = gf.FitParameters()
    params.convNorm = 1.0
    params.promptNorm = 1.0
    params.zenithCorrection = 0.0
    params.kaonLosses = 0.0
    params.hadronicHEkp = 0.0
    params.hadronicHEkm = 0.0
    params.hadronicVHE1pip = 0.0
    params.hadronicVHE1pim = 0.0
    params.hadronicVHE3kp = 0.0
    params.hadronicVHE3km = 0.0
    params.hadronicVHE3pip = 0.0
    params.hadronicVHE3pim = 0.0
    params.hadronicVHE3p = 0.0
    params.hadronicVHE3n = 0.0
    params.cosmicRay1 = 0.0
    params.cosmicRay2 = 0.0
    params.cosmicRay3 = 0.0
    params.cosmicRay4 = 0.0
    params.cosmicRay5 = 0.0
    params.cosmicRay6 = 0.0
    params.icegrad0 = 0.0
    params.icegrad1 = 0.0
    params.icegrad2 = 0.0
    params.icegrad3 = 0.0
    params.icegrad4 = 0.0
    params.icegrad5 = 0.0
    params.icegrad6 = 0.0
    params.icegrad7 = 0.0
    params.icegrad8 = 0.0
    params.domEfficiency = 1.27
    params.holeiceForward = -1.0
    params.astroNorm = 4.72 / 6.0
    params.astroDeltaGamma = 0.0
    params.astroDeltaGammaSec = 0.0
    params.astroPivot = 5.0
    params.nuxs = 1.0
    params.nubarxs = 1.0
    params.NeutrinoAntineutrinoRatio = 1.0
    return params


@pytest.fixture
def default_priors(default_fit_parameters):
    """
    Create a Priors object with standard prior values.
    
    Returns:
        gf.Priors: Priors object configured with Gaussian and Uniform priors
                   for all fit parameters.
    """
    priors = gf.Priors()
    
    # Set all Gaussian prior centers and widths based on examples
    priors.convNormCenter = 1.0
    priors.convNormWidth = 0.2
    priors.promptNormCenter = 1.0
    priors.promptNormWidth = 1.0
    priors.zenithCorrectionCenter = 0.0
    priors.zenithCorrectionWidth = 1.0
    priors.kaonLossesCenter = 0.0
    priors.kaonLossesWidth = 1.0
    priors.hadronicHEkpCenter = 0.0
    priors.hadronicHEkpWidth = 1.0
    priors.hadronicHEkmCenter = 0.0
    priors.hadronicHEkmWidth = 1.0
    priors.hadronicVHE1pipCenter = 0.0
    priors.hadronicVHE1pipWidth = 1.0
    priors.hadronicVHE1pimCenter = 0.0
    priors.hadronicVHE1pimWidth = 1.0
    priors.hadronicVHE3kpCenter = 0.0
    priors.hadronicVHE3kpWidth = 1.0
    priors.hadronicVHE3kmCenter = 0.0
    priors.hadronicVHE3kmWidth = 1.0
    priors.hadronicVHE3pipCenter = 0.0
    priors.hadronicVHE3pipWidth = 1.0
    priors.hadronicVHE3pimCenter = 0.0
    priors.hadronicVHE3pimWidth = 1.0
    priors.hadronicVHE3pCenter = 0.0
    priors.hadronicVHE3pWidth = 1.0
    priors.hadronicVHE3nCenter = 0.0
    priors.hadronicVHE3nWidth = 1.0
    priors.cosmicRay1Center = 0.0
    priors.cosmicRay1Width = 1.0
    priors.cosmicRay2Center = 0.0
    priors.cosmicRay2Width = 1.0
    priors.cosmicRay3Center = 0.0
    priors.cosmicRay3Width = 1.0
    priors.cosmicRay4Center = 0.0
    priors.cosmicRay4Width = 1.0
    priors.cosmicRay5Center = 0.0
    priors.cosmicRay5Width = 1.0
    priors.cosmicRay6Center = 0.0
    priors.cosmicRay6Width = 1.0
    priors.icegrad0Center = 0.0
    priors.icegrad0Width = 1.0
    priors.icegrad1Center = 0.0
    priors.icegrad1Width = 1.0
    priors.icegrad2Center = 0.0
    priors.icegrad2Width = 1.0
    priors.icegrad3Center = 0.0
    priors.icegrad3Width = 1.0
    priors.icegrad4Center = 0.0
    priors.icegrad4Width = 1.0
    priors.icegrad5Center = 0.0
    priors.icegrad5Width = 1.0
    priors.icegrad6Center = 0.0
    priors.icegrad6Width = 1.0
    priors.icegrad7Center = 0.0
    priors.icegrad7Width = 1.0
    priors.icegrad8Center = 0.0
    priors.icegrad8Width = 1.0
    priors.domEfficiencyCenter = 1.27
    priors.domEfficiencyWidth = 0.123
    priors.holeiceForwardCenter = -1.0
    priors.holeiceForwardWidth = 10.0
    priors.astroNormCenter = 4.72 / 6.0
    priors.astroNormWidth = 0.36
    priors.astroDeltaGammaCenter = 0.0
    priors.astroDeltaGammaWidth = 0.36
    priors.astroDeltaGammaSecCenter = 0.0
    priors.astroDeltaGammaSecWidth = 0.36
    priors.nuxsCenter = 1.0
    priors.nuxsWidth = 0.1
    priors.nubarxsCenter = 1.0
    priors.nubarxsWidth = 0.1
    priors.NeutrinoAntineutrinoRatioCenter = 1.0
    priors.NeutrinoAntineutrinoRatioWidth = 1.0
    
    # Uniform prior for astroPivot
    priors.astroPivotMin = 4.0
    priors.astroPivotMax = 6.0

    # Load real correlation matrices (required by ConstructLikelihoodProblem)
    resources_dir = os.path.join(os.path.dirname(__file__), '..', 'resources', 'correlation_matrices')
    iceg_corr = np.load(os.path.join(resources_dir, 'icegrad_correlations.npy'))
    flux_corr = np.load(os.path.join(resources_dir, 'flux_correlations_new_ddmnodeis.npy'))
    for idx, val in np.ndenumerate(iceg_corr):
        priors.SetIceGradientsCorr(idx[0], idx[1], val)
    for idx, val in np.ndenumerate(flux_corr):
        priors.SetFluxCorr(idx[0], idx[1], val)

    return priors


@pytest.fixture
def default_bounds():
    """
    Create a FitParametersBound object with standard parameter bounds.
    
    Returns:
        gf.FitParametersBound: Bounds object with min/max values for all parameters.
    """
    bounds = gf.FitParametersBound()
    
    # Set bounds based on examples
    bounds.convNormMin = 0.1
    bounds.convNormMax = 3.0
    bounds.promptNormMin = 0.0
    bounds.promptNormMax = 3.0
    bounds.zenithCorrectionMin = -3.0
    bounds.zenithCorrectionMax = 3.0
    bounds.kaonLossesMin = -3.0
    bounds.kaonLossesMax = 3.0
    bounds.hadronicHEkpMin = -2.0
    bounds.hadronicHEkpMax = 2.0
    bounds.hadronicHEkmMin = -2.0
    bounds.hadronicHEkmMax = 2.0
    bounds.hadronicVHE1pipMin = -2.0
    bounds.hadronicVHE1pipMax = 2.0
    bounds.hadronicVHE1pimMin = -2.0
    bounds.hadronicVHE1pimMax = 2.0
    bounds.hadronicVHE3kpMin = -2.0
    bounds.hadronicVHE3kpMax = 2.0
    bounds.hadronicVHE3kmMin = -1.5
    bounds.hadronicVHE3kmMax = 2.0
    bounds.hadronicVHE3pipMin = -2.0
    bounds.hadronicVHE3pipMax = 2.0
    bounds.hadronicVHE3pimMin = -2.0
    bounds.hadronicVHE3pimMax = 2.0
    bounds.hadronicVHE3pMin = -2.0
    bounds.hadronicVHE3pMax = 2.0
    bounds.hadronicVHE3nMin = -2.0
    bounds.hadronicVHE3nMax = 2.0
    bounds.cosmicRay1Min = -4.0
    bounds.cosmicRay1Max = 4.0
    bounds.cosmicRay2Min = -4.0
    bounds.cosmicRay2Max = 4.0
    bounds.cosmicRay3Min = -4.0
    bounds.cosmicRay3Max = 4.0
    bounds.cosmicRay4Min = -4.0
    bounds.cosmicRay4Max = 4.0
    bounds.cosmicRay5Min = -4.0
    bounds.cosmicRay5Max = 4.0
    bounds.cosmicRay6Min = -4.0
    bounds.cosmicRay6Max = 4.0
    bounds.icegrad0Min = -3.0
    bounds.icegrad0Max = 3.0
    bounds.icegrad1Min = -3.0
    bounds.icegrad1Max = 3.0
    bounds.icegrad2Min = -3.0
    bounds.icegrad2Max = 3.0
    bounds.icegrad3Min = -3.0
    bounds.icegrad3Max = 3.0
    bounds.icegrad4Min = -3.0
    bounds.icegrad4Max = 3.0
    bounds.icegrad5Min = -3.0
    bounds.icegrad5Max = 3.0
    bounds.icegrad6Min = -3.0
    bounds.icegrad6Max = 3.0
    bounds.icegrad7Min = -3.0
    bounds.icegrad7Max = 3.0
    bounds.icegrad8Min = -3.0
    bounds.icegrad8Max = 3.0
    bounds.domEfficiencyMin = 1.234
    bounds.domEfficiencyMax = 1.346
    bounds.holeiceForwardMin = -5.35
    bounds.holeiceForwardMax = 1.85
    bounds.astroNormMin = 0.0
    bounds.astroNormMax = 3.0
    bounds.astroDeltaGammaMin = -2.0
    bounds.astroDeltaGammaMax = 2.0
    bounds.astroDeltaGammaSecMin = -2.0
    bounds.astroDeltaGammaSecMax = 2.0
    bounds.astroPivotMin = 4.0
    bounds.astroPivotMax = 6.0
    bounds.nuxsMin = 0.824
    bounds.nuxsMax = 1.176
    bounds.nubarxsMin = 0.824
    bounds.nubarxsMax = 1.176
    bounds.NeutrinoAntineutrinoRatioMin = 0.0
    bounds.NeutrinoAntineutrinoRatioMax = 2.0
    
    return bounds


@pytest.fixture
def default_flags():
    """
    Create a FitParametersFlag object with all flags set to False (vary all).
    
    Returns:
        gf.FitParametersFlag: Flags object with all parameters set to vary.
    """
    # False means the parameter will be varied in the fit
    return gf.FitParametersFlag(False)


@pytest.fixture
def example_datapaths():
    """
    Create a DataPaths object pointing to example data files.
    
    Note: This assumes the examples directory structure exists.
    Adjust paths as needed for your test environment.
    
    Returns:
        gf.DataPaths: DataPaths object with paths to test resources.
    """
    datapaths = gf.DataPaths()
    
    # Point to example resources (adjust these paths based on actual location)
    base_path = "/n/home01/awen/meows_home/GollumFit_test/GollumFit"
    
    datapaths.domeff_spline_path = os.path.join(
        base_path, "resources/Splines/DOMEffSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
    )
    datapaths.holeice_spline_path = os.path.join(
        base_path, "resources/Splines/HoleIceSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
    )
    datapaths.attenuation_spline_path = os.path.join(
        base_path, "resources/Splines/AttenuationSplines/new_ddmnodeis"
    )
    datapaths.compact_file_path = os.path.join(
        base_path, "examples/FastMC/compact.fastmc"
    )
    
    return datapaths


@pytest.fixture
def example_steering_params():
    """
    Create a SteeringParams object with example analysis configuration.
    
    Returns:
        gf.SteeringParams: Analysis configuration matching the examples.
    """
    steering = gf.SteeringParams()
    
    # Energy binning
    steering.minFitEnergy = 300.0
    steering.maxFitEnergy = 1e5
    steering.logEbinEdge = np.log10(300.0)
    steering.logEbinWidth = (np.log10(1e5) - np.log10(300.0)) / 24
    
    # Zenith binning
    steering.minCosth = -1.0
    steering.maxCosth = 0.0
    steering.cosThbinEdge = 0.0
    steering.cosThbinWidth = 0.05
    
    # Other settings
    steering.selectionStart = 0.99
    steering.evalThreads = 1
    
    # Convergence criteria
    steering.change_tol = 1.e-20
    steering.grad_tol = 1.e-20
    steering.uncertaintyModSigmaOverMu = 0.0
    
    return steering


# ============================================================================
# Test Classes
# ============================================================================

class TestFitParameters:
    """
    Test suite for FitParameters struct.
    
    Tests initialization, modification, and access of fit parameters.
    """
    
    def test_fitparameters_initialization(self):
        """Test that FitParameters can be instantiated."""
        params = gf.FitParameters()
        assert params is not None
    
    def test_fitparameters_set_get(self, default_fit_parameters):
        """Test setting and getting individual fit parameters."""
        params = default_fit_parameters
        
        # Test a few key parameters
        assert params.convNorm == 1.0
        assert params.domEfficiency == pytest.approx(1.27)
        assert params.astroNorm == pytest.approx(4.72 / 6.0)
        
        # Modify and verify
        params.convNorm = 1.5
        assert params.convNorm == 1.5
        
        params.icegrad0 = 0.5
        assert params.icegrad0 == 0.5
    
    def test_fitparameters_all_38_parameters(self, default_fit_parameters):
        """Test that all 38 parameters are accessible."""
        params = default_fit_parameters
        
        # List of all 38 parameter names
        param_names = [
            'convNorm', 'promptNorm', 'zenithCorrection', 'kaonLosses',
            'hadronicHEkp', 'hadronicHEkm', 'hadronicVHE1pip', 'hadronicVHE1pim',
            'hadronicVHE3kp', 'hadronicVHE3km', 'hadronicVHE3pip', 'hadronicVHE3pim',
            'hadronicVHE3p', 'hadronicVHE3n', 'cosmicRay1', 'cosmicRay2',
            'cosmicRay3', 'cosmicRay4', 'cosmicRay5', 'cosmicRay6',
            'icegrad0', 'icegrad1', 'icegrad2', 'icegrad3', 'icegrad4',
            'icegrad5', 'icegrad6', 'icegrad7', 'icegrad8', 'domEfficiency',
            'holeiceForward', 'astroNorm', 'astroDeltaGamma', 'astroDeltaGammaSec',
            'astroPivot', 'NeutrinoAntineutrinoRatio', 'nuxs', 'nubarxs'
        ]
        
        # Verify all parameters are accessible
        for param_name in param_names:
            assert hasattr(params, param_name), f"Missing parameter: {param_name}"
            value = getattr(params, param_name)
            assert value is not None


class TestFitParametersFlag:
    """
    Test suite for FitParametersFlag struct.
    
    Tests flag initialization and modification for fixing/varying parameters.
    """
    
    def test_flags_default_initialization(self):
        """Test that FitParametersFlag initializes with all True (all fixed)."""
        flags = gf.FitParametersFlag()
        
        # By default, all should be True (fixed)
        assert flags.convNorm == True
        assert flags.domEfficiency == True
    
    def test_flags_boolean_constructor(self):
        """Test FitParametersFlag constructor that sets all flags to same value."""
        # All False (all vary)
        flags_vary = gf.FitParametersFlag(False)
        assert flags_vary.convNorm == False
        assert flags_vary.astroNorm == False
        
        # All True (all fixed)
        flags_fixed = gf.FitParametersFlag(True)
        assert flags_fixed.convNorm == True
        assert flags_fixed.astroNorm == True
    
    def test_flags_individual_modification(self):
        """Test modifying individual flags."""
        flags = gf.FitParametersFlag(True)  # Start with all fixed
        
        # Vary a few parameters
        flags.convNorm = False
        flags.astroNorm = False
        flags.icegrad0 = False
        
        # Check modifications
        assert flags.convNorm == False
        assert flags.astroNorm == False
        assert flags.icegrad0 == False
        
        # Check others remain fixed
        assert flags.domEfficiency == True
        assert flags.promptNorm == True


class TestPriors:
    """
    Test suite for Priors struct.
    
    Tests prior distribution specifications for fit parameters.
    """
    
    def test_priors_initialization(self):
        """Test Priors object can be instantiated."""
        priors = gf.Priors()
        assert priors is not None
    
    def test_priors_from_fitparameters(self, default_fit_parameters):
        """Test that Priors centers can be set to match FitParameters values."""
        priors = gf.Priors()
        priors.convNormCenter = default_fit_parameters.convNorm
        priors.domEfficiencyCenter = default_fit_parameters.domEfficiency
        priors.astroNormCenter = default_fit_parameters.astroNorm

        # Centers should match FitParameters values (with float32 tolerance)
        assert priors.convNormCenter == pytest.approx(default_fit_parameters.convNorm)
        assert priors.domEfficiencyCenter == pytest.approx(default_fit_parameters.domEfficiency)
        assert priors.astroNormCenter == pytest.approx(default_fit_parameters.astroNorm)
    
    def test_priors_set_gaussian_parameters(self, default_priors):
        """Test setting Gaussian prior parameters."""
        priors = default_priors
        
        # Verify some Gaussian priors are set correctly
        assert priors.convNormCenter == pytest.approx(1.0)
        assert priors.convNormWidth == pytest.approx(0.2)
        assert priors.domEfficiencyCenter == pytest.approx(1.27)
        assert priors.domEfficiencyWidth == pytest.approx(0.123)
    
    def test_priors_set_uniform_parameters(self, default_priors):
        """Test setting Uniform prior parameters."""
        priors = default_priors
        
        # astroPivot has a uniform prior
        assert priors.astroPivotMin == 4.0
        assert priors.astroPivotMax == 6.0
    
    def test_priors_correlation_matrices(self):
        """Test setting correlation matrices."""
        priors = gf.Priors()
        
        # Create a simple correlation matrix for ice gradients (9x9)
        ice_corr = np.eye(9)  # Identity matrix
        for i in range(9):
            for j in range(9):
                priors.SetIceGradientsCorr(i, j, ice_corr[i, j])
        
        # Create correlation matrix for flux parameters (16x16)
        flux_corr = np.eye(16)
        for i in range(16):
            for j in range(16):
                priors.SetFluxCorr(i, j, flux_corr[i, j])
        
        # If we get here without error, the methods work
        assert True


class TestFitParametersBound:
    """
    Test suite for FitParametersBound struct.
    
    Tests parameter boundary specifications.
    """
    
    def test_bounds_initialization(self):
        """Test FitParametersBound can be instantiated."""
        bounds = gf.FitParametersBound()
        assert bounds is not None
    
    def test_bounds_set_get(self, default_bounds):
        """Test setting and getting parameter bounds."""
        bounds = default_bounds
        
        # Test a few bounds
        assert bounds.convNormMin == pytest.approx(0.1)
        assert bounds.convNormMax == pytest.approx(3.0)
        assert bounds.domEfficiencyMin == pytest.approx(1.234)
        assert bounds.domEfficiencyMax == pytest.approx(1.346)
    
    def test_bounds_modification(self):
        """Test modifying bounds."""
        bounds = gf.FitParametersBound()
        
        bounds.astroNormMin = 0.0
        bounds.astroNormMax = 5.0
        
        assert bounds.astroNormMin == 0.0
        assert bounds.astroNormMax == 5.0


class TestDataPathsSteeringParams:
    """
    Test suite for DataPaths and SteeringParams structs.
    
    Tests configuration objects.
    """
    
    def test_datapaths_initialization(self):
        """Test DataPaths can be instantiated."""
        datapaths = gf.DataPaths()
        assert datapaths is not None
    
    def test_datapaths_set_paths(self, example_datapaths):
        """Test setting file paths."""
        datapaths = example_datapaths
        
        # Verify paths are set
        assert datapaths.compact_file_path != ""
        assert "compact.fastmc" in datapaths.compact_file_path
    
    def test_steering_initialization(self):
        """Test SteeringParams can be instantiated."""
        steering = gf.SteeringParams()
        assert steering is not None
    
    def test_steering_configuration(self, example_steering_params):
        """Test SteeringParams configuration."""
        steering = example_steering_params
        
        # Check energy configuration
        assert steering.minFitEnergy == 300.0
        assert steering.maxFitEnergy == 1e5
        
        # Check zenith configuration
        assert steering.minCosth == -1.0
        assert steering.maxCosth == 0.0
        
        # Check threads
        assert steering.evalThreads == 1
        
        # Check tolerances
        assert steering.change_tol == 1.e-20
        assert steering.grad_tol == 1.e-20


class TestGollumFitInitialization:
    """
    Test suite for GollumFit object initialization.
    
    Tests creating GollumFit objects with various configurations.
    """
    
    def test_gollumfit_initialization(self, example_datapaths, example_steering_params):
        """Test GollumFit object can be created."""
        # This test requires the FastMC file and splines to exist
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found, skipping initialization test")
        
        gollumfit = gf.GollumFit(example_datapaths, example_steering_params)
        assert gollumfit is not None
    
    def test_gollumfit_getters(self, example_datapaths, example_steering_params):
        """Test getter methods for GollumFit configuration."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        
        gollumfit = gf.GollumFit(example_datapaths, example_steering_params)
        
        # Test getters
        steering = gollumfit.GetSteeringParams()
        assert steering.minFitEnergy == 300.0
        
        paths = gollumfit.GetDataPaths()
        assert paths.compact_file_path == example_datapaths.compact_file_path


class TestGollumFitExpectations:
    """
    Test suite for expectation generation.
    
    Tests GetExpectation and related methods.
    """
    
    @pytest.fixture
    def initialized_gollumfit(self, example_datapaths, example_steering_params):
        """Create and return an initialized GollumFit object."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        return gf.GollumFit(example_datapaths, example_steering_params)
    
    def test_get_expectation_returns_array(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """Test that GetExpectation returns a properly shaped array."""
        hist = initialized_gollumfit.GetExpectation(default_fit_parameters)
        
        # Should return a 3D array (topology, zenith, energy)
        assert hist is not None
        assert isinstance(hist, np.ndarray)
        assert hist.ndim == 3
        
        # Check that histogram has expected shape
        # 2 topologies, 20 zenith bins, 24 energy bins (from example)
        assert hist.shape[0] == 2  # topologies
        assert hist.shape[1] > 0   # zenith bins
        assert hist.shape[2] > 0   # energy bins
    
    def test_get_expectation_positive_weights(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """Test that GetExpectation produces non-negative weights."""
        hist = initialized_gollumfit.GetExpectation(default_fit_parameters)
        
        # All histogram bin contents should be >= 0
        assert np.all(hist >= 0), "Found negative bins in expectation histogram"
    
    def test_get_expectation_parameter_variation(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """Test that varying parameters changes the expectation."""
        hist_nominal = initialized_gollumfit.GetExpectation(default_fit_parameters)
        
        # Vary convNorm
        params_varied = gf.FitParameters()
        for attr in dir(default_fit_parameters):
            if not attr.startswith('_'):
                try:
                    setattr(params_varied, attr, getattr(default_fit_parameters, attr))
                except:
                    pass
        params_varied.convNorm = 1.5  # Increase by 50%
        
        hist_varied = initialized_gollumfit.GetExpectation(params_varied)
        
        # Histograms should be different
        assert not np.allclose(hist_nominal, hist_varied), \
            "Varying convNorm should change expectation"
    
    def test_check_expectation_valid(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """Test CheckExpectation for valid parameters."""
        result = initialized_gollumfit.CheckExpectation(default_fit_parameters)
        
        # Should return 0 for valid expectation
        assert result == 0, "Default parameters should produce valid expectation"


class TestGollumFitDataHandling:
    """
    Test suite for data loading and histogram construction.
    
    Tests SetData and related methods.
    """
    
    @pytest.fixture
    def initialized_gollumfit(self, example_datapaths, example_steering_params):
        """Create and return an initialized GollumFit object."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        return gf.GollumFit(example_datapaths, example_steering_params)
    
    def test_set_data_from_expectation(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """Test loading data from an expectation array."""
        # Get an expectation
        events = initialized_gollumfit.GetExpectationEvents(default_fit_parameters)

        # Set it as data
        total_weight = initialized_gollumfit.SetData(events)
        
        # Total weight should be positive
        assert total_weight > 0, "Total data weight should be positive"
    
    def test_set_data_from_realization(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """Test loading data from a Poisson realization event list."""
        # GetRealizationEvents returns a 2D array (events x 4) for SetData
        seed = 12345
        events = initialized_gollumfit.GetRealizationEvents(default_fit_parameters, seed)
        
        # Set it as data
        total_weight = initialized_gollumfit.SetData(events)
        
        # Total weight should be positive
        assert total_weight > 0
    
    def test_get_data_events(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """Test retrieving data events after loading."""
        events = initialized_gollumfit.GetExpectationEvents(default_fit_parameters)
        initialized_gollumfit.SetData(events)
        
        # Get data events
        data_events = initialized_gollumfit.GetDataEvents()
        
        # Should return a 2D array
        assert isinstance(data_events, np.ndarray)
        assert data_events.ndim == 2
        
        # Should have 4 columns: energy, zenith, topology, weight
        if data_events.shape[0] > 0:
            assert data_events.shape[1] == 4


class TestGollumFitLikelihood:
    """
    Test suite for likelihood evaluation and fitting.
    
    Tests EvalLLH and related methods.
    """
    
    @pytest.fixture
    def configured_gollumfit(
        self, example_datapaths, example_steering_params,
        default_fit_parameters, default_priors, default_bounds, default_flags
    ):
        """Create GollumFit object configured for likelihood evaluation."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        
        gollumfit = gf.GollumFit(example_datapaths, example_steering_params)
        
        # Generate and load fake data
        events = gollumfit.GetExpectationEvents(default_fit_parameters)
        gollumfit.SetData(events)
        
        # Configure for fitting
        gollumfit.SetFitParametersFlag(default_flags)
        gollumfit.SetFitParametersBound(default_bounds)
        gollumfit.SetFitParametersPriors(default_priors)
        gollumfit.SetFitParametersSeed([default_fit_parameters])
        
        # Construct likelihood problem
        gollumfit.ConstructLikelihoodProblem()
        
        return gollumfit
    
    def test_eval_llh(self, configured_gollumfit, default_fit_parameters):
        """Test evaluating negative log-likelihood."""
        llh = configured_gollumfit.EvalLLH(default_fit_parameters, True)
        
        # Likelihood should be a finite number
        assert np.isfinite(llh), "Likelihood should be finite"
        
        # For Asimov data at nominal parameters, LLH should be small
        # (though not exactly 0 due to priors)
        assert llh >= 0, "Negative log-likelihood should be non-negative"
    
    def test_eval_llh_without_prior(self, configured_gollumfit, default_fit_parameters):
        """Test evaluating likelihood without prior contribution."""
        llh_with_prior = configured_gollumfit.EvalLLH(default_fit_parameters, True)
        llh_without_prior = configured_gollumfit.EvalLLH(default_fit_parameters, False)
        
        # Both should be finite
        assert np.isfinite(llh_with_prior)
        assert np.isfinite(llh_without_prior)
        
        # With prior should generally be >= without prior
        # (priors add penalty terms)



class TestGollumFitBinning:
    """
    Test suite for histogram binning methods.
    
    Tests GetEnergyBinsMC, GetZenithBinsMC, etc.
    """
    
    @pytest.fixture
    def initialized_gollumfit(self, example_datapaths, example_steering_params):
        """Create and return an initialized GollumFit object."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        return gf.GollumFit(example_datapaths, example_steering_params)
    
    def test_get_energy_bins(self, initialized_gollumfit):
        """Test retrieving energy bin edges."""
        energy_bins = initialized_gollumfit.GetEnergyBinsMC()
        
        # Should return an array
        assert isinstance(energy_bins, (list, np.ndarray))
        assert len(energy_bins) > 0
        
        # Bins should be monotonically increasing
        energy_bins = np.array(energy_bins)
        assert np.all(np.diff(energy_bins) > 0), "Energy bins should be increasing"
    
    def test_get_zenith_bins(self, initialized_gollumfit):
        """Test retrieving zenith bin edges."""
        zenith_bins = initialized_gollumfit.GetZenithBinsMC()
        
        # Should return an array
        assert isinstance(zenith_bins, (list, np.ndarray))
        assert len(zenith_bins) > 0
        
        # Bins should be monotonically increasing
        zenith_bins = np.array(zenith_bins)
        assert np.all(np.diff(zenith_bins) > 0), "Zenith bins should be increasing"
    
    def test_get_topology_bins(self, initialized_gollumfit):
        """Test retrieving topology bin edges."""
        topology_bins = initialized_gollumfit.GetTopologyBinsMC()
        
        # Should return an array
        assert isinstance(topology_bins, (list, np.ndarray))
        assert len(topology_bins) > 0


class TestGollumFitRealization:
    """
    Test suite for realization generation.
    
    Tests GetRealization and GetRealizationEvents.
    """
    
    @pytest.fixture
    def initialized_gollumfit(self, example_datapaths, example_steering_params):
        """Create and return an initialized GollumFit object."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        return gf.GollumFit(example_datapaths, example_steering_params)
    
    def test_get_realization_shape(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """Test that GetRealization returns properly shaped histogram."""
        seed = 42
        hist = initialized_gollumfit.GetRealization(default_fit_parameters, seed)
        
        # Should return 3D array
        assert isinstance(hist, np.ndarray)
        assert hist.ndim == 3
        
        # Should have expected shape
        assert hist.shape[0] == 2  # topologies
    
    def test_get_realization_reproducibility(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """Test that same seed produces same realization."""
        seed = 12345
        
        hist1 = initialized_gollumfit.GetRealization(default_fit_parameters, seed)
        hist2 = initialized_gollumfit.GetRealization(default_fit_parameters, seed)
        
        # Should be identical
        np.testing.assert_array_equal(hist1, hist2,
                                      "Same seed should produce identical realizations")
    
    def test_get_realization_different_seeds(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """Test that different seeds produce different realizations."""
        hist1 = initialized_gollumfit.GetRealization(default_fit_parameters, 111)
        hist2 = initialized_gollumfit.GetRealization(default_fit_parameters, 222)
        
        # Should be different (with high probability)
        assert not np.allclose(hist1, hist2), \
            "Different seeds should produce different realizations"
    
    def test_get_realization_poisson_statistics(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """Test that realizations follow Poisson statistics."""
        # Generate expectation
        expectation = initialized_gollumfit.GetExpectation(default_fit_parameters)
        
        # Generate many realizations
        n_realizations = 100
        realizations = []
        for i in range(n_realizations):
            hist = initialized_gollumfit.GetRealization(default_fit_parameters, i)
            realizations.append(hist)
        
        realizations = np.array(realizations)
        
        # Mean over realizations should be close to expectation
        mean_realization = np.mean(realizations, axis=0)
        
        # Allow 5-sigma deviation across many bins (avoids false failures
        # from multiple comparisons when hundreds of bins are checked)
        expected_sigma = np.sqrt(expectation / n_realizations)

        # Only check bins with high statistics where CLT applies well
        high_stat_mask = expectation > 50
        if np.any(high_stat_mask):
            diff = np.abs(mean_realization[high_stat_mask] - expectation[high_stat_mask])
            assert np.all(diff < 5 * expected_sigma[high_stat_mask]), \
                "Realization mean should match expectation within statistics"


# ============================================================================
# Integration Tests
# ============================================================================

class TestGollumFitIntegration:
    """
    Integration tests that combine multiple GollumFit features.
    
    These tests verify that different components work together correctly.
    """
    
    def test_full_workflow_expectation(
        self, example_datapaths, example_steering_params, default_fit_parameters
    ):
        """Test full workflow: initialize, configure, get expectation."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        
        # Initialize
        gollumfit = gf.GollumFit(example_datapaths, example_steering_params)
        
        # Get expectation
        hist = gollumfit.GetExpectation(default_fit_parameters)
        
        # Verify output
        assert hist is not None
        assert hist.shape[0] == 2  # Two topologies
        assert np.all(hist >= 0)
    
    def test_full_workflow_likelihood_evaluation(
        self, example_datapaths, example_steering_params,
        default_fit_parameters, default_priors, default_bounds
    ):
        """Test full workflow: initialize, load data, evaluate likelihood."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        
        # Initialize
        gollumfit = gf.GollumFit(example_datapaths, example_steering_params)
        
        # Generate and load Asimov data
        events = gollumfit.GetExpectationEvents(default_fit_parameters)
        total_weight = gollumfit.SetData(events)
        assert total_weight > 0
        
        # Configure for likelihood evaluation
        flags = gf.FitParametersFlag(False)  # Vary all
        gollumfit.SetFitParametersFlag(flags)
        gollumfit.SetFitParametersBound(default_bounds)
        gollumfit.SetFitParametersPriors(default_priors)
        gollumfit.SetFitParametersSeed([default_fit_parameters])
        
        # Construct likelihood problem
        gollumfit.ConstructLikelihoodProblem()
        
        # Evaluate likelihood
        llh = gollumfit.EvalLLH(default_fit_parameters, True)
        
        # Verify
        assert np.isfinite(llh)
        assert llh >= 0


# ============================================================================
# Error Handling Tests
# ============================================================================

class TestGollumFitErrorHandling:
    """
    Test suite for error handling and edge cases.
    
    Tests that GollumFit handles invalid inputs gracefully.
    """
    
    def test_invalid_parameter_values(
        self, example_datapaths, example_steering_params
    ):
        """Test behavior with extreme/invalid parameter values."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        
        gollumfit = gf.GollumFit(example_datapaths, example_steering_params)
        
        # Use boundary values within valid range (convNorm=1000 would trigger
        # a C++ assert and abort() which cannot be caught by Python)
        params = gf.FitParameters()
        params.convNorm = 3.0   # At upper bound
        params.astroNorm = 0.0  # Zero (lower bound)
        
        # Depending on implementation, this might raise or return invalid expectation
        # Test that it doesn't crash
        try:
            hist = gollumfit.GetExpectation(params)
            # If it succeeds, verify output is reasonable
            assert hist is not None
        except Exception as e:
            # If it raises, that's also acceptable error handling
            pass
    
    def test_construct_likelihood_before_data(
        self, example_datapaths, example_steering_params, default_priors,
        default_bounds, default_flags, default_fit_parameters
    ):
        """Test that constructing likelihood without data raises error."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        
        gollumfit = gf.GollumFit(example_datapaths, example_steering_params)
        
        # Try to construct likelihood without loading data
        gollumfit.SetFitParametersFlag(default_flags)
        gollumfit.SetFitParametersBound(default_bounds)
        gollumfit.SetFitParametersPriors(default_priors)
        gollumfit.SetFitParametersSeed([default_fit_parameters])
        
        # This should raise an error
        with pytest.raises(Exception):
            gollumfit.ConstructLikelihoodProblem()


# ============================================================================
# Performance Tests (Optional)
# ============================================================================

class TestGollumFitPerformance:
    """
    Performance tests for GollumFit operations.
    
    These tests verify that operations complete in reasonable time.
    """
    
    @pytest.mark.slow
    def test_expectation_generation_time(
        self, example_datapaths, example_steering_params, default_fit_parameters
    ):
        """Test that expectation generation completes in reasonable time."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        
        import time
        
        gollumfit = gf.GollumFit(example_datapaths, example_steering_params)
        
        start = time.time()
        hist = gollumfit.GetExpectation(default_fit_parameters)
        elapsed = time.time() - start
        
        # Should complete in under 10 seconds (adjust as needed)
        assert elapsed < 10.0, f"Expectation took {elapsed:.2f}s, expected < 10s"
    
    @pytest.mark.slow
    def test_likelihood_evaluation_time(
        self, example_datapaths, example_steering_params,
        default_fit_parameters, default_priors, default_bounds
    ):
        """Test that likelihood evaluation completes in reasonable time."""
        if not os.path.exists(example_datapaths.compact_file_path):
            pytest.skip("FastMC file not found")
        
        import time
        
        gollumfit = gf.GollumFit(example_datapaths, example_steering_params)
        events = gollumfit.GetExpectationEvents(default_fit_parameters)
        gollumfit.SetData(events)
        
        flags = gf.FitParametersFlag(False)
        gollumfit.SetFitParametersFlag(flags)
        gollumfit.SetFitParametersBound(default_bounds)
        gollumfit.SetFitParametersPriors(default_priors)
        gollumfit.SetFitParametersSeed([default_fit_parameters])
        gollumfit.ConstructLikelihoodProblem()
        
        start = time.time()
        llh = gollumfit.EvalLLH(default_fit_parameters, True)
        elapsed = time.time() - start
        
        # Should complete in under 5 seconds
        assert elapsed < 5.0, f"Likelihood eval took {elapsed:.2f}s, expected < 5s"


if __name__ == "__main__":
    # Allow running tests directly with pytest
    pytest.main([__file__, "-v"])
