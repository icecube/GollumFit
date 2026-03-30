# GollumFit Unit Tests

This directory contains comprehensive unit tests for the GollumFit framework.

## Overview

GollumFit is a binned-likelihood fitting framework for neutrino telescope analyses. These unit tests verify the correctness of all major components:

- **Parameter Structures**: FitParameters, Priors, Bounds, Flags
- **Data Handling**: Loading data/MC, histogram construction  
- **Expectations**: Generating expected event distributions
- **Likelihood**: Evaluating negative log-likelihood
- **Realizations**: Generating Poisson-fluctuated pseudo-data
- **Error Handling**: Testing edge cases and invalid inputs

## Test Structure

### Test Files

- `test_gollumfit.py` - Main unit test suite (1000+ lines)

### Test Classes

1. **TestFitParameters** - Tests for the 38-parameter fit parameter struct
2. **TestFitParametersFlag** - Tests for parameter fixing/varying flags
3. **TestPriors** - Tests for Gaussian and Uniform prior specifications
4. **TestFitParametersBound** - Tests for parameter boundary constraints
5. **TestDataPathsSteeringParams** - Tests for configuration objects
6. **TestGollumFitInitialization** - Tests for GollumFit object creation
7. **TestGollumFitExpectations** - Tests for expectation generation
8. **TestGollumFitDataHandling** - Tests for data loading and access
9. **TestGollumFitLikelihood** - Tests for likelihood evaluation (EvalLLH)
10. **TestGollumFitBinning** - Tests for histogram binning retrieval
11. **TestGollumFitRealization** - Tests for pseudo-data generation
12. **TestGollumFitIntegration** - End-to-end workflow tests
13. **TestGollumFitErrorHandling** - Tests for error conditions
14. **TestGollumFitPerformance** - Performance benchmarks (optional)

## Requirements

### Python Dependencies

```bash
pip install pytest numpy h5py
```

### GollumFit Installation

The tests require GollumFit Python bindings (`GollumFitPy`) to be installed. Follow the GollumFit installation instructions to build and install the Python module.

### Test Data

Some tests require example data files:
- FastMC file: `examples/FastMC/compact.fastmc`
- Spline files in `resources/Splines/`
- Correlation matrices in `resources/correlation_matrices/`

Tests will automatically skip if required files are not found.

## Running the Tests

### Run All Tests

```bash
# From the GollumFit root directory
pytest tests/test_gollumfit.py -v

# Or from the tests directory
cd tests
pytest test_gollumfit.py -v
```

### Run Specific Test Class

```bash
# Test only parameter structures
pytest tests/test_gollumfit.py::TestFitParameters -v

# Test only expectation generation
pytest tests/test_gollumfit.py::TestGollumFitExpectations -v

# Test only likelihood functionality
pytest tests/test_gollumfit.py::TestGollumFitLikelihood -v
```

### Run Specific Test

```bash
# Run a single test function
pytest tests/test_gollumfit.py::TestFitParameters::test_fitparameters_all_38_parameters -v
```

### Run with Coverage

```bash
# Generate coverage report
pytest tests/test_gollumfit.py --cov=GollumFitPy --cov-report=html
```

### Skip Slow Tests

```bash
# Skip performance tests marked as slow
pytest tests/test_gollumfit.py -v -m "not slow"
```

## Test Documentation

Each test includes comprehensive docstrings explaining:
- **Purpose**: What the test verifies
- **Setup**: What fixtures or data it uses
- **Assertions**: What conditions must be true
- **Expected Behavior**: What should happen

### Example Test Documentation

```python
def test_get_expectation_returns_array(
    self, initialized_gollumfit, default_fit_parameters
):
    """
    Test that GetExpectation returns a properly shaped array.
    
    This test verifies that the GetExpectation method:
    1. Returns a numpy array
    2. Has the correct dimensionality (3D: topology, zenith, energy)
    3. Has expected shape matching the binning configuration
    
    Fixtures:
        initialized_gollumfit: GollumFit object with loaded FastMC
        default_fit_parameters: FitParameters with nominal values
    
    Assertions:
        - Result is a numpy ndarray
        - Array has 3 dimensions
        - Shape[0] == 2 (number of event topologies)
        - Shape[1], Shape[2] > 0 (zenith and energy bins)
    """
    hist = initialized_gollumfit.GetExpectation(default_fit_parameters)
    
    assert hist is not None
    assert isinstance(hist, np.ndarray)
    assert hist.ndim == 3
    assert hist.shape[0] == 2  # topologies
    assert hist.shape[1] > 0   # zenith bins
    assert hist.shape[2] > 0   # energy bins
```

## Pytest Fixtures

Fixtures provide reusable test data and configurations:

### Parameter Fixtures
- `default_fit_parameters` - FitParameters with nominal values
- `default_priors` - Priors with standard Gaussian/Uniform settings
- `default_bounds` - FitParametersBound with typical constraints
- `default_flags` - FitParametersFlag with all parameters varying

### Configuration Fixtures
- `example_datapaths` - DataPaths pointing to test resources
- `example_steering_params` - SteeringParams with test binning

### GollumFit Fixtures
- `initialized_gollumfit` - Basic GollumFit object
- `configured_gollumfit` - GollumFit ready for likelihood evaluation

## Expected Test Results

### Quick Tests (< 1 second each)
- Parameter structure tests
- Initialization tests
- Data loading tests

### Medium Tests (1-5 seconds each)
- Expectation generation
- Histogram binning
- Likelihood evaluation (single point)

### Slow Tests (> 5 seconds)
- Extensive realizations
- Performance benchmarks

## Common Test Failures

### Missing Data Files

```
SKIPPED [1] test_gollumfit.py:xxx: FastMC file not found
```
**Solution**: Generate FastMC file using `examples/FastMC/generate_fastMC.py`

### Import Errors

```
ModuleNotFoundError: No module named 'GollumFitPy'
```
**Solution**: Build and install GollumFit Python bindings

### Assertion Failures

If tests fail, check:
1. GollumFit version compatibility
2. File paths in fixtures match your installation
3. Correlationmatrices and splines are present

## Extending the Tests

### Adding New Tests

1. **Identify Component**: Determine which test class is appropriate
2. **Write Test Function**: Start with `test_` prefix
3. **Add Docstring**: Clearly document purpose and assertions
4. **Use Fixtures**: Reuse existing fixtures where possible
5. **Add Assertions**: Use descriptive assertion messages

### Example New Test

```python
class TestGollumFitExpectations:
    def test_expectation_normalization(
        self, initialized_gollumfit, default_fit_parameters
    ):
        """
        Test that expectation normalizations are consistent.
        
        Verifies that varying normalization parameters (convNorm, 
        astroNorm) scales the expectation proportionally.
        """
        # Get baseline
        hist_base = initialized_gollumfit.GetExpectation(default_fit_parameters)
        total_base = np.sum(hist_base)
        
        # Scale convNorm by 2x
        params_scaled = default_fit_parameters
        params_scaled.convNorm = 2.0
        hist_scaled = initialized_gollumfit.GetExpectation(params_scaled)
        total_scaled = np.sum(hist_scaled)
        
        # Should scale proportionally
        assert np.isclose(total_scaled, 2.0 * total_base, rtol=0.01)
```

### Adding Test Fixtures

```python
@pytest.fixture
def custom_fit_parameters():
    """Create FitParameters with custom values for specific test scenarios."""
    params = gf.FitParameters()
    # Set custom values
    params.convNorm = 2.0
    params.astroNorm = 1.5
    return params
```

## Continuous Integration

These tests can be integrated into CI/CD pipelines:

### GitHub Actions Example

```yaml
name: Run GollumFit Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: '3.9'
      - name: Install dependencies
        run: |
          pip install pytest numpy h5py
          # Build and install GollumFit...
      - name: Run tests
        run: |
          pytest tests/test_gollumfit.py -v
```

## Test Coverage Goals

Target coverage by component:
- **Parameter Structures**: 100% (simple data structures)
- **Initialization**: 90%+ (various configuration paths)
- **Expectation Generation**: 85%+ (core functionality)
- **Likelihood Evaluation**: 80%+ (complex numerical operations)
- **Error Handling**: 75%+ (edge cases)

Check coverage with:
```bash
pytest tests/test_gollumfit.py --cov=GollumFitPy --cov-report=term-missing
```

## Debugging Failed Tests

### Enable Verbose Output

```bash
pytest tests/test_gollumfit.py::TestName::test_function -vv -s
```

### Print GollumFit Internal State

```python
def test_debug_example(self, initialized_gollumfit):
    """Debug test with internal state printing."""
    steering = initialized_gollumfit.GetSteeringParams()
    print(f"\nEnergy bins: {steering.minFitEnergy} - {steering.maxFitEnergy}")
    print(f"Zenith bins: {steering.minCosth} - {steering.maxCosth}")
    # ... rest of test
```

### Use PDB Debugger

```bash
pytest tests/test_gollumfit.py::test_name --pdb
```

## Contributing

When contributing tests:
1. Follow existing naming conventions
2. Add comprehensive docstrings
3. Use descriptive variable names
4. Group related tests in classes
5. Use fixtures for shared setup
6. Test both success and failure paths
7. Include edge cases

## Contact

For questions about the tests or GollumFit:
- See GollumFit documentation: https://docs.icecube.aq/gollumfit/
- Open an issue on the GollumFit repository
- Contact the GollumFit development team

## License

These tests are distributed under the same license as GollumFit (see LICENSE file in repository root).
