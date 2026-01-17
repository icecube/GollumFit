# Code Structure {#codestructure}

## Overview

The GollumFit program is organized in a straightforward, modular architecture:

### Core Components

1. `GollumFit` Class: The central class that handles the main functionality of the code, accepting different inputs depending on the desired functionality.

2. Helper Structures: A collection of structures that organize variables passed throughout the program:
   - `FitParameters`: Contains nuisance parameter values
   - `FitResult`: Stores minimization results
   - `SteeringParams`: Defines analysis configuration (binning, energy ranges, etc.)
   - `DataPaths`: Specifies locations of input files and resources

3. `WeighterMaker` Class: Handles event re-weighting by constructing individual weighters that calculate weight adjustments for each event based on nuisance parameter values.

### Component Diagram

<img src="GFdiagram.png" alt="component diagram" style="width: 80%; max-width: 900px;">
*Figure: High-level graph showing the relationship between main components.* 

## Nuisance Parameters

GollumFit includes handling of common nuisance parameters for all neutrino telescopes, as well as IceCube-specific parameters. For detailed physics descriptions and prior values, refer to the GollumFit publication.

### Common Neutrino Telescope Parameters
**Flux Parameters:**
- `convNorm` - Conventional atmospheric flux normalization
- `promptNorm` - Prompt atmospheric flux normalization
- `astroNorm` - Astrophysical flux normalization
- `astroDeltaGamma` - Astrophysical spectral index variation (first break)
- `astroDeltaGammaSec` - Astrophysical spectral index variation (second break)
- `astroPivot` - Astrophysical spectrum break energy

**Atmospheric Flux Systematics:**
- `zenithCorrection` - Atmospheric density correction
- `kaonLosses` - Kaon interaction losses
- `cosmicRay1` through `cosmicRay6` - Global Spline Fit (GSF) cosmic ray parameters

**Hadronic Interaction Parameters:**
- `hadronicHEkp`, `hadronicHEkm` - High-energy kaon production ($K^+$, $K^-$ at 158 GeV)
- `hadronicVHE1pip`, `hadronicVHE1pim` - Very-high-energy pion production ($\pi^+$, $\pi^-$ at 20 TeV)
- `hadronicVHE3kp`, `hadronicVHE3km` - Very-high-energy kaon production ($K^+$, $K^-$ at 2 PeV)
- `hadronicVHE3pip`, `hadronicVHE3pim` - Very-high-energy pion production ($\pi^+$, $\pi^-$ at 2 PeV)
- `hadronicVHE3p`, `hadronicVHE3n` - Very-high-energy nucleon production (p, n at 2 PeV)

**Cross Section and Interaction:**
- `NeutrinoAntineutrinoRatio` - Ratio of neutrino to antineutrino flux
- `nuxs` - Neutrino cross section attenuation
- `nubarxs` - Antineutrino cross section attenuation

### IceCube-Specific Detector Parameters

**Detector Efficiency:**
- `domEfficiency` - Digital Optical Module efficiency
- `holeiceForward` - Hole ice forward scattering

**Ice Model Parameters:**
- `icegrad0` through `icegrad4` - Ice absorption gradients
- `icegrad5` through `icegrad8` - Ice phase (scattering) gradients

<!-- | Parameter                                                                                         | `GollumFit` Variable        | Prior Value ± Width       | Weighting Method    |
|:--------------------------------------------------------------------------------------------------|:----------------------------|:--------------------------|:--------------------|
| **Common to all neutrino telescopes**                                                             |                             |                           |                     |
| \( \textrm{convNorm}  \)                                                                             | `convNorm`                  | \( 1.0 \pm 0.2 \)           | scale factor        |
| \( \rho_{\textrm{atm}} \)                                                                           | `zenithCorrection`          | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \sigma_{\textrm{K-Air}} \)                                                                        | `kaonLosses`                | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{K}_{158G}^{+} \)                                                                          | `hadronicHEkp`              | \(  0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{K}_{158G}^{-} \)                                                                          | `hadronicHEkm`              | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \pi_{20T}^{+} \)                                                                                  | `hadronicVHE1pip`           | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \pi_{20T}^{-} \)                                                                                  | `hadronicVHE1pim`           | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{K}_{2P}^{+} \)                                                                            | `hadronicVHE3kp`            | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{K}_{2P}^{-} \)                                                                            | `hadronicVHE3km`            | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \pi_{2P}^{+} \)                                                                                   | `hadronicVHE3pip`           | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \pi_{2P}^{-} \)                                                                                   | `hadronicVHE3pim`           | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{p}_{2P} \)                                                                                | `hadronicVHE3p`             | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{n}_{2P} \)                                                                                | `hadronicVHE3n`             | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{GSF}_1 \)                                                                                 | `cosmicRay1`                | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{GSF}_2 \)                                                                                 | `cosmicRay2`                | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{GSF}_3 \)                                                                                 | `cosmicRay3`                | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{GSF}_4 \)                                                                                 | `cosmicRay4`                | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{GSF}_5 \)                                                                                 | `cosmicRay5`                | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{GSF}_6 \)                                                                                 | `cosmicRay6`                | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \Phi^{\textrm{HE}}/10^{-18}\,\textrm{GeV}^{-1}\,\textrm{sr}^{-1}\,\textrm{s}^{-1}\,\textrm{cm}^{-2} \) | `astroNorm`                 | \( 0.787 \pm 0.36 \)        | scale factor        |
| \( \Delta\gamma_{1}^{\textrm{HE}} \)                                                                 | `astroDeltaGamma`           | \( 0.0 \pm 0.36 \)          | power law formula   |
| \( \Delta\gamma_{2}^{\textrm{HE}} \)                                                                 | `astroDeltaGammaSec`        | \( 0.0 \pm 0.36 \)          | power law formula   |
| \( \log_{10}\left(\textrm{E}_{\textrm{break}}^{\textrm{HE}}/\textrm{GeV}\right) \)                    | `astroPivot`                | \( 5.0 \pm 1.0 \)           | power law formula   |
| \( \textrm{promptNorm} \)                                                                             | `promptNorm`                | \( 1.0 \pm 1.0 \)           | scale factor        |
| \( \nu/\bar{\nu} \)                                                                                  | `NeutrinoAntineutrinoRatio` | \( 1.0 \pm 1.0 \)           | scale factor        |
| \( \nu\textrm{ Att} \)                                                                               | `nuxs`                      | \( 1.0 \pm 0.1 \)           | spline              |
| \( \bar{\nu}\textrm{ Att} \)                                                                         | `nubarxs`                   | \( 1.0 \pm 0.1 \)           | spline              |
| **IceCube-specific Monte Carlo parameters**                                                        |                             |                           |                     |
| \( \textrm{DOM} \: \textrm{eff} \)                                                                  | `domEfficiency`             | \( 1.27 \pm 0.123 \)        | spline              |
| \( \textrm{Hole Ice} \)                                                                              | `holeiceForward`            | \( -1.0 \pm 10.0 \)         | spline              |
| \( \textrm{Ice A}_{0} \)                                                                             | `icegrad0`                  | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{Ice A}_{1} \)                                                                             | `icegrad1`                  | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{Ice A}_{2} \)                                                                             | `icegrad2`                  | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{Ice A}_{3} \)                                                                             | `icegrad3`                  | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{Ice A}_{4} \)                                                                             | `icegrad4`                  | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{Ice Phs}_{1} \)                                                                           | `icegrad5`                  | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{Ice Phs}_{2} \)                                                                           | `icegrad6`                  | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{Ice Phs}_{3} \)                                                                           | `icegrad7`                  | \( 0.0 \pm 1.0 \)           | gradient            |
| \( \textrm{Ice Phs}_{4} \)                                                                           | `icegrad8`                  | \( 0.0 \pm 1.0 \)           | gradient            | -->


<!-- <table border="1" cellspacing="0" cellpadding="4">
  <tr>
    <th>Parameter</th>
    <th><code>GollumFit</code> Variable</th>
    <th>Prior Value ± Width</th>
    <th>Weighting Method</th>
  </tr>
  <tr>
    <td colspan="4" align="center"><strong>Common to all neutrino telescopes</strong></td>
  </tr>
  <tr>
    <td>\f$\textrm{convNorm}\f$</td>
    <td><code>convNorm</code></td>
    <td>\f$1.0 \pm 0.2\f$</td>
    <td>scale factor</td>
  </tr>
  <tr>
    <td>\f$\rho_{\textrm{atm}}\f$</td>
    <td><code>zenithCorrection</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\sigma_{\textrm{K-Air}}\f$</td>
    <td><code>kaonLosses</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{K}_{158G}^{+}\f$</td>
    <td><code>hadronicHEkp</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{K}_{158G}^{-}\f$</td>
    <td><code>hadronicHEkm</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\pi_{20T}^{+}\f$</td>
    <td><code>hadronicVHE1pip</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\pi_{20T}^{-}\f$</td>
    <td><code>hadronicVHE1pim</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{K}_{2P}^{+}\f$</td>
    <td><code>hadronicVHE3kp</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{K}_{2P}^{-}\f$</td>
    <td><code>hadronicVHE3km</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\pi_{2P}^{+}\f$</td>
    <td><code>hadronicVHE3pip</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\pi_{2P}^{-}\f$</td>
    <td><code>hadronicVHE3pim</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{p}_{2P}\f$</td>
    <td><code>hadronicVHE3p</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{n}_{2P}\f$</td>
    <td><code>hadronicVHE3n</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{GSF}_1\f$</td>
    <td><code>cosmicRay1</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{GSF}_2\f$</td>
    <td><code>cosmicRay2</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{GSF}_3\f$</td>
    <td><code>cosmicRay3</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{GSF}_4\f$</td>
    <td><code>cosmicRay4</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{GSF}_5\f$</td>
    <td><code>cosmicRay5</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{GSF}_6\f$</td>
    <td><code>cosmicRay6</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\Phi^{\textrm{HE}}/10^{-18}\,\textrm{GeV}^{-1}\,\textrm{sr}^{-1}\,\textrm{s}^{-1}\,\textrm{cm}^{-2}\f$</td>
    <td><code>astroNorm</code></td>
    <td>\f$0.787 \pm 0.36\f$</td>
    <td>scale factor</td>
  </tr>
  <tr>
    <td>\f$\Delta\gamma_{1}^{\textrm{HE}}\f$</td>
    <td><code>astroDeltaGamma</code></td>
    <td>\f$0.0 \pm 0.36\f$</td>
    <td>power law formula</td>
  </tr>
  <tr>
    <td>\f$\Delta\gamma_{2}^{\textrm{HE}}\f$</td>
    <td><code>astroDeltaGammaSec</code></td>
    <td>\f$0.0 \pm 0.36\f$</td>
    <td>power law formula</td>
  </tr>
  <tr>
    <td>\f$\log_{10}\left(\textrm{E}_{\textrm{break}}^{\textrm{HE}}/\textrm{GeV}\right)\f$</td>
    <td><code>astroPivot</code></td>
    <td>\f$5.0 \pm 1.0\f$</td>
    <td>power law formula</td>
  </tr>
  <tr>
    <td>\f$\textrm{promptNorm}\f$</td>
    <td><code>promptNorm</code></td>
    <td>\f$1.0 \pm 1.0\f$</td>
    <td>scale factor</td>
  </tr>
  <tr>
    <td>\f$\nu/ \bar{\nu}\f$</td>
    <td><code>NeutrinoAntineutrinoRatio</code></td>
    <td>\f$1.0 \pm 1.0\f$</td>
    <td>scale factor</td>
  </tr>
  <tr>
    <td>\f$\nu\textrm{ Att}\f$</td>
    <td><code>nuxs</code></td>
    <td>\f$1.0 \pm 0.1\f$</td>
    <td>spline</td>
  </tr>
  <tr>
    <td>\f$\bar{\nu}\textrm{ Att}\f$</td>
    <td><code>nubarxs</code></td>
    <td>\f$1.0 \pm 0.1\f$</td>
    <td>spline</td>
  </tr>
  <tr>
    <td colspan="4" align="center"><strong>IceCube-specific Monte Carlo parameters</strong></td>
  </tr>
  <tr>
    <td>\f$\textrm{DOM}\:{\textrm{eff}}\f$</td>
    <td><code>domEfficiency</code></td>
    <td>\f$1.27 \pm 0.123\f$</td>
    <td>spline</td>
  </tr>
  <tr>
    <td>\f$\textrm{Hole Ice}\f$</td>
    <td><code>holeiceForward</code></td>
    <td>\f$-1.0 \pm 10.0\f$</td>
    <td>spline</td>
  </tr>
  <tr>
    <td>\f$\textrm{Ice A}_{0}\f$</td>
    <td><code>icegrad0</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{Ice A}_{1}\f$</td>
    <td><code>icegrad1</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{Ice A}_{2}\f$</td>
    <td><code>icegrad2</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{Ice A}_{3}\f$</td>
    <td><code>icegrad3</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{Ice A}_{4}\f$</td>
    <td><code>icegrad4</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{Ice Phs}_{1}\f$</td>
    <td><code>icegrad5</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{Ice Phs}_{2}\f$</td>
    <td><code>icegrad6</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{Ice Phs}_{3}\f$</td>
    <td><code>icegrad7</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>\f$\textrm{Ice Phs}_{4}\f$</td>
    <td><code>icegrad8</code></td>
    <td>\f$0.0 \pm 1.0\f$</td>
    <td>gradient</td>
  </tr>
</table> -->


