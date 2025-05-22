# Code Structure

## Overview

The GollumFit program is organized in a straightforward way. There is a central `GollumFit` class which handles the main functionality of the code, 
accepting different necessary inputs depending on the desired functionality. 

There are a handful of helper structures like `FitParameters`, `FitResult` or `SteeringParams` which organize the various variables
that are passed around. 

Re-weighting of events is handled by the `WeighterMaker` class, which constructs the various individual weighters
that returns a weight for each event. 

## Nuisance parameters

By default, GollumFit includes handling of common nuisance parameters to all neutrino telescopes, as well as IceCube-specific neutrino telescopes. The following table lists the nuisance parameters that are included.

| Parameter                                                                                         | `GollumFit` Variable       | Prior Value ± Width      | Weighting Method    |
|:--------------------------------------------------------------------------------------------------|:---------------------------|:-------------------------|:--------------------|
| **Common to all neutrino telescopes**                                                            |                              |                         |                      |
| \f$\textrm{convNorm}\f$                                                                           | `convNorm`                 | \f$1.0 \pm 0.2\f$        | scale factor        |
| \f$\rho_{\textrm{atm}}\f$                                                                         | `zenithCorrection`         | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\sigma_{\textrm{K-Air}}\f$                                                                      | `kaonLosses`               | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{K}_{158G}^{+}\f$                                                                        | `hadronicHEkp`             | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{K}_{158G}^{-}\f$                                                                        | `hadronicHEkm`             | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\pi_{20T}^{+}\f$                                                                                | `hadronicVHE1pip`          | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\pi_{20T}^{-}\f$                                                                                | `hadronicVHE1pim`          | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{K}_{2P}^{+}\f$                                                                          | `hadronicVHE3kp`           | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{K}_{2P}^{-}\f$                                                                          | `hadronicVHE3km`           | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\pi_{2P}^{+}\f$                                                                                 | `hadronicVHE3pip`          | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\pi_{2P}^{-}\f$                                                                                 | `hadronicVHE3pim`          | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{p}_{2P}\f$                                                                              | `hadronicVHE3p`            | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{n}_{2P}\f$                                                                              | `hadronicVHE3n`            | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{GSF}_1\f$                                                                               | `cosmicRay1`               | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{GSF}_2\f$                                                                               | `cosmicRay2`               | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{GSF}_3\f$                                                                               | `cosmicRay3`               | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{GSF}_4\f$                                                                               | `cosmicRay4`               | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{GSF}_5\f$                                                                               | `cosmicRay5`               | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{GSF}_6\f$                                                                               | `cosmicRay6`               | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\Phi^{\textrm{HE}}/10^{-18}\,\textrm{GeV}^{-1}\,\textrm{sr}^{-1}\,\textrm{s}^{-1}\,\textrm{cm}^{-2}\f$ | `astroNorm`                | \f$0.787 \pm 0.36\f$     | scale factor        |
| \f$\Delta\gamma_{1}^{\textrm{HE}}\f$                                                               | `astroDeltaGamma`          | \f$0.0 \pm 0.36\f$       | power law formula   |
| \f$\Delta\gamma_{2}^{\textrm{HE}}\f$                                                               | `astroDeltaGammaSec`       | \f$0.0 \pm 0.36\f$       | power law formula   |
| \f$\log_{10}\left(\textrm{E}_{\textrm{break}}^{\textrm{HE}}/\textrm{GeV}\right)\f$               | `astroPivot`               | \f$5.0 \pm 1.0\f$        | power law formula   |
| \f$\textrm{promptNorm}\f$                                                                          | `promptNorm`               | \f$1.0 \pm 1.0\f$        | scale factor        |
| \f$\nu/ \bar{\nu}\f$                                                                              | `NeutrinoAntineutrinoRatio`| \f$1.0 \pm 1.0\f$        | scale factor        |
| \f$\nu\textrm{ Att}\f$                                                                             | `nuxs`                     | \f$1.0 \pm 0.1\f$        | spline              |
| \f$\bar{\nu}\textrm{ Att}\f$                                                                       | `nubarxs`                  | \f$1.0 \pm 0.1\f$        | spline              |
| **IceCube-specific Monte Carlo parameters**                                                                                                                                   |
| \f$\textrm{DOM}\:{\textrm{eff}}\f$                                                                 | `domEfficiency`            | \f$1.27 \pm 0.123\f$     | spline              |
| \f$\textrm{Hole Ice}\f$                                                                            | `holeiceForward`           | \f$-1.0 \pm 10.0\f$      | spline              |
| \f$\textrm{Ice A}_{0}\f$                                                                            | `icegrad0`                 | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{Ice A}_{1}\f$                                                                            | `icegrad1`                 | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{Ice A}_{2}\f$                                                                            | `icegrad2`                 | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{Ice A}_{3}\f$                                                                            | `icegrad3`                 | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{Ice A}_{4}\f$                                                                            | `icegrad4`                 | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{Ice Phs}_{1}\f$                                                                          | `icegrad5`                 | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{Ice Phs}_{2}\f$                                                                          | `icegrad6`                 | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{Ice Phs}_{3}\f$                                                                          | `icegrad7`                 | \f$0.0 \pm 1.0\f$        | gradient            |
| \f$\textrm{Ice Phs}_{4}\f$                                                                          | `icegrad8`                 | \f$0.0 \pm 1.0\f$        | gradient            |

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


