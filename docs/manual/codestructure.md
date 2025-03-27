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

<table border="1" cellspacing="0" cellpadding="4">
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
    <td>$\textrm{convNorm}$</td>
    <td><code>convNorm</code></td>
    <td>$1.0 \pm 0.2$</td>
    <td>scale factor</td>
  </tr>
  <tr>
    <td>$\rho_{\textrm{atm}}$</td>
    <td><code>zenithCorrection</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\sigma_{\textrm{K-Air}}$</td>
    <td><code>kaonLosses</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{K}_{158G}^{+}$</td>
    <td><code>hadronicHEkp</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{K}_{158G}^{-}$</td>
    <td><code>hadronicHEkm</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\pi_{20T}^{+}$</td>
    <td><code>hadronicVHE1pip</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\pi_{20T}^{-}$</td>
    <td><code>hadronicVHE1pim</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{K}_{2P}^{+}$</td>
    <td><code>hadronicVHE3kp</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{K}_{2P}^{-}$</td>
    <td><code>hadronicVHE3km</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\pi_{2P}^{+}$</td>
    <td><code>hadronicVHE3pip</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\pi_{2P}^{-}$</td>
    <td><code>hadronicVHE3pim</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{p}_{2P}$</td>
    <td><code>hadronicVHE3p</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{n}_{2P}$</td>
    <td><code>hadronicVHE3n</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{GSF}_1$</td>
    <td><code>cosmicRay1</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{GSF}_2$</td>
    <td><code>cosmicRay2</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{GSF}_3$</td>
    <td><code>cosmicRay3</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{GSF}_4$</td>
    <td><code>cosmicRay4</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{GSF}_5$</td>
    <td><code>cosmicRay5</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{GSF}_6$</td>
    <td><code>cosmicRay6</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\Phi^{\textrm{HE}}/10^{-18}\,\textrm{GeV}^{-1}\,\textrm{sr}^{-1}\,\textrm{s}^{-1}\,\textrm{cm}^{-2}$</td>
    <td><code>astroNorm</code></td>
    <td>$0.787 \pm 0.36$</td>
    <td>scale factor</td>
  </tr>
  <tr>
    <td>$\Delta\gamma_{1}^{\textrm{HE}}$</td>
    <td><code>astroDeltaGamma</code></td>
    <td>$0.0 \pm 0.36$</td>
    <td>power law formula</td>
  </tr>
  <tr>
    <td>$\Delta\gamma_{2}^{\textrm{HE}}$</td>
    <td><code>astroDeltaGammaSec</code></td>
    <td>$0.0 \pm 0.36$</td>
    <td>power law formula</td>
  </tr>
  <tr>
    <td>$\log_{10}\left(\textrm{E}_{\textrm{break}}^{\textrm{HE}}/\textrm{GeV}\right)$</td>
    <td><code>astroPivot</code></td>
    <td>$5.0 \pm 1.0$</td>
    <td>power law formula</td>
  </tr>
  <tr>
    <td>$\textrm{promptNorm}$</td>
    <td><code>promptNorm</code></td>
    <td>$1.0 \pm 1.0$</td>
    <td>scale factor</td>
  </tr>
  <tr>
    <td>$\nu/ \bar{\nu}$</td>
    <td><code>NeutrinoAntineutrinoRatio</code></td>
    <td>$1.0 \pm 1.0$</td>
    <td>scale factor</td>
  </tr>
  <tr>
    <td>$\nu\textrm{ Att}$</td>
    <td><code>nuxs</code></td>
    <td>$1.0 \pm 0.1$</td>
    <td>spline</td>
  </tr>
  <tr>
    <td>$\bar{\nu}\textrm{ Att}$</td>
    <td><code>nubarxs</code></td>
    <td>$1.0 \pm 0.1$</td>
    <td>spline</td>
  </tr>
  <tr>
    <td colspan="4" align="center"><strong>IceCube-specific Monte Carlo parameters</strong></td>
  </tr>
  <tr>
    <td>$\textrm{DOM}\:{\textrm{eff}}$</td>
    <td><code>domEfficiency</code></td>
    <td>$1.27 \pm 0.123$</td>
    <td>spline</td>
  </tr>
  <tr>
    <td>$\textrm{Hole Ice}$</td>
    <td><code>holeiceForward</code></td>
    <td>$-1.0 \pm 10.0$</td>
    <td>spline</td>
  </tr>
  <tr>
    <td>$\textrm{Ice A}_{0}$</td>
    <td><code>icegrad0</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{Ice A}_{1}$</td>
    <td><code>icegrad1</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{Ice A}_{2}$</td>
    <td><code>icegrad2</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{Ice A}_{3}$</td>
    <td><code>icegrad3</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{Ice A}_{4}$</td>
    <td><code>icegrad4</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{Ice Phs}_{1}$</td>
    <td><code>icegrad5</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{Ice Phs}_{2}$</td>
    <td><code>icegrad6</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{Ice Phs}_{3}$</td>
    <td><code>icegrad7</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
  <tr>
    <td>$\textrm{Ice Phs}_{4}$</td>
    <td><code>icegrad8</code></td>
    <td>$0.0 \pm 1.0$</td>
    <td>gradient</td>
  </tr>
</table>


