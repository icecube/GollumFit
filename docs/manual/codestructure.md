# Code Structure

## Overview

The GollumFit program is organized in a straightforward way. There is a central `GollumFit` class which handles the main functionality of the code, 
accepting different necessary inputs depending on the desired functionality. 

There are a handful of helper structures like `FitParameters`, `FitResult` or `SteeringParams` which organize the various variables
that are passed around. 

Re-weighting of events is handled by the `WeighterMaker` class, which constructs the various individual weighters
that returns a weight for each event. 

