ThermoModel
===========

Creating input files
--------------------
We are including the code we used to generate our own inputs, but the folder also contains all of the input files generated so you can write your own code to generate them as needed.
We are including 62 input files (38 that we fit on and 59 of which we generated predictions for).  Supplementary table T2 specifies which input file corresponds to each construct.  The code below also generates 62 files, which can be changed on line 911 of ``MASTcreateInputs.cpp``.

In order to run this code, you will need to have the mast executable in the create_inputs directory.
To do this, download the MEME suite, and follow the instructions to compile mast.


In the folder create_inputs:

``make MAST``

``./create_inputs [flags]``

``./create_inputs -h`` will give a list of flags to set thresholds, PWMs, interactions, type of cooperativity, type of quenching

For the modeling used in the runs from which data is presented, the following settings were used:
	``Dthresh 0.001``

	``Tthresh 0.001``

	``Sthresh 0.001``

	``dlPWM 1``

	``twPWM 3``

	``snPWM 3``
  
	``inter 1 (nearest-neighbor interactions)``

In order to run the model, move the input files into the folder containing the rest of the code for the model.

Fitting parameters
------------------
For models with fewer parameters, the can be done locally in a reasonable amount time. However, we recommend the use of a dedicated server or compute cluster for larger models.

- Files to be changed depending on runs:

  ``initials.par``
	line 12: currently reads "N 23"
	update this line to reflect the total number of parameters in the model you are running (scaling factors + cooperativity + quenching)
	
  ``ThermoExample.c``
	line 18: typeC = type of cooperativity (same numbers used for create_inputs flags), numC = number of cooperativity parameters
	
  line 19: typeQ = type of quenching (same numbers used for create_inputs flags), numQ = number of quenching parameters
	
- Optional changes:

  ``Thermomodel.cpp``
	line 970 - 971: the 38 refers to the number of constructs the model is fitting on (line 970 also contains a list of which constructs they are).
	If you want to change how many constructs the model is fitting on or which ones, adjust these two lines accordingly

  Make sure your makefile uses the correct compiler for the system you will be compiling and running this code on (currently set to use the GNU g++ compiler; we used the Intel C++ compiler for our runs).
	
- Running the code:
  
  ``make read``
  
  ``./read``
  
  ``make``

  ``./evo``
	
Predicting expression patterns
------------------------------
NOTE: this code currently generates predicted expression patterns for 59 constructs. To change this, in file ``getResults_and_Predict.cpp``, change the loop on line 114 to go to the desired number of constructs (where it currently says ``i < 59``). Similarly, in the loop on line 272, make the same change where it says ``con < 59``, changing to the number of constructs you wish to make predictions for.

Within the directory containing your code:

``make predict``

``./predict``

This generates a series of files (one for each construct) named Prediction_[construct#].txt which have a value for the predicted expression at each of the 17 points along the DV axis of the embryo. We used MATLAB to visualize these values compared to known expression.  
