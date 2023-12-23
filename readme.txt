This is the readme for the model file associated with the paper:
Hoshino O, Zheng M, Fukuoka Y. (2022). Effect of cortical extracellular GABA on motor response. J Comput Neurosci [PubMed]

This c-program was contributed by O Hoshino. It was originally built using Microsoft Visual C++ and also might work in Microsoft Visual Studio 2012 (create a new project, add the c file to it, and build and run).

More usage instructions:

1. Set input current to the sensory network by giving a proper value to "int_inp0_3". "onset_0" and "period_0" define its onset time and duration.
2. Set times for output data by giving values "OUT" (starting time) and "PERIOD" (recording time period).
3. Run. The default (as provided) setting of the model is for Figure 3A (left).
4. Output data files (vPY*_*2.dat, vPY*_*1.dat, vm*_*2.dat,) provide the rasters of N_S (top panel) and N_M (middle panel) P cells and spinal motoneurons (Mn) (bottom panel). In these data files, value -10 was assigned to no spike emission and should be discarded when plotting (Igor Pro used by us). GABA_V2 gives the basal ambient GABA concentration in N_S: [GABA]_0^S.