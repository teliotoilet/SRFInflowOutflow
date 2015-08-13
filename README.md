# SRFInflowOutflow
inletOutlet-based boundary condition for SRF farfield boundaries

This folder should be placed in 
  src/finiteVolume/cfdTools/general/SRF/derivedFvPatchFields

Add the following line to src/finiteVolume/Make/files in the SRF section:
  $(SRF)/derivedFvPatchFields/SRFInflowOutflow/SRFInflowOutflowFvPatchVectorField.C
(tested for v2.3.1)
