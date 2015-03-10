# TIMBIR

Directory structure
> Source Code                     #Contains all the source code
  >> MBIR_4D                      #Contains the source code for 4D model-based iterative reconstruction (4D MBIR)
  >> include                      #Contains header files with function prototypes for the reconstruction routines
  >> lib                          #Contains compiled libraries of the reconstruction routines
  >> reconstruct                  #Code to read the data, reconstruct, and write the data.
    >>> basic                     #Software front-end to reconstruct data in basic binary format
    >>> real_data                 #Software front-end to reconstruct data in HDF format
        >>>> APS_Data_Format      #Reconstruct data in HDF format primarily used at APS (in Argonne national lab)
        >>>> Standard_Data_Format #Reconstruct data in a standard HDF format 
    >>> sim_data                  #Reconstruct simulated data (NOT FUNCTIONAL YET)

Compiling the software
