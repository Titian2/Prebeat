# Prebeat

Current Version 2.0

This function takes ringdown data from mechanical loss systems and allows
for logical inspection of the data. A simple 'pass' / 'reject' can be
invoked allowing data to be kept for further analysis or flat out
rejected. If there are artifcacts which would cause influence an
incorrect fit such as long portions of data where the ringdown has
finished, these can be cropped out. Only the cropped portions of these
regions are saved for further analysis

On the event of a crash or error the outputs of the program will be saved
to the working directory. Upon re-running, the function will look
in the current directory and subfolders for the saved files and
automatically resume at the crashed point.

          prebeat(verbose,'samplename')

### INPUTS

  **verbose**   : Turn on prints to command line 1 = on , 0 =off

 **samplename**                   : Name of Sample / Set of data - used for save files and
                                   plots. should be entered as a string
## OUTPUTS

**PyTotalAnalysis.mat**           :matlab results file which contains the
                                   final analysed data

**Processed_Mode_X.txt**          :file  written to a .txt file in the
                                  selected directory containing a breakdown
                                  of measurements for each mechanical  mode
                                  format = 'freq' 'offset' 'amplitude' 'fitted tau' ' Mechanical loss' 'Origin filename'

**Log_samplename_.txt**     :Log file containing the results of each
                                  ringdown after it has been presented to the user.
