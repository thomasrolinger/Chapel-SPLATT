import sys
import os

def getSubDirs(myDir):
    return [name for name in os.listdir(myDir) if os.path.isdir(os.path.join(myDir, name))]

if __name__ == "__main__":
    # First arg is top-level directory
    topDir = sys.argv[1]
    if topDir[-1] == '/':
        topDir = topDir[:-1]
    # Second arg is output file
    outFile = sys.argv[2]

    # Third arg is tool name
    toolName = sys.argv[3]

    dataSet = topDir.split("/")[-1].split("/")[0].strip()

    # Create the output file.
    fOut = open(outFile, "w")
    
    # Write out header. 
    fOut.write("TOOL,DATASET,THREADS,CPD,MTTKRP,INVERSE,MAT MULT,MAT A^TA,MAT NORM,CPD FIT,IO,SORT,TOTAL\n")
    # Get the subdirs in the top level dir. These should be the 1THD, 2THD, etc dirs
    subDirs = getSubDirs(topDir)

    # Sort the subdirs
    subDirs.sort(key=lambda item: (int(item.partition('THD')[0]) if item[0].isdigit() else float('inf'), item))
    subDirs = [os.path.join(topDir, subDir) for subDir in subDirs]

    # For each subdir...
    for subDir in subDirs:
        # Get the number of threads for this dir
        numThreads = subDir.split("THD")[0].split("/")[-1].strip()
        # This subdir will have some number of files, each corresponding to a trial
        trialFiles = ["%s/%s" %(subDir,f) for f in os.listdir(subDir) if os.path.isfile(os.path.join(subDir, f))] 
        
        # If no files, skip
        if len(trialFiles) != 0:
            # Get num files and num threads
            numFiles = len(trialFiles)
            CPD = 0.0
            MTTKRP = 0.0
            INVERSE = 0.0
            MAT_MULT = 0.0
            MAT_ATA = 0.0
            MAT_NORM = 0.0
            CPD_FIT = 0.0
            IO = 0.0
            SORT = 0.0
            TOTAL= 0.0
            for trialFile in trialFiles:
                try:
                    fTrial = open(trialFile, "r")
                    lines = fTrial.readlines()
                    fTrial.close()
                    # Get TOTAL time
                    TOTAL += [float(line.split("TOTAL")[1].strip().split("s")[0].strip()) for line in lines if "TOTAL" in line][0]
                    # Get CPD time
                    CPD += [float(line.split("CPD")[1].strip().split("s")[0].strip()) for line in lines if "CPD" in line][0]            
                    # Get ATA time
                    MAT_ATA += [float(line.split("MAT A^TA")[1].strip().split("s")[0].strip()) for line in lines if "MAT A^TA" in line][0]
                    # Get MAT NORM time
                    MAT_NORM += [float(line.split("MAT NORM")[1].strip().split("s")[0].strip()) for line in lines if "MAT NORM" in line][0]
                    # Get MTTKRP time
                    MTTKRP += [float(line.split("MTTKRP")[1].strip().split("s")[0].strip()) for line in lines if "MTTKRP" in line][0]
                    # Get INVERSE time
                    INVERSE += [float(line.split("INVERSE")[1].strip().split("s")[0].strip()) for line in lines if "INVERSE" in line][0]
                    # Get MAT_MULT time
                    MAT_MULT += [float(line.split("MAT MULT")[1].strip().split("s")[0].strip()) for line in lines if "MAT MULT" in line][0]
                    # Get CPD_FIT time
                    CPD_FIT += [float(line.split("FIT")[1].strip().split("s")[0].strip()) for line in lines if "FIT" in line][0]
                    # Get IO time
                    IO += [float(line.split("IO")[1].strip().split("s")[0].strip()) for line in lines if "IO " in line][0]
                    # Get SORT time
                    SORT += [float(line.split("SORT")[1].strip().split("s")[0].strip()) for line in lines if "SORT" in line][0]
                except Exception, e:
                    print("****Failed to process file: %s" %trialFile)
                    print("\t%s" %str(e))
                    pass
                # Now that we processed all the files, take averages and write them out
            fOut.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(toolName, dataSet, str(numThreads), str(CPD/numFiles), str(MTTKRP/numFiles), \
                                                                         str(INVERSE/numFiles), str(MAT_MULT/numFiles), str(MAT_ATA/numFiles), str(MAT_NORM/numFiles), \
                                                                         str(CPD_FIT/numFiles), str(IO/numFiles), str(SORT/numFiles), str(TOTAL/numFiles) ) )
    
    fOut.write("\n")
    fOut.close()
