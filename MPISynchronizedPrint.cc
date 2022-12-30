            // write a MPI output file in .txt format
            //
            void
                writeMPIOutput(const std::vector<std::vector<double> > & data,
                        const std::string fileName,
                        const int numberLocalLines,
                        const int numCols)

                {

                    Utils::MPIController::mpi_communicator_type mpiCommunicator;
#if defined(HAVE_MPI)
                    //
                    // get MPI controller
                    //
                    const Utils::MPIController & mpiController = 
                        Utils::MPIControllerSingleton::getInstance();

                    //
                    // get MPI communicator
                    // 
                    mpiCommunicator = mpiController.getCommunicator();
                    const int numberProcessors = mpiController.getNumberProcesses();
                    const int taskId = mpiController.getId();
                    const int rootTaskId = mpiController.getRootId();

#else  
                    mpiCommunicator = MPI_COMM_SELF;     
                    const int numberProcessors = 1;
                    const int taskId = 0;
                    const int rootTaskId = 0;

#endif // HAVE_MPI

                    int ierr, rank, size;
                    MPI_Offset offset;
                    MPI_File   file;
                    MPI_Status status;
                    MPI_Datatype num_as_string;
                    MPI_Datatype localarray;
                    const int charspernum=15;

                    char *const fmt="%8.5f\t";
                    char *const endfmt="%8.5f\n";	

                    std::vector<int> numberLinesInProcs(numberProcessors);
                    MPI_Allgather(&numberLocalLines, 1, MPI_INT, &numberLinesInProcs[0], 1, MPI_FLOAT,
                            MPI_COMM_WORLD);

                    //
                    // get the global number lines
                    //
                    int numberGlobalLines = 0;
                    for(unsigned int iProc = 0; iProc < numberProcessors; ++iProc)
                        numberGlobalLines += numberLinesInProcs[iProc];

                    //
                    // get the start row for current proc
                    //
                    int startRow = 0;
                    for(unsigned int iProc = 0; iProc < taskId; ++iProc)
                        startRow += numberLinesInProcs[iProc];

                    int endRow = startRow + numberLinesInProcs[taskId] - 1;


                    /* each number is represented by charspernum chars */
                    MPI_Type_contiguous(charspernum, MPI_CHAR, &num_as_string); 
                    MPI_Type_commit(&num_as_string);

                    /* convert our data into txt */
                    char * data_as_txt = (char*) malloc(numberLocalLines*numCols*charspernum*sizeof(char));

                    int count = 0;
                    for (int iLine = 0; iLine < numberLocalLines; iLine++) 
                    {
                        for (int iCol = 0; iCol < numCols; iCol++) 
                        {

                            std::ostringstream dataOString;
                            dataOString << data[iLine][iCol];
                            std::string dataString = dataOString.str();
                            int stringLength = dataString.length();
                            for(unsigned int iChar = stringLength; iChar < charspernum-1; ++iChar)
                                dataString += " ";

                            if(iCol == numCols - 1)
                                dataString += "\n";
                            else
                                dataString += "\t";

                            const char * dataStringPtr = dataString.c_str();
                            int iStart = count*charspernum;
                            int iStop = iStart + charspernum;
                            for(unsigned int iChar = 0; iChar < charspernum; ++iChar)
                                data_as_txt[iChar + iStart] = dataStringPtr[iChar];
                            //sprintf(&data_as_txt[count*charspernum], fmt, data[iLine][iCol]);
                            count++;
                        }

                        // last col is added separetely to have '\n' formatting at the end
                        //sprintf(&data_as_txt[count*charspernum], endfmt, data[iLine][numCols-1]);
                        //count++;
                    }
                    //
                    //	std::ostringstream dataOString;
                    //
                    //	for(unsigned int iLine = 0; iLine < numberLocalLines; ++iLine)
                    //	{
                    //
                    //	  for(unsigned int iCol = 0; iCol < numCols; ++iCol)
                    //	  {
                    //
                    //	    dataOString << data[iLine][iCol] << "\t";
                    //	  
                    //	  }
                    //
                    //	  dataOString << "\n";
                    //
                    //	}
                    //
                    //	std::string dataString = dataOString.str();
                    //	const char * data_as_txt = dataString.c_str();

                    /* create a type describing our piece of the array */
                    int globalsizes[2] = {numberGlobalLines, numCols};
                    int localsizes [2] = {numberLocalLines, numCols};
                    int starts[2]      = {startRow, 0};
                    int order          = MPI_ORDER_C;

                    MPI_Type_create_subarray(2, globalsizes, localsizes, starts, order, num_as_string, &localarray);
                    MPI_Type_commit(&localarray);

                    const char * outFile = fileName.c_str();

                    /* open the file, and set the view */
                    MPI_File_open(MPI_COMM_WORLD, outFile, 
                            MPI_MODE_CREATE|MPI_MODE_WRONLY,
                            MPI_INFO_NULL, &file);

                    MPI_File_set_view(file, 0,  MPI_CHAR, localarray, 
                            "native", MPI_INFO_NULL);

                    MPI_File_write_all(file, data_as_txt, numberLocalLines*numCols, num_as_string, &status);
                    MPI_File_close(&file);

                    MPI_Type_free(&localarray);
                    MPI_Type_free(&num_as_string);

                }
