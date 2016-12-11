// Sample execution: ./a.out input/ibuprofen_result_1_compounds.txt input/ibuprofen_result_2_compounds.txt input/ibuprofen_result_1_proteins.txt input/ibuprofen_result_2_proteins.txt input/para.txt ibuprofen 8

//======================================= file = counting_multi_th.cpp ======
//=  Program to implement multi threaded version of the DDI algorithm       =
//===========================================================================
//=  Notes:                                                                 =
//=    1) This file requires compiler and library support for the ISO C++   =
//=       2011 standard. This support is currently experimental, and must   =
//=       be enabled with the -std=c++11 or -std=gnu++11 or -std=c++0x      =
//=       compiler options. Requires g++ version 4.7 or higher              =
//=    2) Put the initial definition as _LINUX or _WIN32 depending on       =
//=       wheather your OS is windows or linux (mac should also work)       =
//=    3) The program takes parameters from the command line in the         =
//=       following order                                                   =
//=                  argv[1]: <DrugName>_result_1_compounds.txt             =
//=                  argv[2]: <DrugName>_result_2_compounds.txt             =
//=                  argv[3]: <DrugName>_result_1_proteins.txt              =
//=                  argv[4]: <DrugName>_result_2_proteins.txt              =
//=                  argv[5]: para.txt                                      =
//=                  argv[6]: <DrugName>                                    =
//=                  argv[7]: # of threads                                  =
//=    4) Program writes to a folder named output. This folder must be      =
//=       created before hand and the user executing this program must      =
//=       have the permission to write to this folder. Following files      =
//=       are writen to this folder                                         =
//=                  output/<DrugName>_compounds_output_p_value.txt         =
//=                  output/<DrugName>_proteins_output_p_value.txt          =
//=                  output/<DrugName>_compounds_output.txt                 =
//=                  output/<DrugName>_proteins_output.txt                  =
//=    5) Uses OpenMP parallel programming model that is readily available  =
//=       in GCC compliers. Number of threads is defined in NUM_THREADS     =
//=       For GCC compliers -fopenmp flag must be used while compiling      =
//=-------------------------------------------------------------------------=
//=  Build: GCC   - g++ -std=c++0x -fopenmp counting_multi_th.cpp           =
//=         INTEL - icpc -std=c++0x -qopenmp counting_multi_th.cpp          =
//=-------------------------------------------------------------------------=
//=  Execute: ./a.out input/ibuprofen_result_1_compounds.txt                =
//=                   input/ibuprofen_result_2_compounds.txt                =
//=                   input/ibuprofen_result_1_proteins.txt                 =
//=                   input/ibuprofen_result_2_proteins.txt                 =
//=                   input/para.txt   ibuprofen                            =
//=-------------------------------------------------------------------------=
//=  Author: Sulav Malla                                                    =
//=          University of South Florida                                    =
//=          WWW: http://eng.usf.edu/~sulavmalla/                           =
//=          Email: sulavmalla@mail.usf.edu                                 =
//===========================================================================
#define _LINUX           // _LINUX or _WIN32

#include <unordered_map>
#include <vector>
#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <omp.h>

#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

using namespace std;

// Remove if already defined
typedef long long int64; typedef unsigned long long uint64;
int  NUM_THREADS;    // Number of threads to create

void randomize(vector<std::size_t>& indexes, int num_random);
void calculate_p_z_value(double x, vector<int>& sampleCount,double& p_value,double& z_value);
uint64 GetTimeMs64(void);

int main(int argc, char** argv) {
    // variables to measure time taken for calculations
    // time is measured in milliseconds
    uint64 startTime, endTime, calculationTime, samplingBTime, totalTime;
    calculationTime = 0;
    samplingBTime = 0;
    totalTime = 0;

    // check for correct number of parameter from command line
    if(argc != 8)
    {
        cout << "Invalid number of parameters (7 required)" << endl;
        exit(1);
    }

    startTime = GetTimeMs64();  // measure time for this section

    string   drugName(argv[6]);
    NUM_THREADS = atoi(argv[7]);

    // read in parameters from the input parameter file
    string   parameter, value;
    int      N_SAMPLE, CUTOFF;
    double   P_VALUE, Z_SCORE;
    string   END_DATE;

    ifstream inputParameterFile(argv[5]);
    while (inputParameterFile >> parameter >> value)
    {
        if(parameter == "sampleTimes")
            N_SAMPLE = stoi(value);
        else if(parameter == "cutoff")
            CUTOFF = stoi(value);
        else if(parameter == "p_value")
            P_VALUE = stod(value);
        else if(parameter == "z_score")
            Z_SCORE = stod(value);
        else if(parameter == "endDate")
            END_DATE = value;
    }

    inputParameterFile.close();



    // hash table to store the count of occurance of substance in group A
    unordered_map<string, int> compoundsHashTableA;
    unordered_map<string, int> proteinsHashTableA;    
    // also create a hash table to store the count of terms (that
    // are also in group A) in each sampling
    unordered_map<string, vector<int>> compoundsSampleCount;
    unordered_map<string, vector<int>> proteinsSampleCount;
    // vector of string to store each record of group B
    // each record in these files is a medical paper (with unique PMID)
    vector<string>             compoundsRecordsB;
    vector<string>             proteinsRecordsB;
    // vector of strings to store all unique substance in Group A
    vector<string>             compoundsInA;
    vector<string>             proteinsInA;
    // vector of double to store p-value and z-value for substances
    vector<double>             compoundsPvalue, proteinsPvalue;
    vector<double>             compoundsZvalue, proteinsZvalue;

    // input stream of the compound and protein file from the two groups
    // create a dictionary of term and count from group A
    // create a string array to store each line from group B
    // we randomly sample from this string array multiple times and
    // create the dictionary of term and count for each sample.
    // term can either be compound or protein
    ifstream compoundsFileA(argv[1]);
    ifstream compoundsFileB(argv[2]);
    ifstream proteinsFileA(argv[3]);
    ifstream proteinsFileB(argv[4]);
    string   line, term;
    int      numberOfCompoundsRecordsA, numberOfCompoundsRecordsB;
    int      numberOfProteinsRecordsA, numberOfProteinsRecordsB;
    int      sampleSize;
    stringstream line_stream;  // to convert string to stream


//=-------------------------------------------------------------------------=
//=  Processing of Group A (records with "drug interaction" keyword)        =
//=-------------------------------------------------------------------------=

    // create the dictionary of term and count from group A
    // Read the file line by line till the end
    numberOfCompoundsRecordsA = 0;
    while (getline(compoundsFileA, line))
    {
        // Make a stream out of the line
        line_stream.clear();    // Need to clear the status
        line_stream.str(line);
        // Read the line, term by term
        // terms are separated using ~
        while (getline(line_stream, term, '~'))
        {
            // Check to see if the term is alread in the
            // hash table. If present, increment the count by 1
            // else, initialize the count to 1 also
            // initialize sampling count of that term to zeros
            if (compoundsHashTableA.count(term) > 0)
                compoundsHashTableA[term]++;
            else
            {
                compoundsHashTableA[term] = 1;
                compoundsInA.push_back(term);
                for(int i = 0; i < N_SAMPLE; ++i)
                    compoundsSampleCount[term].push_back(0);
            }
        }
        numberOfCompoundsRecordsA++;
    }
    // repeat the same for proteins
    numberOfProteinsRecordsA = 0;
    while (getline(proteinsFileA, line))
    {
        // Make a stream out of the line
        line_stream.clear();    // Need to clear the status
        line_stream.str(line);
        while (getline(line_stream, term, '~'))
        {
            if (proteinsHashTableA.count(term) > 0)
                proteinsHashTableA[term]++;
            else
            {
                proteinsHashTableA[term] = 1;
                proteinsInA.push_back(term);
                for(int i = 0; i < N_SAMPLE; ++i)
                    proteinsSampleCount[term].push_back(0);
            }
        }
        numberOfProteinsRecordsA++;
    }

//=-------------------------------------------------------------------------=
//=  Processing of Group B (records without keyword that are not in A)      =
//=-------------------------------------------------------------------------=

    // create a string array that contains each line from group B
    // Read the file line by line till the end
    // also get the count of lines (records) in group B
    while (getline(compoundsFileB, line))
    {
        compoundsRecordsB.push_back(line);
    }
    while (getline(proteinsFileB, line))
    {
        proteinsRecordsB.push_back(line);
    }
    numberOfCompoundsRecordsB = compoundsRecordsB.size();
    numberOfProteinsRecordsB = proteinsRecordsB.size();


    // We are done with the files, so close them
    compoundsFileA.close();
    compoundsFileB.close();
    proteinsFileA.close();
    proteinsFileB.close();

    endTime = GetTimeMs64();
    totalTime += (endTime - startTime);


//=-------------------------------------------------------------------------=
//=  Sampling of Group B (this is done in parallel)                         =
//=-------------------------------------------------------------------------=
    sampleSize = numberOfCompoundsRecordsA;
    omp_set_num_threads(NUM_THREADS); // set the number of threads
    startTime = GetTimeMs64();  // measure time for this calculation

//**********************START OF PARALLEL SECTION****************************

    #pragma omp parallel private(term)
    {
        stringstream line_stream;  // to convert string to stream

        // vector to store randomized permutation
        vector<std::size_t> compoundsRandomIndex(compoundsRecordsB.size());
        vector<std::size_t> proteinsRandomIndex(proteinsRecordsB.size());
        std::iota(compoundsRandomIndex.begin(), compoundsRandomIndex.end(), 0);
        std::iota(proteinsRandomIndex.begin(), proteinsRandomIndex.end(), 0);
        // hash table to store the count of occurance of substance in
        // sampled records in group B
        unordered_map<string, int> compoundsHashTableB;
        unordered_map<string, int> proteinsHashTableB;
        // initialize random seed
        srand(time(NULL));
        // We do sampling for N_SAMPLE number of times
        // each time, we select SampleSize number of records
        // from entire records in Group B and count the occurance
        // of terms. If these terms are also present in Group A,
        // we record the term count in SampleCount table
        // parallelize this for loop by breaking apart iterations into chunks
        #pragma omp for
        for(int j=0; j<N_SAMPLE; ++j)
        {
            // randomize the records such that the first
            // SampleSize elements represents unique
            // sampling from the entire records
            randomize(compoundsRandomIndex, sampleSize);
            randomize(proteinsRandomIndex, sampleSize);
            compoundsHashTableB.clear(); // initialize dictionary to be empty
            proteinsHashTableB.clear(); // initialize dictionary to be empty

            // create the dictionary of term and count from sampled records
            // Read line by line
            for(int i=0; i<sampleSize; ++i)
            {
                // Make a stream out of each record
                line_stream.clear();    // Need to clear the status
                line_stream.str(compoundsRecordsB[compoundsRandomIndex[i]]);
                // Read the line, term by term
                // terms are separated using ~
                while (getline(line_stream, term, '~'))
                {
                    // Check to see if the term is alread in the
                    // hash table. If it is not present, initialize
                    // the count to 1 and record the new term found
                    // else just increment the count by 1
                    if (compoundsHashTableB.count(term) > 0)
                        compoundsHashTableB[term]++;
                    else
                        compoundsHashTableB[term] = 1;
                }
            }
            // Repeat the same for proteins
            for(int i=0; i<sampleSize; ++i)
            {
                line_stream.clear();    // Need to clear the status
                line_stream.str(proteinsRecordsB[proteinsRandomIndex[i]]);
                while (getline(line_stream, term, '~'))
                {
                    if (proteinsHashTableB.count(term) > 0)
                        proteinsHashTableB[term]++;
                    else
                        proteinsHashTableB[term] = 1;
                }
            }

            // if the terms in Group A are also present in the sample, we record
            // the corresponding count from sample. We can safely write to the
            // shared dictionary compoundsSampleCount because we are writing to
            // different memory location within it. This is guaranteed by a
            // different value of variable j for different threads
            for (vector<string>::iterator it = compoundsInA.begin() ; it != compoundsInA.end(); ++it)
            {
                term = *it;
                if (compoundsHashTableB.count(term) > 0)
                    compoundsSampleCount[term][j] = compoundsHashTableB[term];
            }
            for (vector<string>::iterator it = proteinsInA.begin() ; it != proteinsInA.end(); ++it)
            {
                term = *it;
                if (proteinsHashTableB.count(term) > 0)
                    proteinsSampleCount[term][j] = proteinsHashTableB[term];
            }
        }
    }
//**********************END OF PARALLEL SECTION******************************

    endTime = GetTimeMs64();
    samplingBTime = (endTime - startTime);
    calculationTime += samplingBTime;
    totalTime += samplingBTime;
    cout << "Time for sampling B: " << samplingBTime << " ms" << endl;


//=-------------------------------------------------------------------------=
//=  p-value and z-value calculations                                       =
//=-------------------------------------------------------------------------=

    startTime = GetTimeMs64();  // measure time for this calculation
    // iterate through the substances calculating p-value and z-value for each
    double p_value, z_value;
    for (vector<string>::iterator it = compoundsInA.begin() ; it != compoundsInA.end(); ++it)
    {
        term = *it;
        calculate_p_z_value((double)compoundsHashTableA[term],compoundsSampleCount[term], p_value, z_value);
        compoundsPvalue.push_back(p_value);
        compoundsZvalue.push_back(z_value);
    }
    for (vector<string>::iterator it = proteinsInA.begin() ; it != proteinsInA.end(); ++it)
    {
        term = *it;
        calculate_p_z_value((double)proteinsHashTableA[term],proteinsSampleCount[term], p_value, z_value);
        proteinsPvalue.push_back(p_value);
        proteinsZvalue.push_back(z_value);
    }

    endTime = GetTimeMs64();
    calculationTime += (endTime - startTime);
    totalTime += (endTime - startTime);
    cout << "Total time for all calculation: " << calculationTime << " ms" << endl;


//=-------------------------------------------------------------------------=
//=  Find the sort permutation according to z-value (descending order)      =
//=  this is later used to output the result in sorted order                =
//=-------------------------------------------------------------------------=

    startTime = GetTimeMs64();  // measure time for this section

    // vector to store sort permutation (index order)
    vector<std::size_t> compoundsSortedIndex(compoundsZvalue.size());
    vector<std::size_t> proteinsSortedIndex(proteinsZvalue.size());
    std::iota(compoundsSortedIndex.begin(), compoundsSortedIndex.end(), 0);
    std::iota(proteinsSortedIndex.begin(), proteinsSortedIndex.end(), 0);
    std::sort(compoundsSortedIndex.begin(), compoundsSortedIndex.end(),
        [&](std::size_t i, std::size_t j){ return (compoundsZvalue[i] > compoundsZvalue[j]); });
    std::sort(proteinsSortedIndex.begin(), proteinsSortedIndex.end(),
        [&](std::size_t i, std::size_t j){ return (proteinsZvalue[i] > proteinsZvalue[j]); });


//=-------------------------------------------------------------------------=
//=  output the sorted order to a file                                      =
//=-------------------------------------------------------------------------=
    // create output file and write header to it
    ofstream compoundsOutputFile("output/"+drugName+"_compounds_output.txt");
    ofstream compoundsOutputFileWithPValue("output/"+drugName+"_compounds_output_p_value.txt");
    ofstream proteinsOutputFile("output/"+drugName+"_proteins_output.txt");
    ofstream proteinsOutputFileWithPValue("output/"+drugName+"_proteins_output_p_value.txt");
    compoundsOutputFile << "Compound,Count in A,Distribution in B,z-score,p-value\n";
    compoundsOutputFileWithPValue << "Compound,Count in A,Distribution in B,z-score,p-value\n";
    proteinsOutputFile << "Protein,Count in A,Distribution in B,z-score,p-value\n";
    proteinsOutputFileWithPValue << "Protein,Count in A,Distribution in B,z-score,p-value\n";
    // iterate through the sorted index and print the 
    // corresponding elements to the file
    std::size_t index;
    for (vector<std::size_t>::iterator it = compoundsSortedIndex.begin() ; it != compoundsSortedIndex.end(); ++it)
    {
        index = *it;
        term = compoundsInA[index];

        if(compoundsHashTableA[term] >= CUTOFF)
        {
            compoundsOutputFile << drugName << ";" << term << "," << compoundsHashTableA[term] << ",[";
            for (vector<int>::iterator it2 = compoundsSampleCount[term].begin() ; it2 != compoundsSampleCount[term].end(); ++it2)
                compoundsOutputFile << *it2 << "~";
            compoundsOutputFile << "]," << compoundsZvalue[index] << "," << compoundsPvalue[index] << "\n";
        }
        if((compoundsHashTableA[term] >= CUTOFF) && (compoundsPvalue[index] <= P_VALUE))
        {
            compoundsOutputFileWithPValue << drugName << ";" << term << "," << compoundsHashTableA[term] << ",[";
            for (vector<int>::iterator it2 = compoundsSampleCount[term].begin() ; it2 != compoundsSampleCount[term].end(); ++it2)
                compoundsOutputFileWithPValue << *it2 << "~";
            compoundsOutputFileWithPValue << "]," << compoundsZvalue[index] << "," << compoundsPvalue[index] << "\n";
        }
    }
    // Repeat the same for proteins
    for (vector<std::size_t>::iterator it = proteinsSortedIndex.begin() ; it != proteinsSortedIndex.end(); ++it)
    {
        index = *it;
        term = proteinsInA[index];

        if(proteinsHashTableA[term] >= CUTOFF)
        {
            proteinsOutputFile << drugName << ";" << term << "," << proteinsHashTableA[term] << ",[";
            for (vector<int>::iterator it2 = proteinsSampleCount[term].begin() ; it2 != proteinsSampleCount[term].end(); ++it2)
                proteinsOutputFile << *it2 << "~";
            proteinsOutputFile << "]," << proteinsZvalue[index] << "," << proteinsPvalue[index] << "\n";
        }
        if((proteinsHashTableA[term] >= CUTOFF) && (proteinsPvalue[index] <= P_VALUE))
        {
            proteinsOutputFileWithPValue << drugName << ";" << term << "," << proteinsHashTableA[term] << ",[";
            for (vector<int>::iterator it2 = proteinsSampleCount[term].begin() ; it2 != proteinsSampleCount[term].end(); ++it2)
                proteinsOutputFileWithPValue << *it2 << "~";
            proteinsOutputFileWithPValue << "]," << proteinsZvalue[index] << "," << proteinsPvalue[index] << "\n";
        }
    }

    // close the output file
    compoundsOutputFile.close();
    compoundsOutputFileWithPValue.close();
    proteinsOutputFile.close();
    proteinsOutputFileWithPValue.close();

    endTime = GetTimeMs64();
    totalTime += (endTime - startTime);
    cout << "Grand total time: " << totalTime << " ms" << endl;
}


//=-------------------------------------------------------------------------=
//=  function to randomize the first "num_random" elements of the           =
//=  given vector using Fisher-Yates shuffle this is equivalent to          =
//=  sampling m unique records from n                                       =
//=-------------------------------------------------------------------------=
void randomize(vector<std::size_t>& indexes, int num_random)
{
    int left = indexes.size();
    if(num_random <= left)
    {
        for(int i = 0; i < num_random; ++i) {
            swap(indexes[i], indexes[(rand()%left) + i]);
            --left;
        }
    }
}


//=-------------------------------------------------------------------------=
//=  function to calculate the p-value and z-value. parameters are          =
//=  modified in place                                                      =
//=-------------------------------------------------------------------------=
void calculate_p_z_value(double x, vector<int>& sampleCount,double& p_value,double& z_value)
{
    // calculate mean and standard deviation of the distribution
    double sum = std::accumulate(std::begin(sampleCount), std::end(sampleCount), 0.0);
    double mean =  sum / sampleCount.size();
    double accum = 0.0;
    std::for_each (std::begin(sampleCount), std::end(sampleCount), [&](const double d) {
        accum += (d - mean) * (d - mean);
    });
    double stDeviation = sqrt(accum / (sampleCount.size()-1));

    // calculate z-value
    if(stDeviation == 0)
    {
        if(x == mean)
            z_value = -100.0;
        else
            z_value = x * 100.0;
    }
    else
        z_value = (x - mean) / stDeviation;
    
    // calculate p-value as 1 minus CDF of normal distribution
    // some constants for CDF calculation
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign
    int sign = 1;
    if (z_value < 0)
        sign = -1;
    double temp_z_value = fabs(z_value)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*temp_z_value);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-temp_z_value*temp_z_value);

    p_value = 1 - (0.5*(1.0 + sign*y));
}


//=-------------------------------------------------------------------------=
//=  function to return the amount of milliseconds elapsed since the        = 
//=  UNIX epoch. Works on both windows and linux                            =
//=-------------------------------------------------------------------------=
uint64 GetTimeMs64()
{
    #ifdef _WIN32
    /* Windows */
    FILETIME ft;
    LARGE_INTEGER li;
    /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
    * to a LARGE_INTEGER structure. */
    GetSystemTimeAsFileTime(&ft);
    li.LowPart = ft.dwLowDateTime;
    li.HighPart = ft.dwHighDateTime;
    uint64 ret = li.QuadPart;
    ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
    ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */
    return ret;
    #else
    /* Linux */
    struct timeval tv;
    gettimeofday(&tv, NULL);
    uint64 ret = tv.tv_usec;
    /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
    ret /= 1000;
    /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
    ret += (tv.tv_sec * 1000);
    return ret;
    #endif
}