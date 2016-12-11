************************************************************
*******Download data from MEDLINE using python script*******
************************************************************
1. The script "step_1_extract_substance.py" is in "input" folder
2. This script downloads drugs in the file "drugNameList_original.txt" (each drugname must be in a newline)
3. Run the python script using python from command line as follows (requires bio package to be installed beforehand)

python step_1_extract_substance.py drugNameList_original.txt para.txt

4. This will download and give you files like <drugname>_result_1_compounds.txt etc.
5. result_1 are records in Group A (with DDI keyword) and result_2 are records in Group B (without DDI keyword that are not in Group A).
6. Each line in these files represent a record (substance field of a medical paper)
7. Term (compunds or proteins) in a record are separated using ~ character.


************************************************************
*******   Run the analysis algorithm on these files  *******
************************************************************
1. The C++ program is "counting_multi_th.cpp"
2. Compile this file using g++ complier as follows (or any other c++ compiler). If you are using windows modify line 48 to "#define _WIN32"

g++ -std=c++0x -fopenmp -o counting counting_multi_th.cpp

3. This shoud give you an executable program "counting"
4. Input files are in the "input" folder.
5. Outputs are written to "output" folder.
6. To run the program, go to command line.
7. Type the following command

ON LINUX:    ./counting input/ibuprofen_result_1_compounds.txt input/ibuprofen_result_2_compounds.txt input/ibuprofen_result_1_proteins.txt input/ibuprofen_result_2_proteins.txt input/para.txt ibuprofen 8
ON WINDOWS:  counting input/ibuprofen_result_1_compounds.txt input/ibuprofen_result_2_compounds.txt input/ibuprofen_result_1_proteins.txt input/ibuprofen_result_2_proteins.txt input/para.txt ibuprofen 8

where the last number, 8 is the number of threads.

8. Try varying this number and run multiple times.
9. Samping frequency can be changed by changing "sampleTimes" in "para.txt" inside "input" folder