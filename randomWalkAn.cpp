#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <string>
#include <vector>
using namespace std;


int main(int argc, char const* argv[])
{
    int numWalkers = 1000;
    int numSteps = 500;

    //Read binary values from binary file and store all values in an array except first column.
    //TASK: read txt file instead of the binary file
    //http://www.cplusplus.com/doc/tutorial/files/
    //http://courses.cs.vt.edu/cs2604/fall01/binio.html
    ifstream positions("positions2.bin", ios::in | ios::binary);


    //Offstream files
    ofstream posMeanCSV;
    ofstream posMean;

    posMeanCSV.open("quantitiesVStime2.csv");
    posMean.open("quantitiesVStime2.txt");

    //Variables for Squared Distance
    double q1 = 1 / sqrt(numSteps), q2 = 1, q3 = sqrt(numSteps);

    cout << "q1 " << q1 << endl;
    cout << "q2 " << q2 << endl;
    cout << "q3 " << q3 << endl;

    double aStepSquared;
    double aStepSquared2;
    double cosRq1;
    double cosRq2;
    double cosRq3;
    double cosSqrdRq1;
    double cosSqrdRq2;
    double cosSqrdRq3;

    //=====================================================//

    for (int step = 0; step < numSteps; step++) {

        aStepSquared = 0.0;
        aStepSquared2 = 0.0;
        cosRq1 = 0.0;
        cosRq2 = 0.0;
        cosRq3 = 0.0;
        cosSqrdRq1 = 0.0;
        cosSqrdRq2 = 0.0;
        cosSqrdRq3 = 0.0;

        double val;
        positions.read((char*)&val, sizeof(double));


        for (int walker = 0; walker < numWalkers; walker++) {

            positions.read((char*)&val, sizeof(double));   //extracted from positions then store it to val
                                                        //&val : Initial byte of an object stored in file. http://www.tutorialdost.com/Cpp-Programming-Tutorial/62-Cpp-File-Handling-IO-Read-Write-Object.aspx
                                                        //sizeof(int) : size of object represents the total number of bytes to be read from initial byte.
            double currentPosition = val;           //use the currentposition immediately
            //cout << "Step Walker currentPosition val " << step << " " << walker << " " << currentPosition << " " << val << endl;

            //COMPUTE aStepSquared and cosRq1, cosRq2, cosRq3
            aStepSquared += (currentPosition * currentPosition);
            aStepSquared2 += pow(currentPosition, 4);

            // cos (q*R)
            cosRq1 += cos(q1 * currentPosition);
            cosRq2 += cos(q2 * currentPosition);
            cosRq3 += cos(q3 * currentPosition);

            //cos^2 (q*R)
            cosSqrdRq1 += cos(q1 * currentPosition) * cos(q1 * currentPosition);
            cosSqrdRq2 += cos(q2 * currentPosition) * cos(q2 * currentPosition);
            cosSqrdRq3 += cos(q3 * currentPosition) * cos(q3 * currentPosition);

        }
        //Mean Square
        double meanSquare = (double)aStepSquared / (double)numWalkers;
        double meanSquare2 = (double)aStepSquared2 / (double)numWalkers;
        //<cos(q*R)>
        cosRq1 = (double)cosRq1 / (double)numWalkers;
        cosRq2 = (double)cosRq2 / (double)numWalkers;
        cosRq3 = (double)cosRq3 / (double)numWalkers;
        //<cos^2(q*R)>
        cosSqrdRq1 = (double)cosSqrdRq1 / (double)numWalkers;
        cosSqrdRq2 = (double)cosSqrdRq2 / (double)numWalkers;
        cosSqrdRq3 = (double)cosSqrdRq3 / (double)numWalkers;
        //standard deviation of the mean
        double stdevMeanSquare = (double)sqrt(meanSquare2 - meanSquare * meanSquare);
        double stdErrorOfMean = (double)(stdevMeanSquare / sqrt(numWalkers));
        double stdevCosSqrdRq1 = (double)sqrt(cosSqrdRq1 - (cosRq1 * cosRq1));
        double stdErrorCosSqrdRq1 = (double)(stdevCosSqrdRq1 / sqrt(numWalkers));
        double stdevCosSqrdRq2 = (double)sqrt(cosSqrdRq2 - (cosRq2 * cosRq2));
        double stdErrorCosSqrdRq2 = (double)(stdevCosSqrdRq2 / sqrt(numWalkers));
        double stdevCosSqrdRq3 = (double)sqrt(cosSqrdRq3 - (cosRq3 * cosRq3));
        double stdErrorCosSqrdRq3 = (double)(stdevCosSqrdRq3 / sqrt(numWalkers));

        double exp1 = exp(-((step + 1.0) / 2.0) * q1 * q1);
        double exp2 = exp(-((step + 1.0) / 2.0) * q2 * q2);
        double exp3 = exp(-((step + 1.0) / 2.0) * q3 * q3);

        // .txt file
        posMean << step << " " << meanSquare << " " << step + 1 << " " << stdErrorOfMean << " " << cosRq1 << " " << cosRq2 << " " << cosRq3 << " " << stdErrorCosSqrdRq1 << " " << stdErrorCosSqrdRq2 << " " << stdErrorCosSqrdRq3 << " " << setprecision(9) << exp1 << " " << exp2 << " " << exp3 << endl;
        //posMean << step << " " << meanSquare << " " << step + 1 << " " << stdErrorOfMean << " " << cosRq1 << " " << cosRq2 << " " << cosRq3 << " " << stdErrorCosSqrdRq1 << " " << stdErrorCosSqrdRq2 << " " << stdErrorCosSqrdRq3 << endl;

        posMeanCSV << step << "," << meanSquare << "," << step + 1 << "," << stdErrorOfMean << "," << cosRq1 << "," << cosRq2 << "," << cosRq3 << "," << stdErrorCosSqrdRq1 << "," << stdErrorCosSqrdRq2 << "," << stdErrorCosSqrdRq3 << ", " << setprecision(9) << exp1 << ", " << exp2 << ", " << exp3 << endl;

    }

    cout << "q1 " << q1 << endl;
    cout << "q2 " << q2 << endl;
    cout << "q3 " << q3 << endl;


    posMeanCSV.close();
    posMean.close();

    return 0;
}