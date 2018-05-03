/* 4-digit airfoil generator for arbitrary lift coefficient specification */

#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <map>
#include <utility>
#include <cmath>
#include <stdio.h>
#include <fstream>

#include <stdlib.h>

void linspace(double* vec, double lowVal, double highVal, uint32_t elements);
void cosspace(double* vec, double lowVal, double highVal, const uint32_t elements);

int main(int argc, char** argv) {

    /* constants*/
    const uint32_t numPoints = 100;
    const uint32_t numCambers = 9;
    const double a[5] = {0.2969, -0.126, -0.3516, 0.2843, -0.1015}; // These values are for a open trailing edge; a[4] = -0.1036 for closed TE


    /* variables */
    double camber;
    double posMaxCamber;
    double thickness;
    double x[numPoints];
    double yc[numPoints];
    double yt[numPoints];
    double dycdx[numPoints];
    double theta[numPoints];
    double xLower[numPoints];
    double xUpper[numPoints];
    double yLower[numPoints];
    double yUpper[numPoints];

    /* Parse command line arguments */
    if(argc != 4) {

        std::cout << "ERROR: four-digit requires exactly three arguments of the form 'four-digit 0.### 0.### ##.###" << std::endl;
        std::cout << "       The first argument is camber as a fraction of the chord" << std::endl;
        std::cout << "       The second argument is max camber as a fraction of the chord" << std::endl;
        std::cout << "       The third argument is the thickness of the airfoil, as a fraction of the chord" << std::endl;
        return -1;

    } else {

        try {
            std::string camberString = argv[1]; // camber
            std::string maxCamberString = argv[2]; // max camber position
            std::string thicknessString = argv[3]; // thickness
            std::string::size_type sz;

            camber = std::stod(camberString, &sz);
            posMaxCamber = std::stod(maxCamberString, &sz);
            thickness = std::stod(thicknessString, &sz);

        } catch(...) {
            std::cout << "ERROR: Input arguments formatted incorrectly (stod() failed)" << std::endl;
            return -2;
        }

    }


    /* Do the work */

    // setup the x vector
    double store = thickness;
    cosspace(x, 0.0, 1.0, numPoints);
    thickness = store; // is this a bug in clang?

    // compute the values
    for(uint32_t i = 0; i < numPoints; i++) {
        if(x[i] < posMaxCamber) {
            yc[i] = (camber / pow(posMaxCamber,2.0)) * (2 * posMaxCamber * x[i] - pow(x[i], 2.0));
            dycdx[i] = (2 * camber / pow(posMaxCamber,2.0)) * (posMaxCamber - x[i]);
        } else {
            yc[i] = (camber / pow(1.0 - posMaxCamber, 2.0)) * ((1 - 2 * posMaxCamber) + 2 * posMaxCamber * x[i] - pow(x[i], 2.0));
            dycdx[i] = (2 * camber / pow(1.0 - posMaxCamber, 2.0)) * (posMaxCamber - x[i]);
        }

        yt[i] = (thickness / 0.2) * (a[0]*pow(x[i], 0.5) + a[1]*x[i] + a[2]*pow(x[i],2.0) + a[3]*pow(x[i], 3.0) + a[4]*pow(x[i], 4.0));

        theta[i] = atan(dycdx[i]);

        // assign coordinates
        xUpper[i] = x[i] - yt[i]*sin(theta[i]);
        xLower[i] = x[i] + yt[i]*sin(theta[i]);
        yUpper[i] = yc[i] + yt[i]*cos(theta[i]);
        yLower[i] = yc[i] - yt[i]*cos(theta[i]);

    }


    /* Output to file */
    char filename[numPoints];
    FILE* fp;
    sprintf(filename, "./NACA(%2.2f)(%2.2f)(%2.2f).dat", camber*100.0, posMaxCamber*10.0, thickness*100);

    fp = fopen(filename, "w");

    fprintf(fp, "NACA(%2.2f)(%1.2f)(%2.2f) Airfoil M=%1.2f%% P=%1.2f%% T=%2.2f%%\n", camber*100.0, posMaxCamber*10.0, thickness*100, camber, posMaxCamber*100.0, thickness*100);


    for(int i = numPoints - 1; i >= 0; i--) {
        fprintf(fp, "  %1.6f  %1.6f\n", xUpper[i], yUpper[i]);
    }

    // start this loop at 1 because the (0, 0) point has already been placed in the file
    for(uint32_t i = 1; i < numPoints; i++) {
        fprintf(fp, "  %1.6f %1.6f\n", xLower[i], yLower[i]);
    }

    fclose(fp);

    /* Exit program */

    std::cout << "Output file saved as " << filename << std::endl;
    return 0;

}


// Simple function for creating linearly spaced vectors, inspired by MATLAB
void linspace(double* vec, double lowVal, double highVal, const uint32_t elements) {

    double dv = (highVal - lowVal) / static_cast<double>(elements - 1);

    for(uint32_t i = 0; i < elements; i++) {
        vec[i] = lowVal + static_cast<double>(i) * dv;
    }

}

// simple function for creating cosine-spaced vectors
void cosspace(double* vec, double lowVal, double highVal, const uint32_t elements) {

    double xLin[elements];
    linspace(xLin, lowVal, highVal, elements);

    for(uint32_t i = 0; i <= elements; i++) {
        double frac = (xLin[i] - lowVal) / (highVal - lowVal);
        vec[i] = lowVal + (highVal - lowVal) * (1 - cos(frac * M_PI / 2.0));
    }

}
