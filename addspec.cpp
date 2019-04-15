#include <iostream>
#include <string.h>
#include "fitsio.h"
#include <fstream>
#include <vector>


using namespace std;
void PrintError(int status);
void ReadSpecFile(char* filename, double counts, double exposure, int nRows, char telescop, char instrume, int detchans);

void help()
{
        cout<<endl<<"USAGE"<<endl
        <<"addspec: just use for adding the PHA spectrum of same detID!"<<endl
        <<"addspec inputPHA.txt inputBKG.txt inputRSP.txt spec.fits bkg.fits totalRSP.rsp"<<endl;
}


int main(int argc, char* argv[])
{
    std::ifstream infilelist;
    infilelist.open(argv[1]);
    char filename[200];

    int ifile = 0;
    while(infilelist.getline(filename, sizeof(filename)))
    {
        std::cout << filename << endl;
        double counts;
        double exposure;
        int nRows;
        char telescop;
        char instrume;
        int detchans;
        ReadSpecFile(filename, counts, exposure, nRows, telescop, instrume, detchans);
    }

    return 0;
}

void ReadSpecFile(char* filename, double counts, double exposure, int nRows, char telescop, char instrume, int detchans)
{
    fitsfile* fptr;
    char card[FLEN_CARD];
    int status = 0, datatype, anynull;
    int colnum_COUNTS

    if(fits_open_file(&fptr, filename, READONLY, &status)) PrintError(status);
    if(ffmahd(fptr, 2, &datatype, &status)) PrintError(status);
    if(ffgky(fptr, TDOUBLE, "EXPOSURE", &exposure, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "NAXIS2", &nRows, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TSTRING, "TELESCOP", &telescop, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TSTRING, "INSTRUME", &instrume, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "DETCHANS", &detchans, NULL, &status)) PrintError(status);
    if(fits_get_colnum(fptr, CASEINSEN, "COUNTS", &colnum_COUNTS, &status)) PrintError(status);

    for (int currentRow = 1; currentRow <= nRows; currentRow++)
    {
        if(fits_read_col_dbl(fptr, colnum_COUNTS, currentRow, 1, 1, 0.0, &counts, &anynull, &status)) PrintError(status);

    }

    if(fits_close_file(fptr, &status)) PrintError(status);
}

void WriteSpecFile()
{
//TODO:
}

void ReadRspFile(char* filename, double* energ_lo, double* energ_hi, int N_GRP, int F_CHAN, int detchans, double* matrix, int matrix_nRows)
{
    fitsfile* fptr;
    char card[FLEN_CARD];
    int status = 0, datatype, anynull;

    if(fits_open_file(&fptr, filename, READONLY, &status)) PrintError(status);
    if(ffmahd(fptr, 2, &datatype, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "NAXIS2", &matrix_nRows, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "DETCHANS", &detchans, NULL, &status)) PrintError(status);
    for (int currentRow = 1; currentRow <= matrix_nRows; currentRow++)
    {
    
    }

}

void WriteRspFile()
{
//TODO
}

void PrintError(int status)
{
    if (status)
    {
        fits_report_error(stderr, status);
        exit(status);
    }
    return;
}





