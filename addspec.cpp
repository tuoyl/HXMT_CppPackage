#include <iostream>
#include <string.h>
#include "fitsio.h"
#include <fstream>
#include <vector>


using namespace std;
void PrintError(int status);
void ReadSpecFile(char* filename, double* &counts, double exposure, int nRows, char telescop, char instrume, int detchans);
int InitDataSize(char* filename);

void help()
{
        cout<<endl<<"USAGE"<<endl
        <<"addspec: just use for adding the PHA spectrum of same detID!"<<endl
        <<"addspec inputPHA.txt inputBKG.txt inputRSP.txt spec.fits bkg.fits totalRSP.rsp"<<endl;
}



int main(int argc, char* argv[])
{
    /* Pass Spec file list */
    std::ifstream infilelist;
    infilelist.open(argv[1]);
    char filename[200];
    int ifile = 0;
    double* total_counts;
    double total_exposure = 0;
    int counts_size;

    while(infilelist.getline(filename, sizeof(filename)))
    {
        std::cout << filename << endl;
        double* counts;
        double exposure;
        int nRows;
        char telescop;
        char instrume;
        int detchans;

        /*initial data size*/
        if (total_counts == NULL)
        {
            counts_size = InitDataSize(filename);
            total_counts = new double[counts_size];
        }

        cout << "size of counts" << sizeof(counts) << endl;
        ReadSpecFile(filename, counts, exposure, nRows, telescop, instrume, detchans);
        total_exposure = total_exposure + exposure;
        cout << "size of counts" << sizeof(counts) << endl;
        for (int i=0; i <= nRows; i++)
        {
//            total_counts[i] = total_counts[i] + counts[i];
        }

    }

    return 0;
}

int InitDataSize(char* filename)
{
    int sizenum; 

    fitsfile* fptr;
    char card[FLEN_CARD];
    int status = 0, datatype, anynull;
    if(fits_open_file(&fptr, filename, READONLY, &status)) PrintError(status);
    if(ffmahd(fptr, 2, &datatype, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "NAXIS2", &sizenum, NULL, &status)) PrintError(status);
    if(fits_close_file(fptr, &status)) PrintError(status);
    return sizenum;
}

void ReadSpecFile(char* filename, double* &counts, double exposure, int nRows, char telescop, char instrume, int detchans)
{
    fitsfile* fptr;
    char card[FLEN_CARD];
    int status = 0, datatype, anynull;
    double doublenull = 0;
    int colnum_COUNTS;

    if(fits_open_file(&fptr, filename, READONLY, &status)) PrintError(status);
    if(ffmahd(fptr, 2, &datatype, &status)) PrintError(status);
    if(ffgky(fptr, TDOUBLE, "EXPOSURE", &exposure, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "NAXIS2", &nRows, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TSTRING, "TELESCOP", &telescop, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TSTRING, "INSTRUME", &instrume, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "DETCHANS", &detchans, NULL, &status)) PrintError(status);
    if(fits_get_colnum(fptr, CASEINSEN, "COUNTS", &colnum_COUNTS, &status)) PrintError(status);

    counts = new double[nRows];

    if(fits_read_col(fptr, TDOUBLE, colnum_COUNTS, 1, 1, nRows, &doublenull, &counts, &anynull, &status)) PrintError(status);

    if(fits_close_file(fptr, &status)) PrintError(status);
}

void WriteSpecFile()
{
//TODO:
}

void ReadRspFile(char* filename, double* energ_lo, double* energ_hi, int N_GRP, int F_CHAN, int detchans, double* matrix, int matrix_nRows)
{
    cout << "RSP file: " << filename << endl;
    fitsfile* fptr;
    char card[FLEN_CARD];
    int status = 0, datatype, anynull;
    double doublenull = 0;
    int colnum_ENERG_LO, colnum_ENERG_HI, colnum_N_GRP, colnum_F_CHAN, colnum_MATRIX;

    if(fits_open_file(&fptr, filename, READONLY, &status)) PrintError(status);
    if(ffmahd(fptr, 2, &datatype, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "NAXIS2", &matrix_nRows, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "DETCHANS", &detchans, NULL, &status)) PrintError(status);

    if(fits_get_colnum(fptr, CASEINSEN, "ENERG_LO", &colnum_ENERG_LO, &status)) PrintError(status);
    if(fits_get_colnum(fptr, CASEINSEN, "ENERG_HI", &colnum_ENERG_HI, &status)) PrintError(status); 
    if(fits_get_colnum(fptr, CASEINSEN, "N_GRP", &colnum_N_GRP, &status)) PrintError(status);
    if(fits_get_colnum(fptr, CASEINSEN, "F_CHAN", &colnum_F_CHAN, &status)) PrintError(status);
    if(fits_get_colnum(fptr, CASEINSEN, "MATRIX", &colnum_MATRIX, &status)) PrintError(status);

    matrix = new double[detchans*matrix_nRows];
    energ_lo = new double[matrix_nRows];
    energ_hi = new double[matrix_nRows];

    if(fits_read_col(fptr, TDOUBLE, colnum_ENERG_LO, 1, 1, matrix_nRows, &doublenull, energ_lo, &anynull, &status)) PrintError(status);
    if(fits_read_col(fptr, TDOUBLE, colnum_ENERG_HI, 1, 1, matrix_nRows, &doublenull, energ_hi, &anynull, &status)) PrintError(status);
    if(fits_read_col(fptr, TDOUBLE, colnum_MATRIX, 1, 1, matrix_nRows*detchans, &doublenull, matrix, &anynull, &status)) PrintError(status);

    if(fits_close_file(fptr, &status)) PrintError(status);
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





