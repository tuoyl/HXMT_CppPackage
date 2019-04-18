#include <iostream>
#include <string.h>
#include "fitsio.h"
#include <fstream>
#include <vector>
#include <stdlib.h>


using namespace std;
void PrintError(int status);
int InitDataSize(char* filename, int ext_num);
void ReadSpecFile(char* filename, double* &out_counts, double &out_exposure, int &out_nRows, int &out_detchans);
void ReadRspFile(char* filename, double* &energ_lo, double* &energ_hi, int &N_GRP, int &F_CHAN, int &detchans, double* &matrix, int &matrix_nRows, double* &e_min, double* &e_max);
void WriteSpecFile(char* filename, double exposure, double* counts, int detchans, const char* telescop,const char* instrume, int quality, int grouping);
void WriteRspFile(char* filename, double* matrix, int matrix_nRows, double* energ_lo, double* energ_hi, int N_GRP, int F_CHAN, double* E_MIN, double* E_MAX, int detchans, const char* telescop, const char* instrume);

void help()
{
        cout<<endl<<"USAGE"<<endl
        <<"addspec: just use for adding the PHA spectrum of same detID!"<<endl
        <<"addspec inputPHA.txt inputBKG.txt inputRSP.txt spec.fits bkg.fits totalRSP.rsp"<<endl;
}

string out_telescop;
string out_instrume;
std::vector<double> spec_exposure;
std::vector<double> bkg_exposure;
double total_spec_exposure=0;
double* out_matrix;

int main(int argc, char* argv[])
{
    for (int argvnum=1; argvnum <= 2; argvnum++)
    {
        /* Pass Spec file list */
        std::ifstream infilelist;
        infilelist.open(argv[argvnum]);
        char filename[200];
        int ifile = 0;
        double* total_counts;
        double total_exposure = 0;
        int counts_size;
        int detchans=0;

        while(infilelist.getline(filename, sizeof(filename)))
        {
            std::cout << filename << endl;
            double* counts;
            double exposure;
            int nRows;

            /*initial data size*/
            if (total_counts == NULL)
            {
                counts_size = InitDataSize(filename, 2);
                total_counts = new double[counts_size];
                for(int i=0;i<counts_size;i++) total_counts[i]=0.0;
            }
            /*initial end*/

            ReadSpecFile(filename, counts, exposure, nRows,detchans);

            if (argvnum == 1)
            {
                spec_exposure.push_back(exposure);
                total_exposure = total_exposure + exposure;
                total_spec_exposure = total_spec_exposure + exposure;
                for (int i=0; i <= nRows; i++)
                {
                    total_counts[i] = total_counts[i] + counts[i];
                }
            }
            else if (argvnum == 2)
            {
                bkg_exposure.push_back(exposure);
                total_exposure = total_exposure + exposure;
                total_spec_exposure = total_spec_exposure + exposure;
                for (int i=0; i <= nRows; i++)
                {
                    total_counts[i] = total_counts[i] + counts[i]*spec_exposure[ifile]/bkg_exposure[ifile];
                }
            }

            ifile++;
        }
        infilelist.close();
        if (strcmp(argv[argvnum+3], "None") && strcmp(argv[argvnum+3], "none") != 0)
        {
            WriteSpecFile(argv[argvnum+3], total_exposure, total_counts, detchans, out_telescop.c_str(), out_instrume.c_str(), 0, 1);
        }
    }

    /* read and write rsp */
    std::ifstream infilelist;
    infilelist.open(argv[3]);
    char filename[200];
    int ifile = 0;
    
    double* energ_lo;
    double* energ_hi;
    double* e_min;
    double* e_max;
    int N_GRP=0;
    int F_CHAN=0;
    double* out_matrix;
    std::vector<double*> matrix_tmp;
    int detchans=0;
    int matrix_nRows;

    while(infilelist.getline(filename, sizeof(filename)))
    {
        double* matrix;

        /*initial data size (expect matrix)*/
        if (energ_lo == NULL)
        {
            int energ_size=0;
            energ_size = InitDataSize(filename, 2);
            energ_lo = new double[energ_size];
            energ_hi = new double[energ_size];
        }
        if (e_min == NULL)
        {
            int e_min_size=0;
            e_min_size=InitDataSize(filename, 3);
            e_min = new double[e_min_size];
            e_max = new double[e_min_size];
        }

        ReadRspFile(filename, energ_lo, energ_hi, N_GRP, F_CHAN, detchans, matrix, matrix_nRows, e_min, e_max);
        matrix_tmp.push_back(matrix);

        ifile++;
    }
    
    /*initial matrix size*/
    out_matrix = new double[detchans*matrix_nRows];
    for (int i=0;i<detchans*matrix_nRows;i++) out_matrix[i] = 0;

    /*merge matrix together*/
    for (int i=0; i<matrix_tmp.size(); i++)
    {
        for (int j=0; j<detchans*matrix_nRows; j++)
        {
            out_matrix[j] = out_matrix[j] + matrix_tmp[i][j]*spec_exposure[i]/total_spec_exposure;
        }
    }

    /* write RSP file */
    //TODO
    WriteRspFile(argv[6], out_matrix, matrix_nRows, energ_lo, energ_hi, N_GRP, F_CHAN, e_min, e_max, detchans, out_telescop.c_str(), out_instrume.c_str());



    return 0;
}

int InitDataSize(char* filename, int ext_num)
{
    int sizenum; 

    fitsfile* fptr;
    int status = 0, datatype, anynull;
    if(fits_open_file(&fptr, filename, READONLY, &status)) PrintError(status);
    if(ffmahd(fptr, ext_num, &datatype, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "NAXIS2", &sizenum, NULL, &status)) PrintError(status);
    if(fits_close_file(fptr, &status)) PrintError(status);
    return sizenum;
}

void ReadSpecFile(char* filename, double* &out_counts, double &out_exposure, int &out_nRows, int &out_detchans)
{
    fitsfile* fptr;
    int status = 0, datatype, anynull;
    double doublenull = 0;
    int colnum_COUNTS;

    double counts;
    double exposure;
    int nRows=0;
    char telescop[20];
    char instrume[20];
    int detchans;

    if(fits_open_file(&fptr, filename, READONLY, &status)) PrintError(status);
    if(ffmahd(fptr, 2, &datatype, &status)) PrintError(status);
    if(ffgky(fptr, TDOUBLE, "EXPOSURE", &exposure, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "NAXIS2", &nRows, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TSTRING, "TELESCOP", &telescop, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TSTRING, "INSTRUME", &instrume, NULL, &status)) PrintError(status);
    if(ffgky(fptr, TINT, "DETCHANS", &detchans, NULL, &status)) PrintError(status);
    if(fits_get_colnum(fptr, CASEINSEN, "COUNTS", &colnum_COUNTS, &status)) PrintError(status);


    out_counts = new double[nRows];
    int channel;
    for (int currentRow = 1; currentRow <= nRows; currentRow++)
    {
        if(fits_read_col_dbl(fptr, colnum_COUNTS, currentRow, 1, 1, 0, &counts, &anynull, &status)) PrintError(status);
        out_counts[currentRow-1] = counts;
    }

    out_nRows = nRows;
    out_exposure = exposure;
    out_telescop = telescop;
    out_instrume = instrume;
    out_detchans = detchans;

    if(fits_close_file(fptr, &status)) PrintError(status);
}

void WriteSpecFile(char* filename, double exposure, double* counts, int detchans, const char* telescop,const char* instrume, int quality, int grouping)
{
    string tmp = filename;
    string outfilename = "!" + tmp;

    fitsfile *fptr;
    static char* ttype[] = {"CHANNEL", "COUNTS", "QUALITY", "GROUPING"};
    static char* tform[] = {"1J", "1J", "1I", "1I"};
    static char* tunit[] = {"chan", "counts", " ", " "};
    int status = 0;
    if(fits_create_file(&fptr, outfilename.c_str(), &status)) PrintError(status);
    if(fits_create_tbl(fptr, BINARY_TBL, 0, 4, ttype, tform, tunit, "SPECTRUM", &status)); PrintError(status);

    if(fits_write_key_lng(fptr, "DETCHANS", detchans, "Total no. detector channels available", &status)) PrintError(status);
    if(fits_write_key_str(fptr, "BACKFILE", "NONE","background FITS file for this object", &status)) PrintError(status);
    if(fits_write_key_flt(fptr, "BACKSCAL", 1.0, 10, "background scaling factor", &status)) PrintError(status);
    if(fits_write_key_str(fptr, "CORRFILE", "NONE","correlation FITS file for this object", &status)) PrintError(status);
    if(fits_write_key_flt(fptr, "CORRSCAL", 1.0, 10, "correlation scaling factor", &status)) PrintError(status);
    if(fits_write_key_str(fptr, "RESPFILE", "NONE"," ", &status)) PrintError(status);
    if(fits_write_key_str(fptr, "ANCRFILE", "NONE"," ", &status)) PrintError(status);
    if(fits_write_key_str(fptr, "FILTER",   "NONE", "redistribution matrix file(RMF)", &status)) PrintError(status);
    if(fits_write_key_str(fptr, "PHAVERSN", "1992a", " ", &status)) PrintError(status);
    if(fits_write_key_log(fptr, "STATERR" , 0, "no statisical error specified", &status)) PrintError(status);
    if(fits_write_key_log(fptr, "SYSERR"  , 0, "no systematic error", &status)) PrintError(status);
    if(fits_write_key_log(fptr, "POISSERR", 1, "Poissonian statistical errors to be assumed", &status)) PrintError(status);//XXX:check the written content is correct
    if(fits_write_key_lng(fptr, "GROUPING", 1, "grouping of the data has been defined", &status)) PrintError(status);
    if(fits_write_key_lng(fptr, "QUALITY",  1, "data quality information specified", &status)) PrintError(status);
    if(fits_write_key_lng(fptr, "AREASCAL", 1,"area scaling factor",&status)) PrintError(status);
    if(fits_write_key_dbl(fptr, "EXPOSURE", exposure, 10, "exposure time",&status)) PrintError(status);
    if(fits_write_key_dbl(fptr, "LIVETIME", 1, 5,   "Total spectrum accumulation time", &status)) PrintError(status);
    if(fits_write_key_dbl(fptr, "DEADC",    0, 5, "Deadtime correction factor", &status)) PrintError(status);
    if(fits_write_key_dbl(fptr, "DETID",    0, 5, " ", &status)) PrintError(status);
    if(fits_write_key_str(fptr, "CHANTYPE", "PI", " ", &status)) PrintError(status);
    if(fits_write_key_lng(fptr, "TLMIN2", 0, " ", &status)) PrintError(status);
    if(fits_write_key_lng(fptr, "TLMAX2", detchans-1, " ", &status)) PrintError(status);
    if(fits_write_key_str(fptr, "TELESCOP", telescop, "Telescope (mission) name", &status)) PrintError(status);
    if(fits_write_key_str(fptr, "INSTRUME", instrume, "Instrument name", &status)) PrintError(status);

    /*write colum*/
    for (int i = 0; i<detchans; i++)
    {
        if(fits_write_col_int(fptr, 1, i+1, 1, 1, &i, &status)) PrintError(status);
        if(fits_write_col_dbl(fptr, 2, i+1, 1, 1, &counts[i], &status)) PrintError(status);
        if(fits_write_col_int(fptr, 3, i+1, 1, 1, &quality, &status)) PrintError(status);
        if(fits_write_col_int(fptr, 4, i+1, 1, 1, &grouping, &status)) PrintError(status);
    }
    if(fits_close_file(fptr, &status)) PrintError(status);

}

void ReadRspFile(char* filename, double* &energ_lo, double* &energ_hi, int &N_GRP, int &F_CHAN, int &detchans, double* &matrix, int &matrix_nRows, double* &e_min, double* &e_max)
{
    cout << "RSP file: " << filename << endl;
    fitsfile* fptr;
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

    /* Read third Extension */
    int colnum_e_min, colnum_e_max;
    if(ffmahd(fptr, 3, &datatype, &status)) PrintError(status);
    if(fits_get_colnum(fptr, CASEINSEN, "E_MIN", &colnum_e_min, &status)) PrintError(status);
    if(fits_get_colnum(fptr, CASEINSEN, "E_MAX", &colnum_e_max, &status)) PrintError(status); 

    if(fits_read_col(fptr, TDOUBLE, colnum_e_min, 1, 1, detchans, &doublenull, e_min, &anynull, &status)) PrintError(status);
    if(fits_read_col(fptr, TDOUBLE, colnum_e_max, 1, 1, detchans, &doublenull, e_max, &anynull, &status)) PrintError(status);

    if(fits_close_file(fptr, &status)) PrintError(status);
}

void WriteRspFile(char* filename, double* matrix, int matrix_nRows, double* energ_lo, double* energ_hi, int N_GRP, int F_CHAN, double* E_MIN, double* E_MAX, int detchans, const char* telescop, const char* instrume)
{
    string tmp = filename;
    string outfilename = "!" + tmp;

    fitsfile* fptr;
    int status = 0;
    char matrix_form[10];
    sprintf(matrix_form,"%dE", detchans);

    static char* matrix_colname[] = {"ENERG_LO", "ENERG_HI", "N_GRP", "F_CHAN", "N_CHAN", "MATRIX"};
    static char* matrix_colform[] = {"1E", "1E", "1I","1I","1I",matrix_form};
    static char* matrix_colunit[] = {"keV", "keV", " ", " ", " ", " "};

    if(fits_create_file(&fptr, outfilename.c_str(), &status)) PrintError(status);
    if(fits_create_tbl(fptr, BINARY_TBL, 0, 6, matrix_colname, matrix_colform, matrix_colunit, "MATRIX", &status)); PrintError(status);

    fits_write_key_str(fptr, "TUNIT1", "keV", "", &status);
    fits_write_key_str(fptr, "TUNIT2", "keV", "", &status);
    fits_write_key_str(fptr, "EXTNAME", "MATRIX", "", &status);
    fits_write_key_str(fptr, "TELESCOP", telescop, "Telescope (mission) name", &status);
    fits_write_key_str(fptr, "INSTRUME", instrume, "Instrument name", &status);
    fits_write_key_str(fptr, "DETNAM", instrume, "", &status);
    fits_write_key_lng(fptr, "DETCHANS", detchans, "", &status);
    fits_write_key_str(fptr, "CHANTYPE", "PI", "", &status);
    fits_write_key_str(fptr, "HDUVERS", "1.3.0", "", &status);
    fits_write_key_str(fptr, "CCLS0001", "BCF", "", &status);
    fits_write_key_str(fptr, "HDUCLASS", "OGIP", "", &status);
    fits_write_key_str(fptr, "HDUCLAS1", "RESPONSE", "", &status);
    fits_write_key_str(fptr, "HDUCLAS2", "RSP_MATRIX", "", &status);
    fits_write_key_str(fptr, "CDTP0001", "DATA", "", &status);
    fits_write_key_str(fptr, "TLMIN4", "0", "", &status);

    
    fits_write_col(fptr, TINT, 1, 1, 1, matrix_nRows, energ_lo, &status);
    fits_write_col(fptr, TINT, 2, 1, 1, matrix_nRows, energ_hi, &status);
    for (int i=0; i<matrix_nRows; i++)
    {
        fits_write_col_int(fptr, 3, i+1, 1, 1, &N_GRP, &status);
        fits_write_col_int(fptr, 4, i+1, 1, 1, &F_CHAN, &status);
        fits_write_col_int(fptr, 5, i+1, 1, 1, &detchans, &status);
    }
    fits_write_col(fptr, TDOUBLE, 6, 1, 1, matrix_nRows*detchans, matrix, &status);

    static char* ebounds_colname[] = {"CHANNEL", "E_MIN", "E_MAX"};
    static char* ebounds_colform[] = {"1I",      "1E",    "1E"};
    static char* ebounds_colunit[] = {" ",       "keV",   "keV"};
    fits_create_tbl(fptr, BINARY_TBL, 0, 3, ebounds_colname, ebounds_colform, ebounds_colunit, "EBOUNDS", &status);

    fits_write_key_str(fptr, "TUNIT2", "keV", "", &status);
    fits_write_key_str(fptr, "TUNIT3", "keV", "", &status);
    fits_write_key_str(fptr, "EXTNAME", "EBOUNDS", "", &status);
    fits_write_key_str(fptr, "TELESCOP", telescop, "Telescope (mission) name", &status);
    fits_write_key_str(fptr, "INSTRUME", instrume, "Instrument name", &status);
    fits_write_key_str(fptr, "DETNAM", instrume, "", &status);
    fits_write_key_lng(fptr, "DETCHANS", detchans, "", &status);
    fits_write_key_str(fptr, "CHANTYPE", "PI", "", &status);
    fits_write_key_str(fptr, "HDUVERS", "1.3.0", "", &status);
    fits_write_key_str(fptr, "CCLS0001", "BCF", "", &status);
    fits_write_key_str(fptr, "HDUCLASS", "OGIP", "", &status);
    fits_write_key_str(fptr, "HDUCLAS1", "RESPONSE", "", &status);
    fits_write_key_str(fptr, "HDUCLAS2", "RSP_MATRIX", "", &status);
    fits_write_key_str(fptr, "CDTP0001", "DATA", "", &status);
    fits_write_key_str(fptr, "TLMIN4", "0", "", &status);

    for (int i=0;i<detchans;i++)
    {
        fits_write_col_int(fptr, 1, i+1, 1, 1, &i, &status);
        fits_write_col_dbl(fptr, 2, i+1, 1, 1, &E_MIN[i], &status);
        fits_write_col_dbl(fptr, 3, i+1, 1, 1, &E_MAX[i], &status);
    }

    fits_close_file(fptr, &status);
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





