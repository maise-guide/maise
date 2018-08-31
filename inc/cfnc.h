#ifndef CFNC
#define CFNC

extern const double Pi;
extern const int D3;

int    FIND_NDIM(Cell *C);
void   KILL_DBL (Cell *C, double tol);
void   APPL_SG  (Cell *C, double tol);
int    FIND_WYC (Cell *C, Cell *D, double tol, char *ISO, int J);
void   READ_CIF (Cell *C, char file[], double tol, int NM);
int    FIND_MTY (Cell *C, double tol);
void   FIND_PRS (Cell *C, Cell *D, double tol);
void   FIND_CXC (Cell *C, Cell *D, int argc, char **argv);
void   COMP_STR (Cell *C, Cell *D, int argc, char **argv);
void   INIT_CELL(Cell *C, char filename[], int N, int NM, int J);
void   CELL_EXAM(Cell *C, Cell *D, int argc, char **argv);

#endif
