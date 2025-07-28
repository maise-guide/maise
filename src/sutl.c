#include "sutl.h"

extern const double Pi;

//========================================================
//Read space group operations
//========================================================
void READ_SG(Cell *C)
{
  int    rotation[192][3][3];
  double translation[192][3];
  int    j,SGN,HN=0;

  SGN=C->SGN;

  if(SGN==   1) HN=   1;
  if(SGN==   2) HN=   2;
  if(SGN==   3) HN=   3;
  if(SGN==   4) HN=   6;
  if(SGN==   5) HN=   9;
  if(SGN==   6) HN=  18;
  if(SGN==   7) HN=  21;
  if(SGN==   8) HN=  30;
  if(SGN==   9) HN=  39;
  if(SGN==  10) HN=  57;
  if(SGN==  11) HN=  60;
  if(SGN==  12) HN=  63;
  if(SGN==  13) HN=  72;
  if(SGN==  14) HN=  81;
  if(SGN==  15) HN=  90;
  if(SGN==  16) HN= 108;
  if(SGN==  17) HN= 109;
  if(SGN==  18) HN= 112;
  if(SGN==  19) HN= 115;
  if(SGN==  20) HN= 116;
  if(SGN==  21) HN= 119;
  if(SGN==  22) HN= 122;
  if(SGN==  23) HN= 123;
  if(SGN==  24) HN= 124;
  if(SGN==  25) HN= 125;
  if(SGN==  26) HN= 128;
  if(SGN==  27) HN= 134;
  if(SGN==  28) HN= 137;
  if(SGN==  29) HN= 143;
  if(SGN==  30) HN= 149;
  if(SGN==  31) HN= 155;
  if(SGN==  32) HN= 161;
  if(SGN==  33) HN= 164;
  if(SGN==  34) HN= 170;
  if(SGN==  35) HN= 173;
  if(SGN==  36) HN= 176;
  if(SGN==  37) HN= 182;
  if(SGN==  38) HN= 185;
  if(SGN==  39) HN= 191;
  if(SGN==  40) HN= 197;
  if(SGN==  41) HN= 203;
  if(SGN==  42) HN= 209;
  if(SGN==  43) HN= 212;
  if(SGN==  44) HN= 215;
  if(SGN==  45) HN= 218;
  if(SGN==  46) HN= 221;
  if(SGN==  47) HN= 227;
  if(SGN==  48) HN= 228;
  if(SGN==  49) HN= 230;
  if(SGN==  50) HN= 233;
  if(SGN==  51) HN= 239;
  if(SGN==  52) HN= 245;
  if(SGN==  53) HN= 251;
  if(SGN==  54) HN= 257;
  if(SGN==  55) HN= 263;
  if(SGN==  56) HN= 266;
  if(SGN==  57) HN= 269;
  if(SGN==  58) HN= 275;
  if(SGN==  59) HN= 278;
  if(SGN==  60) HN= 284;
  if(SGN==  61) HN= 290;
  if(SGN==  62) HN= 292;
  if(SGN==  63) HN= 298;
  if(SGN==  64) HN= 304;
  if(SGN==  65) HN= 310;
  if(SGN==  66) HN= 313;
  if(SGN==  67) HN= 316;
  if(SGN==  68) HN= 322;
  if(SGN==  69) HN= 334;
  if(SGN==  70) HN= 335;
  if(SGN==  71) HN= 337;
  if(SGN==  72) HN= 338;
  if(SGN==  73) HN= 341;
  if(SGN==  74) HN= 343;
  if(SGN==  75) HN= 349;
  if(SGN==  76) HN= 350;
  if(SGN==  77) HN= 351;
  if(SGN==  78) HN= 352;
  if(SGN==  79) HN= 353;
  if(SGN==  80) HN= 354;
  if(SGN==  81) HN= 355;
  if(SGN==  82) HN= 356;
  if(SGN==  83) HN= 357;
  if(SGN==  84) HN= 358;
  if(SGN==  85) HN= 359;
  if(SGN==  86) HN= 361;
  if(SGN==  87) HN= 363;
  if(SGN==  88) HN= 364;
  if(SGN==  89) HN= 366;
  if(SGN==  90) HN= 367;
  if(SGN==  91) HN= 368;
  if(SGN==  92) HN= 369;
  if(SGN==  93) HN= 370;
  if(SGN==  94) HN= 371;
  if(SGN==  95) HN= 372;
  if(SGN==  96) HN= 373;
  if(SGN==  97) HN= 374;
  if(SGN==  98) HN= 375;
  if(SGN==  99) HN= 376;
  if(SGN== 100) HN= 377;
  if(SGN== 101) HN= 378;
  if(SGN== 102) HN= 379;
  if(SGN== 103) HN= 380;
  if(SGN== 104) HN= 381;
  if(SGN== 105) HN= 382;
  if(SGN== 106) HN= 383;
  if(SGN== 107) HN= 384;
  if(SGN== 108) HN= 385;
  if(SGN== 109) HN= 386;
  if(SGN== 110) HN= 387;
  if(SGN== 111) HN= 388;
  if(SGN== 112) HN= 389;
  if(SGN== 113) HN= 390;
  if(SGN== 114) HN= 391;
  if(SGN== 115) HN= 392;
  if(SGN== 116) HN= 393;
  if(SGN== 117) HN= 394;
  if(SGN== 118) HN= 395;
  if(SGN== 119) HN= 396;
  if(SGN== 120) HN= 397;
  if(SGN== 121) HN= 398;
  if(SGN== 122) HN= 399;
  if(SGN== 123) HN= 400;
  if(SGN== 124) HN= 401;
  if(SGN== 125) HN= 402;
  if(SGN== 126) HN= 404;
  if(SGN== 127) HN= 406;
  if(SGN== 128) HN= 407;
  if(SGN== 129) HN= 408;
  if(SGN== 130) HN= 410;
  if(SGN== 131) HN= 412;
  if(SGN== 132) HN= 413;
  if(SGN== 133) HN= 414;
  if(SGN== 134) HN= 416;
  if(SGN== 135) HN= 418;
  if(SGN== 136) HN= 419;
  if(SGN== 137) HN= 420;
  if(SGN== 138) HN= 422;
  if(SGN== 139) HN= 424;
  if(SGN== 140) HN= 425;
  if(SGN== 141) HN= 426;
  if(SGN== 142) HN= 428;
  if(SGN== 143) HN= 430;
  if(SGN== 144) HN= 431;
  if(SGN== 145) HN= 432;
  if(SGN== 146) HN= 433;
  if(SGN== 147) HN= 435;
  if(SGN== 148) HN= 436;
  if(SGN== 149) HN= 438;
  if(SGN== 150) HN= 439;
  if(SGN== 151) HN= 440;
  if(SGN== 152) HN= 441;
  if(SGN== 153) HN= 442;
  if(SGN== 154) HN= 443;
  if(SGN== 155) HN= 444;
  if(SGN== 156) HN= 446;
  if(SGN== 157) HN= 447;
  if(SGN== 158) HN= 448;
  if(SGN== 159) HN= 449;
  if(SGN== 160) HN= 450;
  if(SGN== 161) HN= 452;
  if(SGN== 162) HN= 454;
  if(SGN== 163) HN= 455;
  if(SGN== 164) HN= 456;
  if(SGN== 165) HN= 457;
  if(SGN== 166) HN= 458;
  if(SGN== 167) HN= 460;
  if(SGN== 168) HN= 462;
  if(SGN== 169) HN= 463;
  if(SGN== 170) HN= 464;
  if(SGN== 171) HN= 465;
  if(SGN== 172) HN= 466;
  if(SGN== 173) HN= 467;
  if(SGN== 174) HN= 468;
  if(SGN== 175) HN= 469;
  if(SGN== 176) HN= 470;
  if(SGN== 177) HN= 471;
  if(SGN== 178) HN= 472;
  if(SGN== 179) HN= 473;
  if(SGN== 180) HN= 474;
  if(SGN== 181) HN= 475;
  if(SGN== 182) HN= 476;
  if(SGN== 183) HN= 477;
  if(SGN== 184) HN= 478;
  if(SGN== 185) HN= 479;
  if(SGN== 186) HN= 480;
  if(SGN== 187) HN= 481;
  if(SGN== 188) HN= 482;
  if(SGN== 189) HN= 483;
  if(SGN== 190) HN= 484;
  if(SGN== 191) HN= 485;
  if(SGN== 192) HN= 486;
  if(SGN== 193) HN= 487;
  if(SGN== 194) HN= 488;
  if(SGN== 195) HN= 489;
  if(SGN== 196) HN= 490;
  if(SGN== 197) HN= 491;
  if(SGN== 198) HN= 492;
  if(SGN== 199) HN= 493;
  if(SGN== 200) HN= 494;
  if(SGN== 201) HN= 495;
  if(SGN== 202) HN= 497;
  if(SGN== 203) HN= 498;
  if(SGN== 204) HN= 500;
  if(SGN== 205) HN= 501;
  if(SGN== 206) HN= 502;
  if(SGN== 207) HN= 503;
  if(SGN== 208) HN= 504;
  if(SGN== 209) HN= 505;
  if(SGN== 210) HN= 506;
  if(SGN== 211) HN= 507;
  if(SGN== 212) HN= 508;
  if(SGN== 213) HN= 509;
  if(SGN== 214) HN= 510;
  if(SGN== 215) HN= 511;
  if(SGN== 216) HN= 512;
  if(SGN== 217) HN= 513;
  if(SGN== 218) HN= 514;
  if(SGN== 219) HN= 515;
  if(SGN== 220) HN= 516;
  if(SGN== 221) HN= 517;
  if(SGN== 222) HN= 518;
  if(SGN== 223) HN= 520;
  if(SGN== 224) HN= 521;
  if(SGN== 225) HN= 523;
  if(SGN== 226) HN= 524;
  if(SGN== 227) HN= 525;
  if(SGN== 228) HN= 527;
  if(SGN== 229) HN= 529;
  if(SGN== 230) HN= 530;
  
  C->NSG = spg_get_symmetry_from_database(rotation,translation,HN);

  for(j = 0; j < C->NSG; j++) 
  {
    C->SG[j][ 0] = (double) rotation[j][0][0];
    C->SG[j][ 1] = (double) rotation[j][0][1];
    C->SG[j][ 2] = (double) rotation[j][0][2];
    C->SG[j][ 3] = (double) rotation[j][1][0];
    C->SG[j][ 4] = (double) rotation[j][1][1];
    C->SG[j][ 5] = (double) rotation[j][1][2];
    C->SG[j][ 6] = (double) rotation[j][2][0];
    C->SG[j][ 7] = (double) rotation[j][2][1];
    C->SG[j][ 8] = (double) rotation[j][2][2];

    C->SG[j][ 9] = translation[j][0];
    C->SG[j][10] = translation[j][1];
    C->SG[j][11] = translation[j][2];
  }

  return;
}
//========================================================
//Interface function for spglib
//========================================================
int sym_dataset(Cell *C, const double origin_shift[3], char w[], int equivalent_atoms[], double tol, int J) 
{
  SpglibDataset *dataset;
  const char    *wl = "abcdefghijklmnopqrstuvwxyz";
  char          symbol[11];
  int           i, q, N;
  int           types[C->A];
  double        lattice[3][3], position[C->A][3];

  for(i=0;i<C->N;i++)
    types[i] = C->ATMN[i]+1;

  for(i=0;i<D3;i++)
    for(q=0;q<D3;q++)
      lattice[i][q] = C->L[q][i];

  for(i=0;i<C->N;i++)
    for(q=0;q<D3;q++)
      position[i][q] = C->X[i][q];

  for(i=0;i<C->N;i++)
    for(q=0;q<D3;q++) 
      position[i][q] += origin_shift[q];


  if(J==0)
  {
    N       = spg_standardize_cell(lattice, position, types, C->N, 1, 0, tol);
    dataset = spg_get_dataset(lattice, position, types, N, tol); 
  }
  else
  {
    N       = spg_standardize_cell(lattice, position, types, C->N, 0, 0, tol);
    dataset = spg_get_dataset(lattice, position, types, N, tol); 

    C->N = dataset->n_std_atoms;

    for(i=0;i<C->NSPC;i++)
      C->SPCN[i] = 0;

    for(i=0;i<C->N;i++)
    {
      C->ATMN[i] = dataset->std_types[i]-1;
      C->ATMZ[i] = C->SPCZ[dataset->std_types[i]-1];
      C->SPCN[C->ATMN[i]]++;

      for(q=0;q<D3;q++)
        C->X[i][q] = dataset->std_positions[i][q];
    }

    for(i=0;i<D3;i++)
      for(q=0;q<D3;q++)
        C->L[q][i] = dataset->std_lattice[i][q];

    ORDER(C);

    for(i=0;i<D3;i++)
      for(q=0;q<D3;q++)
        lattice[i][q] = C->L[q][i];

    for(i=0;i<C->N;i++)
      for(q=0;q<D3;q++)
        position[i][q] = C->X[i][q];

    for(i=0;i<C->N;i++)
      types[i] = C->ATMN[i]+1;

    //===== call again to make Wyckoff positions properly ordered by species =====
    dataset = spg_get_dataset(lattice, position, types, C->N, tol);

    for(i=0;i<C->N;i++)
    {
      w[i]                = wl[dataset->wyckoffs[i]];
      equivalent_atoms[i] = dataset->equivalent_atoms[i];
    }

  }

  spg_get_international(symbol, dataset->std_lattice, dataset->std_positions, dataset->std_types, dataset->n_std_atoms, tol);
  strncpy(C->SGS,symbol,9);

  C->SGN = dataset->spacegroup_number;

  spg_free_dataset(dataset);

  return N;
}
//========================================================
// Neural Network Job
//========================================================
int FIND_WYC(Cell *C, Cell *D, double tol, int J)
{
  int i,ii,j,q,m,n,SGN;

  char   s[200],s1[200];
  char   w[C->A];
  int    equivalent_atoms[C->A];
  double a=180.0/Pi;
  double origin_shift[]={0.0,0.0,0.0};

  FILE *out;

  Relative(C);
  sym_dataset(C,origin_shift,w,equivalent_atoms,tol,J);
  SGN = C->SGN;

  if(J==0)
  {
    Real(C);
    return SGN;
  }

  Relative(D);

  for(i=0;i<D->N;i++)
    for(q=0;q<D3;q++)
    {
      if(D->X[i][q] <0.0)
        D->X[i][q] += 1.0;

      if(D->X[i][q]>=1.0)
        D->X[i][q] -= 1.0;

      for(m=0;m<12;m++)
        if(fabs(D->X[i][q]-(double)m/12.0)< (double)0.0002)
          D->X[i][q] = (double)m/12.0;
    }

  abc(C);

  out = fopen("str.cif","w");

  fprintf(out,"_symmetry_Int_Tables_number %d\n",C->SGN);
  fprintf(out,"_cell_length_a   % .16lf\n",C->LAT[0]);
  fprintf(out,"_cell_length_b   % .16lf\n",C->LAT[1]);
  fprintf(out,"_cell_length_c   % .16lf\n",C->LAT[2]);
  fprintf(out,"_cell_angle_alpha% .16lf\n",C->ANG[0]*a);
  fprintf(out,"_cell_angle_beta % .16lf\n",C->ANG[1]*a);
  fprintf(out,"_cell_angle_gamma% .16lf\n",C->ANG[2]*a);
  fprintf(out,"_symmetry_Int_Tables_number %d\n",C->SGN);
  fprintf(out,"_chemical_formula_sum\n");
  fprintf(out,"'");
  for(i=0;i<C->ATMN[C->N-1]+1;i++)
  {
    if(C->SPCZ[0]>0)
    {
      atom_symb(C->SPCZ[i],s1);    
      fprintf(out,"%s ",s1);
    }
    else
    {
      fprintf(out,"%c ",65+i);
    }
  }
  fprintf(out,"'\n");
  fprintf(out,"loop_\n");
  fprintf(out,"_atom_site_label\n");
  fprintf(out,"_atom_site_type_symbol\n");
  fprintf(out,"_atom_site_site_symmetry_multiplicity\n");
  fprintf(out,"_atom_site_occupancy wyckoff\n");
  fprintf(out,"_atom_site_fract_x\n");
  fprintf(out,"_atom_site_fract_y\n");
  fprintf(out,"_atom_site_fract_z\n");

  n=0;
  j=0;

  for(i=0;i<C->N;i++)
  {
    if(i>0 && C->ATMZ[i]!=C->ATMZ[i-1])
      j=0;

    if(C->SPCZ[0]>0)
      atom_symb(C->ATMZ[i],s);

    for(q=0;q<D3;q++) {
      if(C->X[i][q] <0.0)
        C->X[i][q] += 1.0;

      if(C->X[i][q]>=1.0)
        C->X[i][q] -= 1.0;

      for(m=0;m<12;m++)
        if(fabs(C->X[i][q]-(double)m/12.0)< (double)0.0002)
          C->X[i][q] = (double)m/12.0;
    }

    if( equivalent_atoms[i]==i )
    {
      for(ii=0;ii<C->N;ii++)
        if(equivalent_atoms[ii]==equivalent_atoms[i])
          n++;

      fprintf(out,"%s%d %s  %3d %c % 8.16lf % 8.16lf % 8.16lf\n",s,j+1,s,n,w[i],C->X[i][0],C->X[i][1],C->X[i][2]);
      j++;
      n=0;
    }
  }
  fprintf(out,"#End\n");
  fclose(out);

  for(i=0;i<C->ATMN[C->N-1]+1;i++)
  {
    if(C->SPCZ[i]>0) 
      atom_symb(C->SPCZ[i],s1); 
    else 
      strcpy(s1,"A");
    sprintf(C->TAG+2*i,"%s ",s1);
  }

  Real(C);
  Real(D);

  return SGN;
}
//========================================================
