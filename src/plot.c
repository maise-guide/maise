#include "plot.h"

//=============================================================================
void PLOT(ANN *R, int type, int N1, int N2, int M, char *test, char *CMP)
{
  double xmin = 2.0;
  double xmax = 4.9;
  double ymin =-2.3;
  double ymax =-0.0;

  int smart_range = 1;  //set to 0 for manually adjusting ranges, 1 for finding them from CONTCARS

  int     color[]={0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9};

  FILE   *in[101],*ot[101];
  FILE   *out;
  int    i,j,t0;
  char   s[101][200],u[101][200],buf[101][200],nntp[100][2];
  double tmp,data[101][101],t1,t2;
  char   name[400];
  char   compound[400];
  char   ctmp[200];

  printf("%s\n",CMP);

  sprintf(nntp[ 1],"M");
  sprintf(nntp[ 2],"G");
  sprintf(nntp[ 3],"S");
  sprintf(nntp[ 4],"L");

  sprintf(s[ 0],"%s/%s%s00bcc.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[ 1],"%s/%s%s01fcc.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[ 2],"%s/%s%s02hcp.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[ 3],"%s/%s%s03sc_.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[ 4],"%s/%s%s04dia.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[ 5],"%s/%s%s05bsn.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[ 6],"%s/%s%s06dim.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[ 7],"%s/%s%s07tri.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[ 8],"%s/%s%s08sqr.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[ 9],"%s/%s%s09nps.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[20],"%s/%s%s20000.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[21],"%s/%s%s21001.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[22],"%s/%s%s22002.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[23],"%s/%s%s23003.dat",R->otpt,nntp[R->MODT],CMP);
  sprintf(s[24],"%s/%s%s24004.dat",R->otpt,nntp[R->MODT],CMP);

  sprintf(u[ 0],"bcc");
  sprintf(u[ 1],"fcc");
  sprintf(u[ 2],"hcp");
  sprintf(u[ 3],"sc_");
  sprintf(u[ 4],"dia");
  sprintf(u[ 5],"bsn");
  sprintf(u[ 6],"dim");
  sprintf(u[ 7],"tri");
  sprintf(u[ 8],"sqr");
  sprintf(u[ 9],"nps");
  sprintf(u[20],"000");
  sprintf(u[21],"001");
  sprintf(u[22],"002");
  sprintf(u[23],"003");
  sprintf(u[24],"004");

  //===== EOS =====
  if( type == 0 )
  {  
    for(i=N1; i < N2;i++)
    {
      in[i]=fopen(s[i],"r");
      ot[i]=fopen(u[i],"w");
    }

    for(i=N1; i < N2;i++)
      for(j=0; j < M;j++)
      {
	fscanf( in[i],"%lf %lf %lf %lf\n",&data[j][i*3],&data[j][i*3+1],&data[j][i*3+2],&tmp);
	fprintf(ot[i],"%lf %lf %lf\n",data[j][i*3],data[j][i*3+1],data[j][i*3+2]);
      }

    for(i=N1; i < N2;i++) 
    {
      fclose(in[i]);
      fclose(ot[i]);
    } 
    //Finding min and max range automatically
    if(smart_range)
    {
      xmax = xmin = data[0][N1*3+0];
      for(i=N1; i < N2;i++) 
	for(j=0; j < M;j++)
	{
	  if( data[j][i*3] < xmin ) xmin=data[j][i*3];
	  if( data[j][i*3] > xmax ) xmax=data[j][i*3];
	}

      ymax = ymin = data[0][N1*3+1];
      for(i=N1; i < N2;i++) 
	for(j=0; j < M;j++)
	{
	  if( data[j][i*3+1] < ymin ) ymin=data[j][i*3+1];
          if( data[j][i*3+2] < ymin ) ymin=data[j][i*3+2];
	  if( data[j][i*3+1] > ymax ) ymax=data[j][i*3+1];
          if( data[j][i*3+2] > ymax ) ymax=data[j][i*3+2];
	}

      if(ymax>3.0)
	ymax = 3.0;

      out=popen("gnuplot","w");
      fprintf(out,"set terminal pngcairo size 800, 600\n");
      fprintf(out,"set title \"%s (%s)\"\n",CMP,test);
      fprintf(out,"set xrange [%lf:%lf]\n",xmin-abs(xmin*0.05),xmax+abs(xmax*0.05));
      fprintf(out,"set yrange [%lf:%lf]\n",ymin-abs(ymin*0.05),ymax+abs(ymax*0.05));
      fprintf(out,"set autoscale fix\n");
      fprintf(out,"set size 1.0,1.0\n");
      fprintf(out,"set ylabel \"energy (eV/atom)\"\n");
      fprintf(out,"set xlabel \"N-N distance (A)\"\n");
      fprintf(out,"set output \'%s/%s-%s.png\'\n",R->otpt,CMP,test);
      fprintf(out,"set multiplot\n");
      fprintf(out,"set key top\n");
      fprintf(out,"plot ");

      for(i=N1; i < N2;i++) 
      {
	fprintf(out,"  \"%s\" using 1:3 with lines t \"%s\" lw 2 lc %d",u[i],u[i],color[i]); 
	if( i != N2-1 ) 
	  fprintf(out,",");
      }
      fprintf(out,"\nset key off\nplot ");

      for(i=N1; i < N2;i++) 
      {
	fprintf(out," \"%s\" u 1:2 w points pt 7 lc %d",u[i],color[i]); 
	if( i != N2-1 ) 
	  fprintf(out,",");
      }
      fprintf(out,"\n");
      fprintf(out,"unset multiplot\n");
      fclose(out);


      for(i=N1; i < N2;i++) 
      {
	sprintf(name,"rm %s",u[i]);  
	system(name);  
      }
    }
  }
  //===== VAC =====
  if( type == 1 )
  {    
    sprintf(name,"%s/vac-%s.dat",R->otpt,CMP);
    in[0]=fopen(name,"r");
    fscanf(in[0],"%d %lf %lf %s %lf\n",&t0,&t1,&t2,buf[0],&tmp);
    ymin=t1;
    ymax=t1;
    for(i=N1; i < N2;i++)
    {
      if( ymin > t1 ) ymin=t1;
      if( ymax < t1 ) ymax=t1;
      if( ymin > t2 ) ymin=t2;
      if( ymax < t2 ) ymax=t2;
      fscanf(in[0],"%d %lf %lf %s %lf\n",&t0,&t1,&t2,buf[0],&tmp);
    }
    xmin=0.0;
    xmax=(double)N2-1;
    fclose(in[0]);

    strcpy(buf[0],"");
    for(i=N1; i < N2;i++)
    {strcat(buf[0],u[i]);strcat(buf[0],"-");}

    out=popen("gnuplot","w");
    fprintf(out,"set terminal pngcairo size 800, 600\n");
    fprintf(out,"set title \"%s (%s)\"\n",CMP,test);
    fprintf(out,"set xrange [%lf:%lf]\n",xmin-abs(xmin*0.05),xmax+abs(xmax*0.05));
    fprintf(out,"set yrange [%lf:%lf]\n",ymin-abs(ymin*0.05),ymax+abs(ymax*0.05));
    fprintf(out,"set autoscale fix\n");
    fprintf(out,"set size 1.0,1.0\n");
    fprintf(out,"set ylabel \"energy (eV/vac)\"\n");
    fprintf(out,"set xlabel \"structure\"\n");
    fprintf(out,"set output \'%s/%s-vac.png\'\n",R->otpt,CMP);
    fprintf(out,"set multiplot\n");
    fprintf(out,"set style line 1 lc rgb \'#0060ad\' lt 1 lw 2 pt 4  ps 1.5\n");
    fprintf(out,"set style line 2 lc rgb \'#dd181f\' lt 1 lw 2 pt 13 ps 1.5\n");
    fprintf(out,"set key top right\n");
    fprintf(out,"plot \"%s\" using 2:xticlabels(4) with linespoints t \"DFT\" ls 1, \"%s\" using 3 with linespoints t \"ANN\" ls 2\n",name,name);
    fprintf(out,"set key off\n");
    fprintf(out,"plot \"\" using 1:3:5:xticlabels(4) with labels offset 0,char 1 tc rgb \"red\"\n");
    fprintf(out,"unset multiplot\n");
    fclose(out);
  }
  //======================== SUB 
  if( type == 2 )
  {  
    sprintf(compound,"%s/sub-%s.dat",R->otpt,CMP);
    in[0]=fopen(compound,"r");
    fgets(name,200,in[0]);
    sscanf(name,"%d %lf %lf %lf %lf %s %lf %lf\n",&j,&data[0][0],&data[0][1],&data[0][2],&data[0][3],ctmp,&tmp,&tmp);
    printf("%d\n",j);
    xmin=xmax=(double)j;
    if( data[0][0] < data[0][1]) {ymin=data[0][0];ymax=data[0][1];}
    else {ymin=data[0][1];ymax=data[0][0];}
    if( data[0][2] < data[0][3] ) {t1=data[0][2];t2=data[0][3];}
    else {t1=data[0][3];t2=data[0][2];}
    printf("XY %lf %lf\n",xmin,xmax);
    while( fgets(name,200,in[0]) )
    {
      sscanf(name,"%d %lf %lf %lf %lf %s %lf %lf\n",&j,&data[0][0],&data[0][1],&data[0][2],&data[0][3],ctmp,&tmp,&tmp);
      if( data[0][0] < ymin ) ymin=data[0][0];
      if( data[0][1] < ymin ) ymin=data[0][1];
      if( data[0][0] > ymax ) ymax=data[0][0];
      if( data[0][1] > ymax ) ymax=data[0][1];
      if( data[0][2] < t1 ) t1=data[0][2];
      if( data[0][3] < t1 ) t1=data[0][3];
      if( data[0][2] > t2 ) t2=data[0][2];
      if( data[0][3] > t2 ) t2=data[0][3];
      if( j < xmin ) xmin=(double)j;
      if( j > xmax ) xmax=(double)j;
      printf("XY %lf %lf\n",xmin,xmax);
    }

    for(i=1; i < 3;i++)
    {
      out=popen("gnuplot","w");
      fprintf(out,"set terminal pngcairo size 800, 600\n");
      fprintf(out,"set title \"%s (%s-%d in %d)\"\n",CMP,test,((i+2)%2)+1,i);
      fprintf(out,"set xrange [%lf:%lf]\n",xmin-abs(xmin*0.05),xmax+abs(xmax*0.05));
      fprintf(out,"set yrange [%lf:%lf]\n",ymin-abs(ymin*0.05),ymax+abs(ymax*0.05));
      fprintf(out,"set autoscale fix\n");
      fprintf(out,"set size 1.0,1.0\n");
      fprintf(out,"set ylabel \"energy (eV/vac)\"\n");
      fprintf(out,"set xlabel \"structure\"\n");
      fprintf(out,"set output \'%s/%s-sub-%d.png\'\n",R->otpt,CMP,i);
      fprintf(out,"set multiplot\n");
      fprintf(out,"set style line 1 lc rgb \'#0060ad\' lt 1 lw 2 pt 4  ps 1.5\n");
      fprintf(out,"set style line 2 lc rgb \'#dd181f\' lt 1 lw 2 pt 13 ps 1.5\n");
      fprintf(out,"set key top right\n");
      fprintf(out,"plot \"%s\" using %d:xticlabels(6) with linespoints t \"DFT\" ls 1, \"%s\" using %d with linespoints t \"ANN\" ls 2\n",compound,i*2,compound,i*2+1);
      fprintf(out,"set key off\n");
      fprintf(out,"plot \"\" using 1:%d:%d:xticlabels(6) with labels offset 0,char 1 tc rgb \"red\"\n",i*2+1,6+i);
      fprintf(out,"unset multiplot\n");
      fclose(out);
      ymin=t1;
      ymax=t2;
    }
    fclose(in[0]);
  }

  //======================== PHD
  if( type == 3 )
  {  
    sprintf(compound,"%s/phd-%s00.dat",R->otpt,CMP);
    in[0]=fopen(compound,"r");
    fgets(name,200,in[0]);
    sscanf(name,"%lf %lf %lf %d\n",&data[0][0],&data[0][1],&data[0][2],&j);
    xmin=xmax=data[0][0];
    if( data[0][1] < data[0][2] ) {ymin=data[0][1];ymax=data[0][2];}
    else {ymin=data[0][2];ymax=data[0][1];}

    while( fgets(name,200,in[0]) )
    {
      sscanf(name,"%lf %lf %lf %d\n",&data[0][0],&data[0][1],&data[0][2],&j);
      if( data[0][1] < ymin ) ymin=data[0][1];
      if( data[0][2] < ymin ) ymin=data[0][2];
      if( data[0][1] > ymax ) ymax=data[0][1];
      if( data[0][2] > ymax ) ymax=data[0][2];
      if( data[0][0] < xmin ) xmin=data[0][0];
      if( data[0][0] > xmax ) xmax=data[0][0];
    }

    out=popen("gnuplot","w");
    fprintf(out,"set terminal pngcairo size 800, 600\n");
    fprintf(out,"set title \"%s (phd)\"\n",CMP);
    fprintf(out,"set xrange [%lf:%lf]\n",xmin-abs(xmin*0.05),xmax+abs(xmax*0.05));
    fprintf(out,"set yrange [%lf:%lf]\n",ymin-abs(ymin*0.05),ymax+abs(ymax*0.05));
    fprintf(out,"set autoscale fix\n");
    fprintf(out,"set size 1.0,1.0\n");
    fprintf(out,"set ylabel \"energy (eV/vac)\"\n");
    fprintf(out,"set xlabel \"composition ratio\"\n");
    fprintf(out,"set output \'%s/%s-phd.png\'\n",R->otpt,CMP);
    fprintf(out,"set multiplot\n");
    fprintf(out,"set style line 1 lc palette lt 1 lw 2 pt 4  ps 1.5\n");
    fprintf(out,"set style line 2 lc palette lt 1 lw 2 pt 13 ps 1.5\n");
    fprintf(out,"set key bottom right\n");
    fprintf(out,"unset colorbox\n");
    fprintf(out,"plot \"%s\" using 1:2:4:xticlabels(1) t \"DFT\" with points ls 1,\"%s\" using 1:3:4:xticlabels(1) t \"ANN\" with points ls 2 \n",compound,compound);
    fprintf(out,"unset multiplot\n");
    fclose(out);
    fclose(in[0]);
  }
  //===== SRF =====
  if( type == 4 )
  {    
    sprintf(name,"%s/srf-%s.dat",R->otpt,CMP);
    in[0]=fopen(name,"r");
    fscanf(in[0],"%d %lf %lf %s %lf\n",&t0,&t1,&t2,buf[0],&tmp);
    ymin=t1;
    ymax=t1;
    for(i=N1; i < N2;i++)
    {
      if( ymin > t1 ) ymin=t1;
      if( ymax < t1 ) ymax=t1;
      if( ymin > t2 ) ymin=t2;
      if( ymax < t2 ) ymax=t2;
      fscanf(in[0],"%d %lf %lf %s %lf\n",&t0,&t1,&t2,buf[0],&tmp);
    }
    xmin=0.0;
    xmax=(double)N2-1;
    fclose(in[0]);

    strcpy(buf[0],"");
    for(i=N1; i < N2;i++)
    {strcat(buf[0],u[i]);strcat(buf[0],"-");}

    out=popen("gnuplot","w");
    fprintf(out,"set terminal pngcairo size 800, 600\n");
    fprintf(out,"set title \"%s (%s)\"\n",CMP,test);
    fprintf(out,"set xrange [%lf:%lf]\n",xmin-abs(xmin*0.05),xmax+abs(xmax*0.05));
    fprintf(out,"set yrange [%lf:%lf]\n",ymin-abs(ymin*0.05),ymax+abs(ymax*0.05));
    fprintf(out,"set autoscale fix\n");
    fprintf(out,"set size 1.0,1.0\n");
    fprintf(out,"set ylabel \"energy (eV/A^2)\"\n");
    fprintf(out,"set xlabel \"structure\"\n");
    fprintf(out,"set output \'%s/%s-srf.png\'\n",R->otpt,CMP);
    fprintf(out,"set multiplot\n");
    fprintf(out,"set style line 1 lc rgb \'#0060ad\' lt 1 lw 2 pt 4  ps 1.5\n");
    fprintf(out,"set style line 2 lc rgb \'#dd181f\' lt 1 lw 2 pt 13 ps 1.5\n");
    fprintf(out,"set key top right\n");
    fprintf(out,"plot \"%s\" using 2:xticlabels(4) with linespoints t \"DFT\" ls 1, \"%s\" using 3 with linespoints t \"ANN\" ls 2\n",name,name);
    fprintf(out,"set key off\n");
    fprintf(out,"plot \"\" using 1:3:5:xticlabels(4) with labels offset 0,char 1 tc rgb \"red\"\n");
    fprintf(out,"unset multiplot\n");
    fclose(out);
  }
}
