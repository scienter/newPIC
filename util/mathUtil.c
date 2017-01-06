#include <stdio.h>
#include <gsl/gsl_linalg.h>


double rounding(double x, int digit);


void findFactor(double *result,double *a_data,double *b_data,int row,int col,int indexI,int indexJ)
{
  int s,i,j,digit=4;
  double data[16],tmp,unit=0.999999;

  for(i=0; i<col*row; i++)  {
    tmp=a_data[i]*unit;
    data[i]=rounding(tmp,digit);
    a_data[i]=data[i];
  }

  gsl_matrix_view m  =gsl_matrix_view_array(a_data,row,col);
//  gsl_matrix_view inv=gsl_matrix_view_array(invA,row,col);
  gsl_vector_view b  = gsl_vector_view_array(b_data,row);


//  printf("The a_data matrix\n");
//  for(i=0; i<4; i++)
//    for(j=0; j<4; j++)
//      printf(j==3? "%g\n":"%g ",a_data[i*4+j]);

//  printf("The b_data matrix\n");
//  for(i=0; i<4; i++)
//      printf(j==3? "%g\n":"%g ",b_data[i]);


  gsl_vector *x = gsl_vector_alloc(col);
  gsl_vector *work = gsl_vector_alloc(col);
  gsl_matrix *V = gsl_matrix_alloc(col,col);
  gsl_vector *S = gsl_vector_alloc(col);
  gsl_permutation * p = gsl_permutation_alloc(row);

  gsl_linalg_SV_decomp(&m.matrix,V,S,work);
  gsl_linalg_SV_solve(&m.matrix,V,S,&b.vector,x);

//  printf("x=\n");
//  gsl_vector_fprintf(stdout,x,"%g");
  for(i=0; i<4; i++)
    result[i]=gsl_vector_get(x,i);

//  for(i=0; i<4; i++)
//    printf("h[%d]=%g\n",i,result[i]);

  gsl_permutation_free(p);
  gsl_vector_free(x);
  gsl_vector_free(work);
  gsl_vector_free(S);
  gsl_matrix_free(V);

}

double rounding(double x, int digit)
{
  double tmp,result;

  tmp=pow(10.0,digit);
  result=floor(x*tmp+0.5f)/tmp;
  return result;
}

