#include<stdio.h>
#include<math.h>
#include<stdlib.h>

long double invh[10][10];    //array to store the inverse of the Hessian matrix
long double bpr1[10];       //arrays for values to be returned by bounding phase function
long double bpr2[10];
long double ihr[10];        //array for value to be returned by interval hlaving function
FILE *fpo;                  //output file pointer is made a global variable so that any fuction can write to output file

long double objf(long double * x, int n, int q) {             //objective function, returns different values based on value of q and n
  long double y=0;
  int i=0;
  if (q==1) {                   //sum squares function
    for (i=0;i<n;i++) {
      y=y+((i+1)*(x[i]*x[i]));
    }
  }
  else if (q==2) {              //Rosenbrock function
    for (i=0;i<n-1;i++) {
      y=y+(100*pow(x[i+1]-x[i]*x[i],2)+pow(x[i]-1,2));
    }
  }
  else if (q==3) {             //Dixon Price function
    y=pow(x[0]-1,2);
    for (i=1;i<n;i++) {
      y=y+((i+1)*pow(2*x[i]*x[i]-x[i-1],2));
    }
  }
  else if (q==4) {            //Trid function
    long double y1=0;
    long double y2=0;
    for (i=0;i<n;i++) {
      y1=y1+pow(x[i]-1,2);
    }
    for (i=1;i<n;i++) {
      y2=y2+(x[i]*x[i-1]);
    }
    y=y1-y2;
  }
  else if (q==5) {            //Zakharov function
    long double y1=0;
    long double y2=0;
    for (i=0;i<n;i++) {
      y1=y1+(x[i]*x[i]);
    }
    for (i=0;i<n;i++) {
      y2=y2+(0.5*(i+1)*x[i]);
    }
    y=y1+pow(y2,2)+pow(y2,4);
  }
  return y;
}

long double * objfd(long double * x, int n, int q) {        //gradient of the objective function, returns different values based on the values of q and n
  static long double * y;            //array to be returned
  y=(long double *) malloc(n*sizeof(long double));
  int i=0;
  if (q==1) {                       //sum squares function
    for (i=0;i<n;i++) {
      y[i] = 2*(i+1)*x[i];
    }
  }
  else if (q==2) {                 //Rosenbrock function
    y[0] = 2*x[0]-2-400*x[0]*(x[1]-x[0]*x[0]);
    for (i=1;i<n-1;i++) {
      y[i] = -400*x[i]*(x[i+1]-x[i]*x[i])+2*x[i]-2+200*(x[i]-x[i-1]*x[i-1]);
    }
    y[n-1] = 2*x[n-1]-2+200*(x[n-1]-x[n-2]*x[n-2]);
  }
  else if (q==3) {                //Dixon Price function
    y[0]=2*(x[0]-1)-4*(2*x[1]*x[1]-x[0]);
    y[n-1]=8*n*x[n-1]*(2*x[n-1]*x[n-1]-x[n-2]);
    for (i=1;i<n-1;i++) {
      y[i]=8*(i+1)*x[i]*(2*x[i]*x[i]-x[i-1])-2*(i+2)*(2*x[i+1]*x[i+1]-x[i]);
    }
  }
  else if (q==4) {              //Trid function
    y[0]=2*(x[0]-1)-x[1];
    y[n-1]=2*(x[n-1]-1)-x[n-2];
    for (i=1;i<n-1;i++) {
      y[i]=2*(x[i]-1)-x[i-1]-x[i+1];
    }
  }
  else if (q==5) {            //Zakharov function
    long double y1=0;
    for (i=0;i<n;i++) {
      y1=y1+(0.5*(i+1)*x[i]);
    }
    for (i=0;i<n;i++) {
      y[i]=2*x[i]+(i+1)*y1+2*(i+1)*pow(y1,3);
    }
  }
  return y;
}

void objfdd(long double * x, int n, int q) {      //This function calculates the Hessian matrix based on the value of q, n and x. Then it inverts the Hessian matrix and stores the resulting matrix in the global variable invh
  long double y[10][20];
  long double factor=0;
  int i=0;
  int j=0;
  if (q==1) {                     //Sum squares function
    for (i=0;i<n;i++) {
      for (j=0;j<n;j++) {
        if (i==j) {
          y[i][j]=2*(i+1);
        }
        else {
          y[i][j]=0;
        }
      }
    }
  }
  else if (q==2) {                //Rosenbrock function
    y[0][0]=-400*(x[1]-3*x[0]*x[0])+2;
    y[n-1][n-1]=200;
    for (i=0;i<n;i++) {
      for (j=i;j<n;j++) {
        if (i==0 && j==0) {
          continue;
        }
        else if (i==n-1 && j==n-1) {
          continue;
        }
        else if (i==j) {
          y[i][j]=-400*(x[i+1]-3*x[i]*x[i])+202;
        }
        else if (j==i+1) {
          y[i][j]=-400*x[i];
          y[j][i]=y[i][j];
        }
        else if (j>i+1) {
          y[i][j]=0;
          y[j][i]=y[i][j];
        }
      }
    }
  }
  else if (q==3) {                  //Dixon Price function
    y[0][0]=6;
    y[n-1][n-1]=8*n*(6*x[n-1]*x[n-1]-x[n-2]);
    for (i=0;i<n;i++) {
      for (j=i;j<n;j++) {
        if (i==0 && j==0) {
          continue;
        }
        else if (i==n-1 && j==n-1) {
          continue;
        }
        else if (j==i) {
          y[i][j]=8*(i+1)*(6*x[i]*x[i]-x[i-1])+2*(i+2);
        }
        else if (j==i+1) {
          y[i][j]=-8*(i+2)*x[i+1];
          y[j][i]=y[i][j];
        }
        else if (j>i+1) {
          y[i][j]=0;
          y[j][i]=y[i][j];
        }
      }
    }
  }
  else if (q==4) {              //Trid function
    for (i=0;i<n;i++) {
      for (j=i;j<n;j++) {
        if (i==j) {
          y[i][j]=2;
        }
        else if (j==i+1) {
          y[i][j]=-1;
          y[j][i]=y[i][j];
        }
        else if (j>i+1) {
          y[i][j]=0;
          y[j][i]=y[i][j];
        }
      }
    }
  }
  else if (q==5) {            //Zakharov function
    long double y1=0;
    for (i=0;i<n;i++) {
      y1=y1+(0.5*(i+1)*x[i]);
    }
    for (i=0;i<n;i++) {
      for (j=0;j<n;j++) {
        if (i==j) {
          y[i][j]=2+0.5*(i+1)*(i+1)+3*(i+1)*(i+1)*y1*y1;
        }
        else {
          y[i][j]=0;
        }
      }
    }
  }

  int k=0;                  //the augmented matrix is prepared by inserting the nXn identity matrix in front of the Hessian matrix in y
  int l=0;
  for (i=0;i<n;i++) {
    for (j=n;j<2*n;j++) {
      if (j==i+n) {
        y[i][j]=1;
      }
      else {
        y[i][j]=0;
      }
    }
  }

  for (i=0;i<n;i++) {         //the inverse of the Hessian matrix is found using Gaussian row elimination
    for (j=i+n;j>=i;j--) {
      y[i][j]=y[i][j]/y[i][i];
    }
    for (k=i+1;k<n;k++) {
      factor=y[k][i];
      for (l=i;l<=i+n;l++) {
        y[k][l]=y[k][l]-(factor*y[i][l]);
      }
    }
  }

  for (i=n-1;i>=0;i--) {
    for (k=i-1;k>=0;k--) {
      factor = y[k][i];
      for (l=(2*n)-1;l>=i;l--) {
        y[k][l]=y[k][l]-(factor*y[i][l]);
      }
    }
  }

  for (i=0;i<n;i++) {           //the inverse of the Hessian matrix is stored in the global variable invh
    for (j=0;j<n;j++) {
      invh[i][j]=y[i][j+n];
    }
  }
}

void bphase(int q, int var, long double * x0, long double * s0) {       //function for bounding phase method

  long double * x1;
  x1=(long double *) malloc(var*sizeof(long double));
  long double * x2;
  x2=(long double *) malloc(var*sizeof(long double));
  long double * xn;
  xn=(long double *) malloc(var*sizeof(long double));
  long double * xn1;
  xn1=(long double *) malloc(var*sizeof(long double));
  long double delta=0;
  int k=0;
  long double f0 = 0;
  long double f1 = 0;
  long double f2 = 0;
  long double fn = 0;
  int i=0;


  //step 1
  //initial guess is x0
  delta = 0.001;      //increment

  //step 2
  for (i=0;i<var;i++) {
    x1[i] = x0[i]-(delta*s0[i]);        //x1 = x0 - delta
    x2[i] = x0[i]+(delta*s0[i]);        //x2= x0 + delta
  }
  f0=objf(x0,var,q);
  f1=objf(x1,var,q);
  f2=objf(x2,var,q);

  if (f1>=f0 && f0>=f2) {           //if f(x0 - delta) > f(x0) > f(x0 + delta), delta is positive

  }
  else if (f1<=f0 && f0<=f2) {        //if f(x0 - delta) < f(x0) < f(x0 + delta), delta is negative
    delta = -delta;
  }

  //step 3
  label:
  for (i=0;i<var;i++) {               //compute the new point, xn
    xn[i] = x0[i] + (pow(2,k)*delta*s0[i]);
  }

  //step 4
  fn = objf(xn,var,q);
  if (fn<f0) {        //if f(xn) < f(x1), k = k+1 and go to step 3
    k++;
    for (i=0;i<var;i++) {
      xn1[i]=x0[i];
    }
    for (i=0;i<var;i++) {
      x0[i]=xn[i];
    }
    goto label;
  }

  for (i=0;i<var;i++) {       //the optima lies between x(k-1) and x(k+1). These values are stored in the global variables bpr1 and bpr2
    bpr1[i]=xn1[i];
    bpr2[i]=xn[i];
  }
}

void ihalving(long double * a, long double * b, int q, int var) {             //function for interval halving method. a and b are the lower and upper bounds respectively

  //step 1
  long double e = 0.001;        //small e is chosen
  long double * xm;             //xm, l and f(xm) are calculated
  xm = (long double *) malloc(var*sizeof(long double));
  long double * l;
  l = (long double *) malloc(var*sizeof(long double));
  int i=0;
  for (i=0;i<var;i++) {        //xm = (a+b)/2
    xm[i] = 0.5 * (a[i]+b[i]);
  }
  for (i=0;i<var;i++) {        //l = b-a
    l[i] = b[i] - a[i];
  }
  long double fm = objf(xm,var,q);        //f(xm) is computed
  long double * x1;
  x1 = (long double *) malloc(var*sizeof(long double));
  long double * x2;
  x2 = (long double *) malloc(var*sizeof(long double));
  long double f1=0;
  long double f2=0;
  long double l1=0;

  label2: ;
  //step 2
  //new equidistant points are calculated
  for (i=0;i<var;i++) {             //x1 = a + l/4
    x1[i] = a[i] + (0.25*l[i]);
  }
  for (i=0;i<var;i++) {             //x2 = b - l/4
    x2[i] = b[i] - (0.25*l[i]);
  }
  f1 = objf(x1,var,q);              //f(x1), f(x2) and f(xm) are calculated
  f2 = objf(x2,var,q);
  fm = objf(xm,var,q);

  //step 3
  if (f1 < fm) {        //if f(x1) < f(xm), eliminate region xm to b and go to termination condition
    for (i=0;i<var;i++) {
      b[i] = xm[i];
      xm[i] = x1[i];
    }
    goto label1;
  }

  //step 4
  if (f2 < fm) {        //if f(x2) < f(xm), eliminate region a to xm and go to termination condition
    for (i=0;i<var;i++) {
      a[i]=xm[i];
      xm[i]=x2[i];
    }
  }
  else {              //else eliminate regions a to x1 and x2 to b and go to termination condition
    for (i=0;i<var;i++) {
      a[i]=x1[i];
      b[i]=x2[i];
    }
  }

  //step 5
  label1:;
  for (i=0;i<var;i++) {       //new value of l = b-a is calculated
    l[i] = b[i] - a[i];
  }
  l1=0;
  for (i=0;i<var;i++) {       //the modulus of l is calculated and stored in l1
    l1 = l1 + pow(l[i],2);
  }
  l1=sqrt(l1);
  if (l1 < e) {                      //termination condition
    for (i=0;i<var;i++) {
      ihr[i] = 0.5*(a[i] + b[i]);
    }
  }
  else {                //if not terminated, iterate again
    goto label2;
  }
}

long double * nwtnm(int q, int var, long double * x0) {     //function for Newton's method
  //step 1
  static long double ret[10];        //coordinates of the minimum point to be returned
  int i=0;
  int j=0;
  int m=100;                         //maximum number of iterations
  long double e1=0.001;             //termination parameter
  long double e2=0.001;             //termination parameter
  int k=0;
  long double * yd;
  long double ydm=0;
  long double * ydk;
  int itc=0;                        //iteration counter to store the current iteration number
  long double * temp;
  temp=(long double *) malloc(var*sizeof(long double));
  long double temp1=0;
  long double l=0;
  long double x0m=0;

  fprintf(fpo,"*************Newton's Method*************\n\n");
  fprintf(fpo,"The initial point is ");
  for (i=0;i<var-1;i++) {
    fprintf(fpo,"%Lf, ",x0[i]);
  }
  fprintf(fpo,"%Lf.\n\n",x0[var-1]);
  fprintf(fpo,"It#\t");
  for(i=0;i<var;i++) {
    fprintf(fpo,"x%d\t",i+1);
  }
  fprintf(fpo,"f(x)\t||f'(x)||\t|f'(x(k+1)).f'(x(k))|\t||x(k+1)-x(k)||/||x(k)||\n");

  //step 2
  label3: ;
  itc++;                      //iteration counter incremented at the start of each new iteration
  yd = objfd(x0,var,q);       //first derivative is calculated at point x(k)

  //step 3
  for (i=0;i<var;i++) {       //calculate the modulus of first derivative
    ydm=ydm+(yd[i]*yd[i]);
  }
  ydm=sqrt(ydm);
  if (ydm<e1) {               //terminate if modulus of first derivative is less than first termination parameter
    for (i=0;i<var;i++) {
      ret[i]=x0[i];
    }
    return ret;
  }
  else if (itc>m) {           //terminate if the maximum number of iterations exceeded
    printf("Algorithm did not converge\n");
    for (i=0;i<var;i++) {
      ret[i]=x0[i];
    }
    return ret;
  }

  //step 4
  objfdd(x0,var,q);           //inverse of the hessian matrix at x(k) is calculated. The inverse is stored in 2-D global variable matrix invh
  for (i=0;i<var;i++) {        //Inverse of the hessian matrix is then multiplied with gradient of the function and the resulting 1-D vector is stored in temp
    temp[i]=0;
    for (j=0;j<var;j++) {
      temp[i]=temp[i]+(invh[i][j]*yd[j]);
    }
  }
  bphase(q,var,x0,temp);      //using bounding phase and interval halving methods, a point from x(k) along the direction temp is calculated such that the objective function is minimised.
  ihalving(bpr1,bpr2,q,var);    //The bounding phase function stores the range in the global variables bpr1 and bpr2. These variables are then passed to interval halving function, which stores the minima in the global variable ihr

  fprintf(fpo,"%d\t",itc);    //the output is written to the output file. This includes the iteration number, value of minima found in current iteration, the function value at this point and the termination conditions
  for (i=0;i<var;i++) {
    fprintf(fpo,"%Lf\t",ihr[i]);
  }
  fprintf(fpo,"%Lf\t",objf(ihr,var,q));

  ydk=objfd(ihr,var,q);       //first derivative of the objective function is calculated at the new point. If the dot product of the first derivative at the new point and the first derivative at the old point is less than the termination parameter, the function is terminated and the new point is returned
  fprintf(fpo,"%Lf\t",ydm);
  temp1=0;
  for (i=0;i<var;i++) {       //the dot product of the first derivative at the new and old points is calculated and stored in temp1
    temp1=temp1+(ydk[i]*yd[i]);
  }
  if (temp1<0) {              //temp1 is made positive
    temp1=-temp1;
  }
  fprintf(fpo,"%Lf\t",temp1);
  if (temp1<=e2) {            //terminate if the dot product is less than the second termination parameter
    for (i=0;i<var;i++) {
      ret[i]=ihr[i];
    }
    return ret;
  }

  //step 5
  l=0;
  x0m=0;
  for (i=0;i<var;i++) {         //The modulus of the difference between the old and new points is calculated and stored in l. The modulus od the old point is stored in x0m.
    l=l+pow(ihr[i]-x0[i],2);
    x0m=x0m+pow(x0[i],2);
  }
  l=sqrt(l);
  x0m=sqrt(x0m);
  fprintf(fpo,"%Lf\n",l/x0m);
  if ((l/x0m)<e1) {            //terminate if the ratio of the two moduli is less than the first termination parameter
    for (i=0;i<var;i++) {
      ret[i] = ihr[i];
    }
    return ret;
  }
  else {                      //if not terminated, go to next iteration
    k++;
    for (i=0;i<var;i++) {
      x0[i]=ihr[i];
    }
    goto label3;
  }
}

int main() {
  FILE * fpi;                    //the input and output files are opened
  fpi=fopen("input.txt","r");
  fpo=fopen("output.txt","w");
  int q;                         //the question, the number of variables and the initial point are read from the input file
  int var;
  int i=0;
  fscanf(fpi,"%d",&q);
  fscanf(fpi,"%d",&var);
  static long double * x0;
  x0=(long double *) malloc(var*sizeof(long double));
  for (i=0;i<var;i++) {
    fscanf(fpi,"%Lf",&x0[i]);
  }
  long double * y;
  y=(long double *) malloc(var*sizeof(long double));
  y=nwtnm(q,var,x0);             //the function for Newton's method is called and the returned array is stored in y
  printf("The minima is at ");        //output the minuma and the minimum value to the screen and the output file
  fprintf(fpo,"\n\nThe minima is at ");
  for (i=0;i<var-1;i++) {
    printf("%Lf, ",y[i]);
    fprintf(fpo,"%Lf, ",y[i]);
  }
  printf("%Lf.\n",y[var-1]);
  printf("The minimum value of the function is %Lf.\n",objf(y,var,q));
  fprintf(fpo,"%Lf.\n",y[var-1]);
  fprintf(fpo,"The minimum value of the function is %Lf.\n",objf(y,var,q));
  return 0;
}

