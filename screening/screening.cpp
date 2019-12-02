#include <Rcpparmadillo.h>

#include <cmath>
#include <stdio.h>

// Compute gaussian kernel
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat KernelG(arma::mat x,float sigma) {
  
  int m=x.n_rows;
  
arma::mat fil_ig=arma::repmat(x.row(0)*trans(x.row(0)),1,m);
  
  for(int i = 1; i < m; i++){
    fil_ig=join_cols(fil_ig,arma::repmat(x.row(i)*trans(x.row(i)),1,m));
  }
  
  arma::mat col_ig=trans(fil_ig);
  arma::mat cross=x*trans(x);
  
  return arma::exp(-0.5*(fil_ig+col_ig-2*cross)/pow(sigma,2));
  
}
//************************************************************
// sigma1 sigma x sen censurado
// sigma2 sigma y censurado
// sigma3 sigma x censurado
// pesos1 x sen censurado
// pesos2 y censurado
// pesos3 xy censurados
// non se utiliza sigma3 neste programa

// Main program
// Returns the value of Hilbert-Schmidt criterion statistic
// Input: both vectors of observations (x and y), sigma1,2,3 resective kernel bandwidths
// pesos1: weights for random variable x (non-censored)
// pesos2: weights for random variable y (censored)
// pesos3: weights for joint-distributed variable xy   

// [[Rcpp::export]]
Rcpp::List HS1C(const arma::vec x, const arma::vec y,const double sigma1, 
          const double sigma2, const double sigma3,
            arma::vec pesos1, arma::vec pesos2,  arma::vec pesos3) {
  
long int i=0;
long  int j=0;
long  int r=0;
int n= x.n_elem;
  
  arma::mat Kx(n,n,arma::fill::zeros);
  arma::mat Ky(n,n,arma::fill::zeros);
  
// Call kernel subroutine
  
  Kx=KernelG(x,sigma1);
  Ky=KernelG(y,sigma2);

  double sumapesos1= arma::sum(pesos1);
  double sumapesos2= arma::sum(pesos2);
  double sumapesos3= arma::sum(pesos3);
  
  /*printf("%f\n",sumapesos1);
  printf("%f\n",sumapesos2);
  printf("%f\n",sumapesos3);*/
  
  pesos1= pesos1/sumapesos1;
  pesos2= pesos2/sumapesos2;
  pesos3= pesos3/sumapesos3;

 /* printf(pesos1)
  printf(pesos2)
  printf(pesos3)*/
  
  double term1=0;
  double term2=0;
  double term3=0;
  double suma1=0;
  double suma2=0;
  
  double aux=0;
  double aux1=0;
  double aux2=0;
  
  float contar9=0;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      contar9= contar9+1;
      aux1= Kx(i,j);
      aux2= Ky(i,j);
      aux= aux+ aux1*aux2*pesos3(i)*pesos3(j);
      suma1= suma1+aux1*pesos1(i)*pesos1(j);
      suma2= suma2+aux2*pesos2(i)*pesos2(j);
    }
  }
 
  term1= aux;
  
  double media1= suma1;
  double media2= suma2;
 
  //print(n)
  /* printf("%f\n",media1);
   printf("%f\n",media2);*/

  term2= media1*media2;
  
  double suma3=0;
  double auxnuevo;
  double auxnuevo2;
  double auxnuevo3;
  double productopesos;

  float  mult=0;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(r=0;r<n;r++){
       auxnuevo2= Kx(i,j);
       auxnuevo3= Ky(i,r);
       auxnuevo=auxnuevo2*auxnuevo3;
       productopesos= pesos1(j)*pesos3(i);
       productopesos= productopesos*pesos2(r);
        suma3= suma3+auxnuevo*productopesos;
      }
    }
  }

  term3= -suma3-suma3;
 
  double term= term1+term3+term2;
  Rcpp::List ret;
  //ret["pesos1"]=pesos1;
  //ret["pesos2"]=pesos2;
  //ret["pesos3"]=pesos3;
  //ret["media1"]=media1;
  //ret["media2"]=media2;
  //ret["Kx"]=Kx;
 //ret["Ky"]=Ky;
  //ret["x"]=x;
  //ret["contar7"]=mult;
  //ret["sigma1"]=sigma1;
  ret["term"]=term;
  //ret["term1"]=term1;
  //ret["term2"]=term2;
  //ret["term3"]=term3;
 // ret["contar"]=auxnuevo;
 // ret["contar2"]=auxnuevo2;
//  ret["contar3"]=auxnuevo3;
//  ret["contar4"]=contar;
//  ret["contar5"]=cruzados;
//  ret["contar6"]=contar9;
  
return(ret);
//return(term);
}






