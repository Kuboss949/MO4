#include <iostream>
#include <cmath>

using namespace std;


#define ERR0 1.0e-8
#define ERR1 1.0e-8
#define ERR2 1.0e-8
#define FUNERR0 1.0e-8
#define FUNERR1 1.0e-8
#define FUNERR2 1.0e-8


void jacob(double ** A, double * xn);
void function(double * fxn, double * xn);
const double eps = 1e-12;


void gauss(int n, double **AB, double *X) {

    for(int i=0; i<n; i++) {
        int maxRow = i;
        for(int j=i+1; j<n; j++) {
            if(fabs(AB[j][i]) > fabs(AB[maxRow][i])) {
                maxRow = j;
            }
        }

        if(maxRow != i) {
            for(int k=0; k<n+1; k++) {
                double temp = AB[i][k];
                AB[i][k] = AB[maxRow][k];
                AB[maxRow][k] = temp;
            }
        }
        for(int j=i+1; j<n; j++) {
            double factor = AB[j][i] / AB[i][i];
            for(int k=i+1; k<n+1; k++) {
                AB[j][k] -= factor * AB[i][k];
            }
            AB[j][i] = 0.0;
        }
    }

    for(int i=n-1; i>=0; i--) {
        X[i] = AB[i][n] / AB[i][i];
        for(int j=i-1; j>=0; j--) {
            AB[j][n] -= AB[j][i] * X[i];
            AB[j][i] = 0.0;
        }
    }
}

bool sprawdzWynik(double *x){
    if(fabs(x[0]) <= ERR0 && fabs(x[1]) <= ERR1 && fabs(x[1])<=ERR2){
        return true;
    }
    return false;
}

bool sprawdzResiduum(double *x){
    if(fabs(x[0]) <= FUNERR0 && fabs(x[1]) <= FUNERR1 && fabs(x[1]) <= FUNERR2){
        return true;
    }
    return false;
}

double norma(double *x){
    double result = 0;
    for (int i = 0; i < 3; i++) {
        result+=x[i]*x[i];
    }
    return sqrt(result);
}

void uogNewton(double * initialGuess, int maxIterations){

    double ** J = new double*[3];
    double * delta = new double[3];
    double * xn1 = new double[3];
    double * fxn = new double[3];
    double * estymator = new double[3];
    for(int i=0; i<3; i++){
        J[i]=new double [4];
    }


    for(int i=0; i<maxIterations; i++){

        jacob(J, initialGuess);
        gauss(3, J, delta);
        function(fxn, initialGuess);
        for(int j=0; j<3; j++){
            xn1[j]=initialGuess[j]-delta[j];
            estymator[j]=xn1[j]-initialGuess[j];
            initialGuess[j]=xn1[j];
        }
        cout << "Iteracja nr " << i << ". Przyblizenie: [" << initialGuess[0] << ", "<< initialGuess[1] << ", "<< initialGuess[2] << "]. Estymator: " << norma(delta) << " Residuum: " << norma(fxn) << endl;
        if(sprawdzWynik(delta) && sprawdzResiduum(fxn)){
            cout << "Przerwane ze wzgledu na kryterium dokladnosci wyznaczenia xn oraz na kryterium wiarygodnosci xn jako przyblizenia pierwiastka" << endl;
            break;
        }

        if(i+1==maxIterations){
            cout << "Zakończono z powodu wyczerpania liczby iteracji";
            return;
        }

    }
    cout << "Rozwiazanie: [" << initialGuess[0] << ", "<< initialGuess[1] << ", "<< initialGuess[2] << "]" << endl;


}
void function(double * fxn, double * xn){
    fxn[0]=xn[0]*xn[0]+xn[1]*xn[1]+xn[2]*xn[2]-2;
    fxn[1]=xn[0]*xn[0]+xn[1]*xn[1]-1;
    fxn[2]=xn[0]*xn[0]-xn[1];

}


void jacob(double ** A, double * xn){
    A[0][0]=2*xn[0];
    A[0][1]=2*xn[1];
    A[0][2]=2*xn[2];
    A[0][3]=xn[0]*xn[0]+xn[1]*xn[1]+xn[2]*xn[2]-2;

    A[1][0]=2*xn[0];
    A[1][1]=2*xn[1];
    A[1][2]=0;
    A[1][3]=xn[0]*xn[0]+xn[1]*xn[1]-1;

    A[2][0]=2*xn[0];
    A[2][1]=-1;
    A[2][2]=0;
    A[2][3]=xn[0]*xn[0]-xn[1];

}



int main() {

    double * xn = new double[3];
    xn[0]=10.0;
    xn[1]=5.0;
    xn[2]=-3.0;
    uogNewton(xn, 1000);


    return 0;
}