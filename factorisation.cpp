#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

using namespace std;

typedef vector<float> Vector;
typedef vector<Vector> Matrix;

class Factorisation{
    private:
        Matrix A;
        Vector b; //second membre
        Vector D; //pour stocker les éléments diagonaux de la matrice L
        int n;
        
    public:
        Factorisation();
        void displayMatrix(Matrix M);
        void displayVector(Vector V);
        void initializeData();
        Matrix get_matrix();
        Vector get_vector();
        Vector get_solution();
        Vector get_diag();
        Vector solution;
        void decomposition();
        void lower_resolution();
        void upper_resolution();
};

Factorisation::Factorisation(){
    this->initializeData();
}

Matrix Factorisation::get_matrix(){
    return this->A;
}

Vector Factorisation::get_vector(){
    return this->b;
}

Vector Factorisation::get_diag(){
    return this->D;
}

void Factorisation::displayMatrix(Matrix M){
    for(int i=0;i<int(M.size());i++){
        for(int j=0;j<int(M[i].size());j++){
            cout<<M[i][j]<<"    ";
        }
        cout<<endl;
    }
    cout<<endl;
}

void Factorisation::displayVector(Vector V){
    for(int j=0;j<int(V.size());j++){
        cout<<V[j]<<endl;
    }
}

void Factorisation::initializeData(){
    ifstream file("data.txt");
    if(file){
        file >> this->n;
        for(int i=0;i<this->n;i++){
            Vector temp;
            for(int j=0;j<this->n;j++){
                float element;
                file>>element;
                temp.push_back(element);
            }
            A.push_back(temp);
        }
        for(int j=0;j<this->n;j++){
            float element;
            file>>element;
            this->b.push_back(element);
        }
        D.resize(n);
        for(int i=0;i<this->n;i++){
            this->solution.push_back(0.0);
        }
        file.close();
    }
    else{
        cout<<"Erreur d'ouverture du fichier"<<endl;
    }
}

void Factorisation::decomposition(){
    for(int i=0;i<int(A.size());i++){
        float sum = A[i][i];
        for(int k=0;k<i;k++){
            sum -= A[i][k]*A[i][k]*D[k];
        }
        if(sum <= 0){
            cout << "La matrice n est pas definie positive." << endl;
            break;
        }
        D[i] = sum;
        for(int j=i+1;j<int(A.size());j++){
            sum = A[j][i];
            for(int k=0;k<i;k++){
                sum -= A[j][k]*A[i][k]*D[k];
            }
            A[j][i] = sum / D[i];
        }
    }
}

void Factorisation::lower_resolution(){
    for(int i=0;i<this->n;i++){
        float sum = 0;
        for(int j=0;j<i;j++){
            sum += this->A[i][j] * this->solution[j];
        }
        this->solution[i] = (this->b[i] - sum);
    }
}

void Factorisation::upper_resolution(){
    for(int i=this->n-1;i>=0;i--){
        float sum = 0;
        for(int j=i+1;j<this->n;j++){
            sum += this->A[j][i] * this->solution[j];
        }
        this->solution[i] = round((this->solution[i] - sum) / D[i]);
    }
}

int main(){
    Factorisation f;
    cout<<"Voici la matrice donnee :"<<endl;
    f.displayMatrix(f.get_matrix());
    cout<<"Second membre :"<<endl;
    f.displayVector(f.get_vector());
    f.decomposition();
    f.lower_resolution();
    f.upper_resolution();
    cout<<"\nLa solution :"<<endl;
    f.displayVector(f.solution);
    
    return(0);
}
