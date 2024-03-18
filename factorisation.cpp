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
        Matrix D,L;

        int n;
        
    public:
        Factorisation();
        ~Factorisation();
        void displayMatrix(Matrix M);
        void displayVector(Vector V);
        void initializeData();
        Matrix get_matrix();
        Vector get_vector();
        Vector get_solution();
        Vector solution;
        // Decomposition LDL(transpose)
        void decomposition();
        void diagonal();
        // Resolution par Gauss-Seidel
        void lower_resolution();
        void upper_resolution();
};

Factorisation::Factorisation(){
    this->initializeData();
}

Factorisation::~Factorisation(){
    cout<<endl;
}

Matrix Factorisation::get_matrix(){
    return this->A;
}

Vector Factorisation::get_vector(){
    return this->b;
}
void Factorisation::displayMatrix(Matrix M){
    for(int i=0;i<int(M.size());i++){
        for(int j=0;j<int(M.size());j++){
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
        file>>this->n;

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
        // initialiser la variable Solution en 0
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
    float somme;

    for(int i=0;i<int(A.size());i++){
        for(int j=0;j<int(A.size());j++){
            for(int k=0;k<j;k++){
                somme +=  this->A[i][k] * this->A[k][k] * this->A[j][k];
            }
            this->A[i][j] = 1.0/(this->A[j][j]) * (this->A[i][j] - somme);
            somme = 0;
        }
        for(int k=0;k<i;k++){
            somme += this->A[k][k] * pow(this->A[i][k],2);
        }
        this->A[i][i] -= somme;
        somme = 0;
    }
}

void Factorisation::lower_resolution(){
    float somme = 0;
    for(int i=0;i<this->n;i++){
        somme = 0;
        for(int j=0;j<i;j++){
            somme += this->A[i][j]* this->solution[j];
        }
        this->solution[i] = (this->b[i] - somme);
    }
}

void Factorisation::diagonal(){
    for (int i = 0; i < this->n; i++){
        this->solution[i] /= this->A[i][i];
    }
}

void Factorisation::upper_resolution(){
    float somme = 0;
    for (int i = this->n - 1; i >= 0; i--){
        somme = 0;
        for (int j = i + 1; j < this->n; j++){
            somme += this->A[j][i] * this->solution[j];
        }
        this->solution[i] = ((this->solution[i] - somme));
    }
}


int main(){
    Factorisation f;
    cout<<"Voici la matrice donnee :"<<endl;
    f.displayMatrix(f.get_matrix());
    cout<<"Second membre :"<<endl;
    f.displayVector(f.get_vector());
    f.decomposition();
    cout<<"On obtient alors LDL^t:\n";
    f.displayMatrix(f.get_matrix());
    f.lower_resolution();
    f.diagonal();
    f.upper_resolution();
    cout<<"Solutions :"<<endl;
    f.displayVector(f.solution);
    return(0);
}