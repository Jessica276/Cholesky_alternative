#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

using namespace std;

typedef vector<float> Vector;
typedef vector<Vector> Matrix;

class Profil{
    private:
        Matrix A;
        Vector b;       //second membre
        Vector D;       //pour stocker les éléments diagonaux de la matrice L
        Vector APi;     //Valeur de la matrice
        Vector li;      //nombre de valeur avant la valeur du diagonal
        Vector pi;      //i - li
        Vector nDiag;   //numéro du diagonal
        int n;

    public:
        Profil();
        void displayMatrix(Matrix M);
        void displayVector(Vector V);
        void initializeData();
        Matrix get_matrix();
        Vector get_vector();
        Vector get_APi();
        Vector get_solution();
        Vector get_diag();
        Vector get_num_diag();
        Vector get_li();
        Vector get_pi();
        Vector solution;
        void diagonal();
        void lower_resolution();
        void upper_resolution();
        void stockage();
};

Profil::Profil(){
    this->initializeData();
}

Matrix Profil::get_matrix(){
    return this->A;
}

Vector Profil::get_vector(){
    return this->b;
}

Vector Profil::get_diag(){
    return this->D;
}

Vector Profil::get_APi(){
    return this->APi;
}

Vector Profil::get_num_diag(){
    return this->nDiag;
}

Vector Profil::get_li(){
    return this->li;
}

Vector Profil::get_pi(){
    return this->pi;
}

void Profil::displayMatrix(Matrix M){
    for(int i=0;i<int(M.size());i++){
        for(int j=0;j<int(M[i].size());j++){
            cout<<M[i][j]<<"          ";
        }
        cout<<endl;
    }
    cout<<endl;
}

void Profil::displayVector(Vector V){
    for(int j=0;j<int(V.size());j++){
        cout<<V[j]<<"  ";
    }
}

void Profil::initializeData(){
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

void Profil::stockage(){
    int count = 1;
    for(int i=0;i<this->n;i++){
        int indice = 0;
        bool passage = false;
        for(int j=0;j<i+1;j++){ 
            auto element = this->A[i][j]; 
            //Numerotation du diagonal 
            if(i==j){
                this->nDiag.push_back(count);
            }

            if(element != 0 || A[i][i] == 0){
                indice = j;
                passage = true;
                this->APi.push_back(element);
                count++;
            }
            else{
                if(j>indice && passage){
                    this->APi.push_back(element);
                    count++;
                }
            }
        }

        //Calcul de l[i] = nDiag[i] - nDiag[i-1] - 1
        this->li.push_back(this->nDiag[i] - this->nDiag[i-1] - 1);

        //Calcul de p[i] = i - l[i]
        this->pi.push_back(i + 1 - this->li[i]);
    }
    cout<<"dimension :"<<nDiag.size();
}

void Profil::diagonal()
{
    float sum = 0;
    for(int i=0; i<int(this->nDiag.size());i++){
        for(int j=pi[i]; j<i;j++){
            for(int k=pi[j];k<j;k++){
                if(k>= pi[i])
                    sum += this->APi.at(nDiag[i] - i + k) * this->APi.at(nDiag[k]) * this->APi.at(nDiag[j] - j + k);
            }
            APi[nDiag[i] - i + j] = 1.0 / this->APi.at(nDiag[j]) * (this->APi.at(nDiag[i] - i + j) - sum);
            sum = 0;
        }
        for(int k=pi[i];k<i;k++){
            sum += APi[nDiag[k]] * APi[nDiag[i] - i + k] * APi[nDiag[i] - i + k];
        }
        APi[nDiag[i]] -= sum;
        sum = 0;
    }
}

//Matrice triangulaire inférieur

void Profil::lower_resolution(){
    float sum = 0;
    for(int i=0;i<int(this->nDiag.size());i++){
        sum = 0;
        for(int j=this->pi.at(i);j<i;j++)
        {
            sum += this->APi[nDiag[i] - i + j] * this->solution[j];
        }
        solution[i] = (this->b[i] - sum);
    }
}

//Matrice triangulaire supérieur

void Profil::upper_resolution(){
    float sum = 0;
    for(int i=0;i<int(this->nDiag.size());i++){
        solution[i] /= this->APi[nDiag[i]];
    }
    for(int j=nDiag.size()-1;j>=0;j--){
        sum = 0;
        for(int i=j+1;i<int(this->nDiag.size());i++){
            if(pi[i]<=j){
                sum += APi[nDiag[i] - i + j] * solution[i];
            }
        }
        this->solution[j] = ((this->solution[j] - sum));
    }
}

int main(){
    Profil p;
    cout<<"Voici la matrice donnee :"<<endl;
    p.displayMatrix(p.get_matrix());
    p.stockage();
    cout<<"\nAPi : \n";
    p.displayVector(p.get_APi());
    cout<<endl;
    cout<<"nDiag :\n";
    p.displayVector(p.get_num_diag());
    cout<<endl;
    cout<<"li :\n";
    p.displayVector(p.get_li());
    cout<<endl;
    cout<<"pi :\n";
    p.displayVector(p.get_pi());
    cout<<"\n-----------------------\n";
    p.diagonal();
    p.lower_resolution();
    p.upper_resolution();
    cout<<"La solution est :\n";
    p.displayVector(p.solution);
    cout<<endl;
    


    return(0);
}
