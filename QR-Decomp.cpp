#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

int m, n;
vector <double> weights;



// function to print matrix content at any stage
void print(vector<vector<double>> matrix) {
    cout << "The matrix now: " << endl;
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix[0].size(); j++)
            cout << setw(7) << setprecision(4) << matrix[i][j] << " ";
        cout << endl;
    }
}

// function to get the transpose
vector<vector<double>> transpose(vector<vector<double>> a) {
    int i, j;
    vector <vector<double>> B;
    B.resize(a[0].size());
    for (int i = 0; i < a[0].size(); i++)
        B[i].resize(a.size());
    for (i = 0; i < a[0].size(); i++)
        for (j = 0; j < a.size(); j++)
            B[i][j] = a[j][i];
    return B;
}

// RREF to check independence
bool RREF(vector<vector<double>> A) {
    if (n > m)
        return false;
    
    const int nrows = A.size(); // number of rows
    const int ncols = A[0].size(); // number of columns
    
    int lead = 0;
    
    while (lead < nrows) {
        double d, m;
        
        for (int r = 0; r < nrows; r++) { // for each row ...
            /* calculate divisor and multiplier */
            d = A[lead][lead];
            m = A[r][lead] / A[lead][lead];
            
            for (int c = 0; c < ncols; c++) { // for each column ...
                if (r == lead)
                    A[r][c] /= d; // make pivot = 1
                else
                    A[r][c] -= A[lead][c] * m;  // make other = 0
                
                if (A[r][c] == -0) A[r][c] = 0;
            }
        }
        
        lead++;
        
    }
    bool Indep = true;
    for (int i = 0; i < n; i++) {
        int flag = 0;
        for (int j = 0; j < m; j++)
            if (A[i][j] != 0)
                flag = -1;
        if (flag == 0) {
            Indep = false;
            break;
        }
    }
    
    return Indep;
    
}

// function to compute the inner product
double product(vector<double>a, vector<double> b) {
    double result = 0;
    for (int i = 0; i < weights.size(); i++)
        result += a[i] * b[i] * weights[i];
    return result;
}

// function to subtract vectors
vector<double> subtract(vector<double>v1, vector<double>v2) {
    vector<double>result;
    for (int i = 0; i < v2.size(); i++) {
        result.push_back(v1[i] - v2[i]);
    }
    return result;
}

// function to add vectors
vector<double> add(vector<double>v1, vector<double>v2) {
    vector<double>result;
    for (int i = 0; i < v2.size(); i++) {
        result.push_back(v1[i] + v2[i]);
    }
    return result;
}

// function to compute the projection (// proj q2 v1)
vector <double> projection(vector<double> first, vector<double> second) {
    vector<double> result;
    
    double num = product(first, second) / product(first, first);
    for (int i = 0; i < first.size(); i++) {
        result.push_back(num * first[i]);
    }
    return result;
}

// functionn to compute the orthogonal vector-columns without normalizing it
vector<vector<double>> calculate_Q(vector<vector<double>> matrix) {
    vector <vector<double>>q;
    
    vector <vector<double>>result_matrix;
    
    for (int i = 0; i < n; i++) {
        vector <double>temp1(m, 0);
        if (i == 0) {
            result_matrix.push_back(matrix[0]);
            q.push_back(matrix[0]);
        }
        else {
            
            for (int c = 0; c < i; c++) {//3
                temp1 = subtract(temp1, projection(q[c], matrix[i]));
                
            }
            temp1 = add(matrix[i], temp1);
            result_matrix.push_back(temp1);
            q.push_back(temp1);
        }
        
        
    }
    return result_matrix;
    
}
//function to normalize
vector <vector<double>> normalize(vector<vector<double>> matrix) {
    vector<vector<double>> orthonormal(matrix.size());
    for (int i = 0; i < matrix.size(); i++)
    {
        vector <double> row(matrix[0].size(), 0); //3
        double length = sqrt(product(matrix[i], matrix[i]));
        
        for (int j = 0; j < row.size(); j++)
        {
            row[j] = matrix[i][j] / length;
            orthonormal[i] = row;
        }
        
    }
    return transpose(orthonormal);
}
// function to calculate the R vector
vector <vector<double>> R_calculator(vector<vector<double>> matrix, vector<vector<double>> Q) {
    vector <vector<double>> R(n);
    for (int i = 0; i < n; i++)
        R[i].resize(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            if (i < j)
                R[i][j] = 0;
            else
                R[i][j] = product(matrix[i], Q[j]);
        }
    return R;
}

int main() {
    // Enter the matrix
    cout << "Please enter the number of rows (m) & columns (n) of your matrix" << endl;
    cin >> m >> n;
    vector<vector<double>> matrix(m);
    cout << "Enter your matrix cells row by row" << endl;
    for (int i = 0; i < m; i++)
    {
        matrix[i].resize(n);
        for (int j = 0; j < n; j++)
            cin >> matrix[i][j];
    }
    
    
    // Check if it's linearly independent
    /* if (!RREF(matrix))
     {
     cout << "There is no QR-decomosition for this matrix, since its columns are not linearly independent" << endl;
     cout << "Program Terminated ..." << endl;
     return 0;
     }
     else
     {*/
    // Enter the weights of the inner product
    cout << "The inner product will be in the form: ((u1, u2, .. , un), (v1, v2, .. , vn)) = a1u1v1 + a2u2v2 + .... + anvnun" << endl;
    cout << "Please enter the weights: a1, a2, ... , an" << endl;
    weights.resize(m);
    for (int i = 0; i < m; i++)
        cin >> weights[i];
    
    cout << endl;
    // Apply Gram-Schmidt
    vector <vector<double>> trans = transpose(matrix);
    vector <vector<double>> Q = calculate_Q(trans);
    Q = transpose(Q);
    cout << "Vector Q: " << endl;
    print(Q);
    cout << endl;
    
    //Normalize
    cout << "Q after Normalization ";
    Q = normalize(transpose(Q));
    print(Q);
    cout << endl;
    
    // Calculating R
    vector <vector<double>> R = R_calculator(trans, transpose(Q));
    R = transpose(R);
    cout << "Vector R: " << endl;
    print(R);
    //}
    
}



