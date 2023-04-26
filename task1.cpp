#include <iostream>

#include <vector>
#include<iomanip>

template<typename T>
class SquareMatrix;

template <typename T>
class Matrix {
    private:
    std::vector<std::vector<T>> matrix;
    int n, m;
    public:
    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        for (int i = 0; i < n; ++i) {
            matrix.push_back(std::vector<T>{});
            for (int j = 0; j < m; ++j) {
                matrix[i].push_back(0);
            }
        }
    }

    Matrix operator=(Matrix B) {
        int B_n, B_m;
        B_n = B.size_n();
        B_m = B.size_m();
        this->n = B_n;
        this->m = B_m;
        matrix.clear();
        for (int i = 0; i < B_n; ++i) {
            matrix.push_back(B[i]);
        }
        return *this;
    }
    std::vector<T>& operator[](int i) {
        return matrix[i];
    }
    
    int size_n() {
        return n;
    }

    int size_m() {
        return m;
    }

    Matrix operator+(Matrix &B) {
        int size_B[2] = {B.size_n(), B.size_m()};
        if (size_B[0] == n && size_B[1] == m) {
            Matrix result{n, m};
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    result[i][j] = B[i][j] + matrix[i][j];
                }
            }
            
            return result;
        } else {
            
            return Matrix{0, 0};
        }
    }

    Matrix operator-(Matrix& B) {
        int size_B[2] = {B.size_n(), B.size_m()};
        if (size_B[0] == n && size_B[1] == m) {
            Matrix result{n, m};
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    result[i][j] = matrix[i][j] - B[i][j];
                }
            }
            
            return result;
        } else {
            
            return Matrix{0, 0};
        }
    }

    Matrix operator*(Matrix& B) {
        int size_B[2] = {B.size_n(), B.size_m()};
        if (m == size_B[0]) {
            Matrix result{n, size_B[1]};
            Matrix trans = B.transpose();
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < size_B[1]; ++j) {
                    T mult = 0;
                    std::vector<T> v1 = matrix[i];
                    std::vector<T> v2 = trans[j];
                    for (int k = 0; k < m; ++k) {
                        mult += v1[k] * v2[k];
                    }
                    result[i][j] = mult;
                }
            }
            
            return result;
        } else {
            
            return Matrix{0, 0};
        }
    }

    Matrix transpose() {
        Matrix result{m, n};
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result[j][i] = matrix[i][j];
            }
        }
        
        return result;
    }

    friend void operator>>(std::istream& input, Matrix<T>& matrix) {
        for (int i = 0; i < matrix.size_n(); ++i) {
            for (int j = 0; j < matrix.size_m(); ++j) {
                T element;
                input >> element;
                matrix[i][j] = element;
            }
        }
    }

    friend void operator<<(std::ostream& output, Matrix<T>& matrix) {
        for (int i = 0; i < matrix.size_n(); ++i) {
            for (int j = 0; j < matrix.size_m(); ++j) {
                if (std::abs(matrix[i][j]) < 0.00005) {
                    output << "0.0000 ";
                } else {
                    output << matrix[i][j] << " ";
                }
            }
            output << std::endl;
        }
    }
    
    void multiplyRow(int row, double scalar) {
        for (int i = 0; i < m; ++ i) {
            matrix[row][i] *= scalar;
        }
    }

    SquareMatrix<T> operatorSquareMatrix() {
        SquareMatrix<T> result{std::min(n, m)};
        for (int i = 0; i < std::min(n, m); ++i) {
            for (int j = 0; j < std::min(n, m); ++j) {
                result[i][j] = (*this)[i][j];
            }
        }
    }

    
};

template<typename T>
class PermutationMatrix;

template<typename T>
class ElemenationMatrix;

template<typename T>
class IdentityMatrix;

template<typename T>
class SquareMatrix : public Matrix<T> {
    public:
    SquareMatrix(int n) : Matrix<T>{n, n} {
    }

    int size() {return Matrix<T>::size_n();}

    SquareMatrix inverse() {
        SquareMatrix<T> inverse{this->size()};
        inverse = IdentityMatrix<T>{this->size()};
        SquareMatrix<T> matrix = *this;
        for (int i = 0; i < matrix.size(); ++ i) { 
            int next = i; 
            for (int j = i; j < matrix.size(); ++j) {
                if (std::abs(matrix[j][i]) > std::abs(matrix[next][i])) {
                    next = j;
                }
            }
            if (next != i) {
                PermutationMatrix<T> P{next+1, i+1, matrix.size_n()};
                matrix = P*matrix;
                inverse = P*inverse;
                
            }
            for (int j = i+1; j < matrix.size(); ++ j) {
                if (std::abs(matrix[j][i]) <= 10E-10) {
                    continue;
                }
                ElemenationMatrix<T> E{j+1, i+1, matrix};
                matrix = E*matrix;
                inverse = E*inverse;
                
            }
        }

        for (int i = matrix.size()-1; i >= 0; -- i) {

            for (int j = i-1; j >= 0; -- j) {
                ElemenationMatrix<T> E{j+1, i+1, matrix};
                matrix = E*matrix;
                inverse = E*inverse;
                
            }
        }
        
        for(int i = 0; i < matrix.size(); i++) {
            inverse.multiplyRow(i, 1/matrix[i][i]);
            matrix.multiplyRow(i, 1/matrix[i][i]);
        }
        return inverse;
    }

    SquareMatrix<T> operator=(Matrix<T> B) {
        int B_n, B_m;
        B_n = B.size_n();
        B_m = B.size_m();
        *this = SquareMatrix<T>(std::min(B_n, B_m));

        for (int i = 0; i < B_n; ++i) {
            for (int j = 0; j < B_n; ++j) {
                (*this)[i][j] = B[i][j];
            }
        }
        return *this;
    }
    
};

template<typename T>
class IdentityMatrix : public SquareMatrix<T> {
    public:
    IdentityMatrix(int n) : SquareMatrix<T>{n} {
        for (int i = 0; i < n; ++i) {
            (*this)[i][i] = 1;
        }
    }

};


template<typename T>
class ElemenationMatrix : public IdentityMatrix<T> {
    public:
    ElemenationMatrix(int i, int j, Matrix<T>& A) : IdentityMatrix<T>{A.size_n()} {
        (*this)[i-1][j-1] = -1*A[i-1][j-1] / A[j-1][j-1];
    }
};

template<typename T>
class PermutationMatrix : public IdentityMatrix<T> {
    public:
    PermutationMatrix(int i, int j, int size) : IdentityMatrix<T>{size} {
        auto temp = (*this)[i-1];
        (*this)[i-1] = (*this)[j-1];
        (*this)[j-1] = temp;
    }
};

template<typename T>
class ColumnVector : public Matrix<T> {
    public:

    
    ColumnVector(int size) : Matrix<T>{size, 1} {
    }

    
    int size() const { return Matrix<T>::size_n; }

    T& operator[](int i) {
        return Matrix<T>::operator[](i)[0];
    }

    
    T operator*(ColumnVector& other) {
        if (size() == other.size()) {
            int result;
            for (int i = 0; i < size(); ++ i) {
                result += other[i] * Matrix<T>::at(i)[0];
            }
            return result;
        }
    }

    ColumnVector<T> operator=(Matrix<T> B) {
        int B_n;
        B_n = B.size_n();
        *this = ColumnVector<T>{B_n};

        for (int i = 0; i < B_n; ++i) {
            (*this)[i] = B[i][0];
        }
        return *this;
    }
    
};




int main() {
    int m;
    std::cin >> m;

    ColumnVector<double> b{m};

    double t[m];
    for (int i = 0; i < m; i++) {
        std::cin >> t[i] >> b[i];
    }

    int degree;
    std::cin >> degree;
    Matrix<double> A{m, degree+1};
    for (int i = 0; i < m; ++i) {
        double temp = 1;
        for (int j = 0; j <= degree; j++) {
            A[i][j] = temp;
            temp *= t[i];
        }
    }

    std::cout << std::setprecision(4) << std::fixed;

    std::cout << "A:\n" << A;

    SquareMatrix<double> As{degree+1};
    As = A.transpose()*A;

    std::cout << "A_T*A:\n" << As;

    As = As.inverse();
    std::cout << "(A_T*A)^-1:\n" << As;

    b = A.transpose()*b;
    std::cout << "A_T*b:\n" << b;

    b = As*b;
    std::cout << "x~:\n" << b;
    return 0;
}