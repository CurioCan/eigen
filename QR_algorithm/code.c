#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Complex number structure
typedef struct {
    double r;  // Real part
    double i;  // Imaginary part
} Complex;

// Operations to Handle complex values
Complex add(Complex a, Complex b) {
    return (Complex){a.r + b.r, a.i + b.i};
}

Complex sub(Complex a, Complex b) {
    return (Complex){a.r - b.r, a.i - b.i};
}

Complex mul(Complex a, Complex b) {
    return (Complex){a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r};
}

Complex c_div(Complex a, Complex b) {
    double denom = b.r * b.r + b.i * b.i;
    return (Complex){(a.r * b.r + a.i * b.i) / denom, (a.i * b.r - a.r * b.i) / denom};
}

Complex sqrtC(Complex a) {
    double mag = sqrt(a.r * a.r + a.i * a.i);
    return (Complex){sqrt((mag + a.r) / 2), (a.i >= 0 ? 1 : -1) * sqrt((mag - a.r) / 2)};
}

Complex c_conj(Complex a) {
    return (Complex){a.r, -a.i};
}

double absC(Complex a) {
    return sqrt(a.r * a.r + a.i * a.i);
}

// Matrixverse
typedef struct {
    int n;    // Dimension
    Complex **m;  // Matrix elements
} Matrix;

Matrix* createMatrix(int n) {
    Matrix *mat = (Matrix*)malloc(sizeof(Matrix));
    mat->n = n;
    mat->m = (Complex**)malloc(n * sizeof(Complex*));
    for (int i = 0; i < n; i++) {
        mat->m[i] = (Complex*)calloc(n, sizeof(Complex));
    }
    return mat;
}

// Free matrix memory
void freeMatrix(Matrix *mat) {
    for (int i = 0; i < mat->n; i++) {
        free(mat->m[i]);
    }
    free(mat->m);
    free(mat);
}

// Multiply two matrices
Matrix* multiply(int n, Matrix *A, Matrix *B) {
    Matrix *C = createMatrix(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Complex sum = {0, 0};
            for (int k = 0; k < n; k++) {
                sum = add(sum, mul(A->m[i][k], B->m[k][j]));
            }
            C->m[i][j] = sum;
        }
    }
    return C;
}

// Gram-Schmidt orthogonalization
void gramSchmidt(int n, Matrix *A, Matrix *Q, Matrix *R) {
    for (int i = 0; i < n; i++) {
        Q->m[i][i] = A->m[i][i];  // Q's columns from A
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            Complex dot = {0, 0};
            for (int k = 0; k < n; k++) {
                dot = add(dot, mul(Q->m[k][i], c_conj(Q->m[k][j])));
            }
            Complex scale = mul(dot, Q->m[i][j]);
            for (int k = 0; k < n; k++) {
                Q->m[k][i] = sub(Q->m[k][i], scale);
            }
        }
        Complex norm = {0, 0};
        for (int k = 0; k < n; k++) {
            norm = add(norm, mul(Q->m[k][i], c_conj(Q->m[k][i])));
        }
        norm = sqrtC(norm);
        for (int k = 0; k < n; k++) {
            Q->m[k][i] = c_div(Q->m[k][i], norm);
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            Complex dot = {0, 0};
            for (int k = 0; k < n; k++) {
                dot = add(dot, mul(Q->m[k][i], A->m[k][j]));
            }
            R->m[i][j] = dot;
        }
    }
}

// QR algorithm 
Matrix* qrAlgorithm(int n, Matrix *A) {
    Matrix *Q = createMatrix(n);
    Matrix *R = createMatrix(n);
    int maxIterations = 1000;

    for (int iter = 0; iter < maxIterations; iter++) {
        gramSchmidt(n, A, Q, R);
        Matrix *newA = multiply(n, R, Q);
        double offDiagonalSum = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    offDiagonalSum += absC(newA->m[i][j]);
                }
            }
        }

        
        if (offDiagonalSum < 1e-12) { 
            break;
        }

        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A->m[i][j] = newA->m[i][j];
            }
        }

        freeMatrix(newA);
    }
    freeMatrix(Q);
    freeMatrix(R);
    return A;
}


Matrix* readMatrixFromUser() {
    int n;
    printf("Enter the dimension of matrix: ");
    scanf("%d", &n);
    Matrix *mat = createMatrix(n);

    printf("Enter the elements of the matrix ( real+imagi for complex numbers):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("Element [%d][%d]: ", i + 1, j + 1);
            scanf("%lf+%lfi", &mat->m[i][j].r, &mat->m[i][j].i);
        }
    }

    return mat;
}

// Extract eigenvalues for asymmetric matrices
void extract_eigenvalues_asymmetric(Matrix *A, Complex *eigenvalues) {
    for (int i = 0; i < A->n; i++) {
        eigenvalues[i] = A->m[i][i];
    }
}

// Print eigenvalues (real part)
void printEigenvalues(Matrix *mat) {
    for (int i = 0; i < mat->n; i++) {
        printf("Eigenvalue %d: %.2lf + %.2lf i\n", i + 1, mat->m[i][i].r, mat->m[i][i].i);
    }
}

int main() {
    
    Matrix *A = readMatrixFromUser();
    if (A == NULL) return 1;

    
    Matrix *result = qrAlgorithm(A->n, A);

    
    Complex eigenvalues[A->n];
    extract_eigenvalues_asymmetric(result, eigenvalues);
    printEigenvalues(result);

    
    freeMatrix(A);
    return 0;
}

