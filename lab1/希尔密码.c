#include <stdio.h>
#include <stdlib.h>
// #include <math.h>
#include <string.h>
int pow(int x,int y){
    int res = 1;
    while(y){
        if(y&1) res = res*x;
        x = x*x;
        y >>= 1;
    }
    return res;
}
// 计算模逆
int modInverse(int a, int mod) {
    a = a % mod;
    for (int x = 1; x < mod; x++) {
        if ((a * x) % mod == 1) {
            return x; // Found the modular inverse
        }
    }
    fprintf(stderr, "Modular inverse doesn't exist\n");
    exit(EXIT_FAILURE);
}

// 计算模
int mod(int a) {
    return (a % 26 + 26) % 26;
}

// 计算行列式
int getDet(int **arr, int n) {
    if (n == 1) return arr[0][0];
    if (n == 2) return mod(arr[0][0] * arr[1][1] - arr[0][1] * arr[1][0]);

    int det = 0;
    for (int i = 0; i < n; i++) {
        int **subMatrix = malloc((n - 1) * sizeof(int *));
        for (int j = 0; j < n - 1; j++) {
            subMatrix[j] = malloc((n - 1) * sizeof(int));
        }
        for (int j = 1; j < n; j++) {
            int l = 0;
            for (int k = 0; k < n; k++) {
                if (k == i) continue;
                subMatrix[j - 1][l++] = arr[j][k];
            }
        }
        det = mod(det + arr[0][i] * pow(-1, i) * getDet(subMatrix, n - 1));
        for (int j = 0; j < n - 1; j++) {
            free(subMatrix[j]);
        }
        free(subMatrix);
    }
    return det;
}

// 计算伴随矩阵
int **getAdjoint(int **arr, int n) {
    int **adj = malloc(n * sizeof(int *));
    for (int i = 0; i < n; i++) {
        adj[i] = malloc(n * sizeof(int));
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int **subMatrix = malloc((n - 1) * sizeof(int *));
            for (int k = 0; k < n - 1; k++) {
                subMatrix[k] = malloc((n - 1) * sizeof(int));
            }
            int p = 0;
            for (int k = 0; k < n; k++) {
                if (k == i) continue;
                int q = 0;
                for (int l = 0; l < n; l++) {
                    if (l == j) continue;
                    subMatrix[p][q++] = arr[k][l];
                }
                p++;
            }
            adj[j][i] = mod(pow(-1, i + j) * getDet(subMatrix, n - 1));
            for (int k = 0; k < n - 1; k++) {
                free(subMatrix[k]);
            }
            free(subMatrix);
        }
    }
    return adj;
}

// 计算逆矩阵
int **inverse(int **arr, int n) {
    int det = getDet(arr, n);
    int invdet = modInverse(det, 26);
    int **adj = getAdjoint(arr, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            adj[i][j] = mod(adj[i][j] * invdet);
        }
    }
    return adj;
}

// 矩阵乘法
int **multiply(int **A, int **B, int n, int m, int p) {
    // if (m != p) {
    //     fprintf(stderr, "Matrix dimensions do not match for multiplication.\n");
    //     exit(EXIT_FAILURE);
    // }
    int **C = malloc(n * sizeof(int *));
    for (int i = 0; i < n; i++) {
        C[i] = malloc(m * sizeof(int));
        for (int j = 0; j < m; j++) {
            C[i][j] = 0;
            for (int k = 0; k < p; k++) {
                C[i][j] = mod(C[i][j] + A[i][k] * B[k][j]);
            }
        }
    }
    return C;
}

// 打印矩阵
void printMatrix(int **matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main() {
    int m;
    char plain_text[100], cipher_text[100];
    scanf("%d %s %s", &m, plain_text, cipher_text);

    int n = strlen(plain_text);

    int **X = malloc((m + 1) * sizeof(int *));
    int **Y = malloc((m + 1) * sizeof(int *));
    for (int i = 0; i < m + 1; i++) {
        X[i] = malloc((m + 1) * sizeof(int));
        Y[i] = malloc(m * sizeof(int));
    }

    for (int i = 0; i < m * (m + 1); i++) {
        X[i / m][i % m] = plain_text[i] - 'A';
        Y[i / m][i % m] = cipher_text[i] - 'A';
    }
    for(int i = 0; i < m + 1; i++) {
        X[i][m] = 1;
    }
    int **X_1 = inverse(X, m + 1);

    int **K = multiply(X_1, Y, m + 1, m, m + 1);

    printMatrix(K, m + 1, m);

    for (int i = 0; i < m + 1; i++) {
        free(X[i]);
        free(Y[i]);
        free(X_1[i]);
        free(K[i]);
    }
    free(X);
    free(Y);
    free(X_1);
    free(K);

    return 0;
}