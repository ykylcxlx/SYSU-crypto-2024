#include<stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include<iostream>
#include <vector>
#include<math.h>
#include<string>
#include <algorithm>
#include<x86intrin.h>
#include <cstring>
using namespace std;
#define uint8_t unsigned char
#define MGF1_BUF_SIZE 256
void print(uint8_t *data, int len) {
    for (int i = 0; i < len; i++) {
        printf("%02X ", data[i]);
        if((i+1)%32==0)
            printf("\n");
    }
    printf("\n");
}

#define uint8_t unsigned char
#define uint32_t unsigned int
#define uint64_t unsigned long long
#define BLOCK_SIZE 64  // 512 bits = 64 bytes
#define BYTE unsigned char


class BigInt {
public:
    static const uint64_t BASE = 0x100000000;
    static const int SIZE = 8 ;
    uint32_t P[SIZE]={0};
    uint32_t R[SIZE] = {0};
    uint32_t R2[SIZE] = {0};
    uint32_t P_[SIZE] = {0};
    uint32_t ZERO[SIZE] = {0};
    uint32_t ONE[SIZE] = {1};
    uint32_t TWO[SIZE] = {2};
    uint32_t RES[SIZE] = {0};
    uint32_t THREE[SIZE] = {3};

    int P_bits = 0;
    int R_bits = 0;

    uint32_t pow2[32] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                         1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288,
                         1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912,
                         1073741824, 2147483648};
public:
    int getBits(const uint32_t a[SIZE])
    {
        for (int i = SIZE-1; i >= 0; i--)
        {
            if (a[i])
                for (int j = 31; j >= 0; j--)
                    if (a[i] & pow2[j])
                        return i * 32 + j;
        }

        return 0;
    }
    
    void exculid(uint32_t a[SIZE], uint32_t b[SIZE], uint32_t x[SIZE], uint32_t y[SIZE])
    {
        bool flag = false;
        for (int i = 0; i < SIZE; i++)
            if (b[i])
            {
                flag = true;
                break;
            }

        if (!flag)
        {
            x[0] = 1;
            y[0] = 0;
            return;
        }
        uint32_t temp[SIZE] = {0};
        divInternal(temp, a, b);
        mulInternal(temp, b, temp);
        subInternal(temp, a, temp);

        exculid(b, temp, y, x);

        uint32_t temp2[SIZE] = {0};
        divInternal(temp2, a, b);
        mulInternal(temp2, temp2, x);

        if (isBigger(y, temp2))
        {
            subInternal(temp2, y, temp2);
            modR(temp2, temp2);
        }
        else
        {
            subInternal(temp2, temp2, y);
            modR(temp2, temp2);
            subInternal(temp2, R, temp2);
        }

        for (int i = 0; i < SIZE; i++)
            y[i] = temp2[i];
    }
    void gcd(uint32_t r0[SIZE], uint32_t a[SIZE], uint32_t b[SIZE])
    {
        // uint32_t r0[SIZE];
        memcpy(r0, a, SIZE * sizeof(uint32_t));
        uint32_t r1[SIZE];
        memcpy(r1, b, SIZE * sizeof(uint32_t));
        
        uint32_t r2[SIZE]={0};
        uint32_t q[SIZE]={0};

        
        while(!isEqual(r1,ZERO)){
            divInternal(q,r0,r1);
            mulInternal(r2,r1,q);
            subInternal(r2,r0,r2);
            memcpy(r0,r1,SIZE*sizeof(uint32_t));
            memcpy(r1,r2,SIZE*sizeof(uint32_t));
        }

        // memcpy(res,r0,SIZE*sizeof(uint32_t));

    }
    void gcd(__uint128_t a, __uint128_t b, __uint128_t &res)
    {
        while (b != 0)
        {
            __uint128_t temp = a % b;
            a = b;
            b = temp;
        }
        res = a;
    }
    BigInt ex_gcd(BigInt a, BigInt b, BigInt& x, BigInt& y) {
        if (b.isEqual(b.RES,ZERO)) {
            x = ONE;
            y = ZERO;
            // cout<<"return"<<endl;
            return a;
        }

        BigInt d = ex_gcd(b, a % b, x, y);
        BigInt temp = x;
        x = y;
        // cout<<"a="<<a<<endl;
        // cout<<"b="<<b<<endl;
        // cout<<"x="<<x<<endl;
        // cout<<"y="<<y<<endl;
        y = temp - (a / b) * y;
        return d;
    }
    void resInit(const string &s);
    void resInit(uint32_t res[SIZE]);
    void removeLeadingZeros(string &s);
    bool isBigger(uint32_t a[SIZE], uint32_t b[SIZE]);
    bool isEqual(uint32_t a[SIZE], uint32_t b[SIZE]);
    void strDiv2(string &s);
    void strMul2(string &s);
    void strAdd1(string &s);
    void strToBi(uint32_t res[SIZE], string &s);
    string biToStr(uint32_t res[SIZE]);
    void addInternal(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void subInternal(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void mulInternal(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void divInternal(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void modInternal(uint32_t res[SIZE], uint32_t a[SIZE]);
    void modn(uint32_t res[SIZE], uint32_t a[SIZE],uint32_t N[SIZE]);//n为P的因子
    void modAdd(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void modSub(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void modSub(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE],uint32_t n[SIZE]);
    void modMul(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void modPow(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void modR(uint32_t res[SIZE], uint32_t a[SIZE]);
    int mod_int(uint32_t a[SIZE],int n);
    void divR(uint32_t res[SIZE], uint32_t a[SIZE+1]);
    void modReduction(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void modMulMontgomery(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void modPowMontgomery(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void invExculid(uint32_t res[SIZE], uint32_t a[SIZE]);
    void invExculid(uint32_t res[SIZE], uint32_t a[SIZE],uint32_t P[SIZE]);
    void pollardRho(uint32_t n[SIZE], uint32_t alpha[SIZE], uint32_t beta[SIZE]);


    void preCalculate();


        // 运算符重载
    BigInt operator+(BigInt &b){
        BigInt result(this->P);
        result.modAdd(result.RES, this->RES, b.RES); // result = (this + b) mod P
        return result;
    }

    BigInt operator-(const BigInt &b){
        uint32_t temp[SIZE];
        memcpy(temp, b.RES, SIZE * sizeof(uint32_t));
        BigInt result(this->P);
        result.modSub(result.RES, this->RES, temp); // result = (this - b) mod P
        return result;
    }

    BigInt operator*(BigInt &b) {
        BigInt result(this->P);
        // result.modReduction(result.RES, this->RES, b.RES); // result = (this * b) mod P
        // result.modInternal(result.RES,result.RES);
        result.modMul(result.RES, this->RES, b.RES); // result = (this * b) mod P
        return result;
    }

    BigInt operator/(BigInt &b) {
        BigInt result(this->P);
        divInternal(result.RES, this->RES, b.RES);
        return result;
    }

    BigInt operator^(BigInt &b) {
        BigInt result(this->P);
        
        result.preCalculate();
        
        result.modPowMontgomery(result.RES, this->RES, b.RES); // result = (this ^ b) mod P
        // result.modPow(result.RES, this->RES, b.RES); // result = (this ^ b) mod P
        return result;
    }

    BigInt operator%(uint32_t N[SIZE]) {
        BigInt result((const uint32_t*)N);
        // Implement mod operation with N if needed
        result.modn(result.RES, this->RES,N);
        return result;
    }

    BigInt operator%(BigInt &b) {
        BigInt result(this->P);
        // uint32_t N = b.RES;
        // cout<<biToStr(b.RES)<<endl;
        result.modn(result.RES, this->RES,b.RES);
        return result;
    }

    // 求逆运算符重载
    BigInt operator-(){
        //求逆要加括号
        BigInt result(this->P);
        result.invExculid(result.RES,this->RES); // 计算逆元
        return result;
    }
    bool operator==(const BigInt &b) const {
        // return std::equal(std::begin(RES), std::end(RES), std::begin(b.RES));
        for(int i=0;i<SIZE;i++){
            if(RES[i]!=b.RES[i])
                return false;
        }
        return true;
    }

    void update(uint32_t f0[SIZE], uint32_t f1[SIZE], uint32_t f2[SIZE],uint32_t n[SIZE], uint32_t alpha[SIZE], uint32_t beta[SIZE]);
    void update(uint32_t f0[SIZE], __uint128_t &f1,__uint128_t &f2 , __uint128_t &n,uint32_t alpha[SIZE], uint32_t beta[SIZE]);
    // 输出运算符重载
    friend std::ostream &operator<<(std::ostream &os,  BigInt bigInt) {
        // BigInt temp("1");
        // uint32_t temp[SIZE];
        // memcpy(temp, bigInt.RES, SIZE * sizeof(uint32_t));
        std::string s;
        // string a = bigInt.biToStr(bigInt.RES); // 使用 biToStr 方法转换为字符串
        os << bigInt.biToStr(bigInt.RES);
        return os;
    }
    // BigInt(){}
    BigInt(const string &p);
    BigInt(const uint32_t p[SIZE]);
    void performOperations(int n);
};
BigInt::BigInt(const string &p) {
    string reversed = p; // 复制字符串以进行反转
    reverse(reversed.begin(), reversed.end());
    strToBi(P, reversed);
    P_bits = getBits(P); // 确保 getBits 函数已定义

    preCalculate();

    // P_words = P_bits / 32;
    // R_words = P_words + 1;
    // R[R_words] = 1;
    // R[(P_bits + 1) / 32] = pow2[(P_bits + 1) % 32];
    // R_bits = getBits(R);
    // // R_words = R_bits / 32;

    // modInternal(R2, R);
    // mulInternal(R2, R2, R2);
    // modInternal(R2, R2);

    // uint32_t negP[128] = {0};
    // subInternal(negP, R, P);

    // uint32_t x[128] = {0},
    //          y[128] = {0};
    // exculid(negP, R, x, y);
    // for (int i = 0; i < 128; i++)
    //     P_[i] = x[i];

    // cout<<"P="<<biToStr(P)<<endl;
    // cout<<"R="<<biToStr(R)<<endl;
    // cout<<"R2="<<biToStr(R2)<<endl;
    // cout<<"P_="<<biToStr(P_)<<endl;
    // cout<<"P_bits"<<P_bits<<endl;
    // cout<<"P_words"<<P_words<<endl;
    // cout<<"R_bits"<<R_bits<<endl;
    // cout<<"R_words"<<R_words<<endl;
}


void BigInt::resInit(const string &res){
    memset(RES, 0, sizeof(RES));
    string reversed = res; // 复制字符串以进行反转
    reverse(reversed.begin(), reversed.end());
    strToBi(this->RES, reversed);
}

void BigInt::resInit(uint32_t res[SIZE]){
    memcpy(RES, res, SIZE * sizeof(uint32_t));
}
BigInt::BigInt(const uint32_t p[SIZE]) {

    // for(int i = 0;i<SIZE;i++)
    //     P[i] = p[i];
    memcpy(P, p, SIZE*sizeof(uint32_t));
    P_bits = getBits(P); // 确保 getBits 函数已定义

    // P_words = P_bits / 32;
    // R_words = P_words + 1;
    // R[R_words] = 1;
    // R[(P_bits + 1) / 32] = pow2[(P_bits + 1) % 32];
    // R_bits = getBits(R);
    // // R_words = R_bits / 32;

    // modInternal(R2, R);
    // mulInternal(R2, R2, R2);
    // modInternal(R2, R2);

    // uint32_t negP[128] = {0};
    // subInternal(negP, R, P);

    // uint32_t x[128] = {0},
    //          y[128] = {0};
    // exculid(negP, R, x, y);
    // for (int i = 0; i < 128; i++)
    //     P_[i] = x[i];

}

// Implementations of the methods

void BigInt::removeLeadingZeros(string &s) {
    size_t end = s.find_last_not_of('0');
    if (end != string::npos)
        s.erase(end + 1);
    else
        s = "0";
}

bool BigInt::isBigger(uint32_t a[SIZE], uint32_t b[SIZE]) {
    for (int i = SIZE-1; i >= 0; i--) {
        if (a[i] > b[i]) return true;
        else if (a[i] < b[i]) return false;
    }
    return false;
}

bool BigInt::isEqual(uint32_t a[SIZE], uint32_t b[SIZE]) {
    for (int i = 0; i < SIZE; i++)
        if (a[i] != b[i]) return false;
    return true;
}

void BigInt::strDiv2(string &s) {
    int carry = 0;
    for (int i = s.length() - 1; i >= 0; i--) {
        int x = s[i] - '0';
        s[i] = (x + carry * 10) / 2 + '0';
        carry = x % 2;
    }
    removeLeadingZeros(s);
}

void BigInt::strMul2(string &s) {
    int carry = 0;
    for (std::string::size_type i = 0; i < s.length(); i++) {
        int x = s[i] - '0';
        s[i] = (x * 2 + carry) % 10 + '0';
        carry = (x * 2 + carry) / 10;
    }
    if (carry) s = s + "1";
}

void BigInt::strAdd1(string &s) {
    int carry = 1;
    for (std::string::size_type i = 0; i < s.length(); i++) {
        int x = s[i] - '0';
        s[i] = (x + carry) % 10 + '0';
        carry = (x + carry) / 10;
    }
    if (carry) s = s + "1";
}

void BigInt::strToBi(uint32_t res[SIZE], string &s) {
    int count = 0;
    while (s != "0") {
        if ((s[0] - '0') & 1) {
            res[count / 32] += pow2[count % 32];
        }
        strDiv2(s);
        count++;
    }
}

string BigInt::biToStr(uint32_t res[SIZE]) {
    string s = "0";
    for (int i = SIZE-1; i >= 0; i--) {
        for (int j = 31; j >= 0; j--) {
            strMul2(s);
            if (res[i] & pow2[j])
                strAdd1(s);
        }
    }
    reverse(s.begin(), s.end());
    return s;
}

void BigInt::addInternal(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]) {
    uint64_t carry = 0;
    uint32_t temp[SIZE] = {0};
    for (int i = 0; i < SIZE; i++) {
        uint64_t sum = carry + a[i] + b[i];
        temp[i] = sum & 0xffffffff;
        carry = sum >> 32;
    }
    for (int i = 0; i < SIZE; i++)
        res[i] = temp[i];
}


void BigInt::subInternal(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE])
{
    uint64_t borrow = 0;
    uint32_t temp[SIZE] = {0};
    for (int i = 0; i < SIZE; i++)
    {
        uint64_t temp1 = b[i] + borrow;
        if (a[i] < temp1)
        {
            temp[i] = BASE + a[i] - temp1;
            borrow = 1;
        }
        else
        {
            temp[i] = a[i] - temp1;
            borrow = 0;
        }
    }

    for (int i = 0; i < SIZE; i++)
        res[i] = temp[i];
}

void BigInt::mulInternal(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]) {
    uint32_t temp[256] = {0};
    for (int i = 0; i < SIZE; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < SIZE; j++) {
            uint64_t sum = (uint64_t)a[i] * b[j] + temp[i + j] + carry;
            temp[i + j] = sum & 0xffffffff;
            carry = sum >> 32;
        }
        if (carry)
            temp[i + SIZE] += carry;
    }
    for (int i = 0; i < SIZE; i++)
        res[i] = temp[i];
}


void BigInt::divInternal(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE])
{
    if (isEqual(b, ZERO))
    {
        cout << "Error: divide by zero" << endl;
        return;
    }
    uint32_t temp[SIZE] = {0};
    int n = getBits(b) / 32;
    int m = getBits(a) / 32 - n;

    uint32_t d[SIZE] = {0};
    d[0] = BASE / (b[n] + (uint64_t)1);
    uint32_t u_[SIZE+1] = {0}, v_[SIZE] = {0};

    uint64_t carry = 0;
    for (int i = 0; i < SIZE; i++)
    {
        uint64_t temp = (uint64_t)a[i] * d[0] + carry;
        u_[i] = temp & 0xffffffff;
        carry = temp >> 32;
    }
    if (carry)
        u_[SIZE] = carry;

    mulInternal(v_, b, d);

    int j = m;

    while (j >= 0)
    {
        uint32_t tem[SIZE+1] = {0};
        for (int i = 0; i <= n + 1; i++)
            tem[i] = u_[i + j];
        uint64_t tem2 = (tem[n + 1] * BASE + tem[n]) / v_[n];
        tem2 = min(tem2, BASE - 1);

        uint32_t q_hat[SIZE] = {static_cast<uint32_t>(tem2 & 0xffffffff)};
        uint32_t qv[SIZE] = {0};
        mulInternal(qv, v_, q_hat);
        while (isBigger(qv, tem))
        {
            q_hat[0]--;
            subInternal(qv, qv, v_);
        }

        subInternal(tem, tem, qv);
        for (int i = 0; i <= n + 1; i++)
            u_[i + j] = tem[i];

        temp[j] = q_hat[0];

        j--;
    }

    for (int i = 0; i < SIZE; i++)
        res[i] = temp[i];
}

void BigInt::modInternal(uint32_t res[SIZE], uint32_t a[SIZE]) {
    uint32_t temp[SIZE] = {0};
    divInternal(temp, a, P);
    mulInternal(temp, temp, P);
    subInternal(res, a, temp);
}

void BigInt::modn(uint32_t res[SIZE], uint32_t a[SIZE],uint32_t N[SIZE]) {//N是P的因子
    uint32_t temp[SIZE] = {0};
    divInternal(temp, a, N);
    mulInternal(temp, temp, N);
    subInternal(res, a, temp);
}

int BigInt::mod_int(uint32_t a[SIZE],int n){
    __uint128_t d = 0;
    for(int i = SIZE - 1;i >= 0; i--){
        d = (d * BASE + a[i]) % n;
    }
    return static_cast<int> (d);
}



void BigInt::modAdd(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]) {
    addInternal(res, a, b);
    modInternal(res, res);
}

void BigInt::modSub(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]) {
    if (isBigger(a, b)) {
        subInternal(res, a, b);
        modInternal(res, res);
    } else {
        subInternal(res, b, a);
        modInternal(res, res);
        subInternal(res, P, res);
    }
}

void BigInt::modSub(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE],uint32_t n[SIZE]) {
    if (isBigger(a, b)) {
        subInternal(res, a, b);
        modn(res, res,n);
    } else {
        subInternal(res, b, a);
        modn(res, res,n);
        subInternal(res, n, res);
    }
}

void BigInt::modMul(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]) {
    uint32_t temp[SIZE] = {0};
    mulInternal(temp, a, b);
    modInternal(res, temp);
}

void BigInt::modPow(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]) {
    uint32_t temp[SIZE] = {1};
    uint32_t base[SIZE] = {0};
    for (int i = 0; i < SIZE; i++)
        base[i] = a[i];

    for (int i = 0; i < SIZE; i++) {
        if (isBigger(temp, ONE) && b[i] == 0)
            break;
        for (int j = 0; j < 32; j++) {
            if (b[i] & pow2[j]) {
                modMul(temp, temp, base);
            }
            modMul(base, base, base);
        }
    }
    for (int i = 0; i < SIZE; i++)
        res[i] = temp[i];
}

    void BigInt::modR(uint32_t res[SIZE], uint32_t a[SIZE]) {
        uint32_t temp[SIZE] = {0};
        for (int i = 0; i < R_bits / 32; i++)
            temp[i] = a[i];
        for (int i = 0; i < R_bits % 32; i++)
            if (a[R_bits / 32] & pow2[i])
                temp[R_bits / 32] += pow2[i];

        for (int i = 0; i < SIZE; i++)
            res[i] = temp[i];
    }

    // // 新的 divR 函数
    void BigInt::divR(uint32_t res[SIZE], uint32_t a[SIZE+1]) {
        uint32_t temp[SIZE] = {0};

        int word_shift = R_bits / 32;
        int bits_shift = R_bits % 32;
        for (int i = 0; i + word_shift < SIZE+1; i++)
            temp[i] = a[i + word_shift];

        if (bits_shift != 0) {
            for (int i = 0; i < SIZE - 1; i++)
                temp[i] = temp[i] >> bits_shift | (temp[i + 1] << (32 - bits_shift));
            temp[SIZE-1] = temp[SIZE-1] >> bits_shift;
        }

        for (int i = 0; i < SIZE; i++)
            res[i] = temp[i];
        // cout<<"dicR"<<biToStr(res)<<endl;

    }
    //蒙哥马利模乘
    void BigInt::modReduction(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]) {
        uint32_t T[SIZE] = {0};
        
        mulInternal(T, a, b);
        // cout<<"t: "<<biToStr(T)<<endl;
        // exit(0);
        
        uint32_t m[SIZE] = {0};
        uint32_t t[SIZE] = {0};
        modR(m, T); // m=T mod R

        mulInternal(m, m, P_); // m=(T mod R) * P_
        modR(m, m);   // m=((T mod R) * P_) mod R
        mulInternal(t, m, P);  // t=m*P
        
        uint64_t carry = 0;
        uint32_t te[129] = {0};
        for (int i = 0; i < SIZE; i++)
        {
            uint64_t sum = carry + t[i] + T[i];
            te[i] = sum & 0xffffffff;
            carry = sum >> 32;
        }
        if (carry)
            te[SIZE] = carry;

        divR(t, te);
        // cout<<biToStr(t)<<endl;
        if (isBigger(t, P) || isEqual(t, P)){
            subInternal(t, t, P);
        }

        for (int i = 0; i < SIZE; i++)
            res[i] = t[i];


        // cout<<"res: "<<biToStr(res)<<endl;
        // exit(0);
        
    }
    void BigInt::modMulMontgomery(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]){
        uint32_t x_[SIZE] = {0}, y_[SIZE] = {0};
        // cout<<"a: "<<biToStr(a)<<endl;
        // cout<<"b: "<<biToStr(b)<<endl;
        modReduction(x_, a, R2);
        modReduction(y_, b, R2);
        modReduction(x_, x_, y_);
        modReduction(x_, x_, ONE);

        // cout<<"modMulMontgomery: "<<biToStr(x_)<<endl;

        for (int i = 0; i < SIZE; i++)
            res[i] = x_[i];
    }
    // 新的模幂函数
    void BigInt::modPowMontgomery(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]) {

            uint32_t temp[SIZE] = {0};
            uint32_t base[SIZE] = {0};
            
            modReduction(temp, ONE, R2);
            modReduction(base, a, R2);
            
            for (int i = 0; i < SIZE; i++)
            {
                if (isBigger(temp, ONE) && b[i] == 0)
                    break;
                for (int j = 0; j < 32; j++)
                {
                    if (b[i] & pow2[j])
                        modReduction(temp, temp, base);

                    modReduction(base, base, base);
                }
            }

            modReduction(res, temp, ONE);
    }
        
        void BigInt::invExculid(uint32_t res[SIZE], uint32_t a[SIZE]){
            uint32_t x1[SIZE] = {0}, y1[SIZE] = {1}, x2[SIZE] = {0}, y2[SIZE] = {0};
            for (int i = 0; i < SIZE; i++)
                x1[i] = a[i];
            for (int i = 0; i < SIZE; i++)
                x2[i] = P[i];

            while (!(isEqual(x1, ONE) || isEqual(x2, ONE)))
            {
                if ((x1[0] & 1) == 0 && (y1[0] & 1) == 0)
                {
                    divInternal(x1, x1, TWO);
                    divInternal(y1, y1, TWO);
                }
                else if ((x1[0] & 1) == 0 && (y1[0] & 1) == 1)
                {
                    divInternal(x1, x1, TWO);
                    addInternal(y1, y1, P);
                    divInternal(y1, y1, TWO);
                }

                if ((x2[0] & 1) == 0 && (y2[0] & 1) == 0)
                {
                    divInternal(x2, x2, TWO);
                    divInternal(y2, y2, TWO);
                }
                else if ((x2[0] & 1) == 0 && (y2[0] & 1) == 1)
                {
                    divInternal(x2, x2, TWO);
                    addInternal(y2, y2, P);
                    divInternal(y2, y2, TWO);
                }

                if (isBigger(x1, x2))
                {
                    modSub(x1, x1, x2);
                    modSub(y1, y1, y2);
                }
                else
                {
                    modSub(x2, x2, x1);
                    modSub(y2, y2, y1);
                }
        }

        if (isEqual(x1, ONE))
            for (int i = 0; i < SIZE; i++)
                res[i] = y1[i];
        else
            for (int i = 0; i < SIZE; i++)
                res[i] = y2[i];
            
    }

    void BigInt::invExculid(uint32_t res[SIZE], uint32_t a[SIZE],uint32_t P[SIZE]){


            uint32_t x1[SIZE] = {0}, y1[SIZE] = {1}, x2[SIZE] = {0}, y2[SIZE] = {0};
            for (int i = 0; i < SIZE; i++)
                x1[i] = a[i];
            for (int i = 0; i < SIZE; i++)
                x2[i] = P[i];

            while (!(isEqual(x1, ONE) || isEqual(x2, ONE)))
            {
                if ((x1[0] & 1) == 0 && (y1[0] & 1) == 0)
                {
                    divInternal(x1, x1, TWO);
                    divInternal(y1, y1, TWO);
                }
                else if ((x1[0] & 1) == 0 && (y1[0] & 1) == 1)
                {
                    divInternal(x1, x1, TWO);
                    addInternal(y1, y1, P);
                    divInternal(y1, y1, TWO);
                }

                if ((x2[0] & 1) == 0 && (y2[0] & 1) == 0)
                {
                    divInternal(x2, x2, TWO);
                    divInternal(y2, y2, TWO);
                }
                else if ((x2[0] & 1) == 0 && (y2[0] & 1) == 1)
                {
                    divInternal(x2, x2, TWO);
                    addInternal(y2, y2, P);
                    divInternal(y2, y2, TWO);
                }

                if (isBigger(x1, x2))
                {
                    modSub(x1, x1, x2,P);
                    modSub(y1, y1, y2,P);
                }
                else
                {
                    modSub(x2, x2, x1,P);
                    modSub(y2, y2, y1,P);
                }
        }

        if (isEqual(x1, ONE))
            for (int i = 0; i < SIZE; i++)
                res[i] = y1[i];
        else
            for (int i = 0; i < SIZE; i++)
                res[i] = y2[i];
        
    }

void BigInt::preCalculate()
{
    // P_words = P_bits / 32;
    // R_words = P_words + 1;
    // R[R_words] = 1;
    R[(P_bits + 1) / 32] = pow2[(P_bits + 1) % 32];
    R_bits = getBits(R);
    // R_words = R_bits / 32;

    modInternal(R2, R);
    mulInternal(R2, R2, R2);
    modInternal(R2, R2);

    uint32_t negP[SIZE] = {0};
    subInternal(negP, R, P);

    uint32_t x[SIZE] = {0},
             y[SIZE] = {0};
    exculid(negP, R, x, y);
    for (int i = 0; i < SIZE; i++)
        P_[i] = x[i];

    // cout<<"P="<<bi2str(P)<<endl;
    // cout<<"R="<<bi2str(R)<<endl;
    // cout<<"R2="<<bi2str(R2)<<endl;
    // cout<<"P_="<<bi2str(P_)<<endl;
    // cout<<"P_bits"<<P_bits<<endl;
    // cout<<"P_words"<<P_words<<endl;
    // cout<<"R_bits"<<R_bits<<endl;
    // cout<<"R_words"<<R_words<<endl;
}
void BigInt::update(uint32_t f0[SIZE], uint32_t f1[SIZE], uint32_t f2[SIZE],uint32_t n[SIZE], uint32_t alpha[SIZE], uint32_t beta[SIZE])
{

    // uint32_t temp[SIZE]={0};
    // modn(temp,f0,THREE);
    int temp = mod_int(f0,3);

    if(temp==1){
        // S1
        modMul(f0,f0,beta);
        addInternal(f2,f2,ONE);
        modn(f2,f2,n);
    }else if(temp==0){
        modMul(f0,f0,f0);
        mulInternal(f1,f1,TWO);
        mulInternal(f2,f2,TWO);
        // uint32_t f0_[SIZE]={0};
        // for(int i = 0;i<SIZE;i++){
        //     f0_[2*i] = f0[i];
        // }
        // memcpy(f0,f0_,SIZE);
        // for(int i = 0;i<SIZE;i++){
        //     f1[i] = f1[i] << 1;
        //     f2[i] = f2[i] << 1;
        // }
        modn(f1,f1,n);
        modn(f2,f2,n);
    }else{
        modMul(f0,f0,alpha);
        addInternal(f1,f1,ONE);
        modn(f2,f2,n);
        
    }
}
void BigInt::update(uint32_t f0[SIZE], __uint128_t &f1, __uint128_t &f2,__uint128_t& n, uint32_t alpha[SIZE], uint32_t beta[SIZE]){
    int temp = mod_int(f0,3);

    if(temp==1){
        // S1
        modMul(f0,f0,beta);
        f2 = f2+1;
        f2 = f2 % n;
    }else if(temp==0){
        modMul(f0,f0,f0);
        f1 = f1*2%n;;
        f2 = f2*2%n;;
    }
    else{
        modMul(f0,f0,alpha);
        f1 = f1+1;
        f2 = f2 % n;
    }

}

#define outputFile stdout
// 实现RSA加密算法
int main(int argc, char *argv[]) {
    #ifdef ONLINE_JUDGE
    #ifdef _WIN32
        setmode(fileno(stdin), O_BINARY);
        setmode(fileno(stdout), O_BINARY);
    #endif
    #define file stdin
    #define outputFile stdout
    #else
    FILE *file = fopen("dump1.bin", "rb");
    if (!file) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }
    #endif
    BigInt temp("1");
    
    // uint8_t P[16]={0};
    // uint8_t G[16]={0};
    // uint8_t S[16]={0};
    uint32_t p[8]={0};
    uint32_t g[8]={0};
    uint32_t s[8]={0};
    uint32_t p_div2[8]={0};


    uint32_t n;
    uint8_t ans[10000];
    string  str;
    int cnt0=0,cnt1=0;

    fread(p,1,16,file);
    fread(g,1,16,file);
    fread(s,1,16,file);
    fread(&n,4,1,file);
    BigInt cal(p);
    cal.preCalculate();
    cal.divInternal(p_div2,p,cal.TWO);
    // cout<<temp.biToStr(p_div2)<<endl;
    int t = n*8;
    for(int i = 0; i<t;i++){
        cal.modPowMontgomery(s,g,s);
        if(!cal.isBigger(s,p_div2)){

        }else{

            ans[i/8] |= 0x01<<(i%8);
            //1的数量
            cnt1++;
        }
    }
    // for(int i = 0; i<n;i++){
    //     printf("%02x ",ans[i]);
    // }
    fwrite(ans, n,1,outputFile);
    cnt0 = t-cnt1;
    fwrite(&cnt0,4,1,outputFile);
    fwrite(&cnt1,4,1,outputFile);
    return EXIT_SUCCESS;
}
 