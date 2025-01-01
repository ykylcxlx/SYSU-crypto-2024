#include<stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include<iostream>
#include <vector>
#include<math.h>
#include<string>
#include <algorithm>
#include<x86intrin.h>
using namespace std;
#define uint8_t unsigned char
#define MGF1_BUF_SIZE 256



    // uint8_t seed[32]= {
    // 0x84, 0xEE, 0x1D, 0x92, 0xBE, 0xFC, 0x55, 0x96, 0x6F, 0x20,
    // 0xB6, 0xFD, 0x18, 0xFC, 0x45, 0x20, 0xCF, 0x17, 0x0B, 0xD8,
    // 0x92, 0xE2, 0xE0, 0xD3, 0x4F, 
    // 0xE7, 0xA0, 0x10, 0x17, 0xC8,
    // 0x6B, 0x67};
    uint8_t hash_L[32]= {
    0xE3, 0xB0, 0xC4, 0x42, 0x98, 0xFC, 0x1C, 0x14, 0x9A, 0xFB, 0xF4, 0xC8, 0x99, 0x6F, 0xB9, 0x24,
    0x27, 0xAE, 0x41, 0xE4, 0x64, 0x9B, 0x93, 0x4C, 0xA4, 0x95, 0x99, 0x1B, 0x78, 0x52, 0xB8, 0x55};


void print(uint8_t *data, int len) {
    for (int i = 0; i < len; i++) {
        printf("%02X ", data[i]);
        if((i+1)%32==0)
            printf("\n");
    }
    printf("\n");
}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>

#include <memory.h>
#include <stdlib.h>
#include<x86intrin.h>
#include <stddef.h>
// #include <stdint.h>
// #include <Windows.h>

#define uint8_t unsigned char
#define uint32_t unsigned int
#define uint64_t unsigned long long
#define BLOCK_SIZE 64  // 512 bits = 64 bytes
#define BYTE unsigned char

size_t bitlen=0;
size_t datalen=0;
static const uint32_t K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

// 初始哈希值
const uint32_t H0[8] = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
};

// 右旋转
uint32_t rotr(uint32_t x, uint32_t n) {
    return (x >> n) | (x << (32 - n));
}
void process_block(const uint8_t *block) {
    // 这里处理每个512位的块
    // 例如，可以调用SHA-256的核心算法
    // 这里只是一个占位符
    printf("Processing block:\n");
    for (int i = 0; i < BLOCK_SIZE; i++) {
        printf("%02x ", block[i]);
        if(i % 8 ==7){
            printf("\n");
        }
    }
    printf("\n");
}
//SHA-256压缩函数
void sha256_compress(uint32_t *state,const uint8_t *block) {
    uint32_t W[64];
    uint32_t a, b, c, d, e, f, g, h, t1, t2;

    // 消息调度
    for (int t = 0; t < 16; t++) {
        W[t] = (block[t * 4] << 24) | (block[t * 4 + 1] << 16) | (block[t * 4 + 2] << 8) | block[t * 4 + 3];
    }
    for (int t = 16; t < 64; t++) {
        W[t] = W[t - 16] + (rotr(W[t - 15], 7) ^ rotr(W[t - 15], 18) ^ (W[t - 15] >> 3)) + W[t - 7] + (rotr(W[t - 2], 17) ^ rotr(W[t - 2], 19) ^ (W[t - 2] >> 10));
    }

    // 初始化工作变量
    a = state[0];
    b = state[1];
    c = state[2];
    d = state[3];
    e = state[4];
    f = state[5];
    g = state[6];
    h = state[7];

    // 主循环
    for (int t = 0; t < 64; t++) {
        t1 = h + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ ((~e) & g)) + K[t] + W[t];
        t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
        h = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
    }

    // 更新状态
    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    state[5] += f;
    state[6] += g;
    state[7] += h;
}


void sha256_update(BYTE data[], uint32_t state[],  uint8_t buffer[],size_t len)
{
	int i;

	for (i = 0; i < len; ++i) {
		buffer[datalen] = data[i];
		datalen++;
		if (datalen == 64) {
			sha256_compress(state, buffer);
			bitlen += 512;
			datalen = 0;
		}
	}
}

// 9B 11 FF 85 57 C0 FA D9 1B 72 1E 94 A3 51 45 CC B3 28 66 E7 35 40 99 37 5A 01 EE E0 4D 3E 59 F4 C6 1D EB 24 C9 09 2C 9A 5F E0 DA 4A E5 17 98 AA D2 D8 3C 92 CC 8D 9A 15 9D 3F C1 BC E7 89 D5 B7 
// 10 50 60 EC 93 3A 3C EC D8 CB 3D A6 F1 DE 37 1D 88 78 5D 20 82 DF BD 5E 4E 69 D9 95 1B 80 71 4E DD A1 E5 14 89 8A 31 CE 9A D3 EA 82 B1 46 DA E9 DA EE 15 66 E0 AF F6 FF 8D 2B 8A 82 3D 43 7B 46 
// 50 C2 6B 40 51 3A FE 2E 79 62 B6 03 C3 B1 6F BA 9A DF 42 EF 4B 72 58 E7 41 00 B4 16 9B 31 DD E3 01 C1 D5 B5 6E DA 7B 39 A6 00 A0 9E 91 C2 07 23 8E 57 74 57 08 3B B2 8A 0A D8 D4 E8 22 DD 3F EC 
// C3 97 BE 90 53 16 C8 F9 F4 55 3F 32 41 5D FA 64 A4 E9 B1 4B 8B 5A 70 2B 2C 18 3D E0 FB AE AA

void sha256_final(BYTE data[], uint32_t state[], BYTE digest[])
{
	int i;

	i = datalen;
    

	// Pad whatever data is left in the buffer.
	if (datalen < 56) {
		data[i++] = 0x80;
		while (i < 56)
			data[i++] = 0x00;
	}
	else {
		data[i++] = 0x80;
		while (i < 64)
			data[i++] = 0x00;

		sha256_compress(state, data);

        // process_block(data);
		memset(data, 0, 56);
	}

	// Append to the padding the total message's length in bits and transform.
	bitlen += datalen * 8;
    // printf("bitlen = %d\n",bitlen);
    for (int i = 0; i < 8; i++) {
        data[BLOCK_SIZE -8 + i] = (bitlen >> (56 - i * 8)) & 0xFF;
    }

    // process_block(data);
    sha256_compress(state, data);

	// Since this implementation uses little endian byte ordering and SHA uses big endian,
	// reverse all the bytes when copying the final state to the output hash.
	for (i = 0; i < 4; ++i) {
		digest[i]      = (state[0] >> (24 - i * 8)) & 0x000000ff;
		digest[i + 4]  = (state[1] >> (24 - i * 8)) & 0x000000ff;
		digest[i + 8]  = (state[2] >> (24 - i * 8)) & 0x000000ff;
		digest[i + 12] = (state[3] >> (24 - i * 8)) & 0x000000ff;
		digest[i + 16] = (state[4] >> (24 - i * 8)) & 0x000000ff;
		digest[i + 20] = (state[5] >> (24 - i * 8)) & 0x000000ff;
		digest[i + 24] = (state[6] >> (24 - i * 8)) & 0x000000ff;
		digest[i + 28] = (state[7] >> (24 - i * 8)) & 0x000000ff;
	}
    //!!!!!!!!!!!!!!
    memcpy(state, H0, sizeof(H0));
    datalen = bitlen = 0;
}

// void bufferext(uint8_t *buffer, uint32_t counter,unsigned int len)
// {
//     for(int i = len;i<len+4;i++){
//         buffer[i] = (counter >> ((len+3-i)*8)) & 0xff;
//     }   
// }

int MGF1(const unsigned char *mgfSeed, unsigned int mgfSeedLen, unsigned int maskLen, unsigned char *mask)
{
    unsigned char buf[MGF1_BUF_SIZE], *p;
    unsigned char digest[32]; 
    unsigned long digestLen=32;
    unsigned long counter, restLen;

    unsigned char buffer_sha[64];
    unsigned char data[MGF1_BUF_SIZE];
    memset(buffer_sha, 0, 64);

    uint32_t state[8];
    memcpy(state, H0, sizeof(H0));
    if (mgfSeedLen > MGF1_BUF_SIZE - 4)
    {
        printf("MGF1 buffer is not long enough!\n");
        return -1;
    }

    // copy mgfSeed to buffer
    memcpy(buf, mgfSeed, mgfSeedLen);



    // clear rest buffer to 0
    p = buf + mgfSeedLen;
    memset(p, 0, MGF1_BUF_SIZE-mgfSeedLen);


    

    counter = 0;
    restLen = maskLen;

    while (restLen > 0)
    {
        p[0] = (counter >> 24) & 0xff;
        p[1] = (counter >> 16) & 0xff;
        p[2] = (counter >>  8) & 0xff;
        p[3] = counter & 0xff;
        memcpy(data, buf, MGF1_BUF_SIZE);
        
        // if(mgfSeedLen==223){
        // printf("data =");
        // print(data,223);
        // }


        memset(buffer_sha, 0, 64);
        
        if (restLen >= digestLen)
        {
            // HASH(alg, buf, mgfSeedLen+4, (unsigned char *)mask);
            sha256_update(data, state, buffer_sha, mgfSeedLen+4);

            sha256_final(buffer_sha, state, (unsigned char *)mask);

            // for(int i = 0; i < digestLen; i++)
            //     printf("%02x", mask[i]);
            // printf("\n");
        

            restLen -= digestLen;
            mask += digestLen;

            counter ++;
        }
        else // 剩余的不足单次哈希长度的部分
        {
            // HASH(alg, buf, mgfSeedLen+4, (unsigned char *)digest);
            sha256_update(data, state, buffer_sha, mgfSeedLen+4);
            sha256_final(buffer_sha, state, (unsigned char *)digest);
            memcpy(mask, digest, restLen);

            restLen = 0;
        }
    }

    return 0;
}



// 78 A1 3B C7 CF 3C E6 CD 81 89 EA 5C 3A 3E FC E8 94 86 27 03 51 DB 0A 7B FE 94 77 FB 35 6C E1 A1 
// C6 1D EB 24 C9 09 2C 9A 5F E0 DA 4A E5 17 98 AA D2 D8 3C 92 CC 8D 9A 15 9D 3F C1 BC E7 89 D5 B7 
// 10 50 60 EC 93 3A 3C EC D8 CB 3D A6 F1 DE 37 1D 88 78 5D 20 82 DF BD 5E 4E 69 D9 95 1B 80 71 4E 
// DD A1 E5 14 89 8A 31 CE 9A D3 EA 82 B1 46 DA E9 DA EE 15 66 E0 AF F6 FF 8D 2B 8A 82 3D 43 7B 46 
// 50 C2 6B 40 51 3A FE 2E 79 62 B6 03 C3 B1 6F BA 9A DF 42 EF 4B 72 58 E7 41 00 B4 16 9B 31 DD E3 
// 01 C1 D5 B5 6E DA 7B 39 A6 00 A0 9E 91 C2 07 23 8E 57 74 57 08 3B B2 8A 0A D8 D4 E8 22 DD 3F EC 
// C3 97 BE 90 53 16 C8 F9 F4 55 3F 32 41 5D FA 64 A5 A2 DE 26 EE 33 1A 42 0C 53 52 89 88 C6 C3
class BigInt {
private:
    static const uint64_t BASE = 0x100000000;
    static const int SIZE = 128;
    uint32_t P[SIZE]={0};
    uint32_t R[SIZE] = {0};
    uint32_t R2[SIZE] = {0};
    uint32_t P_[SIZE] = {0};
    uint32_t ZERO[SIZE] = {0};
    uint32_t ONE[SIZE] = {1};
    uint32_t TWO[SIZE] = {2};
    int P_bits = 0;
    int R_bits = 0;

    uint32_t pow2[32] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                         1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288,
                         1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912,
                         1073741824, 2147483648};
public:
    int getBits(const uint32_t a[128])
    {
        for (int i = 127; i >= 0; i--)
        {
            if (a[i])
                for (int j = 31; j >= 0; j--)
                    if (a[i] & pow2[j])
                        return i * 32 + j;
        }

        return 0;
    }
    
    void exculid(uint32_t a[128], uint32_t b[128], uint32_t x[128], uint32_t y[128])
    {
    bool flag = false;
    for (int i = 0; i < 128; i++)
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
    }
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
    void modAdd(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void modSub(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void modMul(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);
    void modPow(uint32_t res[SIZE], uint32_t a[SIZE], uint32_t b[SIZE]);

    void preCalculate();




    BigInt(const string &p);
    void performOperations(int n);
};
BigInt::BigInt(const string &p) {
    string reversed = p; // 复制字符串以进行反转
    reverse(reversed.begin(), reversed.end());
    strToBi(P, reversed);
    P_bits = getBits(P); // 确保 getBits 函数已定义

    preCalculate();
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
    for (int i = 127; i >= 0; i--) {
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
    for (int i = 127; i >= 0; i--) {
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



void BigInt::subInternal(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    uint64_t borrow = 0;
    uint32_t temp[128] = {0};
    for (int i = 0; i < 128; i++)
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

    for (int i = 0; i < 128; i++)
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


void BigInt::divInternal(uint32_t res[128], uint32_t a[128], uint32_t b[128])
{
    if (isEqual(b, ZERO))
    {
        cout << biToStr(a) << endl;
        cout << biToStr(b) << endl;
        cout << "Error: divide by zero" << endl;
        return;
    }
    uint32_t temp[128] = {0};
    int n = getBits(b) / 32;
    int m = getBits(a) / 32 - n;

    uint32_t d[128] = {0};
    d[0] = BASE / (b[n] + (uint64_t)1);
    uint32_t u_[129] = {0}, v_[128] = {0};

    uint64_t carry = 0;
    for (int i = 0; i < 128; i++)
    {
        uint64_t temp = (uint64_t)a[i] * d[0] + carry;
        u_[i] = temp & 0xffffffff;
        carry = temp >> 32;
    }
    if (carry)
        u_[128] = carry;

    mulInternal(v_, b, d);

    int j = m;

    while (j >= 0)
    {
        uint32_t tem[129] = {0};
        for (int i = 0; i <= n + 1; i++)
            tem[i] = u_[i + j];
        uint64_t tem2 = (tem[n + 1] * BASE + tem[n]) / v_[n];
        tem2 = min(tem2, BASE - 1);

        uint32_t q_hat[128] = {static_cast<uint32_t>(tem2 & 0xffffffff)};
        uint32_t qv[128] = {0};
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

    for (int i = 0; i < 128; i++)
        res[i] = temp[i];
}

void BigInt::modInternal(uint32_t res[SIZE], uint32_t a[SIZE]) {
    uint32_t temp[SIZE] = {0};
    divInternal(temp, a, P);
    mulInternal(temp, temp, P);
    subInternal(res, a, temp);
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

void BigInt::preCalculate() {
    R[(P_bits + 1) / 32] = pow2[(P_bits + 1) % 32];
    R_bits = getBits(R);
    modInternal(R2, R);
    mulInternal(R2, R2, R2);
    modInternal(R2, R2);

    uint32_t negP[SIZE] = {0};
    subInternal(negP, R, P);

    uint32_t x[SIZE] = {0}, y[SIZE] = {0};
    exculid(negP, R, x, y);

    
    for (int i = 0; i < SIZE; i++)
        P_[i] = x[i];
}

void BigInt::performOperations(int n) {
    while (n--) {
        string a, b;
        cin >> a >> b;
        reverse(a.begin(), a.end());
        reverse(b.begin(), b.end());

        uint32_t A[SIZE] = {0}, B[SIZE] = {0};
        strToBi(A, a);
        strToBi(B, b);

        uint32_t res1[SIZE] = {0};
        modAdd(res1, A, B);
        cout << biToStr(res1) << '\n';

        uint32_t res2[SIZE] = {0};
        modSub(res2, A, B);
        cout << biToStr(res2) << '\n';

        uint32_t res3[SIZE] = {0};
        modMul(res3, A, B);
        cout << biToStr(res3) << '\n';

        // uint32_t x[SIZE] = {0}, y[SIZE] = {0};
        // inv_exculid(A, P, x, y);
        // cout << biToStr(x) << '\n';

        uint32_t res5[SIZE] = {0};
        modPow(res5, A, B);
        cout << biToStr(res5) << '\n';

        if (n)
            cout << '\n';
    }
}
# define outputFile stdout
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
    FILE *file = fopen("dump.bin", "rb");
    if (!file) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }
    #endif
    uint8_t unknown [16];
    fread(unknown ,1,16,file);

    uint8_t n[256];
    uint8_t e[256];
    uint32_t N[128]={0};
    uint32_t E[128]={0};
    uint32_t M[128]={0};

    uint8_t skip[256];
    uint8_t len;

    // 读取n和e
    // fread(n, 1, 256, file);
    // fread(e, 1, 256, file);
    BigInt cal0("1");

    for(int i = 0;i<64;i++){
        fread(n+i*4, 1, 4, file);
         uint32_t   value = ((uint32_t)n[4*i+0] << 24) |
            ((uint32_t)n[4*i+1] << 16) |
            ((uint32_t)n[4*i+2] << 8)  |
            ((uint32_t)n[4*i+3]);

        // printf("n[%d] = %u\n", i, value);
        // N = N*BigIntTiny(2).pow(BigIntTiny(32)) + BigIntTiny(value);
         N[64-i-1]=value;

    }
    // printf("N=");
    // cout<<cal0.biToStr(N)<<endl;
    BigInt cal(cal0.biToStr(N));

    for(int i = 0;i<64;i++){
        fread(e+i*4, 1, 4, file);
         uint32_t   value = ((uint32_t)e[4*i+0] << 24) | 
            ((uint32_t)e[4*i+1] << 16) |
            ((uint32_t)e[4*i+2] << 8)  |
            ((uint32_t)e[4*i+3]);
        E[64-i-1]=value;
    }
    // printf("E=");
    // cout<<cal.biToStr(E)<<endl;

    fread(skip, 1, 256, file);

    
    fread(&len, 1, 1, file);
    // printf("%d\n", len);

    uint8_t m[260];
    
    m[0]=0x00;


    uint8_t seed[32]; // 存储随机字节的数组
    unsigned int random_value; // 用于接收32位随机数的变量
    int i = 0; 


    while (i < 32) {

        if (_rdrand32_step(&random_value)) {

            seed[i++] = (uint8_t)(random_value & 0xFF); // 低8位
            seed[i++] = (uint8_t)((random_value >> 8) & 0xFF); // 次低8位
            seed[i++] = (uint8_t)((random_value >> 16) & 0xFF); // 次高8位
            seed[i++] = (uint8_t)((random_value >> 24) & 0xFF); // 高8位
        }
    }
    for(int i = 1;i < 33; i++){
        m[i]=seed[i-1];
    }

    for(int i = 33;i<65;i++){
        m[i]=hash_L[i-33];
    }

    m[256-len-1] = 0x01;
    fread(m+256-len,1,len,file);

    for(int i = 65;i< 256-len-1; i++){
        m[i]=0x00;
    }

    uint8_t DB_mask[223];
    uint8_t seed_mask[32];
    MGF1(seed,32,223,DB_mask);

    for(int i = 33;i<256; i++){
        m[i] = m[i] ^ DB_mask[i-33];
    }
    // print(m+33,223);

    MGF1(m+33,223,32,seed_mask);
    // printf("MaskedSeed:");
    // print(seed_mask,32);
    for(int i = 1;i<33;i++){
        m[i] = m[i] ^ seed_mask[i-1];
    }
    //EM
    // print(m,256);

    for(int i = 0;i < 64;i++){
        uint32_t   value = ((uint32_t)m[4*i+0] << 24) |
            ((uint32_t)m[4*i+1] << 16) |
            ((uint32_t)m[4*i+2] << 8)  |
            ((uint32_t)m[4*i+3]);
        
        // M = M*BigIntTiny(2).pow(BigIntTiny(32)) + BigIntTiny(value);
        M[64-i-1]=value;

    }
    // printf("M=");
    // cout<<cal.biToStr(M)<<endl;
    // 计算加密后的数据

    // 关闭文件

    // fclose(file);
    uint32_t res[128] = {0};
    cal.modPow(res, M, E);
    // printf("res=");
    // cout<<cal.biToStr(res)<<endl;
    // for(int i = 0;i<64;i++){
    //     printf("%u ",res[i]);
    // }

    //以小端序输出
    for(int i = 0;i<64;i++){
        uint32_t   value = res[64-i-1];
        uint8_t bytes[4];
        bytes[3] = value & 0xFF;
        bytes[2] = (value >> 8) & 0xFF;
        bytes[1] = (value >> 16) & 0xFF;
        bytes[0] = (value >> 24) & 0xFF;
        fwrite(bytes, 1, 4, outputFile);
        // printf("%08x",value);
        // for(int j = 0;j<4;j++)
        //     printf("%02x",bytes[j]);
    }

    return EXIT_SUCCESS;
}
 