#include <iostream>
#include <stdlib.h>
#include <map>
#include<vector>
#include <bitset>
#include <iomanip>
// #include <immintrin.h>
#include <x86intrin.h>
#define M 192
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

// #define gf vector<char>
using namespace std;
#define gf bitset<M>
#define gf2 bitset<2*M>

bitset<2* M> p("111");
// 以x^130 + x^129 + ...格式输出
void gf_print(const gf a);
// c = a + b
void gf_add(gf &c, const gf a, const gf b);
// c = a * b
void gf_mul(gf &c, const gf a, const gf b);
// c = a^2
void gf_pow2(gf &c, const gf a);
// c = a^n
void gf_pow(gf &c, const gf a, size_t n);
// c = a^(-1)
//使用扩展欧几里得算法求逆元
void gf_inv(gf &c, const gf a);
//使用费马小定理求逆元
void gf_inv1(gf &c, const gf a);

void gf_mul1(gf &c, const gf a, const gf b);
// 可以使用不同的方法实现同一种运算
void gf_pow2_naive(gf c, const gf a);
void gf_pow2_mul(gf c, const gf a);
void gf_pow2_bit(gf c, const gf a);
void gf_pow2_simd(gf c, const gf a);
void gf_mul_n(gf &c, const gf a, const gf b,int n);
void gf_mul_karatsuba(gf &c, const gf a, const gf b,int n);
// 然后快速切换
#ifndef gf_pow2
// #define gf_pow2 (gf_pow2_naive)
// #define gf_pow2 (gf_pow2_mul)
// #define gf_pow2 (gf_pow2_bit)
// #define gf_pow2 (gf_pow2_simd)
#endif
std::vector<bool> byte_to_binary(unsigned char byte) {
    std::vector<bool> binary_array(8);
    for (int i = 7; i >= 0; --i) {
        binary_array[7 - i] = (byte & (1 << i)) ? true : false;
    }
    return binary_array;
}


int main() {
#ifdef ONLINE_JUDGE
    // 在在线评测系统中，从 stdin 读取输入，输出到 stdout
    #ifdef _WIN32
    setmode(fileno(stdin), O_BINARY);
    setmode(fileno(stdout), O_BINARY);
    #endif
    #define mode 0
    #define inputFile stdin
    #define outputFile stdout
#else
// 在本地调试时，从文件读取输入和输出
    #define mode 1
    printf("debug mode\n");
    FILE *inputFile = fopen("dump1.bin", "rb");
    FILE *outputFile = fopen("out.bin", "wb");
 
#endif
    p[0] = 1;
    p[13] = 1;
    p[131] = 1;
    // if (inputFile == NULL || outputFile == NULL) {
    //     perror("File opening failed");
    //     return EXIT_FAILURE;
    // }
    int t;//运算次数

    fread(&t,1,4,inputFile);
    // printf("t=%d\n",t);
    
    for(int k = 0;k<t;k++){
        char op = 0;
        fread(&op,1,1,inputFile);
        
        char byte;
        //bitset(M) ax;
        gf ax;
        gf bx;
        gf cx(24);


        // uint64_t ax1[3]={0,0,0}, bx1[3];
        // for(int i = 0;i < 24;i++){
        //     fread(&byte,1,1,inputFile);
        //     // 将字符左移适当的位数
        //     uint64_t shifted_char = (uint64_t)byte << (8 * i)%64;
        //     printf("shift chat = %d\n",byte);

        //     // 使用按位或操作将字符设置到正确的位置
        //     ax1[i/8] |= shifted_char;
                
        // }
        // for(int i = 0;i<3;i++){
        //     printf("%d\n",ax1[i]);
        // }
        // exit(0);

        // for(int i = 0;i < 24;i++){
        //     fread(&byte,1,1,inputFile);
        //     // 将字符左移适当的位数
        //     uint64_t shifted_char = (uint64_t)byte << (8 * i);
        //     // 使用按位或操作将字符设置到正确的位置
        //     bx1[i/8] |= shifted_char;

        // }

        for(int i = 0;i < 24;i++){
            fread(&byte,1,1,inputFile);
            for(int j = 0;j < 8;j++){
                //ax.set(M-(i+1)*8 +7-j, (byte >> j) & 1);
                //ax.set(i*8 +j, (byte >> j) & 1);
                ax[i * 8 + j] = (byte >> j) & 1;
                
            }
            
        }
        
        //printf("op=%d,binary_array_size=%d\n",op,ax.size());
        for(int i = 0;i < 24;i++){
            fread(&byte,1,1,inputFile);
            for(int j = 0;j < 8;j++){
                //bx.set(M-(i+1)*8 +7-j, (byte >> j) & 1);
                //bx.set(i*8 +j, (byte >> j) & 1);
                bx[i * 8 + j] = (byte >> j) & 1;
            }
        }
        
        // cout<<ax<<endl;
        // cout<<bx<<endl;

        if(op == 0){
            gf_add(cx,ax,bx);
        }
        if(op == 1)
            gf_mul(cx,ax,bx);
            //gf_mul_karatsuba(cx,ax,bx,131);
        if(op == 2)
            gf_pow2(cx,ax);
        if(op == 3)
            gf_inv(cx,ax);

        gf_print(cx);
    }

    
    fclose(inputFile);
    fclose(outputFile);

    return 0;
}





void gf_print(gf cx){
    // int dl=8;
    // for (int i = 0; i < cx.size(); i += dl) {
    //     // 打印当前的 8 位
    //     for (int j = 0; j < dl; ++j) {
    //         if (i + j < cx.size()) {
    //             //char byte = cx[i + j] ? 1 : 0;

    //             cout<<cx[i +7- j];
    //         }
    //     }
    //      printf(" "); // 每 8 位后添加一个空格
    // }

    // printf("\n");
    for (size_t i = 0; i < cx.size() / 8; ++i) {
        unsigned byte = 0; // 初始化为 0
        for (size_t j = 0; j < 8; ++j) {
            // 从 bitset 中提取每一位并设置到 byte 中
            byte |= (cx[i * 8 + j] << j); // 从高位到低位构建字节
        }
        // 输出当前的 char 变量
        //std::cout << "Byte " << i << ": " << static_cast<int>(byte) << std::endl; // 输出 ASCII 值
        fwrite(&byte,1,1,stdout);
    }

}

void gf_add(gf &cx, gf ax, gf bx){
    for(int i = 0;i < ax.size();i++){
        cx[i]=ax[i] ^ bx[i];
    }
}



int degree(const bitset<2 * M> &a)
{
    for (int i = 2 * M - 1; i >= 0; i--)
        if (a[i] || i == 0){
            //printf("degree=%d\n",i);
            return i;
        }
    
    
    return 0;
}

bool judge(const gf2 &a)
{
    bitset<M> temp;
    bitset<2 * M> t(temp.set().to_string() + temp.to_string());
    t &= a;
    return t.any();
}



// bitset<M> mod(bitset<2*M> &a, gf2 &r)
// {
//     bitset<2 * M> x, y;
    
//     x = a;

//     while(degree(x) > 130){
//         int d = degree(x);
//         int d2 = degree(r);
//         x ^= (r << (d - d2));
//     }

//     return bitset<M>(x.to_string().substr(M));
// }

bitset<M> mod(bitset<2 * M> a) {
    bitset<M * 2> irreducible=p;

    for (int i = 2 * M - 1; i >= 131; --i) {
        if (a[i]) {
            a = a ^ (irreducible << (i - 131));
        }
    }

    bitset<M> result;
    for (int i = 0; i < 131; ++i) {
        result[i] = a[i];
    }
    return result;
}


void gf_mul_n(gf &c, gf a, gf b, int n){

    for(int i = 0; i<n; i++){
        if(b[i])
            c ^= a<<i;
    }

}
//使用karatsuba算法实现乘法优化
void gf_mul_karatsuba(gf &c, const gf a, const gf b, int n)
{
    bitset<2 * M> temp;
    if (n<=64)
    {
        // c[0] = a[0] & b[0];

        // gf_mul_n(c,a,b,n);
        __m128i a1 = _mm_set_epi64x(0, a.to_ullong());
        __m128i b1 = _mm_set_epi64x(0, b.to_ullong());
        __m128i c1 = _mm_clmulepi64_si128(a1, b1, 0x00);

        char data_array[16];
        _mm_storeu_si128((__m128i*)data_array, c1);


        // 将数组中的每个元素放入 bitset
        for (int i = 0; i < 16; ++i) {

            for(int j = 0; j<8;j++){
                c[i*8+j] = (data_array[i] >> j) & 1;
            }
            
        }
    
    // 输出 bitset 的内容


        return;
    }
    else
    {
        gf a1, a0, b1, b0;
        int n2 = n / 2 ;

        a1 = a >> n2; //取高n/2位
        for(int i = 0;i < n2;i++)
            a0[i] = a[i];
        b1 = b >> n2;
        for(int i = 0;i < n2;i++)
            b0[i] = b[i];
        gf z0;
        gf_mul_karatsuba(z0, a0, b0, n2);
        gf z1;
        gf_mul_karatsuba(z1, a1^a0, b1^b0, n-n2);
        gf z2;
        gf_mul_karatsuba(z2, a1, b1, n-n2);

        z1 = z1 ^ z2 ^ z0;

        bitset <2*M> t1(z0.to_string());
        bitset <2*M> t2(z1.to_string());
        bitset <2*M> t3(z2.to_string());
        temp = t1 ^ (t2 << n2) ^ (t3 << 2 * n2);

        
        c = mod(temp);

    }

}



void gf_mul(gf &c, const gf a, const gf b)
{

    bitset<2 * M> x;
    bitset<2 * M> temp(a.to_string());
    for (int i = 0; i < 131; i++)
    {
        if (b[i])
            x ^= temp << i;
    }
    c = mod(x);

}

// void gf_mul(gf &c, const gf a, const gf b) {
//         gf a2,a1, a0, b2,b1, b0;
//         int n2 = 64;

//         a1 = a >> n2; //取高n/2位
//         for(int i = 0;i < n2;i++){
//             a2[i] = a[i+2*n2];
//             a1[i] = a[i+n2];
//             a0[i] = a[i];
//         }
            

//         b1 = b >> n2;
//         for(int i = 0;i < n2;i++){
//             b2[i] = b[i+2*n2];
//             b1[i] = b[i+n2];
//             b0[i] = b[i];
//         }


//         gf T0 = a1 ^ a2; // T0 = A[1] ⊕ A[2]
//         gf T1 = b1 ^ b2; // T1 = B[1] ⊕ B[2]
//         gf T2 = a0 ^ a2; // T2 = A[0] ⊕ A[2]
//         gf T3 = b0 ^ b2; // T3 = B[0] ⊕ B[2]
//         gf T4 = a0 ^ a1; // T4 = A[0] ⊕ A[1]
//         gf T5 = b0 ^ b1; // T5 = B[0] ⊕ B[1]
    
//     // // 将 bitset 转换为 __m128i 类型
//     __m128i a0_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(a0.to_ullong()));
//    //printf("%s\n",a0_vec);
//    exit(0);
    
//     __m128i a1_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(a1.to_ullong()));
//     __m128i a2_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(a2.to_ullong()));
//     __m128i b0_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(b0.to_ullong()));
//     __m128i b1_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(b1.to_ullong()));
//     __m128i b2_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(b2.to_ullong()));

//     __m128i t0 = _mm_clmulepi64_si128(a0_vec, b0_vec, 0x00);
//     __m128i t1 = _mm_clmulepi64_si128(a1_vec, b1_vec, 0x00);
//     __m128i t2 = _mm_clmulepi64_si128(a2_vec, b2_vec, 0x00);
//     __m128i t3 = _mm_clmulepi64_si128(a0_vec, b1_vec, 0x00);
//     __m128i t4 = _mm_clmulepi64_si128(a1_vec, b0_vec, 0x00);

//     // 组合结果
//     __m128i t5 = _mm_xor_si128(t0, t1);
//     __m128i t6 = _mm_xor_si128(t2, t3);
//     __m128i t7 = _mm_xor_si128(t4, t5);
//     __m128i t8 = _mm_xor_si128(t6, t7);

//    bitset<2*M> t8_bitset;
//     for(int i = 0; i < 2*M; i++){
//         t8_bitset[i] = t8[i/64] & (1 << (i%64)); // 将 __m128i 转换为 bitset
//     }
//     c = mod(t8_bitset, p);
// }

void gf_pow2(gf &cx, gf ax){
    
    bitset<2 * M> x;
    for (int i = 0; i < M; i++)
        x[i * 2] = ax[i];
    cx = mod(x);
}


void gf_inv(gf &cx ,const gf a)
{

    bitset<2 * M> b, c, u, v, temp;
    bitset<M> r(p.to_string().substr(1));
    int j;
    b[0] = 1;
    u = bitset<2 * M>(a.to_string());
    v = p;
    while (degree(u))
    {
        j = degree(u) - degree(v);
        if (j < 0)
        {
            j = -j;
            temp = u;
            u = v;
            v = temp;
            temp = b;
            b = c;
            c = temp;
        }
        u = u ^ (v << j);
        b = b ^ (c << j);
    }
    cx = mod(b);
}

void gf_inv1(gf &cx ,const gf a){
    //使用费马小定理求逆元
    //x1
    gf_pow2(cx,a);
    gf_mul(cx,cx,a);
    gf x1 = cx;
    //x2
    gf_pow2(cx,cx);
    gf_pow2(cx,cx);
    gf_mul(cx,cx,x1);
    gf x2 = cx;
    //x3
    for(int i = 0; i<4;i++){
        gf_pow2(cx,cx);
    }
    gf_mul(cx,cx,x2);
    gf x3 = cx;
    //x4
    for(int i = 0; i<8;i++){
        gf_pow2(cx,cx);
    }
    gf_mul(cx,cx,x3);
    gf x4 = cx;
    //x5
    for(int i = 0; i<16;i++){
        gf_pow2(cx,cx);
    }
    gf_mul(cx,cx,x4);
    gf x5 = cx;
    //x6
    for(int i = 0; i<32;i++){
        gf_pow2(cx,cx);
    }
    gf_mul(cx,cx,x5);
    gf_pow2(cx,cx);
    gf_mul(cx,cx,a);
    gf x6 = cx;

    //x7

    for(int i = 0; i<65;i++){
        gf_pow2(cx,cx);
    }
    gf_mul(cx,cx,x6);
    //x7^2
    gf_pow2(cx,cx);

    return;



}

void gf_pown(gf &cx, gf base, int power){
    gf result;
    result[0]=1;   //用于存储项累乘与返回最终结果，由于要存储累乘所以要初始化为1
    while (power > 0)           //指数大于0说明指数的二进制位并没有被左移舍弃完毕
    {
        if (power & 1)          //指数的当前计算二进制位也就是最末尾的位是非零位也就是1的时候
            gf_mul(result, result, base);         //例如1001的当前计算位就是1， 100*1* 星号中的1就是当前计算使用的位
            //累乘当前项并存储
        gf_pow2(base, base);          //计算下一个项，例如当前是n^2的话计算下一项n^2的值
                                //n^4 = n^2 * n^2;
        power >>= 1;            //指数位右移，为下一次运算做准备
                                //一次的右移将舍弃一个位例如1011(2)一次左移后变成101(2)
    }
    cx = result;              //返回最终结果

}



