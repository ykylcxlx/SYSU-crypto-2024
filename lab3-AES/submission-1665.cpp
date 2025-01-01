
#include<stdio.h>
#include <x86intrin.h>
#include<stdlib.h>
#include <fcntl.h>
#define uint8_t unsigned char
class AES128_NI
{
public:
    // AES128
    static const int Nk = 4;
    static const int Nb = 4;
    static const int Nr = 10;

protected:
    __m128i w128[Nr + 1];
    __m128i dw128[Nr + 1];
    bool decrypting;

public:
    AES128_NI(const unsigned char key[4 * Nk], bool decrypting)
        : decrypting(decrypting)
    {
        KeyExpansion(key);
    }

    void Cipher(const unsigned char in[4 * Nb], unsigned char out[4 * Nb]);
    void EqInvCipher(const unsigned char in[4 * Nb], unsigned char out[4 * Nb]);

protected:
    void KeyExpansion(const unsigned char key[4 * Nk]);
};

// #include "AES128_NI.h"

//Cipher(byte in[4*Nb], byte out[4*Nb], word w[Nb*(Nr+1)])
//begin
void AES128_NI::Cipher(const unsigned char in[4 * Nb], unsigned char out[4 * Nb])
{
    //  ASSERT(!decrypting);

    //state = in
    __m128i state = _mm_loadu_si128(reinterpret_cast<const __m128i *>(in));

    //AddRoundKey(state, w[0, Nb-1])
    // just XOR
    state = _mm_xor_si128(state, w128[0]);

    //for round = 1 step 1 to Nr-1
        //SubBytes(state)
        //ShiftRows(state)
        //MixColumns(state)
        //AddRoundKey(state, w[round*Nb, (round+1)*Nb-1])
    //end for
    state = _mm_aesenc_si128(state, w128[1]);
    state = _mm_aesenc_si128(state, w128[2]);
    state = _mm_aesenc_si128(state, w128[3]);
    state = _mm_aesenc_si128(state, w128[4]);
    state = _mm_aesenc_si128(state, w128[5]);
    state = _mm_aesenc_si128(state, w128[6]);
    state = _mm_aesenc_si128(state, w128[7]);
    state = _mm_aesenc_si128(state, w128[8]);
    state = _mm_aesenc_si128(state, w128[9]);

    // The last round
    //SubBytes(state)
    //ShiftRows(state)
    //AddRoundKey(state, w[Nr*Nb, (Nr+1)*Nb-1])
    state = _mm_aesenclast_si128(state, w128[Nr]);

    //out = state
    _mm_storeu_si128(reinterpret_cast<__m128i *>(out), state);

    //end
}

//EqInvCipher(byte in[4*Nb], byte out[4*Nb], word dw[Nb*(Nr+1)])
//begin
void AES128_NI::EqInvCipher(const unsigned char in[4 * Nb], unsigned char out[4 * Nb])
{
    //  ASSERT(decrypting);

    //byte  state[4,Nb]
    //state = in
    __m128i state = _mm_loadu_si128(reinterpret_cast<const __m128i *>(in));

    //AddRoundKey(state, dw[Nr*Nb, (Nr+1)*Nb-1])
    state = _mm_xor_si128(state, w128[Nr]);

    //for round = Nr-1 step -1 downto 1
        //InvSubBytes(state)
        //InvShiftRows(state)
        //InvMixColumns(state)
        //AddRoundKey(state, dw[round*Nb, (round+1)*Nb-1])
    //end for
    state = _mm_aesdec_si128(state, dw128[9]);
    state = _mm_aesdec_si128(state, dw128[8]);
    state = _mm_aesdec_si128(state, dw128[7]);
    state = _mm_aesdec_si128(state, dw128[6]);
    state = _mm_aesdec_si128(state, dw128[5]);
    state = _mm_aesdec_si128(state, dw128[4]);
    state = _mm_aesdec_si128(state, dw128[3]);
    state = _mm_aesdec_si128(state, dw128[2]);
    state = _mm_aesdec_si128(state, dw128[1]);

    //InvSubBytes(state)
    //InvShiftRows(state)
    //AddRoundKey(state, dw[0, Nb-1])
    state = _mm_aesdeclast_si128(state, dw128[0]);

    //out = state
    _mm_storeu_si128(reinterpret_cast<__m128i *>(out), state);
    //end
}

//static const BYTE Rcon[] = {
//  0x8d, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a,
//};

//KeyExpansion(byte key[4*Nk], word w[Nb*(Nr+1)], Nk)
//begin
void AES128_NI::KeyExpansion(const unsigned char key[4 * Nk])
{
    //word  temp
    //i = 0

    //while (i < Nk)
        //w[i] = word(key[4*i], key[4*i+1], key[4*i+2], key[4*i+3])
        //i = i+1
    //end while
    __m128i work = _mm_loadu_si128(reinterpret_cast<const __m128i*>(key));  // work = w3 : w2 : w1 : w0
    w128[0] = work;

    //i = Nk
    //while (i < Nb * (Nr+1)]
        //temp = w[i-1]
        //if (i mod Nk = 0)
            //temp = SubWord(RotWord(temp)) xor Rcon[i/Nk]
        //else if (Nk > 6 and i mod Nk = 4)
            //temp = SubWord(temp)
        //end if
        //w[i] = w[i-Nk] xor temp
        //i = i + 1
    //end while


    __m128i t, t2;

    // f(w) = SubWord(RotWord(w)) xor Rcon

    // work = w3 : w2 : w1 : w0
    // w4 = w0^f(w3)
    // w5 = w1^w4 = w1^w0^f(w3)
    // w6 = w2^w5 = w2^w1^w0^f(w3)
    // w7 = w3^w6 = w3^w2^w1^w0^f(w3)
#define EXPAND(n, RCON) \
    t  = _mm_slli_si128(work, 4);       /* t    = w2                : w1             : w0          : 0        */\
    t  = _mm_xor_si128(work, t);        /* t    = w3^w2             : w2^w1          : w1^w0       : w0       */\
    t2 = _mm_slli_si128(t, 8);          /* t2   = w1^w0             : w0             : 0           : 0        */\
    t  = _mm_xor_si128(t, t2);          /* t    = w3^w2^w1^w0       : w2^w1^w0       : w1^w0       : w0       */\
    work = _mm_aeskeygenassist_si128(work, RCON);           /* work = f(w3) : - : - : - */                      \
    work = _mm_shuffle_epi32(work, 0xFF);/*work = f(w3)             : f(w3)          : f(w3)       ; f(w3)    */\
    work = _mm_xor_si128(t, work);      /* work = w3^w2^w1^w0^f(w3) : w2^w1^w0^f(w3) : w1^w0^f(w3) : w0^f(w3) */\
    w128[n] = work;                     /* work = w7 : w6 : w5 : w4 */

    EXPAND(1, 0x01);

    // Go on...
    EXPAND(2, 0x02);
    EXPAND(3, 0x04);
    EXPAND(4, 0x08);
    EXPAND(5, 0x10);
    EXPAND(6, 0x20);
    EXPAND(7, 0x40);
    EXPAND(8, 0x80);
    EXPAND(9, 0x1b);
    EXPAND(10, 0x36);

    // Additional process for EqInvCipher
    if (decrypting) {
        //for i = 0 step 1 to (Nr+1)*Nb-1
            //dw[i] = w[i]
        //end for
        dw128[0] = w128[0];

        //for round = 1 step 1 to Nr-1
            //InvMixColumns(dw[round*Nb, (round+1)*Nb-1])
        //end for
        dw128[1] = _mm_aesimc_si128(w128[1]);
        dw128[2] = _mm_aesimc_si128(w128[2]);
        dw128[3] = _mm_aesimc_si128(w128[3]);
        dw128[4] = _mm_aesimc_si128(w128[4]);
        dw128[5] = _mm_aesimc_si128(w128[5]);
        dw128[6] = _mm_aesimc_si128(w128[6]);
        dw128[7] = _mm_aesimc_si128(w128[7]);
        dw128[8] = _mm_aesimc_si128(w128[8]);
        dw128[9] = _mm_aesimc_si128(w128[9]);

        dw128[Nr] = w128[Nr];
    }
    //end
}

void print(uint8_t * text,int len,int mode){
    if(mode==0){
        fwrite(text,1,len,stdout);
    }else{
        for(int i=0;i<len;i++){
        printf("%02X ",text[i]);
    }
    printf("\n");
    }
    
}
int main(){
        int mode=0;
    #ifdef ONLINE_JUDGE
    // 在在线评测系统中，从 stdin 读取输入，输出到 stdout
    #ifdef _WIN32
        setmode(fileno(stdin), O_BINARY);
        setmode(fileno(stdout), O_BINARY);

    #endif
    
    #define inputFile stdin
    #define outputFile stdout

#else
// 在本地调试时，从文件读取输入和输出
    mode=1;
    printf("debug mode\n");
    FILE *inputFile = fopen("dump1.bin", "rb");
    FILE *outputFile = fopen("out.bin", "wb");

#endif
    int Nk=AES128_NI::Nk;
    uint8_t op;
    uint8_t key[4*Nk],iv[4*Nk];
    int len;//plainText的长度

    fread(&op, 1, 1, inputFile);

    for(int i = 0; i < 4*Nk; i++)
            fread(&key[i], 1, 1, inputFile);
    for(int i = 0; i < 4*Nk; i++)
            fread(&iv[i], 1, 1, inputFile);      
    
    // for(int i = 0; i < 4*Nk; i++)
    //         printf("%02x ",key[i]);
    // printf("\n");
    fread(&len, 4, 1, inputFile);
    
    int true_len,pad_len;
    if(len % 16 != 0){
        true_len = (len / 16 + 1) * 16;
        pad_len = true_len - len;
    }
    else{
        true_len = len+16;
        pad_len = 16;
    }
    uint8_t state[4*Nk];



    unsigned char * plain_text= (unsigned char*)malloc(16*sizeof(unsigned char));
    unsigned char * chipher_text = (unsigned char*)malloc(16*sizeof(unsigned char));
    if(op ==0x01){
        AES128_NI aes128(key, 0);
        for(int i=0;i<4*Nk;i++){
            state[i]=iv[i];
        }
       for(int i=0;i<true_len;i+=16){

        // getChars(plain_text,len,pad_len,inputFile);
        if(i==true_len-16){
            if(len%16==0){
                for(int j = 0;j<16;j++){
                    plain_text[j]=0x10;
                }
            }else if(len%16!=0){
                fread(plain_text, 1, len%16, inputFile);
                for(int j=len%16;j<16;j++){
                    plain_text[j]=(unsigned char)pad_len;
            }
            }

        }else{
            fread(plain_text, 1, 16, inputFile);
        }

        for(int i=0;i<16;i++)
            plain_text[i] ^= state[i];
    
        aes128.Cipher(plain_text, state);
        // for(int i=0;i<16;i++){
        //         printf("%02x ",state[i]);
        // }
        // printf("\n");
        fwrite(state, 1, 16, stdout);
    
    }
    }else if(op ==0x81){
        AES128_NI aes128(key, 1);
        unsigned char y[4*Nk];
        for(int i=0;i<4*Nk;i++){
            y[i]=iv[i];
        }

        for(int i=0;i<len;i+=16){
            fread(chipher_text, 1, 16, inputFile);


            // for(int i=0;i<16;i++)
            //     printf("%02x ",chipher_text[i]);
            // printf("\n");
            aes128.EqInvCipher(chipher_text, state);
            for(int i=0;i<16;i++)
                state[i] ^= y[i];
            //fwrite(state, 1, 16, outputFile);
            
            if(i!=len-16){
                fwrite(state, 1, 16, stdout);
                //print(state,16,1);
            }else{
                int un_pad_len = 16-(int)state[15];
                    fwrite(state, 1, un_pad_len, stdout);
                    //print(state,un_pad_len,1);
            }
            for(int i=0;i<4*Nk;i++)
                y[i]=chipher_text[i];
        }
    }
    return 0;
}