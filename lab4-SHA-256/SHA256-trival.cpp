#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>

#define uint8_t unsigned char
#define uint32_t unsigned int
#define uint64_t unsigned long long
#define BLOCK_SIZE 64  // 512 bits = 64 bytes
#define BYTE unsigned char

//为什么用int不行？
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

// SHA-256压缩函数
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


// SHA-256主函数
// void sha256(const uint8_t *message, size_t length, uint8_t *digest) {
//     uint32_t state[8];
//     memcpy(state, H0, sizeof(H0));

//     for (size_t i = 0; i < length; i += BLOCK_SIZE) {
//         sha256_compress(state, &message[i]);
//     }

//     // 输出摘要
//     for (int i = 0; i < 8; i++) {
//         digest[i * 4] = (state[i] >> 24) & 0xFF;
//         digest[i * 4 + 1] = (state[i] >> 16) & 0xFF;
//         digest[i * 4 + 2] = (state[i] >> 8) & 0xFF;
//         digest[i * 4 + 3] = state[i] & 0xFF;
//     }
// }

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

void sha256_padding(int file_size, uint8_t *buffer, size_t *buffer_len) {
    // 计算文件的总长度（以字节为单位）     
    
    // 计算填充后的总长度
    size_t total_bits = file_size * 8;
    size_t padding_bits = 448 - (total_bits % 512);
    if (padding_bits <= 0) {
        padding_bits += 512;
    }
    size_t padding_bytes = padding_bits / 8;

    // 
    size_t bytes_read = file_size%BLOCK_SIZE;
    
    // fread(buffer, 1, bytes_read, file);
    // 在末尾添加一个1
    buffer[bytes_read] = 0x80;
    bytes_read++;

    // 填充0直到达到512位的倍数
    for (size_t i = bytes_read; i < BLOCK_SIZE - 8; i++) {
        buffer[i] = 0x00;
    }

    // 添加原始消息的长度（64位）
    uint64_t original_length = total_bits;
    for (int i = 0; i < 8; i++) {
        buffer[BLOCK_SIZE -8 + i] = (original_length >> (56 - i * 8)) & 0xFF;
    }

    *buffer_len = file_size + padding_bytes;
}


void sha256_update(BYTE data[], uint32_t state[],  uint8_t buffer[],size_t len)
{
	int i;

	for (i = 0; i < len; ++i) {
		buffer[datalen] = data[i];
		datalen++;
		if (datalen == 64) {
			// sha256_transform(ctx, ctx->data);
			sha256_compress(state, buffer);
			bitlen += 512;
			datalen = 0;
		}
	}
}


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
		memset(data, 0, 56);
	}

	// Append to the padding the total message's length in bits and transform.
	bitlen += datalen * 8;
    // printf("bitlen = %d\n",bitlen);

    // data[63] = bitlen;
	// data[62] = bitlen >> 8;
	// data[61] = bitlen >> 16;
	// data[60] = bitlen >> 24;
	// data[59] = bitlen >> 32;
	// data[58] = bitlen >> 40;
	// data[57] = bitlen >> 48;
	// data[56] = bitlen >> 56;
    for (int i = 0; i < 8; i++) {
        data[BLOCK_SIZE -8 + i] = (bitlen >> (56 - i * 8)) & 0xFF;
    }
	// sha256_transform(ctx, ctx->data);
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
}
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

    uint8_t buffer[BLOCK_SIZE];  // 缓冲区大小为64字节
    size_t buffer_len;
    size_t file_size;
    size_t bytes_read;
    // uint8_t* data = (uint32_t* )malloc(128);
    uint8_t* data = (uint8_t* )malloc(13631872);
    uint32_t state[8];
    uint8_t digest[32];
    memcpy(state, H0, sizeof(H0));

    int bytesRead = 0;

    while ((bytesRead = fread(data, 1,13631872 , file)) > 0) 
    {

        sha256_update(data,state, buffer, bytesRead);

    }
    
    sha256_final(buffer, state, digest);




    // 关闭文件
    fclose(file);




    // 输出摘要
    // printf("SHA-256 digest: ");
    // for (int i = 0; i < 32; i++) {
    //     printf("%02x", digest[i]);
    // }

    fwrite(digest, 1, 32, stdout);

    return EXIT_SUCCESS;
}
