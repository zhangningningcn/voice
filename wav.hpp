#ifndef _WAV_HPP
#define _WAV_HPP

#include "cstdio"

typedef unsigned int uint32_t;
typedef unsigned short int uint16_t;

enum ERRORNO{
    EOPENFILE=1,
    EREADFILE,
    ERIFFFLG,
    EWAVEFLG,
    EFMTFLG,
    EFMTLEN,
    ETAG,
    ESAMPFREQ,
    EBITSAMP,
    EDATASTAMP,
    
};

class WAV{
    public:
    WAV();
    ~WAV();
    int Open(const char *filename);
    int ReadBlock(float *data,int channel);
    int getChannels(void);
    int getSampleRate(void);
    private:
    struct wavinfo{
        char riff[4];
        uint32_t fsize;
        char wave[4];
        char fmt[4];
        uint32_t fmt_len; //过滤字节(一般为00000010H)
        uint16_t tag;
        uint16_t channels;
        uint32_t samp_freq;
        uint32_t byte_rate;
        uint16_t block_align;  //块对齐字节数 = channles * bit_samp / 8  
        uint16_t bit_samp;
        char datastamp[4];
        uint32_t datalength;
        void * data;
    }wavinfo;
    int framlen;
    int channels;
    int block_align;
    int byte_samp;
    FILE *fp;
};

#endif