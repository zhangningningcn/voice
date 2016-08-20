#include "cstdio"
#include "cstring"
#include "wav.hpp"

WAV::WAV() {
    fp = NULL;
}
WAV::~WAV() {
    if(fp != NULL) {
        fclose(fp);
    }
}
int WAV::Open(const char *filename) {
    //int rlen;
    fp = fopen(filename,"r");
    if(fp == NULL) {
        return -EOPENFILE;
    }
    //rlen = fread(wavinfo.riff,sizeof(char),4,fp);
    if(4 != fread(wavinfo.riff,sizeof(char),4,fp)) {
        return -EREADFILE;
    }
    if(memcmp(wavinfo.riff,"RIFF",sizeof(char)*4) != 0) {
        return -ERIFFFLG;
    }
    if(1 != fread(&wavinfo.fsize,sizeof(uint32_t),1,fp)) {
        return -EREADFILE;
    }
    if(4 != fread(wavinfo.wave,sizeof(char),4,fp)) {
        return -EREADFILE;;
    }
    if(memcmp(wavinfo.wave,"WAVE",sizeof(char)*4) != 0) {
        return -EWAVEFLG;
    }
    if(4 != fread(wavinfo.fmt,sizeof(char),4,fp)) {
        return -EREADFILE;
    }
    if(memcmp(wavinfo.fmt,"fmt ",sizeof(char)*4) != 0) {
        return -EFMTFLG;
    }
    
    if(1 != fread(&wavinfo.fmt_len,sizeof(uint32_t),1,fp)) {
        return -EREADFILE;
    }
    if(wavinfo.fmt_len != 0x10) {
        return -EFMTLEN;
    }
    if(1 != fread(&wavinfo.tag,sizeof(uint16_t),1,fp)) {
        return -EREADFILE;
    }
    if(wavinfo.fmt_len != 0x10) {
        return -ETAG;
    }
    if(1 != fread(&wavinfo.channels,sizeof(uint16_t),1,fp)) {
        return -EREADFILE;
    }
    if(1 != fread(&wavinfo.samp_freq,sizeof(uint32_t),1,fp)) {
        return -EREADFILE;
    }
    if((wavinfo.samp_freq > 25938)||(wavinfo.samp_freq < 4224)) {
        return -ESAMPFREQ;
    }
    if(wavinfo.samp_freq > 13672) {
        framlen = 512;
    }
    else{
        framlen = 128;
    }
    if(1 != fread(&wavinfo.byte_rate,sizeof(uint32_t),1,fp)) {
        return -EREADFILE;
    }
    if(1 != fread(&wavinfo.block_align,sizeof(uint16_t),1,fp)) {
        return -EREADFILE;
    }
    block_align = wavinfo.block_align;
    if(1 != fread(&wavinfo.bit_samp,sizeof(uint16_t),1,fp)) {
        return -EREADFILE;
    }
    
    byte_samp = wavinfo.bit_samp / 8;
    if(wavinfo.bit_samp % 8 != 0) {
        return -EBITSAMP;
    }
    if(4 != fread(&wavinfo.datastamp,sizeof(char),4,fp)) {
        return -EDATASTAMP;
    }
    if(1 != fread(&wavinfo.datalength,sizeof(uint32_t),1,fp)) {
        return -EREADFILE;
    }
    
    return framlen;
        // uint32_t bits;
        
    
}
int WAV::getChannels(void) {
    return wavinfo.channels;
}
int WAV::getSampleRate(void) {
    return wavinfo.samp_freq;
}
int WAV::ReadBlock(float *data,int channel) {
    char databuffer[block_align];
    char *p_data;
    if(channel >= channels){
        return 0;
    }
    int i;
    int temp;
    //int bytes;
    for(i=0;i<framlen;i++) {
        if( (unsigned int)block_align != fread(databuffer,1,block_align,fp) ) {
            break;
        }
        p_data = databuffer + channel*byte_samp;
        temp = 0;
        memcpy(&temp,p_data,byte_samp);
        if(byte_samp == 2) {
            if(temp & 0x8000) {
                temp |= 0xFFFF0000;
            }
        }
        *data = (float)temp/32768;
        data++;
    }
    if(i == 0) {
        return 0;
    }
    for(;i<framlen;i++) {
        *data = 0.f;
        data++;
    }
    return block_align;
}