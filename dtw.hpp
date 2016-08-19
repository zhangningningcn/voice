#ifndef _DTW_HPP_
#define _DTW_HPP_

typedef struct DataTable_ {
    float *v_mfcc; 
    struct DataTable_ *next;
    struct DataTable_ *prec;
}DataTable;
// typedef struct Distance_ {
    // float *v_mfcc; 
    // float val;
// }Distance;

class DTW {
public:
    DTW(int n_mfcc,const char *filename,float min_distance = 2.f);
    ~DTW();
    int Start();
    float getDistance();
    bool AddData(float *in);
    bool MakeRef(float *inout);
private:
    DataTable *ref;
    DataTable *ref_end;
    //Distance *distance;
    //Distance *distance_end;
    int n_col;
    int n_rows;
    float **distance;
    DataTable *p_ref;
    //���ڴ洢��ǰֵ
    float n_distance;
    //��ǰֵ����
    int cn_distance;
    //�ۼ�ֵ
    float s_distancd;
    int cs_distance;
    float Max_diff;
    float min_sdistance;
    
    float getDistance(float *d1,float *d2);
};

#endif