#include "dtw.hpp"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>

using namespace std;
DTW::DTW(int n_mfcc,const char *filename,float min_distance)
{
    bool ps;
    ps = true;
    ref = NULL;
    n_col = 0;
    n_rows = 0;
    n_distance = 0.f;
    FILE *fp = fopen(filename,"r");
    if(fp != NULL) {
        char *s = new char[10000];
        DataTable *dt1;
        DataTable *dt = new DataTable;
        DataTable *tdt = dt;
        //ref = dt;
        while(!feof(fp)) {
            char *sp = fgets(s,10000,fp);
            if(sp == NULL) {
                break;
            }
            int dlen = n_mfcc;
            char *tokenPtr=strtok(s,",");
            if(tokenPtr != NULL) {
                dt1 = new DataTable;
                dt1->prec = dt;
                dt->next = dt1;
                dt = dt1;
                dt->v_mfcc = new float[n_mfcc];
                n_rows++;
            }
            float *p_dt = dt->v_mfcc;
            while(tokenPtr != NULL) {
                if(dlen == 0) {
                    break;
                }
                dlen--;
                *p_dt = atof(tokenPtr);
                p_dt++;
                tokenPtr = strtok(NULL,",");
            }
            if(dlen > 0) {
                ps = false;
                break;
            }
        }
        delete s;
        ref_end = dt1;
        ref_end->next = NULL;
        ref = tdt->next;
        delete tdt;
        if(ps) {
            // dt = new DataTable;
            // tdt = dt;
            // dt->next = test;
            n_col = n_mfcc;
            // for(int i=0;i<n_rows;i++) {
                // dt1 = new DataTable;
                // dt1->prec = dt;
                // dt->next = dt1;
                // dt = dt1;
                // dt->data = NULL;
            // }
            // test_end = dt1;
            // test_end->next = NULL;
            // dt->prec = NULL;
            // dt = ref->next;
            // delete ref;
            //test->prec = NULL;
            distance = new float *[n_rows];
            for(int i=0;i<n_rows;i++) {
                distance[i] = NULL;
            }
            min_sdistance = min_distance;
        }
        else {
            n_rows = 0;
            distance = NULL;
            ref = NULL;
            // test = NULL;
            ref_end = NULL;
            min_sdistance = 0.f;
            //test_end = NULL;
        }
        fclose(fp);
        // dt1 = ref;
        // while(dt1 != NULL) {
            // for(int i=0;i<12;i++) {
                // printf("%f,",dt1->v_mfcc[i]);
            // }
            // printf("\n");
            // dt1=dt1->next;
        // }
    }
}
DTW::~DTW() {
    DataTable *dt;
    
    while(ref != NULL) {
        dt = ref;
        ref = ref->next;
        delete dt->v_mfcc;
        delete dt;
    }
    // while(test != NULL) {
        // dt = test;
        // test = test->next;
        // delete dt->data;
        // delete dt;
    // }
    if(distance != NULL) {
        for(int i=0;i<n_rows;i++) {
            if(distance[i] != NULL) {
                delete distance[i];
            }
        }
        delete distance;
        distance = NULL;
    }
    n_rows = 0;
    n_col = 0;
}
#if 0
void DTW::AddData(float *in) {
    DataTable *dt = test;
    if(dt == NULL) {return;}
    if(n_col == 0) {return;}
    test = test->next;
    delete dt;
    test_end->next = new DataTable;
    test_end = test_end->next;
    float data = new float[n_col];
    test_end->v_mfcc = data;
    memcpy(data,in,sizeof(float)*n_col);
    float *p_distance = distance[0];
    if(p_distance != NULL) {
        DataTable *dr = ref;
        dt = test;
        distance[0][0] = getDistance(dr,dt);
        dt = dt->next;
        for(int i=1;i<n_rows;i++) {
            distance[0][i] = distance[0][i-1] + getDistance(dr,dt);
            dt = dt->next;
        }
        dt = test;
        dr = dr->next;
        for(int i=1;i<n_rows;i++) {
            distance[i][0] = distance[i-1][0] + getDistance(dr,dt);
            dt = dt->next;
        }
        dr = ref;
        for(int i=1;i<n_rows;i++) {
            dr = dr->next;
            dt = test;
            for(int j=1;j<n_rows;j++) {
                dt = dt->next;
                distance[i][j] = distance[i-1][j-1] + getDistance(dr,dt);
            }
        }
    }
}
#endif

bool DTW::AddData(float *in) {
    float ds_current,ds_next;
    //memcmp(tp,in,sizeof(float)*n_col);
    //if(dt == NULL) {return;}
    if(n_col == 0) {return false;}
    ds_current = getDistance(p_ref->v_mfcc,in);
    DataTable *db_next = p_ref->next;
    if(db_next != NULL) {
        ds_next = getDistance(db_next->v_mfcc,in);
        if((ds_next > min_sdistance)&&(ds_current > min_sdistance)){
            Start();
            return false;
        }
        if(ds_next < ds_current) {
            p_ref = db_next;
            s_distancd += n_distance/cn_distance;
            cs_distance ++;
            cn_distance = 1;
            n_distance = ds_next;
        }
        else {
            //ds_current = ds_current;
            n_distance += ds_current;
            cn_distance ++;
        }
    }
    else{
        if(ds_current > min_sdistance){
            Start();
            return false;
        }
    }
    if(getDistance() < min_sdistance) {
        return true;
    }
    else{
        return false;
    }
}

int DTW::Start() {
    if(ref == NULL) return -1;
    
    p_ref = ref;
    n_distance = 0.f;
    cn_distance = 0;
    s_distancd = 0.f;
    cs_distance = 0;
    return 0;
}

float DTW::getDistance(float *d1,float *d2) {
    float sum = 0;
    for(int i=0;i<n_col;i++) {
        sum += fabs(*d1-*d2);
        d1++;
        d2++;
    }
    return sum;
}

float DTW::getDistance() {
    if(p_ref->next == NULL) {
        if(cn_distance > 0) {
            s_distancd += n_distance/cn_distance;
            cs_distance ++;
        }
        return s_distancd /cs_distance;
    }
    return min_sdistance + 10;
}

bool DTW::MakeRef(float *inout) {
    if(ref == NULL) {
        ref = new DataTable;
        ref->next = NULL;
        ref->prec = NULL;
        ref->v_mfcc = NULL;
        ref_end = ref;
    }
    if(ref_end->v_mfcc == NULL) {
        ref_end->v_mfcc = new float[n_col];
        memcpy(ref_end->v_mfcc,inout,n_col*sizeof(float));
    }
    else{
        DataTable *tref = new DataTable;
        memcpy(tref->v_mfcc,inout,n_col*sizeof(float));
        tref->next = NULL;
        float sub = getDistance(tref->v_mfcc,ref_end->v_mfcc);
        
        if(sub > min_sdistance) {
            int refsize = 0;
            ref_end = ref;
            for(int i=0;i<n_col;i++) {
                inout[i] = 0.f;
            }
            while(ref_end != NULL) {
                for(int i=0;i<n_col;i++) {
                    inout[i] += ref_end->v_mfcc[i];
                }
                ref_end = ref_end->next;
                delete ref;
                ref = ref_end;
                refsize ++;
            }
            for(int i=0;i<n_col;i++) {
                inout[i] = inout[i]/refsize;
            }
            ref = tref;
            ref->prec = NULL;
            ref_end = ref;
            return true;
        }
        else{
            ref_end->next = tref;
            tref->prec = ref_end;
            ref_end = tref;
        }
    }
    return false;
    //memcmp(ref->v_mfcc,in);
}
