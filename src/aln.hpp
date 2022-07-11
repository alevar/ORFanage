#include <iostream>
#include "ksw2.h"
#include <algorithm>
#include <cstring>
#include <vector>
#include <fstream>
#include <map>
#include <math.h>

#include "common.hpp"

#define CMATCH  0
#define CINS    1
#define CDEL    2
#define CREF_SKIP   3
#define CSOFT_CLIP  4
#define CHARD_CLIP  5
#define CPAD    6
#define CEQUAL  7
#define CDIFF   8

class Finder{
public:
    Finder() = default;
    Finder(const int8_t *mat,const int8_t *alphabet,int gapo,int gape):mat(mat),alphabet(alphabet),gapo(gapo),gape(gape){};
    ~Finder() = default;
    void align(const char *tseq, const char *qseq, ksw_extz_t& ez){ // from the ksw2 example
        memset(&ez, 0, sizeof(ksw_extz_t));

        int tl = strlen(tseq), ql = strlen(qseq);
        uint8_t *ts, *qs;

        ts = (uint8_t*)malloc(tl);
        qs = (uint8_t*)malloc(ql);
        for (int i = 0; i < tl; ++i) ts[i] = alphabet[(uint8_t)tseq[i]]; // encode to 0/1/2/3
        for (int i = 0; i < ql; ++i) qs[i] = alphabet[(uint8_t)qseq[i]];
        const int flag = KSW_EZ_EXTZ_ONLY;
        ksw_extz2_sse(0, ql, qs, tl, ts, 5, this->mat, gapo, gape, -1, -1, 0, flag, &ez);
        free(ts);
        free(qs);
    }
    int parse(ksw_extz_t& ez,std::string& q, std::string& t, Score& aln){
        // first find the position of the last matched base on the query
        int last_m_q = ez.max_q;
        for (uint16_t c = ez.n_cigar - 1; c >= 0; --c){
            int opcode =  ez.cigar[c] & 0xf;
            int oplen = ez.cigar[c] >> 4;

            if (opcode == CDEL)
                continue;
            else if (opcode == CINS)
                last_m_q-=oplen;
            else if (opcode == CMATCH)
                break;
            else{
                std::cerr << "unidentified cigar operation\n";
                exit(3);
            }
        }

        // need the template coordinate of the first match
        int cur_pos_q = 0;
        int cur_pos_t = 0;

        for (uint16_t c = 0;c < ez.n_cigar; ++c){
            int opcode = ez.cigar[c] & 0xf;
            int oplen = ez.cigar[c] >> 4;

            if (opcode == CMATCH)
                break;
            else if (opcode == CINS)
                continue;
            else if (opcode == CDEL)
                cur_pos_t += oplen;
            else{
                std::cerr << "unidentified cigar operation\n";
                exit(3);
            }
        }

        int first_match_pos_t = cur_pos_t;
        int last_match_pos_t = ez.max_t;
        if (last_match_pos_t - first_match_pos_t < 1)
            return -1;


        // now can begin parsing
        cur_pos_q = 0;
        cur_pos_t = 0;

        for (uint16_t c = 0; c < ez.n_cigar; ++c){
            int opcode = ez.cigar[c] & 0xf;
            int oplen = ez.cigar[c] >> 4;

            for(int i = 0; i < oplen; ++i){
                if (opcode == CMATCH){
                    if (q[cur_pos_q] != t[cur_pos_t]){
                        aln.snvs.emplace_back(cur_pos_t, q[cur_pos_q], t[cur_pos_t]);
                        aln.cigar.emplace_back(1, 'Z');
                    }
                    else{
                        if (!aln.cigar.empty() && aln.cigar.back().second == CMATCH)
                            ++aln.cigar.back().first;
                        else
                            aln.cigar.emplace_back(1, CMATCH);
                    }
                    ++cur_pos_q;
                    ++cur_pos_t;
                }
                else if (opcode == CINS){
                    if (!aln.cigar.empty() && aln.cigar.back().second == CINS){
                        ++aln.cigar.back().first;
                        aln.inss.back().second += q[cur_pos_q];
                    }
                    else{
                        aln.cigar.emplace_back(1, CINS);
                        aln.inss.emplace_back(cur_pos_t, std::string(1, q[cur_pos_q]));
                    }
                    ++cur_pos_q;
                }
                else if (opcode == CDEL){
                    if (!aln.cigar.empty() && aln.cigar.back().second == CDEL){
                        ++aln.cigar.back().first;
                        aln.dels.back().second += t[cur_pos_t];
                    }
                    else{
                        aln.cigar.emplace_back(1, CDEL);
                        aln.dels.emplace_back(cur_pos_t, std::string(1, t[cur_pos_t]));
                    }
                    ++cur_pos_t;
                }
                else{
                    std::cerr << "unidentified cigar operation\n";
                    exit(3);
                }
            }
        }
        return 0;
    }
    void print_cigar(ksw_extz_t& ez){
        for (int i = 0; i < ez.n_cigar; ++i){
            printf("%d%c", ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
        }
        putchar('\n');
    }

protected:
    const int8_t* mat;
    int gapo = 12, gape = 2;
    const int8_t *alphabet;
};

