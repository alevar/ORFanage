#ifndef COMMON_H
#define COMMON_H

#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <fstream>

#include <string>
#include <unordered_map>
#include <map>
#include <array>

#define OUT_INFO    "\033[1;34m"
#define OUT_ERROR   "\033[1;31m"
#define OUT_RESET   "\033[0m"
#define OUT_DEL     "\x1b[K"
#define OUT_BOLD    "\033[1m"

#define CMATCH  0
#define CINS    1
#define CDEL    2
#define CREF_SKIP   3
#define CSOFT_CLIP  4
#define CHARD_CLIP  5
#define CPAD    6
#define CEQUAL  7
#define CDIFF   8

struct Score{
    // general stats which can be computed from the simplest interval comparison
    int qlen = 0;
    int tlen = 0;
    int ulen = 0; // length of the union of query and tempalte chains

    bool start_match = 0;
    bool stop_match = 0;

    int num_bp_extra = 0;
    int num_bp_missing = 0;
    int num_bp_match = 0;
    int num_bp_outframe = 0;
    int num_bp_inframe = 0;

    // alignment data
    std::vector<std::tuple<int,uint8_t,uint8_t> > snvs;
    std::vector<std::pair<int,std::string> > inss;
    std::vector<std::pair<int,std::string> > dels;
    std::vector<std::pair<int,uint8_t> > cigar;
    int aln_match = 0;

    // summary based on interval arithmetic
    float lpi = 0.0; // percent identity by length between the original and reference transcripts.
    float ilpi = 0.0; // percent identity by length of bases in frame of the reference transcript.
    float mlpi = 0.0; // percent identity by length of bases that are in both query and reference.
    float aln_pi = 0.0; // percent identity

    // PhyloCSF++
    float ppp_score = 0.0;

    // TOTALS
    bool pass = true;

    float num_aln_match(){
        int res = 0;
        for(auto& c : this->cigar){
            if(c.second==CMATCH){
                res+=c.first;
            }
        }
        return res;
    };
    float pi(){ // compute percent identity
        // compute number of matching bases
        return 100.0*((float)this->num_aln_match()/((float)this->ulen/3.0));
    };
    std::string stats_header(){
        std::string res =  "query_len\t"
                           "template_len\t"
                           "union_len\t"
                           "pass\t"
                           "len_match\t"
                           "len_inframe\t"
                           "len_outframe\t"
                           "len_extra\t"
                           "len_missing\t"
                           "length_pi\t"
                           "match_length_pi\t"
                           "inframe_length_pi\t"
                           "alignment_match\t"
                           "start_match\t"
                           "stop_match\t"
                           "pi";
        return res;
    };
    friend std::ostream& operator<<(std::ostream& os, const Score& s){
        os << s.qlen << "\t"
           << s.tlen << "\t"
           << s.ulen << "\t"
           << s.pass << "\t"
           << s.num_bp_match << "\t"
           << s.num_bp_inframe << "\t"
           << s.num_bp_outframe << "\t"
           << s.num_bp_extra << "\t"
           << s.num_bp_missing << "\t"
           << s.lpi << "\t"
           << s.mlpi << "\t"
           << s.ilpi << "\t"
           << s.aln_match << "\t"
           << s.start_match << "\t"
           << s.stop_match << "\t"
           << s.aln_pi;

        return os;
    };
};

static constexpr char    table_nt[4]   = { 'A', 'C', 'G', 'T' };

/* This table is used to transform nucleotide letters into numbers. */
static constexpr int8_t nt_table[128] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static constexpr bool IUPAC[256] =
        { 0, 0, 0, 0, 0, 0, 0, 0,   // 0
          0, 0, 0, 0, 0, 0, 0, 0,   // 8
          0, 0, 0, 0, 0, 0, 0, 0,   // 16
          0, 0, 0, 0, 0, 0, 0, 0,   // 24
          0, 0, 0, 0, 0, 0, 0, 0,   // 32
          0, 0, 0, 0, 0, 0, 0, 0,   // 40
          0, 0, 0, 0, 0, 0, 0, 0,   // 48
          0, 0, 0, 0, 0, 0, 0, 0,   // 56
          0, 1, 1, 1, 1, 0, 0, 1,   // 64 ABCDG
          1, 0, 0, 1, 0, 1, 1, 0,   // 72 HKMN
          0, 0, 1, 1, 1, 1, 1, 1,   // 80 RSTUVW
          0, 1, 0, 0, 0, 0, 0, 0,   // 88 Y
          0, 1, 1, 1, 1, 0, 0, 1,   // 96 ABCDG
          1, 0, 0, 1, 0, 1, 1, 0,   // 104 HKMN
          0, 0, 1, 1, 1, 1, 1, 1,   // 112 RSTUVW
          0, 1, 0, 0, 0, 0, 0, 0,   // 120 Y
          0, 0, 0, 0, 0, 0, 0, 0,   // 128
          0, 0, 0, 0, 0, 0, 0, 0,   // 136
          0, 0, 0, 0, 0, 0, 0, 0,   // 144
          0, 0, 0, 0, 0, 0, 0, 0,   // 152
          0, 0, 0, 0, 0, 0, 0, 0,   // 160
          0, 0, 0, 0, 0, 0, 0, 0,   // 168
          0, 0, 0, 0, 0, 0, 0, 0,   // 176
          0, 0, 0, 0, 0, 0, 0, 0,   // 184
          0, 0, 0, 0, 0, 0, 0, 0,   // 192
          0, 0, 0, 0, 0, 0, 0, 0,   // 200
          0, 0, 0, 0, 0, 0, 0, 0,   // 208
          0, 0, 0, 0, 0, 0, 0, 0,   // 216
          0, 0, 0, 0, 0, 0, 0, 0,   // 224
          0, 0, 0, 0, 0, 0, 0, 0,   // 232
          0, 0, 0, 0, 0, 0, 0, 0,   // 240
          0, 0, 0, 0, 0, 0, 0, 0 }; // 248

static constexpr bool IUPAC_noACGT[256] =
        { 0, 0, 0, 0, 0, 0, 0, 0,   // 0
          0, 0, 0, 0, 0, 0, 0, 0,   // 8
          0, 0, 0, 0, 0, 0, 0, 0,   // 16
          0, 0, 0, 0, 0, 0, 0, 0,   // 24
          0, 0, 0, 0, 0, 0, 0, 0,   // 32
          0, 0, 0, 0, 0, 0, 0, 0,   // 40
          0, 0, 0, 0, 0, 0, 0, 0,   // 48
          0, 0, 0, 0, 0, 0, 0, 0,   // 56
          0, 0, 1, 0, 1, 0, 0, 0,   // 64 BD
          1, 0, 0, 1, 0, 1, 1, 0,   // 72 HKMN
          0, 0, 1, 1, 1, 1, 1, 1,   // 80 RSUVW
          0, 1, 0, 0, 0, 0, 0, 0,   // 88 Y
          0, 1, 1, 1, 1, 0, 0, 1,   // 96 BD
          1, 0, 0, 1, 0, 1, 1, 0,   // 104 HKMN
          0, 0, 1, 1, 0, 1, 1, 1,   // 112 RSUVW
          0, 1, 0, 0, 0, 0, 0, 0,   // 120 Y
          0, 0, 0, 0, 0, 0, 0, 0,   // 128
          0, 0, 0, 0, 0, 0, 0, 0,   // 136
          0, 0, 0, 0, 0, 0, 0, 0,   // 144
          0, 0, 0, 0, 0, 0, 0, 0,   // 152
          0, 0, 0, 0, 0, 0, 0, 0,   // 160
          0, 0, 0, 0, 0, 0, 0, 0,   // 168
          0, 0, 0, 0, 0, 0, 0, 0,   // 176
          0, 0, 0, 0, 0, 0, 0, 0,   // 184
          0, 0, 0, 0, 0, 0, 0, 0,   // 192
          0, 0, 0, 0, 0, 0, 0, 0,   // 200
          0, 0, 0, 0, 0, 0, 0, 0,   // 208
          0, 0, 0, 0, 0, 0, 0, 0,   // 216
          0, 0, 0, 0, 0, 0, 0, 0,   // 224
          0, 0, 0, 0, 0, 0, 0, 0,   // 232
          0, 0, 0, 0, 0, 0, 0, 0,   // 240
          0, 0, 0, 0, 0, 0, 0, 0 }; // 248

static constexpr uint8_t ACGT[256] =
        { 0, 0, 0, 0, 0, 0, 0, 0,   // 0
          0, 0, 0, 0, 0, 0, 0, 0,   // 8
          0, 0, 0, 0, 0, 0, 0, 0,   // 16
          0, 0, 0, 0, 0, 0, 0, 0,   // 24
          0, 0, 0, 0, 0, 0, 0, 0,   // 32
          0, 0, 0, 0, 0, 0, 0, 0,   // 40
          0, 0, 0, 0, 0, 0, 0, 0,   // 48
          0, 0, 0, 0, 0, 0, 0, 0,   // 56
          0, 1, 0, 1, 0, 0, 0, 1,    // 64
          0, 0, 0, 0, 0, 0, 0, 0,   // 72
          0, 0, 0, 0, 1, 0, 0, 0,   // 80
          0, 0, 0, 0, 0, 0, 0, 0,   // 88
          0, 1, 0, 1, 0, 0, 0, 1,    // 96
          0, 0, 0, 0, 0, 0, 0, 0,   // 104
          0, 0, 0, 0, 1, 0, 0, 0,   // 112
          0, 0, 0, 0, 0, 0, 0, 0,   // 120
          0, 0, 0, 0, 0, 0, 0, 0,   // 128
          0, 0, 0, 0, 0, 0, 0, 0,   // 136
          0, 0, 0, 0, 0, 0, 0, 0,   // 144
          0, 0, 0, 0, 0, 0, 0, 0,   // 152
          0, 0, 0, 0, 0, 0, 0, 0,   // 160
          0, 0, 0, 0, 0, 0, 0, 0,   // 168
          0, 0, 0, 0, 0, 0, 0, 0,   // 176
          0, 0, 0, 0, 0, 0, 0, 0,   // 184
          0, 0, 0, 0, 0, 0, 0, 0,   // 192
          0, 0, 0, 0, 0, 0, 0, 0,   // 200
          0, 0, 0, 0, 0, 0, 0, 0,   // 208
          0, 0, 0, 0, 0, 0, 0, 0,   // 216
          0, 0, 0, 0, 0, 0, 0, 0,   // 224
          0, 0, 0, 0, 0, 0, 0, 0,   // 232
          0, 0, 0, 0, 0, 0, 0, 0,   // 240
          0, 0, 0, 0, 0, 0, 0, 0 }; // 248

static const int8_t mat50[] = {
        //  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
        5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5,	// A
        -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5,	// R
        -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5,	// N
        -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5,	// D
        -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5,	// C
        -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5,	// Q
        -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5,	// E
        0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5,	// G
        -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5,	// H
        -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5,	// I
        -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5,	// L
        -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5,	// K
        -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5,	// M
        -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5,	// F
        -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5,	// P
        1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5,	// S
        0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5, 	// T
        -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5, 	// W
        -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5, 	// Y
        0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5, 	// V
        -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5, 	// B
        -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5, 	// Z
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5, 	// X
        -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1 	// *
};

/* This table is used to transform amino acid letters into numbers. */
static constexpr int8_t aa_table[128] = {
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
};

static std::map<std::string,char> codon_map = {{"AAA",'K'},{"AAC",'N'},{"AAG",'K'},{"AAR",'K'},{"AAT",'N'},{"AAY",'N'},{"ACA",'T'},{"ACB",'T'},{"ACC",'T'},{"ACD",'T'},{"ACG",'T'},{"ACH",'T'},{"ACK",'T'},{"ACM",'T'},{"ACN",'T'},{"ACR",'T'},{"ACS",'T'},{"ACT",'T'},{"ACV",'T'},{"ACW",'T'},{"ACY",'T'},{"AGA",'R'},{"AGC",'S'},{"AGG",'R'},{"AGR",'R'},{"AGT",'S'},{"AGY",'S'},{"ATA",'I'},{"ATC",'I'},{"ATG",'M'},{"ATH",'I'},{"ATM",'I'},{"ATT",'I'},{"ATW",'I'},{"ATY",'I'},{"CAA",'Q'},{"CAC",'H'},{"CAG",'Q'},{"CAR",'Q'},{"CAT",'H'},{"CAY",'H'},{"CCA",'P'},{"CCB",'P'},{"CCC",'P'},{"CCD",'P'},{"CCG",'P'},{"CCH",'P'},{"CCK",'P'},{"CCM",'P'},{"CCN",'P'},{"CCR",'P'},{"CCS",'P'},{"CCT",'P'},{"CCV",'P'},{"CCW",'P'},{"CCY",'P'},{"CGA",'R'},{"CGB",'R'},{"CGC",'R'},{"CGD",'R'},{"CGG",'R'},{"CGH",'R'},{"CGK",'R'},{"CGM",'R'},{"CGN",'R'},{"CGR",'R'},{"CGS",'R'},{"CGT",'R'},{"CGV",'R'},{"CGW",'R'},{"CGY",'R'},{"CTA",'L'},{"CTB",'L'},{"CTC",'L'},{"CTD",'L'},{"CTG",'L'},{"CTH",'L'},{"CTK",'L'},{"CTM",'L'},{"CTN",'L'},{"CTR",'L'},{"CTS",'L'},{"CTT",'L'},{"CTV",'L'},{"CTW",'L'},{"CTY",'L'},{"GAA",'E'},{"GAC",'D'},{"GAG",'E'},{"GAR",'E'},{"GAT",'D'},{"GAY",'D'},{"GCA",'A'},{"GCB",'A'},{"GCC",'A'},{"GCD",'A'},{"GCG",'A'},{"GCH",'A'},{"GCK",'A'},{"GCM",'A'},{"GCN",'A'},{"GCR",'A'},{"GCS",'A'},{"GCT",'A'},{"GCV",'A'},{"GCW",'A'},{"GCY",'A'},{"GGA",'G'},{"GGB",'G'},{"GGC",'G'},{"GGD",'G'},{"GGG",'G'},{"GGH",'G'},{"GGK",'G'},{"GGM",'G'},{"GGN",'G'},{"GGR",'G'},{"GGS",'G'},{"GGT",'G'},{"GGV",'G'},{"GGW",'G'},{"GGY",'G'},{"GTA",'V'},{"GTB",'V'},{"GTC",'V'},{"GTD",'V'},{"GTG",'V'},{"GTH",'V'},{"GTK",'V'},{"GTM",'V'},{"GTN",'V'},{"GTR",'V'},{"GTS",'V'},{"GTT",'V'},{"GTV",'V'},{"GTW",'V'},{"GTY",'V'},{"MGA",'R'},{"MGG",'R'},{"MGR",'R'},{"NNN",'X'},{"RAY",'B'},{"SAR",'Z'},{"TAA",'.'},{"TAC",'Y'},{"TAG",'.'},{"TAR",'.'},{"TAT",'Y'},{"TAY",'Y'},{"TCA",'S'},{"TCB",'S'},{"TCC",'S'},{"TCD",'S'},{"TCG",'S'},{"TCH",'S'},{"TCK",'S'},{"TCM",'S'},{"TCN",'S'},{"TCR",'S'},{"TCS",'S'},{"TCT",'S'},{"TCV",'S'},{"TCW",'S'},{"TCY",'S'},{"TGA",'.'},{"TGC",'C'},{"TGG",'W'},{"TGT",'C'},{"TGY",'C'},{"TRA",'.'},{"TTA",'L'},{"TTC",'F'},{"TTG",'L'},{"TTR",'L'},{"TTT",'F'},{"TTY",'F'},{"XXX",'X'},{"YTA",'L'},{"YTG",'L'},{"YTR",'L'}};

inline bool file_exists(std::string& fname){
    std::ifstream if_ss;
    if_ss.open(fname);
    int ret = 1;
    if(!if_ss){
        ret = 0;
    }
    if_ss.close();
    return ret;
}

inline void find_all_codons(const std::string & dna, const std::string & codon, std::array<std::vector<uint32_t>, 3> & hits) noexcept{
    size_t pos = dna.find(codon);
    while (pos != std::string::npos){
        hits[pos % 3].push_back(pos);
        pos = dna.find(codon, pos + 1);
    }
}

struct SEQ{
public:
    SEQ()=default;
    SEQ(std::string& exon_nt_seq, uint tx_start = 0, uint tx_end = 0){
        // The tx_start and tx_end sets where on the sequence the transcript begins and ends.
        // This will allow searching up and down stream of the transcript bounds
        this->exon_nt = exon_nt_seq;
        transform(this->exon_nt.begin(), this->exon_nt.end(), this->exon_nt.begin(), ::toupper);
    }


    void append_seq(std::string seq){
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
        this->exon_nt += seq;
    }

    void set_cds(uint start_pos,uint end_pos){
#ifdef DEBUG
        assert(start_pos<this->exon_nt.size());
        assert(end_pos<this->exon_nt.size());
        assert(((end_pos+1)-start_pos)>=3);
        assert(((end_pos+1)-start_pos)%3==0);
#endif
        this->cds_start = start_pos;
        this->cds_end = end_pos;

        this->aa = translate(this->exon_nt,this->cds_start,this->cds_end); // need to be able to set the farthest start/stop coordinates in frame (or a function to compute). then find codon should search nucleotide only - in frame
    }

    uint get_min_inframe(){ // returns the minimum coordinate plausible within the correct frame
        return this->cds_start%3;
    }
    uint get_max_inframe(){
#ifdef DEBUG
        assert(this->cds_end<this->exon_nt.size());
#endif
        return this->exon_nt.size()-((this->exon_nt.size()+1)-this->cds_end)%3;
    }

    std::string get_cds_nt(){return this->exon_nt.substr(this->cds_start,this->cds_nt_len());}
    std::string& get_exon_nt(){return this->exon_nt;}
    std::string& get_aa(){return this->aa;}
    uint get_exon_len(){
#ifdef DEBUG
#endif
        return this->exon_nt.size();
    }
    uint cds_nt_len(){
#ifdef DEBUG
        assert(this->aa.size()*3==(this->cds_end+1)-this->cds_start);
#endif
        return (this->cds_end+1)-this->cds_start;
    }
    uint cds_aa_len(){
#ifdef DEBUG
        assert(this->aa.size()*3==(this->cds_end+1)-this->cds_start);
#endif
        return this->aa.size();
    }

    uint get_cds_start(){return this->cds_start;}
    uint get_cds_end(){return this->cds_end;}

    bool empty(){return this->exon_nt.size()==0;}
    bool is_translated(){return this->aa.size()>0;}
    char get_aa(uint pos){
#ifdef DEBUG
        assert(this->aa.size()>=pos);
#endif
        return this->aa[pos];
    }
    void clear(bool keep_exonic=false){
        if(!keep_exonic){
            this->exon_nt.clear();
        }
        this->aa.clear();
        this->cds_start=0;
        this->cds_end=0;
        this->min_cds_start=0;
        this->max_cds_end=0;
    }
    bool trim_to_stop(){
        size_t stop_aa_pos = this->aa.find('.');
        if(stop_aa_pos!=std::string::npos){ // stop codon was found
            size_t stop_nt_pos = stop_aa_pos*3;
            this->cds_end = this->cds_start+stop_nt_pos-1;
            this->aa.erase(this->aa.begin()+stop_aa_pos,this->aa.end());
#ifdef DEBUG
            assert(this->cds_nt_len()%3==0);
#endif
            return true;
        }
        return false;
    }
    bool trim_to_start(){
        size_t start_aa_pos = this->aa.find('M');
        if(start_aa_pos==0){
            return false;
        }

        if(start_aa_pos!=std::string::npos){ // start codon was found
            size_t start_nt_pos = (start_aa_pos*3);
            this->cds_start+=start_nt_pos;
            this->aa.erase(0,start_aa_pos);
#ifdef DEBUG
            assert(this->cds_nt_len()%3==0);
#endif
            return true;
        }
        return false;
    }
    // find all instances of a codon within the same frame and saves it into a list. the positions in the list are with respect to the nucleotide sequence
    uint find_inframe_codon(char c,char stop_c,std::vector<uint>& positions, uint start_idx,uint stop_idx, bool forward, bool first){
        positions.clear();
        std::string cur_codon;

        uint i=start_idx;
        char cur_c;
        while(true){
            if(!(i<this->exon_nt.size() && i>=0)){
                break;
            }
            cur_codon = this->exon_nt.substr(i,3);
            cur_c = codon_map[cur_codon];
            if(cur_c==stop_c){
                break;
            }
            if(cur_c==c){
                positions.push_back(i);
                if(first){ // if true - only report the first occurrence
                    break;
                }
            }
            bool terminate = forward ? i>stop_idx : i<stop_idx;
            if(terminate){
                break;
            }

            i = forward ? i+3 : i-3;
        }

        return positions.size();
    }
    void extend_to_pos(uint pos){
        if(pos>=this->cds_start && pos<=this->cds_end){ // nothing to do
            return;
        }

        this->cds_start = std::min(pos,this->cds_start);
        this->cds_end = std::max(pos,this->cds_end);

#ifdef DEBUG
        assert(((this->cds_end+1)-this->cds_start)%3==0);
#endif

        // now need to append the nt and aa sequence for these coordinates
        this->aa = translate(this->exon_nt,cds_start,cds_end);
    }
private:
    std::string exon_nt;

    std::string aa;
    uint cds_start=0; // position within the exon nt at which the CDS begins
    uint cds_end=0;

    // farthest plausible coordinates for the cds start and stop in the same frame as the main cds segment set for the sequence
    uint min_cds_start=0;
    uint max_cds_end=0;

    inline std::string translate(std::string& nts,uint start_pos=0,int end_pos=0) {
        end_pos = end_pos==0 ? nts.size() : end_pos;

#ifdef DEBUG
        if(((end_pos+1)-start_pos)%3!=0){
            std::cerr<<"cannot translate sequence of length%3!=0"<<std::endl;
            exit(-1);
        }
        if(nts.size()==0 || (end_pos+1)-start_pos==0){
            std::cerr<<"cannot translate empty sequence"<<std::endl;
            exit(-1);
        }
        assert(end_pos+1<=nts.size());
#endif

        std::string res;
        for(int i=start_pos;i+2<=end_pos;i+=3) {
            res+=codon_map[nts.substr(i,3)];
        }
        return res;
    }
};

#endif