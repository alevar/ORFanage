//
// Created by sparrow on 5/9/22.
//

#ifndef ORFANAGE_TRANSCRIPTOME_H
#define ORFANAGE_TRANSCRIPTOME_H

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <cassert>

#include <gff.h>
#include <map>
#include "common.hpp"
#include "gff_utils.h"

#include <sstream>
#include <set>

class SEGTP{
public:
    SEGTP() = default;
    SEGTP(uint32_t s,uint32_t e){
        assert(e>=s);
        this->start = s;
        this->end = e;
    }
    SEGTP(uint32_t s,uint32_t e,uint8_t p){
        assert(e>=s);
        this->start = s;
        this->end = e;
        this->phase = p;
    }
    ~SEGTP() = default;

    void set_start(uint32_t s){
        assert(this->end>=s);
        this->start=s;
        this->phase = -1;
    }
    void set_end(uint32_t e){
        assert(this->start<=e);
        this->end=e;
        this->phase = -1;
    }
    void set_phase(int8_t p){
        assert(p>=0 && p<=2);
        this->phase = p;
    }

    bool empty(){
        return this->start==UINT32_MAX && this->end==UINT32_MAX;
    }
    bool start_set(){
        return this->start!=UINT32_MAX;
    }
    bool end_set(){
        return this->end!=UINT32_MAX;
    }
    void clear(){
        this->start=UINT32_MAX;
        this->end=UINT32_MAX;
        this->phase=-1;
    }

    friend std::ostream& operator<<(std::ostream& os, const SEGTP& s){
        os << "[" << s.get_start() << "," << s.get_end() << "]";
        return os;
    };

    bool contains(uint pos) const {return pos>=this->start && pos<=this->end;}
    uint32_t get_start() const {return this->start;}
    uint32_t get_start(char strand){return strand=='+'?this->start:this->end;}
    uint32_t get_end() const {return this->end;}
    uint32_t get_end(char strand){return strand=='+'?this->end:this->start       ;}
    uint32_t get_phase() const {return this->phase;}
    uint32_t get_phase(uint pos, char strand) const {
#ifdef DEBUG
        assert(this->contains(pos));
#endif
        uint res_phase = (this->phase+(strand=='+' ?
                                            3-(pos-this->start)%3 :
                                            3-(this->end-pos)%3
                                      )%3)%3;
        return res_phase;
    }

    uint32_t slen() const{return (this->end + 1) - this->start;}

    int intersect(SEGTP& s2,SEGTP& res){
        int is = this->start_set()?std::max(this->start,s2.get_start()):s2.get_start();
        int ie = std::min(this->end,s2.get_end());
        if(is<=ie){
            res = SEGTP(is,ie);
            return (ie-is)+1;
        }
        res.clear();
        return 0;
    }
    void split(SEGTP& s2, SEGTP& res_left, SEGTP& res_inter, SEGTP& res_right){
        res_left.clear();
        res_inter.clear();
        res_right.clear();

        int il = this->intersect(s2,res_inter);

        if(il>0){
            if(res_inter.get_start()>this->get_start()){
                res_left = SEGTP(this->get_start(),res_inter.get_start()-1,this->get_phase());
            }
            if(res_inter.get_end()<this->get_end()){
                res_right = SEGTP(res_inter.get_end()+1,this->get_end(),this->get_phase());
            }
        }
        else{
            if(this->get_start()<s2.get_start()){
                res_left = SEGTP(this->get_start(),this->get_end(),this->get_phase());
            }
            else{
                res_right = SEGTP(this->get_start(),this->get_end(),this->get_phase());
            }
        }
    }
    int get_union(SEGTP& s2, SEGTP& res){
        int is = std::min(this->start,s2.get_start());
        int ie = this->end_set()?std::max(this->end,s2.get_end()):s2.get_end();
        if(is<=ie){
            res = SEGTP(is,ie);
            return (ie-is)+1;
        }
        res.clear();
        return 0;
    }
    int unite(SEGTP& s2){
        return this->get_union(s2,*this);
    }

    bool overlaps(SEGTP& s2){
        if (this->end < s2.start || s2.end < this->start) {
            return false;
        }
        return true;
    }

    bool operator== (const SEGTP& s) const{
        if(this->start!=s.start){
            return false;
        }
        else if(this->end!=s.end){
            return false;
        }
        else if(this->phase!=s.phase){
            return false;
        }
        return true;
    }

    bool operator< (const SEGTP& s) const{
        return this->get_start()==s.get_start() ? this->get_end()<s.get_end() : this->get_start()<s.get_start();
    }

private:
    uint32_t start=UINT32_MAX;
    uint32_t end=UINT32_MAX;
    int8_t phase=-1;
};

class CHAIN{ // chain type - chain of segment type
public:
    CHAIN() = default;

    CHAIN(CHAIN& ch) {this->chain = ch.chain;}
    CHAIN(const CHAIN& ch) {this->chain = ch.chain;}

    // define method to be able to print chain to stdout
    friend std::ostream& operator<<(std::ostream& os, const CHAIN& ch){
        os << "[";
        for(auto& c : ch.chain){
            os << c << ",";
        }
        os << "]";
        return os;
    };

    bool operator== (const CHAIN& ch) const{
        if(this->chain==ch.chain){
            return true;
        }
        return false;
    }

    ~CHAIN() = default;

    void clear(){this->chain.clear();}

    typedef std::vector<SEGTP>::iterator it;
    typedef std::vector<SEGTP>::const_iterator cit;
    it begin() {return this->chain.begin();}
    cit cbegin() const { return this->chain.cbegin();}
    it end() {return this->chain.end();}
    cit cend() const { return this->chain.cend();}

    typedef std::vector<SEGTP>::reverse_iterator rit;
    typedef std::vector<SEGTP>::const_reverse_iterator crit;
    rit rbegin() {return this->chain.rbegin();}
    crit crbegin() const { return this->chain.crbegin();}
    rit rend() {return this->chain.rend();}
    crit crend() const { return this->chain.crend();}

    std::string chain2str(){
        std::string res;
        for(auto& c : this->chain){
            res+=std::to_string(c.get_start())+"-"+std::to_string(c.get_end())+",";
        }
        if(!res.empty()){
            res.pop_back();
        }
        if(res.empty()){
            res = "-";
        }
        return res;
    }

    int clen() const {
        int len=0;
        for(auto& c : this->chain) {
            len += c.slen();
        }
        return len;
    }

    int get_start() const {
        return this->chain.front().get_start();
    }
    int get_end() const {
        return this->chain.back().get_end();
    }

    uint16_t size() const {return this->chain.size();};
    SEGTP& operator[](int idx){return this->chain[idx];}
    void push_back(SEGTP seg){this->chain.emplace_back(seg);}
    void pop_back(){this->chain.pop_back();}
    SEGTP& back(){return this->chain.back();}
    int get_idx(int val){ // return the index of the segment containing a value or -1 if such is not found
        int i=0;
        for(auto& c : this->chain){
            int es = c.get_start(), ee = c.get_end();
            if(val>=es && val<=ee){
                return i;
            }
            i++;
        }
        return -1;
    }

    bool overlap(int s, int e) const {
#ifdef DEBUG
        if (s>e){
            std::cerr<<"Overlapping failed: "<<s<<"\t"<<e<<std::endl;
            exit(-1);
        }
#endif
        return (this->get_start()<=e && this->get_end()>=s);
    }

    int cut(int start,int end){ // cuts its own exon chain - result contains the new length
        if(start<=this->get_start() && end>=this->get_end()){ // nothing to cut
            return this->clen();
        }
        std::vector<SEGTP> tmp;
        for(auto c : this->chain){ // not reference - the copy is being modified
            if(c.get_start()<=start && c.get_end()>=start){
                c.set_start(start);
            }
            if(c.get_end()>=end){
                c.set_end(end);
                tmp.push_back(c);
                break;
            }
            if(c.get_end()<start || c.get_start()>end){
                continue;
            }
            tmp.push_back(c);
        }
        this->chain = tmp;
        return this->clen();
    }

    int get_union(CHAIN& chain2,CHAIN& res){
        // computes union between two chains
        // returns the length of the union

        int i=0,j=0;
        while(i<this->chain.size() && j<chain2.size()){
            SEGTP c1 = this->chain[i];
            SEGTP c2 = chain2[j];

            if (c1.get_end() < c2.get_start()) {
                res.push_back(c1);
                i++;
            }
            else if (c2.get_end() < c1.get_start()) {
                res.push_back(c2);
                j++;
            }
            else {
                int is = std::min(c1.get_start(),c2.get_start());
                int ie = std::max(c1.get_end(),c2.get_end());
                SEGTP union_seg;
                int union_len = c1.get_union(c2,union_seg);
                res.push_back(union_seg);
                i++;
                j++;
            }
        }
        return res.clen();
    }

    int intersection(CHAIN& c2,CHAIN& res,bool segments=false){
        bool found_yet = false, continue_prev = false;
        int start = 0, full_inter_len = 0;;
        SEGTP inter;
        for(auto& s1 : this->chain){
            found_yet = false;
            for(int j=start;j<c2.size();j++){
                int inter_len = s1.intersect(c2[j],inter);
                if(inter_len>0){
                    full_inter_len+=inter_len;
                    found_yet = true;
                    start=j;
                    if(inter.get_end()==c2[j].get_end()){start++;}
                    if(segments){
                        bool ends_match = inter.get_end()==s1.get_end() && s1.get_end()==c2[j].get_end();
                        bool starts_match = inter.get_start()==s1.get_start() && s1.get_start()==c2[j].get_start();
                        if(!continue_prev || !starts_match) { // if previous ends didn't match or current starts don't match - create new intersection chain
                            res.push_back(SEGTP());
                        }
                        res.back().unite(inter);
                        continue_prev = ends_match ? true : false;
                    }
                    else{
                        res.push_back(inter);
                    }
                }
                else{
                    if(found_yet){ // found the last piece in the intersection for the current i1
                        break;
                    }
                    else{
                        continue_prev=false;
                    }
                }
            }
            if(!found_yet && continue_prev){
                continue_prev = false;
            }
        }
        return full_inter_len;
    }

    int next_position(int cp, bool forward){ // returns coordinates of the next position available to extend the cds within the current transcript // TODO: optimize - should not repeat the search for the first exon
        if(forward){
            for(auto& e : this->chain){
                int es = e.get_start(), ee = e.get_end();
                if(cp<es){return es;} // special case for when "cp" matches the last base of the exon
                if(cp>=es && cp<ee){return std::max(cp+1,es);}
            }
        }
        else{
            for(int eidx=this->chain.size()-1;eidx>=0;eidx--){
                int es = this->chain[eidx].get_start(), ee = this->chain[eidx].get_end();
                if(cp>ee){return ee;}
                if(cp>es && cp<=ee){return std::min(cp-1,ee);}
            }
        }
        return -1; // returns -1 if the function couldn't find the next position
    }

    void extend(CHAIN& c2,int numbp,bool forward){
        int lp = forward ? this->get_end() : this->get_start();
        for(int i=0;i<numbp;i++){
            int np = c2.next_position(lp,forward); // TODO: this next_position thing is extremely inefficient...
            int dist = std::abs(np-lp);
            lp = np;
            if(dist>1){ // create new exon
                if(forward){
                    this->chain.push_back(SEGTP(np,np));
                }
                else{
                    this->chain.insert(this->chain.begin(),SEGTP(np,np)); // TODO: inefficient
                }
            }
            if(forward){
                this->chain.back().set_end(np);
            }
            else{
                this->chain.front().set_start(np);
            }
        }
    }
    void extend_to_pos(CHAIN& c2,int pos,bool forward){
        int lp = forward ? this->get_end() : this->get_start();
        int i=0;
        while(true){
            int np = c2.next_position(lp,forward);
            int dist = std::abs(np-lp);
            lp = np;
            if(dist>1){ // create new exon
                if(forward){
                    this->chain.push_back(SEGTP(np,np));
                }
                else{
                    this->chain.insert(this->chain.begin(),SEGTP(np,np)); // TODO: inefficient
                }
            }
            if(forward){
                this->chain.back().set_end(np);
            }
            else{
                this->chain.front().set_start(np);
            }
            if(np==pos){ // found requested position
                break;
            }
            i++;
        }
    }

    void assign_phase(char strand,int start_phase) {
        int cdsacc=start_phase;
        if (strand=='-') { //reverse strand
            for(int i=this->chain.size()-1;i>=0;i--){
                this->chain[i].set_phase((3-cdsacc%3)%3);
                cdsacc+=this->chain[i].get_end()-this->chain[i].get_start()+1;
            }
        }
        else { //forward strand
            for(auto& cds : this->chain){
                cds.set_phase((3-cdsacc%3)%3);
                cdsacc+=cds.get_end()-cds.get_start()+1;
            }
        }
    }
    int genome2nt(int pos,char strand){ // converts genomic to transcriptomic positions // todo: this could really use an optimization to not iterate through all every time...
        int chain_pos = 0;
        bool found_pos = false;
        for(auto& c : this->chain){
            if(pos>=c.get_start() && pos<=c.get_end()){ // found it
                found_pos=true;
                chain_pos+=pos-c.get_start();
                break;
            }
            chain_pos+=c.slen();
        }
        if(!found_pos){
            return -1;
        }
        return strand=='+' ? chain_pos : (this->clen()-chain_pos)-1;
    }
    int nt2genome(int zero_pos,char strand){ // tells which coordinate on the chain corresponds to the coordinate on the nt
        int chain_pos = -1;
        int left_to_stop = zero_pos;
        bool found_pos = false;
        if(strand=='+'){
            for(int i=0;i<this->chain.size();i++){
                size_t clen = this->chain[i].slen();
                if(left_to_stop<clen){ // found the segment with the stop codon
                    chain_pos = this->chain[i].get_start()+left_to_stop;
                    found_pos = true;
                    break;
                }
                left_to_stop-=clen;
            }
            if(!found_pos){ // return the last position
                chain_pos = this->chain.back().get_end();
            }
        }
        else{
            for(int i=chain.size()-1;i>=0;i--){
                size_t clen = this->chain[i].slen();
                if(left_to_stop<clen){ // found the cds segment with the stop codon
                    chain_pos = this->chain[i].get_end()-left_to_stop;
                    found_pos = true;
                    break;
                }
                left_to_stop-=clen;
            }
            if(!found_pos){ // return the last position
                chain_pos = this->chain[0].get_start();
            }
        }
#ifdef DEBUG
        if(chain_pos<0){
            std::cerr<<"unexpected chain_pos<0"<<std::endl;
            exit(-1);
        }
#endif
        return chain_pos;
    }

    void trim_to_pos(int start_nt_pos, int end_nt_pos,char strand){ // from stop indicates whether we are trimming the end or the start of the sequence
#ifdef DEBUG
        assert(end_nt_pos>=start_nt_pos);
#endif
        int start_pos = this->nt2genome(start_nt_pos,strand);
        int end_pos = this->nt2genome(end_nt_pos,strand);
        int s = std::min(start_pos,end_pos);
        int e = std::max(start_pos,end_pos);
        int len = this->cut(s,e);
        int new_start_phase = start_nt_pos%3;
        this->assign_phase(strand,new_start_phase);
    }

    int get_phase(int pos,char strand){ // returns the phase of coordinate in the provided chain
        int cdsacc = strand=='+' ? this->chain.front().get_phase() : this->chain.back().get_phase();
        int phase = strand=='+' ? this->chain.front().get_phase() : this->chain.back().get_phase();
        bool found_pos = false;
        if(strand=='+'){ // TODO: can be simplified to not have the giant if-else
            for(auto & c : this->chain){
                int cs = c.get_start();
                int ce = c.get_end();
                if(pos>=cs && pos<=ce){ // found the cds segment with the requested position
                    phase = c.get_phase(pos,strand);
                    found_pos = true;
                    break;
                }
            }
        }
        else{
            for(int i=this->chain.size()-1;i>=0;i--){
                int cs = this->chain[i].get_start();
                int ce = this->chain[i].get_end();
                if(pos>=cs && pos<=ce){ // found the cds segment with the stop codon
                    phase = this->chain[i].get_phase(pos,strand);
                    found_pos = true;
                    break;
                }
            }
        }
#ifdef DEBUG
        if(!found_pos){
            std::cerr<<"position not found"<<std::endl;
            exit(-1);
        }
#endif

        return phase;
    }

    void compare(CHAIN&t, CHAIN& res){
        res.clear();

        CHAIN intervals;
        for(auto v : this->chain){
            intervals.push_back(SEGTP(v.get_start(),v.get_end(),-1));
        }
        for(auto v : t.chain){
            intervals.push_back(SEGTP(v.get_start(),v.get_end(),1));
        }
        if(intervals.size()==0){
            return;
        }

        std::sort(intervals.begin(), intervals.end());
        res.push_back(intervals[0]);

        SEGTP left,inter,right,tmp;
        for (size_t i = 1; i < intervals.size(); ++i) {
            res.back().split(intervals[i],left,inter,right);
            if(!right.empty()){
#ifdef DEBUG
                assert(inter.slen()==intervals[i].slen());
#endif
            }
            else{
                intervals[i].split(res.back(),tmp,inter,right);
#ifdef DEBUG
                assert(tmp.empty());
#endif
            }
            res.pop_back();

            if(!left.empty()){
                res.push_back(left);
            }
            if(!inter.empty()){
                inter.set_phase(0);
                res.push_back(inter);
            }
            if(!right.empty()){
                res.push_back(right);
            }
        }
    }

protected:
    std::vector<SEGTP> chain;
};

class Bundle;

class TX{
public:
    TX() = default;
    ~TX() = default;
    TX(uint seqid, GffObj* tx,int idx,bool is_templ);
    TX(const TX& tmpl){
        this->is_templ = tmpl.is_templ;
        this->id = tmpl.id;
        this->tid = tmpl.tid;
        this->gid = tmpl.gid;
        this->seqid = tmpl.seqid;
        this->strand = tmpl.strand;
        this->cds_start = tmpl.cds_start;
        this->cds_end = tmpl.cds_end;
        this->cds_phase = tmpl.cds_phase;
        this->is_coding = tmpl.is_coding;
        this->source = tmpl.source;
        this->cds_source = tmpl.cds_source;

        this->attrs = tmpl.attrs;

        this->exons = tmpl.exons;
        this->cds = tmpl.cds;

        this->seq = tmpl.seq;

        this->bundle = tmpl.bundle;
        this->ref_tx = tmpl.ref_tx;
    }

    void set_tid(std::string new_tid) {this->tid = new_tid;}
    bool has_cds() const {return this->is_coding;}
    void set_cds_source(std::string new_source){this->cds_source = new_source;}
    void set_cds_start(int cs){this->cds_start = cs;}
    void set_cds_end(int ce){this->cds_end = ce;}
    void set_cds_phase(int start_phase);
    int get_cds_start() const {return this->cds_start;}
    int get_cds_end() const {return this->cds_end;}
    int get_cds_phase() const {return this->cds_phase;}
    int exon_count() const {return this->exons.size();}
    int get_end() const {return this->exons.get_end();}
    int get_start() const{return this->exons.get_start();}
    CHAIN get_exons() const{return this->exons;}
    CHAIN get_cds() const{return this->cds;}
    int get_invariant_start() const{return this->is_coding ? this->get_start() : this->get_cds_start();}
    int get_invariant_end() const{return this->is_coding ? this->get_end() : this->get_cds_end();}
    int get_seqid() const{return this->seqid;}
    char get_strand() const{return this->strand;}
    std::string get_source() const{return this->source;}
    std::string get_tid() const{return this->tid;}
    int get_id() const{return this->id;}
    std::string get_geneID() const{return this->gid;}
    bool is_template() const{return this->is_templ;}
    int cds_len() const{return this->cds.clen();}
    int cds_alen() const{return this->cds.clen()/3;}
    void set_cds_as_exons();
    bool overlap(int s,int e) const{return this->exons.overlap(s,e);}
    void build_cds();
    void remove_cds();
    bool add_attribute(std::string k,std::string v);
    std::string get_attributes() const;
    void set_bundle_ref(Bundle* b);
    void remove_bundle_ref(){this->bundle = nullptr;}
    void get_nt_seq(std::string& res){_get_nt_seq(this->exons,res);}
    int adjust_stop(); // adjust the chain coordinates to stop at the first stop codon
    int adjust_start(); // adjust the chain coordinates to start at the first start codon
    void assign_phase(int start_phase);
    int positions(int pos,uint num_positions,std::string& nc, bool forward);
    int next_positions(int last_position,uint num_positions,std::string& nc, bool forward); // return the nts of the next positions - these are returned with respect to the + strand and will need to be reverse complemented separately
    char get_codon_aa(uint pos);
    int get_next_codon_nt(uint pos,std::string& nc, bool towards_end);
    void extend_cds_chain(uint num_pos,bool forward){return this->cds.extend(this->exons,num_pos,forward);}
    void extend_cds_chain_to_pos(uint pos,bool forward){return this->cds.extend_to_pos(this->exons,pos,forward);}
    void extend_seq(int extend_len, bool spliced_extend, bool towards_stop);
    void extend_to_stop(int extend_len, bool spliced_extend); // searches downstream of the CDS for the next stop codon in the same frame
    
    void extend_to_start(int new_start,int extend_len, bool spliced_extend);
    uint inframe_len(TX* t);
    bool extend_to_start(TX* t,bool allow_non_aug);
    int rescue_cds(bool allow_non_aug, int extend_len, bool spliced_extend,TX* t=nullptr);
    void remove_seq();
    void load_seq();
    bool seq_loaded(){return !this->seq.empty();}
    uint aa_len(){return this->seq.empty() ? this->cds_alen() : this->seq.cds_aa_len();}
    void correct_chain_len();
    std::string& get_aa(){return this->seq.get_aa();}
    std::string& get_exon_nt(){return this->seq.get_exon_nt();}

    int get_exon_idx(int val){return this->exons.get_idx(val);}
    int get_cds_idx(int val){return this->cds.get_idx(val);}

    void compare(TX& t,CHAIN& res);
    Score score(TX& t);

    CHAIN* cds_chain(){return &this->cds;}
    CHAIN* exon_chain(){return &this->exons;}

    void store_template(TX* t){this->ref_tx = t;}
    TX* get_template(){assert(this->ref_tx); return this->ref_tx;}
    void remove_template(){this->ref_tx= nullptr;}

    bool operator< (const TX& tx) const{
        if(this->get_seqid()!=tx.get_seqid()){
            return this->get_seqid()<tx.get_seqid();
        }
        else if(this->get_strand()!=tx.get_strand()){
            return this->get_strand()<tx.get_strand();
        }
        else if(this->get_start()!=tx.get_start()){
            return this->get_start()<tx.get_start();
        }
        else if(this->get_end()!=tx.get_end()){
            return this->get_end()<tx.get_end();
        }
        else{ // doesn't matter - they definitely overlap
            return false;
        }
    }
    bool operator> (const TX& tx) const{
        if(this->get_seqid()!=tx.get_seqid()){
            return this->get_seqid()>tx.get_seqid();
        }
        else if(this->get_strand()!=tx.get_strand()){
            return this->get_strand()>tx.get_strand();
        }
        else if(this->get_start()!=tx.get_start()){
            return this->get_start()>tx.get_start();
        }
        else if(this->get_end()!=tx.get_end()){
            return this->get_end()>tx.get_end();
        }
        else{ // doesn't matter - they definitely overlap
            return false;
        }
    }

    bool operator== (const TX& tx) const{
        if(this->get_seqid()!=tx.get_seqid()){
            return false;
        }
        else if(this->get_strand()!=tx.get_strand()){
            return false;
        }
        else if(this->get_start()!=tx.get_start()){
            return false;
        }
        else if(this->get_end()!=tx.get_end()){
            return false;
        }
        else if(this->get_exons()==tx.get_exons()){
            return true;
        }
        else{
            return false;
        }
    }

    int merge(TX& tx){
        this->dup_tids.insert(tx.dup_tids.begin(), tx.dup_tids.end() );
        this->dup_tids.insert(tx.get_tid());
        return this->dup_tids.size();
    }

    int num_dups(){
        return this->dup_tids.size();
    }

    std::string str(std::string& seqid){
        std::stringstream os;
        os << seqid << "\t"
              << this->source << "\t"
              << "transcript" << "\t"
              << this->get_start() << "\t"
              << this->get_end() << "\t"
              << "." << "\t"
              << this->get_strand() << "\t"
              << "." << "\t"
              << "transcript_id \""+this->get_tid()+"\"; "
              << "gene_id \""+this->get_geneID()+"\"; "
              << this->get_attributes()
              << std::endl;

        auto e_it = this->exons.cbegin();
        while(e_it != this->exons.cend()){
            os << seqid << "\t"
               << this->source << "\t"
               << "exon" << "\t"
               << e_it->get_start() << "\t"
               << e_it->get_end() << "\t"
               << "." << "\t"
               << this->get_strand() << "\t"
               << "." << "\t"
               << "transcript_id \""+this->get_tid()+"\";"
               << std::endl;
            ++e_it;
        }
        if(this->cds.size()>0){
            auto c_it = this->cds.cbegin();
            while(c_it != this->cds.cend()){
                os << seqid << "\t"
                   << this->cds_source << "\t"
                   << "CDS" << "\t"
                   << c_it->get_start() << "\t"
                   << c_it->get_end() << "\t"
                   << "." << "\t"
                   << this->get_strand() << "\t"
                   << c_it->get_phase() << "\t"
                   << "transcript_id \""+this->get_tid()+"\";"
                   << std::endl;
                ++c_it;
            }
        }
//        os<<"\b \b";
        std::string s = os.str();
        s.pop_back();
        return s;
    }

    SEGTP& operator[](int idx){return this->exons[idx];}

    typedef std::vector<SEGTP>::iterator it;
    typedef std::vector<SEGTP>::const_iterator cit;
    it e_begin() {return this->exons.begin();}
    cit e_cbegin() const { return this->exons.cbegin();}
    it e_end() {return this->exons.end();}
    cit e_cend() const { return this->exons.cend();}

    typedef std::vector<SEGTP>::reverse_iterator rit;
    typedef std::vector<SEGTP>::const_reverse_iterator crit;
    rit e_rbegin() {return this->exons.rbegin();}
    crit e_crbegin() const { return this->exons.crbegin();}
    rit e_rend() {return this->exons.rend();}
    crit e_crend() const { return this->exons.crend();}

    it c_begin() {return this->cds.begin();}
    cit c_cbegin() const { return this->cds.cbegin();}
    it c_end() {return this->cds.end();}
    cit c_cend() const { return this->cds.cend();}

    rit c_rbegin() {return this->cds.rbegin();}
    crit c_crbegin() const { return this->cds.crbegin();}
    rit c_rend() {return this->cds.rend();}
    crit c_crend() const { return this->cds.crend();}

private:
    void _get_nt_seq(CHAIN& chain,std::string& res);

    bool is_templ = false;
    int id = -1;
    std::string tid;
    std::string gid;
    int seqid = -1;
    char strand = '.';
    int cds_start = 0;
    int cds_end = 0;
    int cds_phase = 0;
    bool is_coding = false;
    std::string source;
    std::string cds_source;

    std::map<std::string,std::string> attrs;

    CHAIN exons;
    CHAIN cds;

    SEQ seq;

    Bundle *bundle = nullptr;
    TX* ref_tx = nullptr;

    std::set<std::string> dup_tids;
};

class Bundle{
public:
    Bundle()=default;
    ~Bundle();

    uint32_t blen() const{return (this->bend + 1) - this->bstart;}
    int get_seqid() const{return this->seqid;}
    std::string get_geneID() const{return this->gid;}
    uint32_t get_start() const{return this->bstart;}
    uint32_t get_end() const{return this->bend;}
    char get_nt(int pos) const;
    std::string get_nts(int start,int end) const;
    bool can_add(TX* t,bool use_id=false) const;
    bool add_tx(TX* t,bool use_id=false);
    void load_seq(GFaSeqGet* seq);
    void clear_seq(){bundle_seq.clear();}

    void extend_to(uint32_t pos); // extends coordinates of the bundle to the specified coordinate

    bool has_template() const;
    bool has_query() const;

    uint size() const{return this->txs.size();}

    void clear(){
        this->txs.clear();
        this->seqid = -1;
        this->strand = 0;
        this->bstart = MAX_INT;
        this->bend = 0;
        this->bundle_seq = nullptr;
        this->gid.clear();
    }

    typedef std::vector<TX*>::iterator it;
    typedef std::vector<TX*>::const_iterator cit;
    it begin() {return this->txs.begin();}
    cit cbegin() const {return this->txs.cbegin();}
    it end() {return this->txs.end();}
    cit cend() const {return this->txs.cend();}

    TX* operator[](int idx){return this->txs[idx];}

private:
    std::vector<TX*> txs;
    int seqid = -1;
    char strand = 0;
    int bstart = MAX_INT;
    int bend = 0;
    std::string bundle_seq;
    std::string gid;
};

// uses gffReader to read-in all transcripts
// and then sorts them using custom rules
class Transcriptome{
public:
    Transcriptome()=default;
    ~Transcriptome()=default;

    void set_ref(const std::string& rff);
    void use_non_aug();
    void load_seqids(GffReader& gffReader);
    void add(const std::string& gtf_fname,bool is_templ,bool coding_only);
    void sort(bool use_id=false){
        if(use_id) {
            std::sort(this->tx_vec.begin(), this->tx_vec.end(),[ ]( const TX& lhs, const TX& rhs ){
                if(lhs.get_seqid()!=rhs.get_seqid()){
                    return lhs.get_seqid() < rhs.get_seqid();
                }
                else if(lhs.get_strand()!=rhs.get_strand()) {
                    return lhs.get_strand() < rhs.get_strand();
                }
                else if(lhs.get_geneID()!=rhs.get_geneID()){
                    return lhs.get_geneID() < rhs.get_geneID();
                }
                else{
                    return lhs<rhs;
                }
            });
        }
        else{
            std::sort(this->tx_vec.begin(), this->tx_vec.end());
        }
    }
    int seqid2name(int seqid,std::string& seqid_name);
    GFaSeqGet* get_fasta_seq(int seqid);
    uint bundleup(bool use_id, uint32_t rescue_len=0,bool spliced_extend=false); // create bundles and return the total number of bundles
    void build_cds_chains();
    uint clean_short_orfs(int minlen);
    uint clean_cds(bool rescue, int extend_len, bool spliced_extend, bool use_id);
    uint deduplicate(bool use_id=false);
    void set_cds_as_exons();
    void remove_non_coding();
    void correct_chain_len();
    int size(){return this->tx_vec.size();}

    // Main TX iterators
    typedef std::vector<TX>::iterator it;
    typedef std::vector<TX>::const_iterator cit;
    it begin() {return this->tx_vec.begin();}
    cit cbegin() const { return this->tx_vec.cbegin();}
    it end() {return this->tx_vec.end();}
    cit cend() const { return this->tx_vec.cend();}

    // Additional Bundle iterators
    typedef std::vector<Bundle>::iterator bit;
    typedef std::vector<Bundle>::const_iterator cbit;
    bit bbegin() {return this->bundles.begin();}
    cbit cbbegin() const { return this->bundles.cbegin();}
    bit bend() {return this->bundles.end();}
    cbit cbend() const { return this->bundles.cend();}

    int bsize(){return this->bundles.size();}


private:
    std::vector<TX> tx_vec;
    std::vector<Bundle> bundles;

    std::string ref_fa_fname;
    bool check_ref = false;
    int loaded_seqid = INT_MIN;
    GFaSeqGet* loaded_seq = nullptr;

    std::vector<std::string> seqid_names;
    std::map<std::string,int> seqnames_ids;
    std::pair<std::map<std::string,int>::iterator,bool> n2i_it;

    bool allow_non_aug = false;
};

class CompMat{
public:
    CompMat() = default;
    ~CompMat() = default;

    void clear(){this->mat.clear();}
    void set_size(uint n){
        this->clear();
        for(int ni=n-1;ni>=0;ni--){
            this->mat.emplace_back(std::vector<Score>());
            this->mat.back().resize(ni); // TODO: test
        }
    };

    void set_score(Score& s,uint q,uint t){
        this->mat[this->map(q)][this->map(t)]=s;
    }
    Score* get_score(uint q,uint t){
        return &this->mat[this->map(q)][this->map(t)];
    }
    uint size(){return this->mat.size()+1;}
private:
    uint map(uint i){
        if(i-1>=this->mat.size()){
            std::cerr<<"Out of bounds index in the compmat"<<std::endl;
            exit(4);
        }
        return i-1;
    }

    std::vector<std::vector<Score>> mat;

    // this matrix is strictly for keeping track of scores and does not perform any computation independently
    // the matrix is initialized with predefined sizeor can be resized
    // and the scores can be added into any one index within the matrix
};

#endif //ORFANAGE_TRANSCRIPTOME_H
