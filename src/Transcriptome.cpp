//
// Created by sparrow on 5/9/22.
//

// ========= TX ========

#include <array>
#include <algorithm>
#include "Transcriptome.h"

TX::TX(uint seqid, GffObj* tx,int idx,bool is_templ){
    if(tx->hasCDS()){
        this->is_coding = true;
        this->cds_start = (int)tx->CDstart;
        this->cds_end = (int)tx->CDend;
        switch(tx->CDphase){
            case '0':
                this->cds_phase=0;
                break;
            case '1':
                this->cds_phase=1;
                break;
            case '2':
                this->cds_phase=2;
                break;
            case '3':
                this->cds_phase=3;
                break;
            case '.':
                this->cds_phase=0;
                break;
            default:
                std::cerr<<"unknown CDS phase for the initial segment"<<std::endl;
                exit(-1);
        }
    }

    this->is_templ=is_templ;
    this->id = idx;
    this->tid = tx->getID();
    char* tmp_gid = tx->getGeneID();
    this->gid = tmp_gid==NULL ? this->tid : tmp_gid;
    this->seqid = seqid;
    this->strand = tx->strand;
    this->source = tx->getTrackName();
    this->cds_source = tx->getTrackName();
    for(int i=0;i<tx->exons.Count();i++){
        this->exons.push_back(SEGTP((int)tx->exons.Get(i)->start,(int)tx->exons.Get(i)->end,0));
    }

    // store attributes
    std::pair<std::map<std::string,std::string>::iterator,bool> ait;
    if (tx->attrs!=nullptr) {
        bool trId=false;
        for (int i=0;i<tx->attrs->Count();i++) {
            const char* attrname=tx->names->attrs.getName(tx->attrs->Get(i)->attr_id);
            const char* attrval=tx->attrs->Get(i)->attr_val;
            if (attrval==nullptr || attrval[0]=='\0') continue;
            if (strcmp(attrname, "transcriptID")==0) {
                if (trId) continue;
                trId=true;
            }
            if (strcmp(attrname, "transcript_id")==0 && !trId) {
                attrname="transcriptID";
                trId=true;
            }
            if (Gstrcmp(attrname, "geneID")==0 || strcmp(attrname, "gene_id")==0){
                continue;
            }
            ait = this->attrs.insert(std::make_pair(attrname,tx->attrs->Get(i)->attr_val));
            if(!ait.second){
                std::cerr<<"Attribute: "<<attrname<<" already exists"<<std::endl;
                exit(-1);
            }
        }
    }
}
void TX::set_cds_as_exons(){
    if(this->is_coding && this->cds.size()>0){
        this->exons.clear();
        for(auto c : this->cds){
            c.set_phase(0);
            this->exons.push_back(c);
        }
    }
}
void TX::build_cds(){
    if(this->cds_start==0 && this->cds_end==0){ // cds has not been set
        return;
    }
    if(this->cds_phase!=0){
        if(this->strand=='+'){
            this->cds_start=this->exons.nt2genome(this->exons.genome2nt(this->cds_start,this->strand)+cds_phase,this->strand);
        }
        else{
            this->cds_end=this->exons.nt2genome(this->exons.genome2nt(this->cds_end,this->strand)+cds_phase,this->strand);
        }
        this->cds_phase = 0;
    }
    for(auto e : this->exons){
        if(this->cds_end<e.get_start() || this->cds_start>e.get_end()){
            continue;
        }
        if(this->cds_start>=e.get_start() && this->cds_start<=e.get_end()){
            e.set_start(cds_start);
        }
        if(this->cds_end>=e.get_start() && this->cds_end<=e.get_end()){
            e.set_end(this->cds_end);
        }
        this->cds.push_back(SEGTP(e.get_start(),e.get_end(),0));
    }
    this->cds.assign_phase(this->strand,0);
}
void TX::correct_chain_len(){ // trims the chain to be of length len%3==0
    int num_trim = this->cds.clen()%3;
    if(num_trim>0){
        this->cds.trim_to_pos(0,this->cds.clen()-num_trim-1,this->strand);
        this->cds_start = this->cds.get_start();
        this->cds_end = this->cds.get_end();
    }
}
void TX::remove_cds(){
    this->cds = CHAIN();
    this->cds_start = 0;
    this->cds_end = 0;
    this->cds_phase=0;
    this->is_coding = false;
    this->seq.clear();
}
bool TX::add_attribute(std::string k,std::string v){
    std::pair<std::map<std::string,std::string>::iterator,bool> ait;
    ait = this->attrs.insert(std::make_pair(k,v));
    return ait.second;
}
std::string TX::get_attributes() const {
    std::string res;
    for(auto& a : this->attrs){
        res+=a.first+" \""+a.second+"\"; ";
    }
    if(!res.empty()){
        res.pop_back();
    }
    return res;
}
void TX::set_bundle_ref(Bundle* b){
    this->bundle = b;
}
void TX::_get_nt_seq(CHAIN& chain,std::string& res){
    res.clear();
    if(this->strand=='+'){
        for(auto& c : chain){
            res+=this->bundle->get_nts(c.get_start(),c.get_end());
        }
    }
    else{ // strand=='-'
        for(int ci= chain.size() - 1; ci >= 0; ci--){
            std::string tmp = this->bundle->get_nts(chain[ci].get_start(),chain[ci].get_end());
            for(int i=tmp.size()-1;i>=0;i--){
                res+=ntComplement(tmp[i]);
            }
        }
    }
    return;
}

// TODO: extension - should we create a wrapper class for sequence data such that handles extension? What if extension is enabled for the bundle (bundle load extended data), but this function forgets the argument?
//    bundle should load sequence with extension and adjust start coordinate.
//    when transcripts load sequence, they should store extension separately, thus only querying it if necessary, but otherwise not going over bound
//    if extension is necessary and transcript coordinates are adjusted - the part of the extension from the sequence needs to be poped from extension and merged into main transcript sequence
//    store two extensions - 3' and 5'
int TX::adjust_stop(){ // adjust the chain coordinates to stop at the first stop codon
    if(this->seq.empty()){return 0;}

    int prev_nt_len = this->seq.cds_nt_len();
    bool trimmed = this->seq.trim_to_stop();
    if(this->seq.cds_nt_len()<3){
        this->remove_cds();
        return -1;
    }
    if(trimmed){
        this->cds.trim_to_pos(0,this->seq.cds_nt_len()-1,this->strand);
        this->cds_start = this->cds.get_start();
        this->cds_end = this->cds.get_end();
    }
    int new_nt_len = this->seq.cds_nt_len();
    return prev_nt_len-new_nt_len;
}
int TX::adjust_start(){ // adjust the chain coordinates to start at the first start codon
    if(this->seq.empty()){return 0;}
#ifdef DEBUG
    assert(this->cds.clen()==this->seq.cds_nt_len());
#endif

    int prev_nt_len = this->seq.cds_nt_len();
    uint old_start = this->seq.get_cds_start();
    bool trimmed = this->seq.trim_to_start();

    if(trimmed){
        uint end_nt_pos = this->seq.get_cds_end()-old_start;
        uint start_nt_pos = (end_nt_pos+1) - this->seq.cds_nt_len();
        this->cds.trim_to_pos(start_nt_pos,end_nt_pos,strand);
#ifdef DEBUG
        assert(this->cds.clen()%3==0);
#endif
        this->cds_start = this->cds.get_start();
        this->cds_end = this->cds.get_end();
    }
    int new_nt_len = this->seq.cds_nt_len();
    return prev_nt_len-new_nt_len;
}
int TX::positions(int pos,uint num_positions,std::string& nc, bool forward){ // return the nts starting at position - these are returned with respect to the + strand and will need to be reverse complemented separately
    nc = std::string(num_positions,'-');
    int np = pos;
    for(uint i=0;i<num_positions;i++){
        int res_idx = forward ? i : (num_positions-1)-i;
        nc[res_idx] = this->bundle->get_nt(np);
        if(i<num_positions-1){
            np = this->exons.next_position(np,forward);
            if(np==-1){return np;} // means the requested number of positions is unavailable
        }
    }
    return np;
}
int TX::next_positions(int last_position,uint num_positions,std::string& nc, bool forward){ // return the nts of the next positions - these are returned with respect to the + strand and will need to be reverse complemented separately
    int np = last_position;
    np = this->exons.next_position(np,forward);
    if(np==-1){return -1;}
    return this->positions(np,num_positions,nc,forward);
}
char TX::get_codon_aa(uint pos){
    if(this->seq.is_translated()){
        return this->seq.get_aa(pos);
    }
    else{ // we did not pre-load the sequence and need to specifically load just the codon
        pos*=3;
        int nt_pos = this->cds.nt2genome(pos,this->strand);
        std::string nc;
        int lp = positions(nt_pos,3,nc,strand=='+');
        if(nc.empty()){return '-';}

        if(this->strand=='-'){
            std::string nts_revcomp;
            for(int i=nc.size()-1;i>=0;i--){
                nts_revcomp+=ntComplement(nc[i]);
            }
            nc = nts_revcomp;
        }
        return codon_map[nc];
    }
}
int TX::get_next_codon_nt(uint pos,std::string& nc, bool towards_end){
    int lp = next_positions(pos,3,nc,towards_end == (this->strand=='+'));

    if(lp==-1){return -1;}

    if(this->strand=='-'){
        std::string nts_revcomp;
        for(int i=nc.size()-1;i>=0;i--){
            nts_revcomp+=ntComplement(nc[i]);
        }
        nc = nts_revcomp;
    }
    return lp;
}
void TX::extend_to_stop(int extend_len){ // searches downstream of the CDS for the next stop codon in the same frame
    if(this->seq.get_aa(this->seq.cds_aa_len()-1)=='.'){ // stop already present
        return;
    }

    // if requested - extend sequence
    if(extend_len>0){ // consider additional sequence
        std::string extension_str = "";
        if(this->strand=='+'){
            uint32_t extension_start = std::min(this->exons.get_end()+1,(int)this->bundle->get_end());
            uint32_t extension_end = std::min(this->exons.get_end()+extend_len,(int)this->bundle->get_end());
            extension_str = this->bundle->get_nts(extension_start,extension_end);
        }
        else{ // strand=='-'
            uint32_t extension_start = std::max(this->exons.get_start()-extend_len,(int)this->bundle->get_start());
            uint32_t extension_end = std::max(this->exons.get_start()-1,(int)this->bundle->get_start());
            std::string tmp = this->bundle->get_nts(extension_start,extension_end);
            for(int i=tmp.size()-1;i>=0;i--){
                extension_str+=ntComplement(tmp[i]);
            }
        }
        this->seq.append_seq(extension_str);
    }

    std::vector<uint> stops;
    uint res = this->seq.find_inframe_codon('.',';',stops,this->seq.get_cds_end()+1,this->seq.get_exon_len(),true,true);    
    if(res!=0){ // stop is found - adjust accordingly
        // if extension requested and stop found within extension - set new transcript coordinates
        if(extend_len>0){
            int extension_len = ((stops[0]+2)-this->exons.clen())+1;
            if(strand=='+'){
                (this->exons.end()-1)->set_end(this->exons.get_end()+extension_len);
            }
            else{
                this->exons.begin()->set_start(this->exons.get_start()-extension_len);
            }
        }
        int chain_extension_len = (stops[0]+2)-this->seq.get_cds_end();
        this->seq.extend_to_pos(stops[0]+2);
        // modify the chain accordinly
        this->extend_cds_chain(chain_extension_len,this->strand=='+');
        this->cds_start = this->cds.get_start();
        this->cds_end = this->cds.get_end();
    }
    
    // reload sequence to adjust for the changed coordinates
    this->load_seq();
    
    return;
}
// this version will search for the first available start codon resulting in the longest possible ORF
void TX::extend_to_start(int new_start){ // searches upstream for the available start codons
    if(new_start==-1){ // start coordinate is not provided - search
        if(this->seq.get_aa(0)=='M'){ // start already present
            return;
        }
        std::vector<uint> starts;
        uint res = this->seq.find_inframe_codon('M','.',starts,this->seq.get_cds_start(),0,false,true);
        if(res==0){
            return;
        }
        new_start = starts[0];
    }

    int chain_extension_len = this->seq.get_cds_start()-new_start;
    this->seq.extend_to_pos(new_start);
    // modify the chain accordinly
    this->extend_cds_chain(chain_extension_len,this->strand=='-');
    this->cds_start = this->cds.get_start();
    this->cds_end = this->cds.get_end();
    return;
}
// this version of the extension function will search for the start codon which maximizes the similarity to the original template ORF
void TX::extend_to_start(TX* ref_tx,bool allow_non_aug){ // extends to start while comparing to the closest reference
    // can skip the ckeck if the start coordinate is the same as the reference and is M
    char start_codon = this->get_codon_aa(0);
    bool match_start = this->strand=='+' ? this->get_cds_start()==ref_tx->get_cds_start() : this->get_cds_end()==ref_tx->get_cds_end();
    if((start_codon=='M' || allow_non_aug) && match_start){return;}

    // check non-aug
    if(allow_non_aug){
        // check if reference start matching is possible (much faster) - if yes - use that, otherwise proceed evaluating other possibilities
        int ref_cds_start = this->strand=='+' ? ref_tx->get_cds_start() : ref_tx->get_cds_end();
        int this_cds_start = this->strand=='+' ? this->get_cds_start() : this->get_cds_end();
        int ref_start_in_this = this->exons.genome2nt(ref_cds_start,this->strand);
        int q_start_in_this = this->exons.genome2nt(this_cds_start,this->strand);

        if(ref_start_in_this >=0 && ref_start_in_this%3==q_start_in_this%3){ // found start - can extend and exit
            // need to make sure no stop codons are being introduced in the process
            std::vector<uint> stops;
            uint res = this->seq.find_inframe_codon('.',';',stops,this->seq.get_cds_start(),ref_start_in_this,false,true);
            if(res==0){
                extend_to_start(ref_start_in_this);
                return;
            }
        }
    }

    // otherwise - search for the best start codon
    std::vector<uint> transcriptomic_starts;
    uint res = this->seq.find_inframe_codon('M','.',transcriptomic_starts,this->seq.get_cds_start(),0,false,false);
    if(res==0){return;}

    int max_chain_extension_len = this->seq.get_cds_start()-transcriptomic_starts.back();

    std::vector<uint> genomic_starts;
    for(auto& s : transcriptomic_starts){
        genomic_starts.push_back(this->exons.nt2genome(s,strand));
    }

    // check if the reference start codon is present
    uint si=0;
    int old_start = this->strand=='+' ? ref_tx->get_cds_start() : ref_tx->get_cds_end();
    for(auto& s : genomic_starts){
        if(s==old_start){
            extend_to_start(transcriptomic_starts[si]);
            return;
        }
        si++;
    }

    // STRATEGY
    // 1. extend to the farthest start
    this->extend_cds_chain(max_chain_extension_len,this->strand=='-');

    // TODO: 2. if longest orf requested - that's it

    // 3. else compute intersection and find the last start codon which adds inframe matches == maximizes similarity
    // compute intersection (q and t only)
    CHAIN intersected_chain;
    this->cds.intersection(ref_tx->cds,intersected_chain);

    uint best_start_idx = genomic_starts.size()-1;
    // iterate over them - check phase in both - continue
    auto get_best_start = [&best_start_idx, &transcriptomic_starts, &ref_tx, this](SEGTP &c) {
        // as soon as we find the first valid (with added inframe) - save and terminate
        // transcriptomic positions with respect to q and t
        uint t_tpos = ref_tx->cds.genome2nt(c.get_start(this->strand),this->strand);
        uint q_tpos = this->cds.genome2nt(c.get_start(this->strand),this->strand);

        for(int n=best_start_idx;n>=0;n--){
            if(transcriptomic_starts[n] <= q_tpos){
                best_start_idx = n;
            }
            else{
                break;
            }
        }
        if((t_tpos%3)==(q_tpos%3)){ // found first occurrence of inframe
            return true;
        }
        return false;
    };
    if (this->strand == '+') {
        std::find_if(intersected_chain.begin(), intersected_chain.end(), get_best_start);
    } else {
        std::find_if(intersected_chain.rbegin(), intersected_chain.rend(), get_best_start);
    }

#ifdef DEBUG
    assert(best_start_idx<=transcriptomic_starts.size()-1 && best_start_idx>=0);
#endif

    // 4. cut to length
    int start_nt_pos = transcriptomic_starts[best_start_idx]-transcriptomic_starts.back();
    int end_nt_pos = this->cds.clen();
    this->cds.trim_to_pos(start_nt_pos,end_nt_pos,this->strand);

    this->seq.extend_to_pos(transcriptomic_starts[best_start_idx]);
    this->cds_start = this->cds.get_start();
    this->cds_end = this->cds.get_end();
    return;
}
uint TX::inframe_len(TX* t){ // it doesn't matter that we do this stranded - either way the number of inframe codons should be the same
    uint res = 0;

    // compute intersection (q and t only)
    CHAIN intersected_chain;
    this->cds.intersection(t->cds,intersected_chain);

    // iterate over them - check phase in both - continue
    for(auto& c : intersected_chain){ // get phase at position for 1
        // transcriptomic positions with respect to q and t
        uint t_tpos = t->cds.genome2nt(c.get_start(),this->strand);
        uint q_tpos = this->cds.genome2nt(c.get_start(),this->strand);

        res += (t_tpos%3)==(q_tpos%3) ? c.slen() : 0;
    }
    return res;
}

int TX::rescue_cds(bool allow_non_aug,int extend_len, TX* t){
    if(this->seq.cds_nt_len()<3){
        this->remove_cds();
        return 0;
    }

    adjust_stop();
    if(this->cds_len()==0){
        this->remove_cds();
        return 0;
    }
    extend_to_stop(extend_len);
    if(t == nullptr){
        if(!allow_non_aug){
            extend_to_start();
        }
    }
    else{
        extend_to_start(t,allow_non_aug);
    }
    bool run_adjust_start = true;
    if(allow_non_aug && t!= nullptr){
        run_adjust_start = !(this->strand=='+' ? this->cds_start==t->get_cds_start() : this->cds_end==t->get_cds_end());
    }
    if(run_adjust_start){
        adjust_start();
    }

    this->set_cds_phase(0);
    return 1;
}
void TX::set_cds_phase(int start_phase){
    this->cds_phase=start_phase;
    this->cds.assign_phase(this->strand,0);
}
void TX::remove_seq(){
    this->seq.clear();
}
void TX::load_seq(){
    this->seq.clear();
    std::string nt_seq;
    this->get_nt_seq(nt_seq);

    this->seq = SEQ(nt_seq);
#ifdef DEBUG
    assert(this->cds.size()>0);
    assert(this->cds_phase==0);
#endif
    int stranded_start = this->strand=='+' ? this->cds_start : this->cds_end;
    int stranded_end = this->strand=='+' ? this->cds_end : this->cds_start;
    this->seq.set_cds(this->exons.genome2nt(stranded_start,this->strand),
                      this->exons.genome2nt(stranded_end,this->strand));
}
void TX::compare(TX& t,CHAIN& res){
    return this->cds.compare(t.cds,res);
}
Score TX::score(TX& t) {
    Score s;
    s.qlen = this->cds_len();
    s.tlen = t.cds_len();
    CHAIN union_chain;\
    s.ulen = this->cds.get_union(t.cds,union_chain);

    // 1. compute the total number of matching positions between query and template
    // 2. compute the number of matching positions in frame between query and template
    CHAIN mod_chain;
    this->compare(t, mod_chain);

    s.num_bp_extra = 0;
    s.num_bp_missing = 0;
    s.num_bp_inframe = 0;
    s.num_bp_match = 0;
    s.num_bp_outframe = 0;

    int t_frame = 0;
    int q_frame = 0;
    auto extract_mods = [&s, &t_frame, &q_frame](SEGTP &mc) {
        if (mc.get_phase() == -1) { // extra positions in the query
            s.num_bp_extra += mc.slen();
            q_frame += mc.slen();
        } else if (mc.get_phase() == 1) { // template positions missing from the query
            s.num_bp_missing += mc.slen();
            t_frame += mc.slen();
        } else if (mc.get_phase() == 0) { // matching positions between query and template
            s.num_bp_match += mc.slen();
            if (q_frame % 3 == t_frame % 3) {
                s.num_bp_inframe += mc.slen(); // TODO: shouldn't this be stranded?
            } else {
                s.num_bp_outframe += mc.slen();
            }
        } else {
            std::cerr << "wrong code" << std::endl;
            exit(-6);
        }
    };
    if (strand == '+') {
        std::for_each(mod_chain.begin(), mod_chain.end(), extract_mods);
    } else {
        std::for_each(mod_chain.rbegin(), mod_chain.rend(), extract_mods);
    }

    // compute lpi, ilpi, mlpi, etc
    s.lpi = 100.0 * ((float) s.qlen / (float) s.ulen);
    s.ilpi = 100.0 * ((float) s.num_bp_inframe / (float) s.ulen);
    s.mlpi = 100.0 * ((float) s.num_bp_match / (float) s.ulen);

    return s;
}

// ======= Bundle =======
Bundle::~Bundle(){
    for(auto& tx : this->txs){
        tx->remove_bundle_ref();
    }
}
char Bundle::get_nt(int pos) const{
    int np = pos-this->get_start();
    assert(np>=0);
    assert(np<this->bundle_seq.size());
    return this->bundle_seq[np];
}
std::string Bundle::get_nts(int s,int e) const{
    int ns = s - this->get_start();
    int ne = e - this->get_start();
    int nl = (ne+1)-ns;
    assert(ns>=0);
    return this->bundle_seq.substr(ns,nl);
}
bool Bundle::can_add(TX* t,bool use_id) const{
    if(this->size()==0){
        return true;
    }
    else{
        if(use_id){
            return this->gid == t->get_geneID();
        }
        if(this->seqid == t->get_seqid() && this->strand==t->get_strand()){
            if(t->overlap(this->bstart,this->bend)){
                return true;
            }
            else{
                return false;
            }
        }
        return false;
    }
}
bool Bundle::add_tx(TX* t,bool use_id){
#ifdef DEBUG
    if(!this->can_add(t)){
        return false;
    }
#endif
    this->txs.push_back(t);
    this->seqid = t->get_seqid();
    this->strand = t->get_strand();
    this->bstart = std::min(this->bstart,t->get_start());
    this->bend = std::max(this->bend,t->get_end());
    this->gid = use_id ? t->get_geneID() : "";
//    t->set_bundle_ref(this);
    return true;
}

// extends coordinates of the bundle to the specified coordinate
void Bundle::extend_to(uint32_t pos){
    this->bstart = std::min(this->bstart,(int)pos);
    this->bend = std::max(this->bend,(int)pos);
}
void Bundle::load_seq(GFaSeqGet* seq){
    if(this->bstart+this->blen()>=seq->getseqlen() || this->bstart<0){
        std::cerr<<"bundle+stop codon extends past the end of the reference sequence"<<std::endl;
        exit(2);
    }
    int tmp_blen = this->blen();

    this->bundle_seq = std::string(seq->subseq(this->bstart,tmp_blen),tmp_blen);
}
bool Bundle::has_template() const{
    for(auto& tx : this->txs){
        if(tx->is_template()){
            return true;
        }
    }
    return false;
}
bool Bundle::has_query() const{
    for(auto& tx : this->txs){
        if(!tx->is_template()){
            return true;
        }
    }
    return false;
}

// ======== Transcriptome ========
void Transcriptome::set_ref(const std::string& rff){
    this->ref_fa_fname = rff;
    this->check_ref = true;
    gfasta.init(rff.c_str());
}
void Transcriptome::use_non_aug(){
    this->allow_non_aug=true;
}
void Transcriptome::set_aligner(const int8_t *mat, const int8_t *alphabet, int gapo, int gape){
    this->aligner = Finder(mat,alphabet,gapo,gape);
}
void Transcriptome::load_seqids(GffReader& gffReader){
    for(int i=0;i<gffReader.gflst.Count();++i) {
        GffObj *pGffObj = gffReader.gflst.Get(i);
        this->n2i_it = this->seqnames_ids.insert(std::make_pair(pGffObj->getGSeqName(),this->seqid_names.size()));
        if(this->n2i_it.second){ // new seqid
            this->seqid_names.resize(this->seqid_names.size()+1);
            this->seqid_names[this->n2i_it.first->second] = this->n2i_it.first->first;
        }
    }
}
void Transcriptome::add(const std::string& gtf_fname,bool is_templ,bool coding_only) {
    this->bundles.clear(); // any pointers that previously existed will be invalidated when new data is added
    FILE *gff_file = fopen(gtf_fname.c_str(), "r");
    if (gff_file == nullptr) {
        std::cerr << "@ERROR::Couldn't open the GTF: " << gtf_fname << std::endl;
        exit(1);
    }
    GffReader gffReader(gff_file, false, false);
    gffReader.keepAttrs();
    gffReader.readAll();

    this->load_seqids(gffReader);

    for (int i = 0; i < gffReader.gflst.Count(); ++i) {
        GffObj *pGffObj = gffReader.gflst.Get(i);

        if(!pGffObj->isTranscript()){ // skip records which could not be transformed into a transcript
            continue;
        }

        if (!coding_only || pGffObj->hasCDS()) {
            this->n2i_it.first = this->seqnames_ids.find(pGffObj->getGSeqName());
            TX tmp(this->n2i_it.first->second, pGffObj, tx_vec.size(), is_templ);
            tx_vec.push_back(tmp);
        }
    }
}
int Transcriptome::seqid2name(int seqid,std::string& seqid_name){
    if(this->seqid_names.size()<seqid){
        seqid_name = "";
        return -1;
    }
    seqid_name = this->seqid_names[seqid];
    if(seqid_name.empty()){
        return -1;
    }
    return 0;
}

GFaSeqGet* Transcriptome::get_fasta_seq(int seqid){
    if(this->check_ref){
        if(seqid!=this->loaded_seqid){
            std::string seqid_name;
            int res = this->seqid2name(seqid,seqid_name);
            assert(res==0);
            this->loaded_seq = fastaSeqGet(gfasta, seqid_name.c_str());
            this->loaded_seqid = seqid;
        }
    }
    return this->loaded_seq;
}

uint Transcriptome::bundleup(bool use_id, uint32_t overhang){ // create bundles and return the total number of bundles
    // if rescue_len is enabled - will extend bundles by that number of bases to enable rescue lookup up and down stream of the original locus coordinates
    this->sort(use_id);
    this->bundles.clear();
    this->bundles.push_back(Bundle());

    GFaSeqGet* seqid_seq = nullptr;
    int cur_seqid = -1;
    std::string seqid_name;
    if(this->check_ref && !this->tx_vec.empty()){ // preload reference based on the front transcript
        cur_seqid = this->tx_vec.front().get_seqid();
        int res = this->seqid2name(cur_seqid,seqid_name);
        assert(res==0);
        seqid_seq = fastaSeqGet(gfasta, seqid_name.c_str());
    }

    for(auto& t : this->tx_vec){
        if(!this->bundles.back().can_add(&t,use_id)){
            if(this->check_ref){
                // compute overhangs to not overextend past the end of the reference
                uint32_t new_start_pos = std::max(this->bundles.back().get_start()-overhang,(uint32_t)0);
                uint32_t new_end_pos = std::min(this->bundles.back().get_end()+overhang,(uint32_t)seqid_seq->getseqlen()-1);
                // make sure these operations did not shorten the data
                assert(new_start_pos<=this->bundles.back().get_start());
                assert(new_end_pos>=this->bundles.back().get_end());
                // extend
                this->bundles.back().extend_to(new_start_pos);
                this->bundles.back().extend_to(new_end_pos);

                // load sequence information
                this->bundles.back().load_seq(seqid_seq);
                if(t.get_seqid()!=cur_seqid){
                    cur_seqid = t.get_seqid();
                    int res = this->seqid2name(cur_seqid,seqid_name);
                    assert(res==0);
                    seqid_seq = fastaSeqGet(gfasta, seqid_name.c_str());
                }
            }
            this->bundles.push_back(Bundle());
            this->bundles.back().add_tx(&t,use_id);
        }
        else{
            this->bundles.back().add_tx(&t,use_id);
        }
    }
    if(this->check_ref){
        // compute overhangs to not overextend past the end of the reference
        uint32_t new_start_pos = std::max(this->bundles.back().get_start()-overhang,(uint32_t)0);
        uint32_t new_end_pos = std::min(this->bundles.back().get_end()+overhang,(uint32_t)seqid_seq->getseqlen()-1);
        // make sure these operations did not shorten the data
        assert(new_start_pos<=this->bundles.back().get_start());
        assert(new_end_pos>=this->bundles.back().get_end());
        // extend
        this->bundles.back().extend_to(new_start_pos);
        this->bundles.back().extend_to(new_end_pos);

        this->bundles.back().load_seq(seqid_seq);
    }

    for(auto& b : this->bundles){
        for(auto& t : b){
            t->set_bundle_ref(&b);
        }
    }

    return this->bundles.size();
}
void Transcriptome::build_cds_chains(){
    for(auto& tx : this->tx_vec){
        if(tx.has_cds() && tx.is_template()) {
#ifdef DEBUG
            if(std::strcmp(tx.get_tid().c_str(),"ENST00000313599.8")==0){
                std::cout<<"found 2"<<std::endl;
            }
#endif
            tx.build_cds();
        }
    }
}
void Transcriptome::correct_chain_len(){
    for(auto& tx : this->tx_vec){
        if(tx.has_cds() && tx.is_template()) {
            tx.correct_chain_len();
        }
    }
}
uint Transcriptome::clean_short_orfs(int minlen){
    uint res=0;
    for(auto& tx : this->tx_vec){
        if(tx.cds_len()<=minlen){
            tx.remove_cds();
            res++;
        }
    }
    return res;
}
uint Transcriptome::clean_cds(bool rescue,int extend_len,bool use_id){
    uint res = 0;

    this->bundleup(use_id,extend_len);
    for(auto& bundle : this->bundles){
        for(auto& tx : bundle){
            if(!tx->has_cds()){continue;}
            tx->load_seq();

            if(this->check_ref){
                if(rescue){
                    tx->rescue_cds(this->allow_non_aug,extend_len); // TODO: allow non aug - needs to allow storage of non-aug transcripts (only perform cleaning based on stop codons) and no start extension. And needs to handle querrying appropriately
                }
                char start_codon = tx->get_codon_aa(0);
                char stop_codon = tx->get_codon_aa(tx->aa_len()-1);
                if((start_codon!='M' && !this->allow_non_aug) || stop_codon!='.'){
                    res++;
                    tx->remove_cds();
                    continue;
                }
            }

            if(tx->cds_len()%3 != 0){
                res++;
                tx->remove_cds();
                continue;
            }
        }
    }
    return res;
}

// removes duplicate transcripts (same exon chain)
// if requested - also removes transcripts with the same CDS chains
// transcript IDs are merged into the parenting transcript
// if requested - comparisons are only made between transcripts with the same geneID
bool unique_tx_cmp_with_id(TX& t1, TX& t2){
    if(t1.get_geneID() != t2.get_geneID()){
        return false;
    }
    if(t1==t2){
        t1.merge(t2);
        return true;
    }
    return false;
}
bool unique_tx_cmp(TX& t1, TX& t2){
    if(t1==t2){
        t1.merge(t2);
        return true;
    }
    return false;
}
uint Transcriptome::deduplicate(bool use_id){
    // first deduplicate based on transcripts

    uint old_transcriptome_size = this->tx_vec.size();

    if(use_id){
        auto last = std::unique( this->tx_vec.begin(), this->tx_vec.end(), unique_tx_cmp_with_id);
        this->tx_vec.erase( last, this->tx_vec.end() );
    }
    else{
        auto last = std::unique( this->tx_vec.begin(), this->tx_vec.end(), unique_tx_cmp);
        this->tx_vec.erase( last, this->tx_vec.end() );
    }

    return old_transcriptome_size - this->tx_vec.size();
}
void Transcriptome::set_cds_as_exons(){
    for(auto& tx : this->tx_vec){
        tx.set_cds_as_exons();
    }
}
void Transcriptome::remove_non_coding(){
    this->tx_vec.erase(std::remove_if(
            this->tx_vec.begin(),this->tx_vec.end(),
            [](const TX& tx) {
                return !tx.has_cds();
            }), this->tx_vec.end());
}