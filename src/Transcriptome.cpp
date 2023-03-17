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
    this->gid = tx->getGeneID();
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
void TX::extend_to_stop(){ // searches downstream of the CDS for the next stop codon in the same frame
    if(this->seq.get_aa(this->seq.cds_aa_len()-1)=='.'){ // stop already present
        return;
    }
    std::vector<uint> stops;
    uint res = this->seq.find_inframe_codon('.',';',stops,this->seq.get_cds_end()+1,this->seq.get_exon_len(),true,true);
    if(res!=0){ // stop is found - adjust accordingly
        int chain_extension_len = (stops[0]+2)-this->seq.get_cds_end();
        this->seq.extend_to_pos(stops[0]+2);
        // modify the chain accordinly
        this->extend_cds_chain(chain_extension_len,this->strand=='+');
        this->cds_start = this->cds.get_start();
        this->cds_end = this->cds.get_end();
    }
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
int TX::rescue_cds(bool allow_non_aug,TX* t){
    if(this->seq.cds_nt_len()<3){
        this->remove_cds();
        return 0;
    }

    adjust_stop();
    if(this->cds_len()==0){
        this->remove_cds();
        return 0;
    }
    extend_to_stop();
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

    // compute lpd, ilpd, mlpd, etc
    s.lpd = 100.0 * ((float) s.qlen / (float) s.tlen);
    s.ilpd = 100.0 * ((float) s.num_bp_inframe / (float) s.tlen);
    s.mlpd = 100.0 * ((float) s.num_bp_match / (float) s.tlen);

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
void Bundle::load_seq(GFaSeqGet* seq){
#ifdef DEBUG
//    if(this->bstart+this->blen()>=seq->getseqlen() || this->bstart<0){
//        std::cerr<<"bundle+stop codon extends past the end of the reference sequence"<<std::endl;
//        exit(2);
//    }
#endif
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
void Transcriptome::load_ppptracks(){
    // check that reference data has been loaded already
    if(this->seqid_names.size()==0){
        std::cerr<<"PhyloCSF tracks can not be loaded before the reference data"<<std::endl;
        exit(3);
    }

    // load bw_files
    std::string tmp_ppp_fname = this->ppp_track_fname;
    const size_t bw_path_suffix_pos = tmp_ppp_fname.find("+1");
    if (bw_path_suffix_pos == std::string::npos){
        std::cerr<<"Could not find '+1' in tracks file name. Expecting a name like 'PhyloCSF+1.bw'."<<std::endl;
        exit(3);
    }
    for (uint16_t i = 0; i < 7; ++i){
        std::string suffix;
        if (i == 6) {
            suffix = "power";
        }
        else {
            suffix = ((i < 3) ? "+" : "-") + std::to_string((i % 3) + 1);
        }
        tmp_ppp_fname.replace(bw_path_suffix_pos, 2, suffix); // NOTE: length of "+1" is 2

        this->bw_files[i] = bwOpen(const_cast<char * >(tmp_ppp_fname.c_str()), nullptr, "r");
        if (!this->bw_files[i]){
            // check whether the user has used an unindexed wig file, then print a useful hint
            if (access(tmp_ppp_fname.c_str(), F_OK) == 0 &&
                    tmp_ppp_fname.size() >= 4 && tmp_ppp_fname.compare(tmp_ppp_fname.size() - 4, 4, ".wig") == 0
                    ){
                std::cerr<<"An error occurred while opening the PhyloCSF file "<<tmp_ppp_fname.c_str()<<std::endl;
                std::cerr<<"It seems you provided a *.wig file. You need to simply index them first with wigToBigWig and then use the *.bw files."<<std::endl;
                exit(3);
            }
            else{
                std::cerr<<"Could not find PhyloCSF track file "<<tmp_ppp_fname.c_str()<<std::endl;
                exit(3);
            }
        }
    }

    // load chr_list
    this->chrom_sizes.resize(this->seqid_names.size());
    this->chr_list = this->bw_files[0]->cl;

    int cur_offset=0;
    for (int i = 0; i < chr_list->nKeys; ++i){
        std::string tsi(*this->chr_list->chrom+cur_offset);
        cur_offset+=tsi.size()+1;

        this->n2i_it.first = this->seqnames_ids.find(tsi);
        if(this->n2i_it.first==this->seqnames_ids.end()){
            std::cerr<<"Unable ot find matching chromosome from PhyloCSF tracks in the reference: "<<std::endl;
            exit(3);
        }
        chrom_sizes[this->n2i_it.first->second] = chr_list->len[i];
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

uint Transcriptome::bundleup(bool use_id){ // create bundles and return the total number of bundles
    this->sort(use_id);
    this->bundles.clear();
    this->bundles.push_back(Bundle());

    GFaSeqGet* seqid_seq = nullptr;
    int cur_seqid = -1;
    std::string seqid_name;
    if(this->check_ref && !this->tx_vec.empty()){ // preload reference based on the fron transcript
        cur_seqid = this->tx_vec.front().get_seqid();
        int res = this->seqid2name(cur_seqid,seqid_name);
        assert(res==0);
        seqid_seq = fastaSeqGet(gfasta, seqid_name.c_str());
    }

    for(auto& t : this->tx_vec){
        if(!this->bundles.back().can_add(&t,use_id)){
            if(this->check_ref){
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
uint Transcriptome::clean_cds(bool rescue){
    uint res = 0;

    this->bundleup();
    for(auto& bundle : this->bundles){
        for(auto& tx : bundle){
            if(!tx->has_cds()){continue;}
            tx->load_seq();

            if(this->check_ref){
                if(rescue){
                    tx->rescue_cds(this->allow_non_aug); // TODO: allow non aug - needs to allow storage of non-aug transcripts (only perform cleaning based on stop codons) and no start extension. And needs to handle querrying appropriately
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
// phylo_scores[exon_id][coord_in_exon]
void Transcriptome::compute_PhyloCSF_for_transcript(TX& t, std::array<std::vector<std::vector<float> >, 4> & extracted_scores){
    for (uint8_t i = 0; i < 4; ++i) {
        extracted_scores[i].resize(t.exon_count());
    }

    for (uint16_t exon_id = 0; exon_id < t.exon_count(); ++exon_id){
        const auto & exon = t[exon_id];

        for (uint8_t i = 0; i < 4; ++i) {
            extracted_scores[i][exon_id].resize(exon.get_end() - exon.get_start(), -999.0f);
        }

        bwOverlappingIntervals_t *intervals = nullptr;

        // power
        intervals = bwGetValues(this->bw_files[6], (char*)this->seqid_names[t.get_seqid()].c_str(), exon.get_start(), exon.get_end(), 0);
        if (intervals != nullptr){
            for (uint32_t i = 0; i < (*intervals).l; ++i)
                extracted_scores[3][exon_id][(*intervals).start[i] - exon.get_start()] = (*intervals).value[i];
            bwDestroyOverlappingIntervals(intervals);
            intervals = nullptr;
        }
        if (t.get_strand() == '-') {
            std::reverse(extracted_scores[3][exon_id].begin(), extracted_scores[3][exon_id].end());
        }

        if (t.get_strand() == '+'){
            for (uint8_t phase = 0; phase < 3; ++phase){
                // (cds_phase + exon.beginPos) % 3 == phase
                intervals = bwGetValues(this->bw_files[0 + phase], (char*)this->seqid_names[t.get_seqid()].c_str(), exon.get_start(), exon.get_end(), 0);
                if (intervals != nullptr){
                    // weighten the scores with the power/confidence
                    for (uint32_t i = 0; i < (*intervals).l; ++i){
                        if (extracted_scores[3][exon_id][(*intervals).start[i] - exon.get_start()] == -999.0f){
                            extracted_scores[phase][exon_id][(*intervals).start[i] - exon.get_start()] = -999.0f;
                        }
                        else{
                            extracted_scores[phase][exon_id][(*intervals).start[i] - exon.get_start()] = (*intervals).value[i] * extracted_scores[3][exon_id][(*intervals).start[i] - exon.get_start()];
                        }
                    }
                    bwDestroyOverlappingIntervals(intervals);
                    intervals = nullptr;
                }
            }
        }
        else{ // negative strand
            for (uint8_t phase = 0; phase < 3; ++phase){
                // (chrLen - c.endPos - 1 + cds_phase + 1) % 3 == phase
                intervals = bwGetValues(this->bw_files[3 + phase], (char*)this->seqid_names[t.get_seqid()].c_str(), exon.get_start(), exon.get_end(), 0);
                if (intervals != nullptr){
                    // weighten the scores with the power/confidence
                    for (int32_t i = (*intervals).l - 1; i >= 0; --i){
                        if (extracted_scores[3][exon_id][(exon.get_end() - 1) - (*intervals).start[i]] == -999.0f){
                            extracted_scores[phase][exon_id][(exon.get_end() - 1) - (*intervals).start[i]] = -999.0f;
                        }
                        else{
                            extracted_scores[phase][exon_id][(exon.get_end() - 1) - (*intervals).start[i]] = (*intervals).value[i] * extracted_scores[3][exon_id][(exon.get_end() - 1) - (*intervals).start[i]];
                        }
                    }
                    bwDestroyOverlappingIntervals(intervals);
                    intervals = nullptr;
                }
            }
        }
    }
}
std::pair<SEGTP,float> Transcriptome::compute_tx_cds_ppp(TX& tx){
    const uint64_t chr_len = this->chrom_sizes[tx.get_seqid()];

    // get transcript sequence
    std::string spliced_transcript_seq;
    tx.get_nt_seq(spliced_transcript_seq);

    std::vector<SEGTP> orfs;
    std::array<std::vector<uint32_t>, 3> startpos;
    std::array<std::vector<uint32_t>, 3> stoppos;

    find_all_codons(spliced_transcript_seq, "ATG", startpos);
    find_all_codons(spliced_transcript_seq, "TAA", stoppos);
    find_all_codons(spliced_transcript_seq, "TAG", stoppos);
    find_all_codons(spliced_transcript_seq, "TGA", stoppos);

    for (uint8_t i = 0; i < 3; ++i){
        std::sort(stoppos[i].begin(), stoppos[i].end());
        for (auto start : startpos[i]){
            // translation cannot continue after the first stop-codon
            auto stop_it = std::find_if(stoppos[i].begin(), stoppos[i].end(), [start] (uint32_t element) { return (start < element); });
            if (stop_it != stoppos[i].end()){
                uint32_t stop = *stop_it;
                if (tx.get_strand() == '+'){
                    stop += 2; // include stop codon in sequence

                    uint32_t orf_length = stop - start + 1;
                    assert(orf_length % 3 == 0);
                    if (3 * this->ppp_mincodons <= orf_length){
                        orfs.emplace_back(start, stop);
                    }
                }
                else if (tx.get_strand() == '-'){
                    uint32_t stop_rev = spliced_transcript_seq.size() - start - 3;
                    stop_rev += 2; // include stop codon in sequence
                    uint32_t start_rev = spliced_transcript_seq.size() - stop - 3;

                    uint32_t orf_length = stop_rev - start_rev + 1;
                    assert(orf_length % 3 == 0);

                    if (3 * this->ppp_mincodons <= orf_length){
                        orfs.emplace_back(start_rev, stop_rev);
                    }
                }
            }
        }
    }

    std::array<std::vector<std::vector<float> >, 4> extracted_scores; // phase 0, phase 1, phase 2, power
    compute_PhyloCSF_for_transcript(tx, extracted_scores);

    std::string annotated_cds_seq, computed_cds_seq;
    std::vector<std::tuple<float, std::string, std::vector<CHAIN> > > hits;

    hits.clear();

    std::string longest_or_best_seq;
    bool CDS_outputted_or_found = false;
    bool found_an_orf_satisfying_criteria = false;
    float best_CDS_score = -999.0f;
    SEGTP longest_or_best_CDS; // really only need to store start and stop - can be rebuild afeterthefact - the same for all other evaluators
    std::tuple<float, float> longest_or_best_phylo_stats;
    for (auto const & orf : orfs){
        // get the CDS for the currect ORF
        TX ppp_tx = tx;
        ppp_tx.set_cds_start(orf.get_start());
        ppp_tx.set_cds_end(orf.get_end());
        ppp_tx.build_cds();

        std::tuple<float, float> phylo_stats = this->compute_PhyloCSF(ppp_tx, extracted_scores, ppp_tx.get_strand()=='+'?0:chr_len);

        // output gtf file (with annotation)
        const float phylo_score = std::get<0>(phylo_stats);
        if (phylo_score >= this->ppp_minscore){ // ORF meets criteria
            found_an_orf_satisfying_criteria = true;

            if (this->ppp_mode == "BEST_SCORE"){ // remember what ORF has to best score
                CDS_outputted_or_found = true;
                if (phylo_score > best_CDS_score){
                    best_CDS_score = phylo_score;
                    longest_or_best_CDS = SEGTP(ppp_tx.get_cds_start(),ppp_tx.get_cds_end());
                    longest_or_best_phylo_stats = phylo_stats;
                }
            }
            else if (this->ppp_mode == "LONGEST" && !CDS_outputted_or_found){
                CDS_outputted_or_found = true;
                longest_or_best_CDS = SEGTP(ppp_tx.get_cds_start(),ppp_tx.get_cds_end());
                longest_or_best_phylo_stats = phylo_stats;

                if (ppp_mode == "LONGEST") {
                    break; // we are done (ORFs are sorted by length in descending order)
                }
            }
        }
    }
    if(!CDS_outputted_or_found){
        return std::make_pair(SEGTP(),-1.0);
    }
    return std::make_pair(longest_or_best_CDS,best_CDS_score);
}
std::tuple<float, float> Transcriptome::compute_PhyloCSF(TX& tx,const std::array<std::vector<std::vector<float> >, 4> & extracted_scores, const uint32_t chrLen){
    float total_phylo_sum = 0.0;
    float total_power_sum = 0.0;
    uint32_t total_phylo_count = 0;
    uint32_t total_power_count = 0;

    uint32_t cds_id = 0;

    float phylo_score = NAN; // if score cannot be computed/looked up, it will automatically output NAN
    float phylo_power = NAN;

    int first_exon_id_in_CDS = tx.get_exon_idx(tx.get_cds_start());
    int last_exon_id_in_CDS = tx.get_exon_idx(tx.get_cds_end());

    for (auto it = tx.c_begin(); it != tx.c_end(); ++it, ++cds_id){
        SEGTP& c = *it;

        // new_phase = last CDS has `new_phase` bases extra
        // new_phase2 = new CDS has `new_phase2` bases that we skip (because we exclude them since we already retrieved score from CDS entry before)
        // uint32_t new_phase2 = (3 - new_phase) % 3;
        assert(c.get_phase() < 3);

        uint32_t exon_id;
        uint32_t phylo_start, phylo_end;
        const std::vector<float> * phased_score = nullptr;

        // compute phylo score
        if (tx.get_strand() == '+'){
            exon_id = first_exon_id_in_CDS + cds_id;
            phased_score = &extracted_scores[(c.get_phase() + c.get_start()) % 3][exon_id];
            phylo_start = c.get_start() - tx[exon_id].get_start();
            phylo_end = tx[exon_id].get_end() - c.get_end();
        }
        else{
            // formula derived experimentally from wig files. `start` and `end` are 0-based
            // (chrLen - start) % 3
            // use end-coordinates instead because we start with the last CDS (and we use new_phase2 to deal with codons spanning two CDS entries)
            // (chrLen - (end - 3 + 1)) % 3
            // end-coordinates are stored as half-open interval, hence substract one more
            // (chrLen - (end - 3 + 2)) % 3
            // skip `new_phase2` bases since we considered new_phase = 3 - new_phase2 bases of the same codon in the previous CDS entry already
            // (chrLen - (end - 3 + 2 - new_phase2)) % 3
            exon_id = last_exon_id_in_CDS - cds_id;
            phased_score = &extracted_scores[(chrLen - c.get_end() - 1 + c.get_phase() + 1) % 3][exon_id];
            phylo_start = tx[exon_id].get_end() - c.get_end();
            phylo_end = c.get_start() - tx[exon_id].get_start();
        }

        float phylo_sum = 0.0;
        uint32_t phylo_count = 0;
        for (auto phylo_it = (*phased_score).begin() + phylo_start; phylo_it != (*phased_score).end() - phylo_end; ++phylo_it){
            if (*phylo_it != -999.0f){
                phylo_sum += *phylo_it;
                ++phylo_count;
            }
        }
        total_phylo_count += phylo_count;

        auto & power_scores = extracted_scores[3][exon_id];
        float power_sum = 0.0;
        uint32_t power_count = power_scores.size() - phylo_end - phylo_start;
        total_power_count += power_count;
        for (auto phylo_it = power_scores.begin() + phylo_start; phylo_it != power_scores.end() - phylo_end; ++phylo_it){
            if (*phylo_it != -999.0f)
                power_sum += *phylo_it;
        }

        phylo_score = (phylo_count > 0) ? phylo_sum / power_sum : NAN; // weighten the scores with the power/confidence (or not)
        phylo_power = (power_count > 0) ? power_sum / power_count : NAN;

        total_phylo_sum += phylo_sum;
        total_power_sum += power_sum;
    }

    return std::make_tuple(
            (total_phylo_count > 0) ? total_phylo_sum / total_power_sum : NAN, // weighten the scores with the power/confidence (or not)
            (total_power_count > 0) ? total_power_sum / total_power_count : NAN // total phylo_power
    );
}