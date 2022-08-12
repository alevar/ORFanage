#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
//
// Created by Ales Varabyou
//

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <gff_utils.h>
#include <unordered_set>

#ifndef DEBUG
    #include <omp.h>
#endif

#include "arg_parse.hpp"
#include "Transcriptome.h"

enum find_mode_t{
    LONGEST, // longest complete orf;
    LONGEST_MATCH, // default (highest number of inframe bases shared above thresholds) - if alignment specified - will be decided by alignment instead;
    BEST, // closest to the reference (even if not as long
    VALIDATE // instructs orfanage to use phylocsf tracks for validating novel segments only (does not search for actual ORFs)
};

std::string mode_to_str(const find_mode_t mode) noexcept{
    std::string mode_str;
    switch (mode){
        case LONGEST:
            mode_str = "LONGEST";
            break;
        case LONGEST_MATCH:
            mode_str = "LONGEST_MATCH";
            break;
        case BEST:
            mode_str = "BEST";
            break;
        case VALIDATE:
            mode_str = "VALIDATE";
            break;
        default:
            std::cerr<<"invalid mode selected: "<<mode<<std::endl;
            exit(-7);
    }
    return mode_str;
}

struct Parameters{
    bool clean_query = false;
    bool clean_templ = false;
    bool rescue = false;
    int len_perc_diff = -1; // percent difference by length between the original and reference transcripts. If -1 (default) is set - the check will not be performed.
    int len_frame_perc_diff = -1; // percent difference by length of bases in frame of the reference transcript. If -1 (default) is set - the check will not be performed.
    int len_match_perc_diff = -1; // percent difference by length of bases that are in both query and reference. If -1 (default) is set - the check will not be performed.
    int cds_minlen = -1; // minimum length
    std::vector<find_mode_t> mode_array{LONGEST_MATCH,BEST,LONGEST}; // priority order of modes // TODO: we could also add new mode - match start-end
    int num_threads = 1;

    // Alignment
    int percent_ident = -1; // percent identity
    int gapo = 4; // gap open
    int gape = 2; // gap extend

    std::string reference_fasta_fname;
    std::string query_fname;
    std::vector<std::string> template_fnames = {};

    std::string output_fname;
    std::string stats_fname;

    std::ofstream stats_fp;
    std::ofstream out_gtf_fp;

    // PhyloCSF
    std::string ppp_track_fname;
    find_mode_t ppp_mode = LONGEST;
    int ppp_mincodons = 25;
    float ppp_minscore = 0.0;

    bool use_id = false; // if enabled will use query gene ids to form bundles. the same reference id may be evaluated in multiple bundles if overlaps queries with different gene ids
    bool non_aug = false; // If enabled, non-AUG start codons in reference transcripts will not be discarded and will be considered in overlapping query transcripts on equal grounds with the AUG start codon.
    // TODO: seleno
    bool seleno = false; // If enabled, all premature stop codons in the reference will be treated as selenocysteine instead and if matched in query - will be retained, but only if the exact position is matched in the same frame.
    bool keep_cds = false; // If enabled, any CDS already presernt in the query will be kept unmodified

} global_params;

typedef std::vector<std::vector<std::pair<TX,Score>>> Scores;

struct RunStats{
    uint num_too_short = 0; // number of template isoforms with ORF that is too short
    uint num_dirty = 0; // number of template isoforms with incorrect/unsuitable ORFs (other than short)
    uint template_duplicates = 0; // number of duplcate ORFs discarded
} rstats;

int run(){
    // load the reference GFF
    std::cerr<<"loading reference genome"<<std::endl;
    Transcriptome transcriptome;
    if(!global_params.reference_fasta_fname.empty()){
        transcriptome.set_ref(global_params.reference_fasta_fname);
    }
    if(global_params.non_aug){
        transcriptome.use_non_aug();
    }
    if(global_params.percent_ident>-1){
        assert(!global_params.reference_fasta_fname.empty());
        transcriptome.set_aligner(mat50,aa_table,global_params.gapo,global_params.gape);
    }

    std::cerr<<"loading reference transcriptomes"<<std::endl;
    for(auto& k : global_params.template_fnames){
        transcriptome.add(k,true,true);
    }

    transcriptome.build_cds_chains();

    if(global_params.clean_templ || !global_params.reference_fasta_fname.empty()){ // under these conditions we need to make sure that the chains are len%3==0
        transcriptome.correct_chain_len();
    }

    if(global_params.clean_templ) {
        rstats.num_too_short = transcriptome.clean_short_orfs(global_params.cds_minlen);
    }

    std::cerr<<"sorting reference transcriptome"<<std::endl;
    transcriptome.sort();

    if(global_params.clean_templ) {
        rstats.num_dirty = transcriptome.clean_cds(global_params.rescue); // TODO: if we load sequence data into bundles - this is not cost effective...
    }

    transcriptome.set_cds_as_exons();
    transcriptome.remove_non_coding();

    rstats.template_duplicates = transcriptome.deduplicate(global_params.use_id);

    std::cerr<<"loading query transcriptome"<<std::endl;
    transcriptome.add(global_params.query_fname,false,false);

    if(!global_params.ppp_track_fname.empty()){
        std::cerr<<"loading phylocsf tracks"<<std::endl;
        transcriptome.set_ppp(global_params.ppp_track_fname, mode_to_str(global_params.ppp_mode),global_params.ppp_mincodons,global_params.ppp_minscore);
        transcriptome.load_ppptracks();
    }

    std::cerr<<"bundling transcriptome"<<std::endl;
    transcriptome.bundleup(global_params.use_id);

    std::cerr<<"starting main evaluation"<<std::endl;
#ifndef DEBUG
    #pragma omp parallel for
#endif
    for(auto bundle_it=transcriptome.bbegin(); bundle_it!=transcriptome.bend(); bundle_it++){ // iterate over bundles

        std::vector<std::vector<std::tuple<SEGTP,TX,TX,Score,std::string>>> stats; // segment,query,template,score,notes
        TX *q,*t;
        Finder* fndr = global_params.percent_ident>-1 ? transcriptome.get_aligner() : nullptr;
        CHAIN segments;

        if(!bundle_it->has_template() || !bundle_it->has_query()){
            // todo: still need to output into the gtf and stats
            if(global_params.stats_fp.is_open()){
                for(int qi=0;qi<bundle_it->size();qi++) {
                    q = bundle_it->operator[](qi);
//                    std::cout<<q->get_tid()<<std::endl;
#ifdef DEBUG
                    if(std::strcmp(q->get_tid().c_str(),"CHS.11626.7")==0){ // rna-XM_011520617.2
                        std::cout<<"found"<<std::endl;
                    }
#endif
                    if(q->is_template()){continue;}
                    std::string cur_seqid;
                    transcriptome.seqid2name(q->get_seqid(),cur_seqid);
#ifndef DEBUG
#pragma omp critical
#endif
                    {
                        if(global_params.keep_cds && q->has_cds()) {
                            q->build_cds();
                        }
                        global_params.out_gtf_fp << q->str(cur_seqid) << std::endl;
                        global_params.stats_fp << q->get_tid() << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << "\t"
                                               << "-" << std::endl;
                    }
                }
            }
            continue;
        }

        for(int qi=0;qi<bundle_it->size();qi++){
            stats.clear();
            q=bundle_it->operator[](qi);
#ifdef DEBUG
            if(std::strcmp(q->get_tid().c_str(),"CHS.11626.7")==0){ // rna-XM_011520617.2
                std::cout<<"found"<<std::endl;
            }
#endif
            if(q->is_template()){continue;}
            if(global_params.keep_cds && q->has_cds()){
                q->build_cds();
                std::string cur_seqid;
                transcriptome.seqid2name(q->get_seqid(),cur_seqid);
#ifndef DEBUG
#pragma omp critical
#endif
                {
                    global_params.out_gtf_fp << q->str(cur_seqid) << std::endl;
                    global_params.stats_fp << q->get_tid() << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "keep_cds" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << std::endl;
                }
                continue;
            }
            q->remove_cds();

            for(int ti=0;ti<bundle_it->size();ti++){
                t=bundle_it->operator[](ti);
                if(!t->is_template()){continue;}
                int intlen = q->exon_chain()->intersection(*t->cds_chain(),segments,true);
                stats.emplace_back(std::vector<std::tuple<SEGTP,TX,TX,Score,std::string>>{});
                if(intlen==0){
                    stats.back().emplace_back(std::make_tuple(SEGTP(),*q,*t,Score(),"no-overlap")); // empty score
                    std::get<3>(stats.back().back()).pass = false;
                    continue;
                }

                std::unordered_set<int> orf_starts;
                std::pair<std::unordered_set<int>::iterator,bool> os_rit;

                // now we need to reconstruct a transcript within the query exon chain for each segment of the template ORF
                for(auto& s : segments){
                    stats.back().emplace_back(std::make_tuple(s,*q,*t,Score(),"-"));
                    // TODO: does it work without sequence avaialble?

#ifdef DEBUG
                    if(std::strcmp(q->get_tid().c_str(),"CHS.11626.7")==0){ // rna-XM_011520617.2
                        std::cout<<"found"<<std::endl;
                    }
#endif

                    int seg_phase = t->cds_chain()->get_phase(t->get_strand()=='+'?s.get_start():s.get_end(),t->get_strand());
#ifdef DEBUG
                    assert(seg_phase>=0 && seg_phase<3);
#endif
                    if(s.slen()<3){
                        std::get<4>(stats.back().back())="segment_len<3";
                        std::get<3>(stats.back().back()).pass = false;
                        continue;
                    }
                    s.set_phase(seg_phase);
                    TX qseg = *q;
                    qseg.set_cds_start(s.get_start());
                    qseg.set_cds_end(s.get_end());
                    qseg.set_cds_phase(s.get_phase());
                    qseg.build_cds();
                    if(qseg.cds_len()<3){
                        std::get<1>(stats.back().back())=qseg;
                        std::get<4>(stats.back().back())="cds_len<3";
                        std::get<3>(stats.back().back()).pass = false;
                        continue;
                    }
                    qseg.correct_chain_len();
                    if(qseg.cds_len()<3){
                        std::get<1>(stats.back().back())=qseg;
                        std::get<4>(stats.back().back())="cds_len<3";
                        std::get<3>(stats.back().back()).pass = false;
                        continue;
                    }

                    // if reference is provided
                    if(!global_params.reference_fasta_fname.empty()) {
                        if(qseg.cds_len()<3){
                            std::get<1>(stats.back().back())=qseg;
                            std::get<4>(stats.back().back())="cds_len<3";
                            std::get<3>(stats.back().back()).pass = false;
                            qseg.remove_cds();
                            continue;
                        }
                        qseg.load_seq();

                        int ret = qseg.rescue_cds(global_params.non_aug,t);
                        std::get<1>(stats.back().back())=qseg;
                        if(!ret){
                            std::get<4>(stats.back().back())="not_rescued";
                            std::get<3>(stats.back().back()).pass = false;
                            continue;
                        }

                        if(qseg.get_aa().back()!='.' ||
                           (!global_params.non_aug && qseg.get_aa().front()!='M') ||
                                (global_params.non_aug && qseg.get_aa().front()!='M' && !(qseg.get_strand()=='+' ? qseg.get_cds_start()==t->get_cds_start() : qseg.get_cds_end()==t->get_cds_end()))){
                            std::get<3>(stats.back().back()).pass = false;
                            if(qseg.get_aa().back()!='.'){
                                std::get<4>(stats.back().back())="missing_stop";
                            }
                            else{
                                std::get<4>(stats.back().back())="missing_start";
                            }
                            qseg.remove_cds();
                            continue;
                        }
                    }

#ifdef DEBUG
                    assert(qseg.get_aa().find(".")==qseg.get_aa().size()-1); // assert that no premature stop-codons are found
#endif

                    if(qseg.cds_chain()->clen()<global_params.cds_minlen){
                        std::get<4>(stats.back().back())="minlen";
                        std::get<3>(stats.back().back()).pass = false;
                        continue;
                    }

                    Score qseg_score = qseg.score(*t);
                    std::get<3>(stats.back().back())=qseg_score;
                    if(qseg_score.lpd<global_params.len_perc_diff ||
                       qseg_score.ilpd<global_params.len_frame_perc_diff ||
                       qseg_score.mlpd<global_params.len_match_perc_diff){
                        std::get<3>(stats.back().back()).pass = false;
                        continue;
                    }

                    os_rit = orf_starts.insert(qseg.get_cds_start());
                    if(!os_rit.second){ // the ORF previously existed
                        std::get<4>(stats.back().back())="dup";
                        std::get<3>(stats.back().back()).pass = false;
                        continue;
                    }

                    // otherwise - we can safely add to the curent evaluation
                    qseg.store_template(t);
                    std::get<1>(stats.back().back())=qseg;
                }
#ifdef DEBUG
                // can we select a single segment to proceed with at this point?
                if(stats.back().size()>1){ // TODO: superflous since we are keeping track of all segments - even duplicates here...
                    std::cerr<<"found more than one valid segment"<<std::endl;
                }
#endif
                segments.clear();
            }

            // check percent identity if requested
            if(global_params.percent_ident>-1) {
                for(auto& sc : stats){
                    for(auto& seg : sc){
                        if(!std::get<3>(seg).pass){continue;}

                        if (!std::get<1>(seg).get_template()->seq_loaded()) {std::get<1>(seg).get_template()->load_seq();} // only load sequence if not previously loaded
                        ksw_extz_t ez;
                        fndr->align(std::get<1>(seg).get_template()->get_aa().c_str(), std::get<1>(seg).get_aa().c_str(), ez);
                        if (ez.n_cigar == 0) {
                            std::get<3>(seg).pass =false;
                        }
                        else{
                            int ret = fndr->parse(ez, std::get<1>(seg).get_aa(), std::get<1>(seg).get_template()->get_aa(), std::get<3>(seg));
                            std::get<3>(seg).aln_match = std::get<3>(seg).num_aln_match();
                            std::get<3>(seg).aln_pi = std::get<3>(seg).pi();
                            if (ret < 0 || std::get<3>(seg).aln_pi < global_params.percent_ident) {
                                std::get<3>(seg).pass = false;
                            }
                        }
                        free(ez.cigar);
                    }
                }
            }

            // TODO: PPP should realistically only run if we don't find anything else
            //   but also when validation is enabled
            //   or as a validation for any novel segments of the ORF
            //   PPP does not require alignment (plus we don't know which template to use for alignment in the first place)
            //   We only need to run if no good/suitable orfanage result is available
            //   if phylocsf is requested - compute for the current qseg and add all to the current list for evaluation
            if(!global_params.ppp_track_fname.empty()){
                if(global_params.ppp_mode==LONGEST || global_params.ppp_mode==BEST){ // todo: add to the stats
                    std::pair<SEGTP,float> ppp_chain = transcriptome.compute_tx_cds_ppp(*q);

                    TX qseg = *q;
                    qseg.set_cds_start(std::get<0>(ppp_chain).get_start());
                    qseg.set_cds_end(std::get<0>(ppp_chain).get_end());
                    qseg.set_cds_phase(std::get<0>(ppp_chain).get_phase());
                    qseg.build_cds();

                    if(qseg.cds_chain()->clen()<global_params.cds_minlen){continue;}

                    Score qseg_score = qseg.score(*t);
                    qseg_score.ppp_score = std::get<1>(ppp_chain);
                    stats.back().emplace_back(std::make_tuple(std::get<0>(ppp_chain),qseg,*t,qseg_score,"ppp"));
                    if(qseg_score.lpd<global_params.len_perc_diff ||
                       qseg_score.ilpd<global_params.len_frame_perc_diff ||
                       qseg_score.mlpd<global_params.len_match_perc_diff){
                        std::get<3>(stats.back().back()).pass = false;
                        continue;
                    }
                }
                else if(global_params.ppp_mode==VALIDATE){ // only perform validation of novel transcripts
                    const uint64_t chr_len = transcriptome.get_chrom_len(q->get_seqid());
                    std::array<std::vector<std::vector<float> >, 4> extracted_scores; // phase 0, phase 1, phase 2, power
                    transcriptome.compute_PhyloCSF_for_transcript(*q, extracted_scores);(*q, extracted_scores);

                    for(auto& sc : stats){
                        for(auto& seg : sc){
                            if(!std::get<3>(seg).pass){continue;}

                            bool ppp_pass = true;
                            std::tuple<float, float> ppp_scores = transcriptome.compute_PhyloCSF(std::get<1>(seg), extracted_scores,q->get_strand()=='+'?0:chr_len);
                            std::get<3>(seg).ppp_score = std::get<0>(ppp_scores);
                            if(std::get<0>(ppp_scores) < global_params.ppp_minscore){
                                ppp_pass = false;
                            }

                            if (!ppp_pass) {
                                std::get<3>(seg).pass=false;
                            } // remove
                        }
                    }
                }
                else{
                    std::cerr<<"unknown mode selected for PPP: "<<mode_to_str(global_params.ppp_mode)<<std::endl;
                    exit(-7);
                }
            }

            std::string cur_seqid;
            transcriptome.seqid2name(q->get_seqid(),cur_seqid);

            int template_comp_id = -1; // index of the best template in the list
            int segment_comp_id = -1; // index of the best segment for the best template
            Score best_score;

            int tci = 0;
            for(auto& t_sc : stats){
                int sci = 0;
                for(auto& seg : t_sc){
                    // write stats regardless of whether it is chosen or not

                    if(std::get<3>(seg).pass){
                        // evaluate and pick the best choice based on a strategy
                        for(auto& m : global_params.mode_array){
                            bool found = false;
                            switch(m){
                                case LONGEST: // pick the longest
                                    if(best_score.qlen < std::get<3>(seg).qlen){
                                        template_comp_id = tci;
                                        segment_comp_id = sci;
                                        best_score = std::get<3>(seg);
                                        found = true;
                                    }
                                    if(best_score.qlen > std::get<3>(seg).qlen){
                                        found = true;
                                    }
                                    break;
                                case LONGEST_MATCH: //
                                    if(best_score.num_bp_inframe < std::get<3>(seg).num_bp_inframe){
                                        template_comp_id = tci;
                                        segment_comp_id = sci;
                                        best_score = std::get<3>(seg);
                                        found = true;
                                    }
                                    if(best_score.num_bp_inframe > std::get<3>(seg).num_bp_inframe){
                                        found = true;
                                    }
                                    break;
                                case BEST: // PI if alignment enabled - otherwise ilpd
                                    if(global_params.percent_ident==-1){ // no alignment was being performed - use ilpd
                                        if(best_score.ilpd < std::get<3>(seg).ilpd){
                                            template_comp_id = tci;
                                            segment_comp_id = sci;
                                            best_score = std::get<3>(seg);
                                            found = true;
                                        }
                                        if(best_score.ilpd > std::get<3>(seg).ilpd){
                                            found = true;
                                        }
                                    }
                                    else{
                                        if(best_score.aln_pi < std::get<3>(seg).aln_pi){
                                            template_comp_id = tci;
                                            segment_comp_id = sci;
                                            best_score = std::get<3>(seg);
                                        }
                                        if(best_score.aln_pi > std::get<3>(seg).aln_pi){
                                            found = true;
                                        }
                                    }

                                    break;
                                default:
                                    std::cerr<<"wrong mode selected"<<std::endl;
                                    exit(-1);
                            }
                            if(found){ // found a better candidate - no need to search further down the mode priority array
                                break;
                            }
                        }
                    }

                    sci++;
                }
                tci++;
            }
#ifndef DEBUG
#pragma omp critical
#endif
            {
                if(template_comp_id>=0 && segment_comp_id>=0){ // found something
                    std::get<1>(stats[template_comp_id][segment_comp_id]).add_attribute("orfanage_status","1");
                    std::get<1>(stats[template_comp_id][segment_comp_id]).add_attribute("orfanage_template",std::get<1>(stats[template_comp_id][segment_comp_id]).get_template()->get_tid());
                    global_params.out_gtf_fp<<std::get<1>(stats[template_comp_id][segment_comp_id]).str(cur_seqid)<<std::endl;
                    std::get<4>(stats[template_comp_id][segment_comp_id]) = "gtf";
                }
                else{
                    // nothing found - exit
                    q->add_attribute("orfanage_status","0");
                    global_params.out_gtf_fp<<q->str(cur_seqid)<<std::endl;
                }
            };
            // TODO: need an option to write out discarded transcripts (those that overlap CDS but without a valid CDS)

            // Lastly, write stats about each query
#ifndef DEBUG
#pragma omp critical
#endif
            {
                if(stats.empty()){
                    global_params.stats_fp<<q->get_tid()<<"\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" << "\t"
                                          << "-" <<std::endl;
                }
                else{
                    for(auto& grp : stats){
                        for(auto& seg : grp){
                            global_params.stats_fp<<std::get<1>(seg).get_tid()<<"\t" // query
                                                  <<std::get<2>(seg).get_tid()<<"\t" // template
                                                  <<std::get<0>(seg)<<"\t" // segment
                                                  <<std::get<4>(seg)<<"\t" // notes
                                                  <<std::get<3>(seg)<<std::endl; // score
                        }
                    }
                }
            }
        }
    }

    return 0;
}

int main(int argc, char** argv) {
    ArgParse args("orfanage",
                  "Annotating Open Reading Frames based on reference, phylogeny and sequence similarity.");

    args.add_option("query",ArgParse::Type::STRING, "Path to a GTF query file with transcripts to which CDSs are to be ported", ArgParse::Level::GENERAL,true);
    args.add_option("output", ArgParse::Type::STRING, "Basename for all output files generated by this software", ArgParse::Level::GENERAL,true);
    args.add_option("reference",ArgParse::Type::STRING,"Path to the reference genome file in FASTA format. This parameter is required when the following parameters are used: 1. cleanq; 2. cleant; 3. pd.",ArgParse::Level::GENERAL,false);
    args.add_option("cleanq",ArgParse::Type::FLAG,"If enabled - will ensure all transcripts in the output file will have a valid start and end codons. This option requires the use of --reference parameter",ArgParse::Level::GENERAL,false);
    args.add_option("cleant",ArgParse::Type::FLAG,"If enabled - will ensure all ORFs in the reference annotations start with a valid start codon and end with the first available stop codon. This option requires the use of --reference parameter",ArgParse::Level::GENERAL,false);
    args.add_option("rescue",ArgParse::Type::FLAG,"If enabled - will attempt rescuing the broken ORFs in the reference annotations. This option requires the use of --reference parameter",ArgParse::Level::GENERAL,false);
    args.add_option("lpd",ArgParse::Type::INT,"Percent difference by length between the original and reference transcripts. If -1 (default) is set - the check will not be performed.",ArgParse::Level::GENERAL,false);
    args.add_option("ilpd",ArgParse::Type::INT,"Percent difference by length of bases in frame of the reference transcript. If -1 (default) is set - the check will not be performed.",ArgParse::Level::GENERAL,false);
    args.add_option("mlpd",ArgParse::Type::INT,"Percent difference by length of bases that are in both query and reference. If -1 (default) is set - the check will not be performed.",ArgParse::Level::GENERAL,false);
    args.add_option("minlen",ArgParse::Type::INT,"Minimum length of an open reading frame to consider for the analysis",ArgParse::Level::GENERAL,false);
    args.add_option("mode", ArgParse::Type::STRING, "Which CDS to report: LONGEST, LONGEST_MATCH, BEST. Default: " + mode_to_str(global_params.mode_array.front()), ArgParse::Level::GENERAL, false);
    args.add_option("stats",ArgParse::Type::STRING,"Output a separate file with stats for each query/template pair",ArgParse::Level::GENERAL,false);
    args.add_option("threads",ArgParse::Type::INT,"Number of threads to run in parallel",ArgParse::Level::GENERAL,false);
    args.add_option("use_id",ArgParse::Type::FLAG,"If enabled, only transcripts with the same gene ID from the query file will be used to form a bundle. In this mode the same template transcript may be used in several bundles, if overlaps transcripts with different gene_ids.",ArgParse::Level::GENERAL,false);
    args.add_option("non_aug",ArgParse::Type::FLAG,"If enabled, non-AUG start codons in reference transcripts will not be discarded and will be considered in overlapping query transcripts on equal grounds with the AUG start codon.",ArgParse::Level::GENERAL,false);
    args.add_option("keep_cds",ArgParse::Type::FLAG,"If enabled, any CDS already presernt in the query will be kept unmodified.",ArgParse::Level::GENERAL,false);

    // Alignment
    args.add_option("pi",ArgParse::Type::INT,"Percent identity between the query and template sequences. This option requires --reference parameter to be set. If enabled - will run alignment between passing pairs.", ArgParse::Level::GENERAL,false);
    args.add_option("gapo",ArgParse::Type::INT,"Gap-open penalty", ArgParse::Level::GENERAL,false);
    args.add_option("gape",ArgParse::Type::INT,"Gap-extension penalty", ArgParse::Level::GENERAL,false);

    // PhyloCSF++
    args.add_option("ppp_mode", ArgParse::Type::STRING, "Which CDS to report: LONGEST, BEST. Default: " + mode_to_str(global_params.ppp_mode), ArgParse::Level::GENERAL, false);
    args.add_option("min-score", ArgParse::Type::FLOAT, "Only consider ORFs with a minimum weighted PhyloCSF mean score (range from -15 to +15, >0 more likely to be protein-coding). Default: " + std::to_string(global_params.ppp_minscore), ArgParse::Level::GENERAL, false);
    args.add_option("min-codons", ArgParse::Type::INT, "Only consider ORFs with a minimum codon length. Default: " + std::to_string(global_params.ppp_mincodons), ArgParse::Level::GENERAL, false);
    args.add_option("tracks", ArgParse::Type::STRING, "Path to the bigWig file PhyloCSF+1.bw (expects the other 5 frames to be in the same directory, optionally the power track).",ArgParse::Level::GENERAL, false);

    args.add_positional_argument("templates", ArgParse::Type::STRING, "One or more GFF/GTF files with coding exons to be used as templates.", true, true);

    args.add_option("help", ArgParse::Type::FLAG, "Prints this help message.", ArgParse::Level::HELP, false);

    args.parse_args(argc, argv);
    args.check_args();

    // first create the execution string
    std::string cl = "orfanage ";
    for (int i = 0; i < argc; i++) {
        if (i == 0) {
            cl += argv[i];
        } else {
            cl += " ";
            cl += argv[i];
        }
    }

    // check that all files exist
    global_params.query_fname = args.get_string("query");
    if(!file_exists(global_params.query_fname)){
        std::cerr << "Input file does not exist! "<<args.get_string("query")<<std::endl;
        exit(2);
    }

    // run for every gff file
    for (size_t i = 0; i < args.positional_argument_size(); ++i) {
        std::string templ_fname = args.get_positional_argument(i);
        global_params.template_fnames.push_back(templ_fname);
        if (!file_exists(templ_fname)) {
            std::cerr << "Template file does not exist! " << templ_fname << std::endl;
            exit(2);
        }
    }

    if(args.is_set("reference")){ // much from gpertea bamcons.cpp
        global_params.reference_fasta_fname = args.get_string("reference");
        if(!file_exists(global_params.reference_fasta_fname)){
            std::cerr << "Reference FASTA file does not exist!\n";
            exit(2);
        }

        // get potential fasta index file name
        std::string fa_idx_fname = args.get_string("reference")+".fai";
        GFastaIndex faIdx(args.get_string("reference").c_str(), fa_idx_fname.c_str());
        if (!faIdx.hasIndex()){
            std::cerr<<"No fasta index found for "<<fa_idx_fname<<". Building now"<<std::endl;
            faIdx.buildIndex();
            if (faIdx.getCount() == 0){
                std::cerr<<"Error: no fasta records found!"<<std::endl;
                exit(2);
            }
            FILE* fcreate = fopen(fa_idx_fname.c_str(), "w");
            if (fcreate == nullptr){
                std::cerr<<"Error: cannot create fasta index: "<<fa_idx_fname<<std::endl;
                exit(2);
            }
            if (faIdx.storeIndex(fcreate) < faIdx.getCount()){
                std::cerr<<"Warning: error writing the index file!"<<std::endl;
                exit(2);
            }
            std::cerr<<"FASTA index rebuilt."<<std::endl;
        }
    }

    // Now check for PhyloCSF parameters
    // This is an adaptation of code by Christopher Pockrandt from PhyloCSF++
    if(args.is_set("tracks")){
        std::string bw_path = args.get_string("tracks");
        const size_t bw_path_suffix_pos = bw_path.find("+1");
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
            bw_path.replace(bw_path_suffix_pos, 2, suffix); // NOTE: length of "+1" is 2
            if(!file_exists(bw_path)){
                std::cerr << "PhyloCSF track is not found: "<<bw_path<<std::endl;
                exit(3);
            }
        }

        global_params.ppp_track_fname = args.get_string("tracks");
    }

    global_params.ppp_mincodons = args.is_set("min-codons") ? args.get_int("min-codons") : 25;
    global_params.ppp_minscore = args.is_set("min-score") ? args.get_float("min-score") : 0.0;

    if (args.is_set("mode")){
        global_params.mode_array.clear();

        std::string mode = args.get_string("mode");
        std::transform(mode.begin(), mode.end(), mode.begin(), toupper);

        if (mode == "LONGEST") {
            global_params.mode_array.push_back(LONGEST);
            global_params.mode_array.push_back(LONGEST_MATCH);
            global_params.mode_array.push_back(BEST);
        }
        else if (mode == "LONGEST_MATCH") {
            global_params.mode_array.push_back(LONGEST_MATCH);
            global_params.mode_array.push_back(BEST);
            global_params.mode_array.push_back(LONGEST);
        }
        else if (mode == "BEST") {
            global_params.mode_array.push_back(BEST);
            global_params.mode_array.push_back(LONGEST_MATCH);
            global_params.mode_array.push_back(LONGEST);
        }
        else{
            printf(OUT_ERROR "Please choose a valid mode (LONGEST, LONGEST_MATCH, BEST)!\n" OUT_RESET);
            return -1;
        }
    }

    global_params.ppp_mode = LONGEST;
    if (args.is_set("ppp_mode")){
        std::string ppp_mode = args.get_string("ppp_mode");
        std::transform(ppp_mode.begin(), ppp_mode.end(), ppp_mode.begin(), toupper);

        if (ppp_mode == "LONGEST")
            global_params.ppp_mode = LONGEST;
        else if (ppp_mode == "BEST")
            global_params.ppp_mode = BEST;
        else if (ppp_mode == "VALIDATE")
            global_params.ppp_mode = VALIDATE;
        else{
            printf(OUT_ERROR "Please choose a valid ppp_mode (LONGEST, BEST, VALIDATE)!\n" OUT_RESET);
            return -1;
        }
    }

    global_params.clean_query = args.is_set("cleanq");
    global_params.clean_templ = args.is_set("cleant");
    global_params.rescue = args.is_set("rescue");
    if(global_params.clean_templ || global_params.clean_query || global_params.rescue){
        assert(args.is_set("reference"));
    }

    global_params.cds_minlen = args.is_set("minlen") ? args.get_int("minlen") : 36;
    global_params.len_perc_diff = args.is_set("lpd") ? args.get_int("lpd") : 50;
    global_params.len_frame_perc_diff = args.is_set("ilpd") ? args.get_int("ilpd") : 50;
    global_params.len_match_perc_diff = args.is_set("mlpd") ? args.get_int("mlpd") : 50;
    global_params.percent_ident = args.is_set("pi") ? args.get_int("pi") : -1;

    global_params.gapo = args.is_set("gapo") ? args.get_int("gapo") : 4;
    global_params.gape = args.is_set("gape") ? args.get_int("gape") : 2;

    global_params.num_threads = args.is_set("threads") ? args.get_int("threads") : 1;
#ifndef DEBUG
    omp_set_num_threads(global_params.num_threads);
#endif

    if(args.is_set("use_id")){
        global_params.use_id = true;
    }

    if(args.is_set("non_aug")){
        global_params.non_aug = true;
    }

    if(args.is_set("keep_cds")){
        global_params.keep_cds = true;
    }

    global_params.output_fname = args.get_string("output");
    global_params.out_gtf_fp.open(global_params.output_fname);

    global_params.stats_fname = args.get_string("stats");
    if(args.is_set("stats")){
        global_params.stats_fp.open(global_params.stats_fname);
        Score s;
        global_params.stats_fp<<"query_id\t"
                                "template_id\t"
                                "segment\t"
                                "notes\t"
                                <<s.stats_header()<<std::endl;
    }

    run();

    global_params.stats_fp.close();
    global_params.out_gtf_fp.close();

    return 0;
}
#pragma clang diagnostic pop