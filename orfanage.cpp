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
    #ifdef OPENMP_AVAILABLE
        #include <omp.h>
    #endif
#endif

#include "arg_parse.hpp"
#include "Transcriptome.h"

enum find_mode_t{
    START_MATCH, // if an ORF with a mathcing start codon is available - report it
    ALL, // report all complete orfs
    LONGEST, // longest complete orf;
    LONGEST_MATCH, // default (highest number of inframe bases shared above thresholds) - if alignment specified - will be decided by alignment instead;
    BEST, // closest to the reference (even if not as long
    FIRST, // most upstream ORF
};

std::string mode_to_str(const find_mode_t mode) noexcept{
    std::string mode_str;
    switch (mode){
        case START_MATCH:
            mode_str = "MATCH_START";
            break;
        case ALL:
            mode_str = "ALL";
            break;
        case LONGEST:
            mode_str = "LONGEST";
            break;
        case LONGEST_MATCH:
            mode_str = "LONGEST_MATCH";
            break;
        case BEST:
            mode_str = "BEST";
            break;
        case FIRST:
            mode_str = "FIRST";
            break;
        default:
            std::cerr<<"invalid mode selected: "<<mode<<std::endl;
            exit(-7);
    }
    return mode_str;
}

struct DEFAULTS{
    int len_perc_ident = -1; // percent identity by length between the original and reference transcripts. If -1 (default) is set - the check will not be performed.
    int len_frame_perc_ident = -1; // percent identity by length of bases in frame of the reference transcript. If -1 (default) is set - the check will not be performed.
    int len_match_perc_ident = -1; // percent identity by length of bases that are in both query and reference. If -1 (default) is set - the check will not be performed.
    int cds_minlen = -1; // minimum length
    int num_threads = 1;

    int overhang = 0;
} def_params;

struct Parameters{
    bool clean_query = false;
    bool clean_templ = false;
    bool rescue = false;
    int len_perc_ident = -1; // percent identity by length between the original and reference transcripts. If -1 (default) is set - the check will not be performed.
    int len_frame_perc_ident = -1; // percent identity by length of bases in frame of the reference transcript. If -1 (default) is set - the check will not be performed.
    int len_match_perc_ident = -1; // percent identity by length of bases that are in both query and reference. If -1 (default) is set - the check will not be performed.
    int cds_minlen = -1; // minimum length
    std::vector<find_mode_t> mode_array{LONGEST_MATCH,BEST,START_MATCH,LONGEST,FIRST,ALL}; // priority order of modes
    int num_threads = 1;

    int overhang = 0;
    bool spliced_overhang = false;

    std::string reference_fasta_fname;
    std::string query_fname;
    std::vector<std::string> template_fnames = {};

    std::string output_fname;
    std::string stats_fname;

    std::ofstream stats_fp;
    std::ofstream out_gtf_fp;

    bool use_id = false; // if enabled will use query gene ids to form bundles. the same reference id may be evaluated in multiple bundles if overlaps queries with different gene ids
    bool non_aug = false; // If enabled, non-AUG start codons in reference transcripts will not be discarded and will be considered in overlapping query transcripts on equal grounds with the AUG start codon.
    // TODO: seleno
    bool seleno = false; // If enabled, all premature stop codons in the reference will be treated as selenocysteine instead and if matched in query - will be retained, but only if the exact position is matched in the same frame.
    bool keep_all_cds = false; // If enabled, any CDS already present in the query will be kept unmodified
    bool keep_cds_if_not_found = false; // if enabled CDS present in query transcripts when no ORF can be detected will be preserved.
} global_params;

struct RunStats{
    uint num_too_short = 0; // number of template isoforms with ORF that is too short
    uint num_dirty = 0; // number of template isoforms with incorrect/unsuitable ORFs (other than short)
    uint template_duplicates = 0; // number of duplcate ORFs discarded
} rstats;

bool score_lt(const Score& lhs, const Score& rhs){
    if(lhs.pass!=rhs.pass){
        return lhs.pass < rhs.pass;
    }
    for(auto& m : global_params.mode_array){
        switch(m){
            case START_MATCH:
                if(lhs.start_match != rhs.start_match){
                    return lhs.start_match < rhs.start_match; // return true if left start does not match
                }
                break;
            case ALL: // skip - only relevant later
                break;
            case LONGEST: // pick the longest
                if(lhs.qlen != rhs.qlen){
                    return lhs.qlen < rhs.qlen;
                }
                break;
            case LONGEST_MATCH:
                if(lhs.num_bp_inframe != rhs.num_bp_inframe){
                    return lhs.num_bp_inframe < rhs.num_bp_inframe;
                }
                break;
            case BEST: // PI if alignment enabled - otherwise ilpi
                if(lhs.ilpi != rhs.ilpi){
                    return lhs.ilpi < rhs.ilpi;
                }
                break;
            case FIRST:
                if(lhs.query_start != rhs.query_start){
                    return lhs.query_start > rhs.query_start; // here the gt is used since smaller query start is better (more upstream)
                }
                break;
            default:
                std::cerr<<"wrong mode selected"<<std::endl;
                exit(-1);
        }
    }
    return lhs.tlen < rhs.tlen;
}

bool score_gt(const Score& lhs, const Score& rhs){
    if(lhs.pass!=rhs.pass){
        return lhs.pass > rhs.pass;
    }
    for(auto& m : global_params.mode_array){
        switch(m){
            case START_MATCH:
                if(lhs.start_match != rhs.start_match){
                    return lhs.start_match > rhs.start_match; // return true if left start does match
                }
                break;
            case ALL: // skip - only relevant later
                break;
            case LONGEST: // pick the longest
                if(lhs.qlen != rhs.qlen){
                    return lhs.qlen > rhs.qlen;
                }
                break;
            case LONGEST_MATCH:
                if(lhs.num_bp_inframe != rhs.num_bp_inframe){
                    return lhs.num_bp_inframe > rhs.num_bp_inframe;
                }
                break;
            case BEST: // PI if alignment enabled - otherwise ilpi
                if(lhs.ilpi != rhs.ilpi){
                    return lhs.ilpi > rhs.ilpi;
                }
                break;
            case FIRST:
                if(lhs.query_start != rhs.query_start){
                    return lhs.query_start < rhs.query_start; // here the lt is used since smaller query start is better (more upstream)
                }
                break;
            default:
                std::cerr<<"wrong mode selected"<<std::endl;
                exit(-1);
        }
    }
    return lhs.tlen > rhs.tlen;
}

// groups duplicate ORFs and selects representative template (maximizing score within each group)
// returns map of start/end coordinates as keys to a vector of query/template pairs with scores
std::map<std::pair<int,int>,std::vector<std::tuple<SEGTP,TX,TX,Score,std::string>>> flatten(std::vector<std::vector<std::tuple<SEGTP,TX,TX,Score,std::string>>>& stats){
    std::map<std::pair<int,int>,std::vector<std::tuple<SEGTP,TX,TX,Score,std::string>>> flat;
    std::pair<std::map<std::pair<int,int>,std::vector<std::tuple<SEGTP,TX,TX,Score,std::string>>>::iterator,bool> fit;
    for(auto& grp : stats){
        for(auto& sub : grp){
            int start = std::get<1>(sub).get_cds_start();
            int end = std::get<1>(sub).get_cds_end();
            std::pair<int,int> key = std::make_pair(start,end);

            fit = flat.insert(std::make_pair(key,std::vector<std::tuple<SEGTP,TX,TX,Score,std::string>>{}));
            fit.first->second.emplace_back(sub);
        }
    }
    // sort each vector of transcripts by score
    for(auto& kv : flat){
        std::sort(kv.second.begin(),kv.second.end(), [ ]( const std::tuple<SEGTP,TX,TX,Score,std::string>& lhs, const std::tuple<SEGTP,TX,TX,Score,std::string>& rhs )
        {
            return score_gt(std::get<3>(lhs),std::get<3>(rhs));
        });
    }
    return flat;
}

// returns -1,-1 if nothing passing has been found
std::pair<int,int> find_best_stat(std::map<std::pair<int,int>,std::vector<std::tuple<SEGTP,TX,TX,Score,std::string>>>& stats_flat){
    std::pair<int,int> best_se = stats_flat.begin()->first;
    Score best_score = std::get<3>(stats_flat.begin()->second.front());
    for(auto& kv : stats_flat){
        if (score_gt(std::get<3>(kv.second.front()),best_score)){
            best_score = std::get<3>(kv.second.front());
            best_se = kv.first;
        }
    }
    if(!best_score.pass){ // nothing found
        return std::make_pair(-1,-1);
    }
    return best_se;
}

int  run(){
    #ifdef DEBUG
        std::cout<<"Running in DEBUG"<<std::endl;
    #endif
    // load the reference GFF
    std::cerr<<"loading reference genome"<<std::endl;
    Transcriptome transcriptome;
    if(!global_params.reference_fasta_fname.empty()){
        transcriptome.set_ref(global_params.reference_fasta_fname);
    }
    if(global_params.non_aug){
        transcriptome.use_non_aug();
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
        rstats.num_dirty = transcriptome.clean_cds(global_params.rescue,global_params.overhang, global_params.spliced_overhang, global_params.use_id); // TODO: if we load sequence data into bundles - this is not cost effective...
    }

    transcriptome.set_cds_as_exons();
    transcriptome.remove_non_coding();

    // sort again since now cds is set as exons
    transcriptome.sort();

    std::cerr<<"start removing duplicates"<<std::endl;
    rstats.template_duplicates = transcriptome.deduplicate(global_params.use_id);

    std::cerr<<"loading query transcriptome"<<std::endl;

    transcriptome.add(global_params.query_fname,false,false);

    std::cerr<<"bundling transcriptome"<<std::endl;
    transcriptome.bundleup(global_params.use_id,global_params.overhang,global_params.spliced_overhang);

    std::cerr<<"starting main evaluation"<<std::endl;

#ifndef DEBUG
    #ifdef OPENMP_AVAILABLE
        #pragma omp parallel for schedule(dynamic)
    #endif
#endif
    for(int bi=0;bi<transcriptome.bsize();bi++){ // iterate over bundles
        std::vector<Bundle>::iterator bundle_it = transcriptome.bbegin();
        bundle_it+=bi;
        TX *q,*t;
        TX q_copy;

        // REFERENCELESS BUNDLE
        if(!bundle_it->has_template() || !bundle_it->has_query()){
            for(int qi=0;qi<bundle_it->size();qi++) {
                q = bundle_it->operator[](qi);
#ifdef DEBUG
                if(std::strcmp(q->get_tid().c_str(),"CHS.11626.7")==0){ // rna-XM_011520617.2
                    std::cout<<"found"<<std::endl;
                }
#endif
                if(q->is_template()){continue;}
                std::string cur_seqid;
                transcriptome.seqid2name(q->get_seqid(),cur_seqid);
#ifndef DEBUG
    #ifdef OPENMP_AVAILABLE
        #pragma omp critical
    #endif
#endif
                {
                    if((global_params.keep_all_cds || global_params.keep_cds_if_not_found) && q->has_cds()) {
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
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << std::endl;
                }
            }
            continue;
        }

        // MAIN LOOP

        std::vector<std::vector<std::tuple<SEGTP,TX,TX,Score,std::string>>> stats; // segment,query,template,score,notes
        CHAIN segments;

        for(int qi=0;qi<bundle_it->size();qi++){
            stats.clear();
            q=bundle_it->operator[](qi);
            TX og_q = *q; // original query with the original CDS
#ifdef DEBUG
            if(std::strcmp(q->get_tid().c_str(),"XLOC_000145-mRNA-2")==0){
                std::cout<<"found"<<std::endl;
            }
#endif
            if(q->is_template()){continue;}
            if(global_params.keep_cds_if_not_found && q->has_cds()){ // if we are keeping CDS if not found - we need to build it here and store in the original query transcript
                og_q.build_cds();
            }
            if(global_params.keep_all_cds && q->has_cds()){
                q->build_cds();
                std::string cur_seqid;
                transcriptome.seqid2name(q->get_seqid(),cur_seqid);
#ifndef DEBUG
    #ifdef OPENMP_AVAILABLE
        #pragma omp critical
    #endif
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
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << std::endl;
                }
                continue;
            }
            q->remove_cds();
            q->set_cds_source("ORFanage");

            for(int ti=0;ti<bundle_it->size();ti++){
                t=bundle_it->operator[](ti);
                if(!t->is_template()){continue;}
                int intlen = q->exon_chain()->intersection(*t->cds_chain(),segments,true);
                stats.push_back(std::vector<std::tuple<SEGTP,TX,TX,Score,std::string>>{});
                if(intlen==0){
                    stats.back().push_back(std::make_tuple(SEGTP(),*q,*t,Score(),"no-overlap")); // empty score
                    std::get<3>(stats.back().back()).pass = false;
                    continue;
                }

                std::unordered_set<int> orf_starts;
                std::pair<std::unordered_set<int>::iterator,bool> os_rit;

                // now we need to reconstruct a transcript within the query exon chain for each segment of the template ORF
                for(auto& s : segments){
                    stats.back().push_back(std::make_tuple(s,*q,*t,Score(),"-"));
                    // TODO: does it work without sequence avaialble?

#ifdef DEBUG
                    if(std::strcmp(q->get_tid().c_str(),"CHS.34256.1")==0){ // rna-XM_011520617.2
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

                        int ret = qseg.rescue_cds(global_params.non_aug,global_params.overhang,global_params.spliced_overhang,t);
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
                    if(qseg_score.lpi<global_params.len_perc_ident ||
                       qseg_score.ilpi<global_params.len_frame_perc_ident ||
                       qseg_score.mlpi<global_params.len_match_perc_ident){
                        std::get<3>(stats.back().back()).pass = false;
                        continue;
                    }

                    os_rit = orf_starts.insert(qseg.get_cds_start());
                    if(!os_rit.second){ // the ORF previously existed
                        std::get<4>(stats.back().back())="dup";
                        std::get<3>(stats.back().back()).pass = false;
                        continue;
                    }

                    // otherwise - we can safely add to the current evaluation
                    qseg.store_template(t);
                    qseg.add_attribute("orfanage_status","1");
                    qseg.add_attribute("orfanage_template_source",t->get_source());
                    qseg.add_attribute("orfanage_template",t->get_tid());
                    std::get<1>(stats.back().back())=qseg;
                }
                segments.clear();
            }

#ifdef DEBUG
            if(std::strcmp(q->get_tid().c_str(),"CHS.34256.1")==0){ // rna-XM_011520617.2
                std::cout<<"found"<<std::endl;
            }
#endif
            // flatten out the results
            std::map<std::pair<int,int>,std::vector<std::tuple<SEGTP,TX,TX,Score,std::string>>> stats_flat = flatten(stats); // key: cds start/end; value: segment, query, templates, score, notes

            // find best CDS
            std::pair<int,int> best_se = find_best_stat(stats_flat);

            std::string cur_seqid;
            transcriptome.seqid2name(q->get_seqid(),cur_seqid);

#ifndef DEBUG
    #ifdef OPENMP_AVAILABLE
        #pragma omp critical
    #endif
#endif
            if(stats_flat.empty() || best_se.first<0){
                // if requested keep_cds_if_not_found - handle here
                std::string keep_cds_result = "-";
                if (global_params.keep_cds_if_not_found && og_q.has_cds()){
                    q = &og_q;
                    q->build_cds();
                    keep_cds_result = "keep_cds_if_not_found";
                }
                q->add_attribute("orfanage_status","0");
                global_params.out_gtf_fp<<q->str(cur_seqid)<<std::endl;
                global_params.stats_fp<<q->get_tid()<<"\t"
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
                                           << "-" << "\t"
                                           << "-" << "\t"
                                           << "-" << std::endl;
            }
            else{
                if(global_params.mode_array[0]==ALL){ // output all
#ifdef DEBUG
                    auto sfit = stats_flat.find(best_se);
                    assert(sfit!=stats_flat.end());
#endif

                    // next - iterate over each, merge duplicates and write the front
                    for(auto& kv: stats_flat){
                        auto& first_kv = std::get<1>(kv.second.front());
                        if(!std::get<3>(kv.second.front()).pass){
                            continue;
                        }
                        if(kv.second.size()>1){
                            for(auto v=kv.second.begin();v!=kv.second.end();v++){
                                first_kv.merge(std::get<2>(*v)); // merge into best
                            }
                        }

                        // assign new transcript name
                        std::string og_tid = first_kv.get_tid();
                        std::string new_tid = og_tid+"."+std::to_string(kv.first.first)+"-"+std::to_string(kv.first.second); // new name shall consist of the tid+start/end cds coordinates
                        first_kv.set_tid(new_tid);

                        first_kv.add_attribute("orfanage_status","1");
                        first_kv.add_attribute("orfanage_template_source",first_kv.get_template()->get_source());
                        first_kv.add_attribute("orfanage_template",first_kv.get_template()->get_tid());
                        first_kv.add_attribute("orfanage_duplicity",std::to_string(first_kv.num_dups()));
                        std::get<4>(kv.second.front()) = "gtf";
                        if(kv.first==best_se){
                            // even though we write ALL - mark the best one
                            std::get<4>(kv.second.front()) = "best_gtf";
                        }
                        global_params.out_gtf_fp<<first_kv.str(cur_seqid)<<std::endl;

                        // reset_tid
                        first_kv.set_tid(og_tid);
                    }
                }
                else{ // just output the best
                    // merge all transcripts into the best representation
                    auto sfit = stats_flat.find(best_se);
#ifdef DEBUG
                    assert(sfit!=stats_flat.end());
#endif
                    auto& frontq = std::get<1>(sfit->second.front());

                    if(sfit->second.size()>1){
                        for(auto v=sfit->second.begin()+1;v!=sfit->second.end();v++){
                            frontq.merge(std::get<2>(*v)); // merge into best
                        }
                    }

                    frontq.add_attribute("orfanage_status","1");
                    frontq.add_attribute("orfanage_template_source",frontq.get_template()->get_source());
                    frontq.add_attribute("orfanage_template",frontq.get_template()->get_tid());
                    frontq.add_attribute("orfanage_duplicity",std::to_string(frontq.num_dups()));
                    std::get<4>(sfit->second.front()) = "best_gtf";
                    global_params.out_gtf_fp<<frontq.str(cur_seqid)<<std::endl;
                }
                // write stats
                for(auto& kv : stats_flat){
                    // assign new transcript name
                    auto& first_kv = std::get<1>(kv.second.front());
                    std::string tid = first_kv.get_tid();
                    if(global_params.mode_array[0]==ALL){
                        tid = tid+"."+std::to_string(kv.first.first)+"-"+std::to_string(kv.first.second); // new name shall consist of the tid+start/end cds coordinates
                    }

                    for(auto& v : kv.second){ // for every segment
                        global_params.stats_fp<<tid<<"\t" // query
                                              <<std::get<2>(v).get_tid()<<"\t" // template
                                              <<std::get<0>(v)<<"\t" // segment
                                              <<std::get<4>(v)<<"\t" // notes
                                              <<std::get<1>(v).num_dups()+1<<"\t" // number of duplciates represented by the current ORF
                                              <<std::get<3>(v)<<std::endl; // score
                    }
                }
            }
        };
    }

    std::cerr<<"Done"<<std::endl;

    return 0;
}

// TODO:
//   - ILPI currently divides by the length of the template. Incorrect. Suppose query is 10bp and reference is 5.
//     and suppose all of query matches the reference - this will return 100% ilpi which is incorrect.
//     Instead should be divided by matches+mismatches = 10

int main(int argc, char** argv) {
    ArgParse args("orfanage",
                  "Annotating Open Reading Frames based on reference, phylogeny and sequence similarity.");

    args.add_option("query",ArgParse::Type::STRING, "Path to a GTF query file with transcripts to which CDSs are to be ported", ArgParse::Level::GENERAL,true);
    args.add_option("output", ArgParse::Type::STRING, "Basename for all output files generated by this software", ArgParse::Level::GENERAL,true);
    args.add_option("reference",ArgParse::Type::STRING,"Path to the reference genome file in FASTA format. This parameter is required when the following parameters are used: 1. cleanq; 2. cleant; 3. pi.",ArgParse::Level::GENERAL,false);
    args.add_option("cleanq",ArgParse::Type::FLAG,"If enabled - will ensure all transcripts in the output file will have a valid start and end codons. This option requires the use of --reference parameter",ArgParse::Level::GENERAL,false);
    args.add_option("cleant",ArgParse::Type::FLAG,"If enabled - will ensure all ORFs in the reference annotations start with a valid start codon and end with the first available stop codon. This option requires the use of --reference parameter",ArgParse::Level::GENERAL,false);
    args.add_option("rescue",ArgParse::Type::FLAG,"If enabled - will attempt rescuing the broken ORFs in the reference annotations. This option requires the use of --reference parameter",ArgParse::Level::GENERAL,false);
    args.add_option("lpi",ArgParse::Type::INT,"Percent identity by length between the original and reference transcripts. If -1 (default) is set - the check will not be performed.",ArgParse::Level::GENERAL,false);
    args.add_option("ilpi",ArgParse::Type::INT,"Percent identity by length of bases in frame of the reference transcript. If -1 (default) is set - the check will not be performed.",ArgParse::Level::GENERAL,false);
    args.add_option("mlpi",ArgParse::Type::INT,"Percent identity by length of bases that are in both query and reference. If -1 (default) is set - the check will not be performed.",ArgParse::Level::GENERAL,false);
    args.add_option("minlen",ArgParse::Type::INT,"Minimum length of an open reading frame to consider for the analysis",ArgParse::Level::GENERAL,false);
    args.add_option("mode", ArgParse::Type::STRING, "Which CDS to report: ALL, LONGEST, LONGEST_MATCH, FIRST, BEST, START_MATCH. Default: " + mode_to_str(global_params.mode_array.front()), ArgParse::Level::GENERAL, false);
    args.add_option("stats",ArgParse::Type::STRING,"Output a separate file with stats for each query/template pair",ArgParse::Level::GENERAL,false);
    args.add_option("threads",ArgParse::Type::INT,"Number of threads to run in parallel",ArgParse::Level::GENERAL,false);
    args.add_option("use_id",ArgParse::Type::FLAG,"If enabled, only transcripts with the same gene ID from the query file will be used to form a bundle. In this mode the same template transcript may be used in several bundles, if overlaps transcripts with different gene_ids.",ArgParse::Level::GENERAL,false);
    args.add_option("non_aug",ArgParse::Type::FLAG,"If enabled, non-AUG start codons in reference transcripts will not be discarded and will be considered in overlapping query transcripts on equal grounds with the AUG start codon.",ArgParse::Level::GENERAL,false);

    args.add_option("keep_all_cds",ArgParse::Type::FLAG,"Mutually exclusive with '--keep_cds_if_not_found'. If enabled, any CDS already present in the query will be kept unmodified.",ArgParse::Level::GENERAL,false);
    args.add_option("keep_cds_if_not_found",ArgParse::Type::FLAG,"Mutually exclusive with '--keep_all_cds'. If enabled, will still search for new ORF in each query transcript. If query transcript has CDS annotated, and no ORF can be identified by the method, the original will be kept. Original CDS will be replaced if a valid ORF can be found. Use '--keep_all_cds' to retain all unmodified CDS in query.",ArgParse::Level::GENERAL,false);

    args.add_option("overhang",ArgParse::Type::INT,"If enabled, will also evaluate nucleotide sequence up and downstream up to N bases as set for the argument.",ArgParse::Level::GENERAL,false);
    args.add_option("spliced_overhang",ArgParse::Type::FLAG,"Only in effect when combined with the '--overhang' parameter. If enabled, this option will extend the sequence up to the '--overhang' number of bases, but terminate prematurely if either a splice donor (if extending towards 3') or splice acceptor (if extending towards 5') is detected.",ArgParse::GENERAL,false);

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

    if (args.is_set("mode")){
        global_params.mode_array.clear();

        std::string mode = args.get_string("mode");
        std::transform(mode.begin(), mode.end(), mode.begin(), toupper);

        if (mode == "LONGEST") {
            global_params.mode_array.push_back(LONGEST);
            global_params.mode_array.push_back(LONGEST_MATCH);
            global_params.mode_array.push_back(BEST);
            global_params.mode_array.push_back(START_MATCH);
        }
        else if (mode == "START_MATCH") {
            global_params.mode_array.push_back(START_MATCH);
            global_params.mode_array.push_back(LONGEST);
            global_params.mode_array.push_back(LONGEST_MATCH);
            global_params.mode_array.push_back(BEST);
        }
        else if (mode == "ALL"){
            global_params.mode_array.push_back(ALL);
            global_params.mode_array.push_back(LONGEST_MATCH);
            global_params.mode_array.push_back(BEST);
            global_params.mode_array.push_back(LONGEST);
            global_params.mode_array.push_back(START_MATCH);
        }
        else if (mode == "LONGEST_MATCH") {
            global_params.mode_array.push_back(LONGEST_MATCH);
            global_params.mode_array.push_back(BEST);
            global_params.mode_array.push_back(START_MATCH);
            global_params.mode_array.push_back(LONGEST);
        }
        else if (mode == "BEST") {
            global_params.mode_array.push_back(BEST);
            global_params.mode_array.push_back(LONGEST_MATCH);
            global_params.mode_array.push_back(START_MATCH);
            global_params.mode_array.push_back(LONGEST);
        }
        else if (mode == "FIRST") {
            global_params.mode_array.push_back(FIRST);
            global_params.mode_array.push_back(BEST);
            global_params.mode_array.push_back(LONGEST_MATCH);
            global_params.mode_array.push_back(START_MATCH);
            global_params.mode_array.push_back(LONGEST);
        }
        else{
            printf(OUT_ERROR "Please choose a valid mode (START_MATCH, ALL, LONGEST, LONGEST_MATCH, FIRST, BEST)!\n" OUT_RESET);
            return -1;
        }
    }

    global_params.clean_query = args.is_set("cleanq");
    global_params.clean_templ = args.is_set("cleant");
    global_params.rescue = args.is_set("rescue");
    if(global_params.clean_templ || global_params.clean_query || global_params.rescue){
        assert(args.is_set("reference"));
    }

    global_params.cds_minlen = args.is_set("minlen") ? args.get_int("minlen") : def_params.cds_minlen;
    global_params.len_perc_ident = args.is_set("lpi") ? args.get_int("lpi") : def_params.len_perc_ident;
    global_params.len_frame_perc_ident = args.is_set("ilpi") ? args.get_int("ilpi") : def_params.len_frame_perc_ident;
    global_params.len_match_perc_ident = args.is_set("mlpi") ? args.get_int("mlpi") : def_params.len_match_perc_ident;

    global_params.overhang = args.is_set("overhang") ? args.get_int("overhang") : def_params.overhang;
    global_params.spliced_overhang = args.is_set("spliced_overhang") ? true : false;

    global_params.num_threads = args.is_set("threads") ? args.get_int("threads") : def_params.num_threads;
#ifndef DEBUG
    #ifdef OPENMP_AVAILABLE
        omp_set_num_threads(global_params.num_threads);
    #endif
#endif

    if(args.is_set("use_id")){
        global_params.use_id = true;
    }

    if(args.is_set("non_aug")){
        global_params.non_aug = true;
    }

    if(args.is_set("keep_all_cds")){
        global_params.keep_all_cds = true;
    }
    if(args.is_set("keep_cds_if_not_found")){
        global_params.keep_cds_if_not_found = true;
    }

    if (global_params.keep_all_cds && global_params.keep_cds_if_not_found){
        std::cerr<<"'--keep_all_cds' and '--keep_cds_if_not_found' are mutually exclusive flags. Please choose one or the other"<<std::endl;
        exit(1);
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
                                "num_templates\t"
                                <<s.stats_header()<<std::endl;
    }

    run();

    global_params.stats_fp.close();
    global_params.out_gtf_fp.close();

    return 0;
}
#pragma clang diagnostic pop
