#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <htslib/sam.h>
#include <htslib/hts.h>

// 函数声明
int get_snp_count(bam1_t *aln);
int get_indel_count(bam1_t *aln);
bool is_proper_pair(bam1_t *aln1, bam1_t *aln2, int max_insert_size);
std::pair<bam1_t*, bam1_t*> choose_best_pair(const std::vector<std::pair<bam1_t*, bam1_t*>>& candidate_pairs);
void process_read_group(std::vector<bam1_t*>& current_group, int mapQ_threshold, int max_insert_size, int max_snp, int max_indel,
                        samFile *out_pass, samFile *out_fail, bam_hdr_t *header, int &processed_pairs);
void clean_up_resources(bam1_t *aln, samFile *in, samFile *out_pass, samFile *out_fail, bam_hdr_t *header);

// 计算 SNP 数量
int get_snp_count(bam1_t *aln) {
    uint8_t *nm = bam_aux_get(aln, "NM");
    if (nm != nullptr) {
        return bam_aux2i(nm);
    }
    return 0;
}

// 计算 Indel 数量
int get_indel_count(bam1_t *aln) {
    uint32_t *cigar = bam_get_cigar(aln);
    int indel_count = 0;
    for (int i = 0; i < aln->core.n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CINS || op == BAM_CDEL) {
            ++indel_count;
        }
    }
    return indel_count;
}

// 判断是否是合适的配对
bool is_proper_pair(bam1_t *aln1, bam1_t *aln2, int max_insert_size) {
    int32_t chr1 = aln1->core.tid;
    int32_t chr2 = aln2->core.tid;

    if (chr1 != chr2) {
        return false;
    }

    int32_t pos1 = aln1->core.pos;
    int32_t pos2 = aln2->core.pos;
    int32_t insert_size1 = aln1->core.isize;
    int32_t insert_size2 = aln2->core.isize;

    if (insert_size1 == -insert_size2) {
        int32_t coord_diff = abs(pos2 - pos1);
        if (coord_diff <= max_insert_size && abs(insert_size1) <= max_insert_size) {
            bool strand1 = aln1->core.flag & BAM_FREVERSE;
            bool strand2 = aln2->core.flag & BAM_FREVERSE;
            if (strand1 != strand2) {
                std::cout << "Proper pair detected: " << bam_get_qname(aln1) << " and " << bam_get_qname(aln2) << "\n";
                return true;
            }
        }
    }

    return false;
}

// 选择得分最高的配对
std::pair<bam1_t*, bam1_t*> choose_best_pair(const std::vector<std::pair<bam1_t*, bam1_t*>>& candidate_pairs) {
    if (candidate_pairs.empty()) {
        return {nullptr, nullptr};  // 没有候选配对
    }

    std::pair<bam1_t*, bam1_t*> best_pair = candidate_pairs[0];
    int best_score = -1;

    for (const auto& pair : candidate_pairs) {
        bam1_t* aln1 = pair.first;
        bam1_t* aln2 = pair.second;

        // 计算配对的得分
        int mapq1 = aln1->core.qual;
        int mapq2 = aln2->core.qual;
        int snp1 = get_snp_count(aln1);
        int snp2 = get_snp_count(aln2);
        int indel1 = get_indel_count(aln1);
        int indel2 = get_indel_count(aln2);

        // 打分规则: MAPQ 越高得分越高，SNP 和 Indel 越少得分越高
        int score = mapq1 + mapq2 - (snp1 + snp2 + indel1 + indel2);

        // 找到得分最高的配对
        if (score > best_score) {
            best_score = score;
            best_pair = pair;
        }
    }

    return best_pair;  // 返回得分最高的配对
}


// 处理一组 reads 并选择得分最高的位置
void process_read_group(std::vector<bam1_t*>& current_group, int mapQ_threshold, int max_insert_size, int max_snp, int max_indel,
                        samFile *out_pass, samFile *out_fail, bam_hdr_t *header, int &processed_pairs) {
    std::unordered_set<int> processed_set;  // 使用索引存储已处理的 reads
    std::vector<std::pair<bam1_t*, bam1_t*>> candidate_pairs;

    if (current_group.size() < 2) {
        std::cout << "----- Current group size < 2, no pairs to process.\n";
        
        // 对于 size < 2 的情况，直接输出到 fail 文件
        for (bam1_t* aln : current_group) {
            int snp = get_snp_count(aln);
            int indel = get_indel_count(aln);

            // 记录日志
            std::cout << "Fail: Unpaired read: " << bam_get_qname(aln) << "\n";
            std::cout << "Chromosome position: " << aln->core.tid << ":" << aln->core.pos 
                      << ", SNPs: " << snp << ", Indels: " << indel << "\n";

            // 输出到 fail 文件
            if (sam_write1(out_fail, header, aln) < 0) {
                std::cerr << "Error: could not write unpaired read to fail BAM\n";
                exit(1);
            }
        }

        // 清理当前组
        for (auto aln_ptr : current_group) {
            bam_destroy1(aln_ptr);
        }
        current_group.clear();
        return;
    }

    std::cout << "----- Processing group of size: " << current_group.size() << "\n";

    // 收集符合阈值要求的配对
    for (size_t i = 0; i < current_group.size() - 1; ++i) {
        if (processed_set.find(i) != processed_set.end()) {
            continue;  // 如果该 read 已经被处理，则跳过
        }
        bam1_t *aln1 = current_group[i];
        for (size_t j = i + 1; j < current_group.size(); ++j) {
            if (processed_set.find(j) != processed_set.end()) {
                continue;  // 如果该 read 已经被处理，则跳过
            }
            bam1_t *aln2 = current_group[j];
            if (is_proper_pair(aln1, aln2, max_insert_size)) {
                int snp1 = get_snp_count(aln1);
                int snp2 = get_snp_count(aln2);
                int indel1 = get_indel_count(aln1);
                int indel2 = get_indel_count(aln2);

                // 记录配对的详细信息
                std::cout << "Checking pair: " << bam_get_qname(aln1) << " and " << bam_get_qname(aln2) << "\n";
                std::cout << "Chromosome positions: " << aln1->core.tid << ":" << aln1->core.pos 
                          << " and " << aln2->core.tid << ":" << aln2->core.pos << "\n";
                std::cout << "SNPs: " << snp1 << " + " << snp2 << ", Indels: " << indel1 << " + " << indel2 << "\n";

                // 过滤条件
                if (snp1 + snp2 <= max_snp && indel1 + indel2 <= max_indel) {
                    candidate_pairs.push_back({aln1, aln2});
                } else {
                    std::cout << "Fail: SNP or indel count exceeds threshold.\n";
                }
            }
        }
    }

    // 选择得分最高的配对
    if (!candidate_pairs.empty()) {
        auto best_pair = choose_best_pair(candidate_pairs);
        bam1_t *aln1 = best_pair.first;
        bam1_t *aln2 = best_pair.second;

        // 输出到 pass 文件
        if (sam_write1(out_pass, header, aln1) < 0 || sam_write1(out_pass, header, aln2) < 0) {
            std::cerr << "Error: could not write alignment to pass BAM\n";
            exit(1);
        }

        // 获取 aln1 和 aln2 在 current_group 中的下标并标记为已处理
        for (size_t i = 0; i < current_group.size(); ++i) {
            if (current_group[i] == aln1) {
                processed_set.insert(i);
            }
            if (current_group[i] == aln2) {
                processed_set.insert(i);
            }
        }

        ++processed_pairs;

        std::cout << "Pass: Pair written to pass BAM. Best scoring proper pair.\n";
    }

    // 未被选择的配对输出到 fail 文件，并详细记录 SNP、Indel 数量及染色体位置信息
    for (size_t i = 0; i < current_group.size(); ++i) {
        if (processed_set.find(i) == processed_set.end()) {
            bam1_t *aln = current_group[i];
            int snp = get_snp_count(aln);
            int indel = get_indel_count(aln);

            // 输出 fail 信息到日志
            std::cout << "Fail: Not paired with any read in the group: " << bam_get_qname(aln) << "\n";
            std::cout << "Chromosome position: " << aln->core.tid << ":" << aln->core.pos 
                      << ", SNPs: " << snp << ", Indels: " << indel << "\n";

            // 输出到 fail 文件
            if (sam_write1(out_fail, header, aln) < 0) {
                std::cerr << "Error: could not write alignment to fail BAM\n";
                exit(1);
            }
        }
    }

    // 清理已处理的 reads
    for (auto aln_ptr : current_group) {
        bam_destroy1(aln_ptr);
    }
    current_group.clear();
}

// 资源清理
void clean_up_resources(bam1_t *aln, samFile *in, samFile *out_pass, samFile *out_fail, bam_hdr_t *header) {
    if (aln != nullptr) bam_destroy1(aln);
    if (in != nullptr) sam_close(in);
    if (out_pass != nullptr) sam_close(out_pass);
    if (out_fail != nullptr) sam_close(out_fail);
    if (header != nullptr) bam_hdr_destroy(header);
}

int main(int argc, char *argv[]) {
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] << " <input.bam> <output_pass.bam> <output_fail.bam> <mapQ_threshold> <max_insert_size> <max_snp> <max_indel>\n";
        return 1;
    }

    const char *input_bam = argv[1];
    const char *output_pass_bam = argv[2];
    const char *output_fail_bam = argv[3];
    int mapQ_threshold = std::stoi(argv[4]);
    int max_insert_size = std::stoi(argv[5]);
    int max_snp = std::stoi(argv[6]);
    int max_indel = std::stoi(argv[7]);

    samFile *in = sam_open(input_bam, "rb");
    if (in == nullptr) {
        std::cerr << "Error: could not open input BAM file " << input_bam << "\n";
        return 1;
    }

    samFile *out_pass = sam_open(output_pass_bam, "wb");
    if (out_pass == nullptr) {
        std::cerr << "Error: could not open output pass BAM file " << output_pass_bam << "\n";
        sam_close(in);
        return 1;
    }

    samFile *out_fail = sam_open(output_fail_bam, "wb");
    if (out_fail == nullptr) {
        std::cerr << "Error: could not open output fail BAM file " << output_fail_bam << "\n";
        sam_close(in);
        sam_close(out_pass);
        return 1;
    }

    bam_hdr_t *header = sam_hdr_read(in);
    if (header == nullptr) {
        std::cerr << "Error: could not read BAM header from " << input_bam << "\n";
        clean_up_resources(nullptr, in, out_pass, out_fail, nullptr);
        return 1;
    }

    if (sam_hdr_write(out_pass, header) < 0 || sam_hdr_write(out_fail, header) < 0) {
        std::cerr << "Error: could not write BAM header to output files\n";
        clean_up_resources(nullptr, in, out_pass, out_fail, header);
        return 1;
    }

    bam1_t *aln = bam_init1();
    if (aln == nullptr) {
        std::cerr << "Error: could not initialize BAM structure\n";
        clean_up_resources(nullptr, in, out_pass, out_fail, header);
        return 1;
    }

    std::vector<bam1_t*> current_group;
    int processed_pairs = 0;
    std::string last_qname;

    // 读取 BAM 文件中的每一行
    while (sam_read1(in, header, aln) >= 0) {
        std::string qname(bam_get_qname(aln));
        if (last_qname.empty()) {
            last_qname = qname;
        }

        // 如果当前 read 的名称与上一个相同，继续添加到当前组中
        if (qname == last_qname) {
            current_group.push_back(bam_dup1(aln));
        } else {
            // 否则，处理当前组并清空，开始新组
            process_read_group(current_group, mapQ_threshold, max_insert_size, max_snp, max_indel, out_pass, out_fail, header, processed_pairs);
            last_qname = qname;
            current_group.push_back(bam_dup1(aln));
        }
    }

    // 处理最后一个组
    if (!current_group.empty()) {
        process_read_group(current_group, mapQ_threshold, max_insert_size, max_snp, max_indel, out_pass, out_fail, header, processed_pairs);
    }

    clean_up_resources(aln, in, out_pass, out_fail, header);
    std::cout << "Done! Processed " << processed_pairs << " read pairs.\n";
    return 0;
}

