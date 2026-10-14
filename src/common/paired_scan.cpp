// Exhaustive ungapped paired-arm search with at most K substitutions per arm.
//
// Input (TSV, one pair per line): key, left_target, right_target, gap.
// Targets are the sense-strand transcript segments that each arm must match
// (the reverse complement of the synthesized targeting arm); the gap is the
// number of transcript bases between the two arm footprints.
//
// Pigeonhole guarantee: the longer arm is split into K+1 disjoint seeds, so
// any alignment with <= K substitutions in that arm contains one exact seed.
// Every seed occurrence is verified across both full arms at the designed
// geometry, so the output is complete for sense-strand, ungapped matches that
// contain no ambiguous reference bases. Insertions/deletions are not searched.
//
// Output TSV columns: key subject l_start l_end r_start r_end mm_l mm_r
// mm_pos_l_5to3 mm_pos_r_5to3 matched_l matched_r (coordinates 1-based inclusive).
// The last line is "#DONE" with run statistics; its absence means truncation.
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <sys/resource.h>
using namespace std;

struct Query { string key, L, R; int gap; };
struct Seed { int q, offset; };
static const int MAX_SEED = 10;

static int base(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; }
}
// -1: ambiguous base inside the window; otherwise mismatch count capped at cap+1.
static int mismatches(const string& text, size_t start, const string& query, int cap) {
    int n = 0;
    for (size_t i = 0; i < query.size(); ++i) {
        char c = text[start + i];
        if (base(c) < 0) return -1;
        if (c != query[i] && ++n > cap) return n;
    }
    return n;
}
static string positions(const string& text, size_t start, const string& query) {
    vector<int> p;
    for (size_t i = 0; i < query.size(); ++i)
        if (text[start + i] != query[i]) p.push_back((int)(query.size() - i));
    sort(p.begin(), p.end());
    string out;
    for (int x : p) { if (!out.empty()) out += ","; out += to_string(x); }
    return out;
}
static int fail(int code, const string& message) { cerr << "paired_scan: " << message << "\n"; return code; }

int main(int argc, char** argv) {
    if (argc != 5) return fail(2, "usage: paired_scan queries.tsv reference.fasta output.tsv max_substitutions_per_arm");
    int K;
    try { K = stoi(argv[4]); } catch (...) { return fail(2, "max_substitutions_per_arm must be an integer"); }
    if (K < 0 || K > 10) return fail(2, "max_substitutions_per_arm must be between 0 and 10");
    vector<Query> qs;
    {
        ifstream qf(argv[1]);
        if (!qf) return fail(3, "cannot read query file");
        string line;
        size_t number = 0;
        while (getline(qf, line)) {
            ++number;
            if (line.empty()) continue;
            istringstream s(line);
            Query q;
            if (!(s >> q.key >> q.L >> q.R >> q.gap) || q.gap < 0)
                return fail(4, "malformed query on line " + to_string(number));
            for (const string* arm : {&q.L, &q.R})
                for (char c : *arm) if (base(c) < 0) return fail(4, "non-ACGT arm on line " + to_string(number));
            if ((int)min(q.L.size(), q.R.size()) < K + 1)
                return fail(4, "arm shorter than max_substitutions+1 on line " + to_string(number));
            qs.push_back(q);
        }
    }
    // Seed index per seed length: code -> seeds (query, offset of seed from the pair start).
    vector<vector<vector<Seed>>> idx(MAX_SEED + 1);
    vector<int> ks;
    for (int j = 0; j < (int)qs.size(); ++j) {
        const Query& q = qs[j];
        bool left = q.L.size() >= q.R.size();
        const string& a = left ? q.L : q.R;
        int off = left ? 0 : (int)q.L.size() + q.gap;
        for (int i = 0; i < K + 1; ++i) {
            int start = (int)(i * a.size() / (K + 1)), end = (int)((i + 1) * a.size() / (K + 1));
            int k = min(end - start, MAX_SEED);   // a prefix of an exact segment is still exact
            if (k < 1) return fail(5, "internal: empty seed");
            if (idx[k].empty()) { idx[k].resize((size_t)1 << (2 * k)); ks.push_back(k); }
            unsigned code = 0;
            for (int p = start; p < start + k; ++p) code = (code << 2) | base(a[p]);
            idx[k][code].push_back({j, off + start});
        }
    }
    ofstream out(argv[3]);
    if (!out) return fail(7, "cannot write output file");
    out << "key\tsubject\tl_start\tl_end\tr_start\tr_end\tmm_l\tmm_r\tmm_pos_l_5to3\tmm_pos_r_5to3\tmatched_l\tmatched_r\n";
    long long nseq = 0, nbases = 0, nhits = 0, nambiguous = 0, checks = 0;
    vector<unordered_set<long long>> emitted(qs.size());
    auto scan = [&](const string& id, const string& seq) {
        if (id.empty()) return;
        ++nseq; nbases += (long long)seq.size();
        for (char c : seq) if (base(c) < 0) ++nambiguous;
        for (auto& e : emitted) e.clear();
        for (int k : ks) {
            unsigned code = 0, mask = (1u << (2 * k)) - 1; int run = 0;
            for (size_t i = 0; i < seq.size(); ++i) {
                int b = base(seq[i]);
                if (b < 0) { run = 0; code = 0; continue; }
                code = ((code << 2) | (unsigned)b) & mask; ++run;
                if (run < k) continue;
                for (const Seed& s : idx[k][code]) {
                    const Query& q = qs[s.q];
                    long long st = (long long)i - k + 1 - s.offset;
                    long long rs = st + (long long)q.L.size() + q.gap;
                    if (st < 0 || rs + (long long)q.R.size() > (long long)seq.size()) continue;
                    ++checks;
                    bool anchor_left = q.L.size() >= q.R.size();
                    int ml, mr;
                    if (anchor_left) { ml = mismatches(seq, (size_t)st, q.L, K); if (ml < 0 || ml > K) continue; mr = mismatches(seq, (size_t)rs, q.R, K); }
                    else { mr = mismatches(seq, (size_t)rs, q.R, K); if (mr < 0 || mr > K) continue; ml = mismatches(seq, (size_t)st, q.L, K); }
                    if (ml < 0 || mr < 0 || ml > K || mr > K || !emitted[s.q].insert(st).second) continue;
                    ++nhits;
                    out << q.key << '\t' << id << '\t' << st + 1 << '\t' << st + (long long)q.L.size() << '\t'
                        << rs + 1 << '\t' << rs + (long long)q.R.size() << '\t' << ml << '\t' << mr << '\t'
                        << positions(seq, (size_t)st, q.L) << '\t' << positions(seq, (size_t)rs, q.R) << '\t'
                        << seq.substr((size_t)st, q.L.size()) << '\t' << seq.substr((size_t)rs, q.R.size()) << '\n';
                }
            }
        }
    };
    ifstream ref(argv[2]);
    if (!ref) return fail(8, "cannot read reference FASTA");
    string line, id, seq;
    while (getline(ref, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            scan(id, seq);
            istringstream s(line.substr(1)); s >> id; seq.clear();
            if (id.empty()) return fail(9, "FASTA header without an identifier");
        } else {
            if (id.empty()) return fail(9, "FASTA sequence precedes its header");
            for (char c : line) if (!isspace((unsigned char)c)) seq += (char)toupper((unsigned char)c);
        }
    }
    scan(id, seq);
    if (nseq == 0) return fail(9, "reference FASTA is empty");
    struct rusage usage; getrusage(RUSAGE_SELF, &usage);
#ifdef __APPLE__
    long long peak_kb = usage.ru_maxrss / 1024;   // bytes on macOS
#else
    long long peak_kb = usage.ru_maxrss;          // kilobytes on Linux
#endif
    out << "#DONE\tqueries=" << qs.size() << "\tsequences=" << nseq << "\tbases=" << nbases
        << "\tambiguous_bases=" << nambiguous << "\tverified_candidates=" << checks << "\tpaired_hits=" << nhits
        << "\tpeak_rss_kb=" << peak_kb << "\n";
    out.flush();
    if (!out) return fail(7, "output write failed");
    return 0;
}
