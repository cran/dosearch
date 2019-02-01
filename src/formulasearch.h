#ifndef FORMULASEARCH_H
#define FORMULASEARCH_H

#include <unordered_map>
#include <vector>
#include <string>
#include <queue>
#include <Rcpp.h>
#include "dcongraph.h"
#include "derivation.h"

using namespace std;

struct p {
    int u, c, j, e;
};

struct rindep {
    int xset, yset, c, j;
};

struct output {
    p to, from, rp;
    string rule_name;
    rindep ri;
    bool valid;
};

struct distr {
    string name, formula, rule_name;
    int rule_num, index, score;
    bool primitive;
    unsigned pa1;
    unsigned pa2;
    p pp;
};

struct comp_distr {
    bool operator()(distr const * d1, distr const * d2) {
        return d1->score < d2->score;
    }
};

class formulasearch {
public:
    formulasearch(const int& n, const bool& dd, const bool& da, const bool& dea, const bool& fa, const bool &imp, const char &ms, const bool &heur, const bool& repl, const bool& verb);
    formulasearch(const formulasearch& orig);
    void add_known(const int& u, const int& c, const int& j, const int& e);
    void set_target(const int& u, const int& c, const int& j, const int& e);
    void set_labels(const Rcpp::StringVector& lab);
    void set_graph(dcongraph* g_);
    void set_derivation(derivation* d_);
    void set_options(const vector<int>& r);
    Rcpp::List search_init();
    virtual ~formulasearch();
private:
    int n, md_s, md_p, md_t, md_u, tr, sb, trsb, index, pool;
    bool md, trivial_id, format_do;
    const bool draw_derivation;
    const bool draw_all;
    const bool derive_all;
    const bool formula;
    const bool improve;
    const char md_sym;
    const bool heuristic;
    const bool replace;
    const bool verbose;
    p target;
    dcongraph *g;
    derivation *deriv;
    vector<distr> target_dist;
    vector<string> labels;
    vector<int> z_sets;
    vector<int> rules;
    // vector<unsigned long> rule_counts;
    // vector<unsigned long> rule_counts_deriv;
    unordered_map<int, distr> L;
    priority_queue<distr*, std::vector<distr*>, comp_distr> Q;
    unordered_map<std::string, int> ps;
    output info;
    bool valid_do_rule(const int& ruleid, const int& u, const int& c, const int& j, const bool& primi) const;
    void apply_do_rule(const int& ruleid, const int& u, const int& c, const int& j, const int& e, const int& z);
    void get_ruleinfo(const int& ruleid, const int& y, const int& x, const int& z, const int& w, const int& d, const int& e);
    string to_string(const p& pp) const;
    string dec_to_text(const int& dec, const int& enabled) const;
    void search_bfs();
    bool is_primitive(const bool& pa1_primitive, const bool& pa2_primitive, const int& ruleid);
    bool required_exists(distr &required);
    bool equal_p(const p& p1, const p& p2) const;
    void draw(const distr& dist, const bool& recursive, derivation& d);
    void derive_formula(distr& dist);
    int compute_score(const p& pp) const;
    int compute_score_md(const p& pp) const;
    int get_previous_vertices(const int &set) const;
};

#endif	/* FORMULASEARCH_H */

