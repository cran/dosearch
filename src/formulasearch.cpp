#include "formulasearch.h"
#include "set.h"
#include "derivation.h"
#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace std;

formulasearch::formulasearch(const int& n_, const bool& dd, const bool& da, const bool& dea, const bool& fa, const bool& imp, const char& ms, const bool& heur, const bool& repl, const bool& verb): 
  n(n_), draw_derivation(dd), draw_all(da), derive_all(dea), formula(fa), improve(imp), md_sym(ms), heuristic(heur), replace(repl), verbose(verb) {
}

formulasearch::~formulasearch() {
}

void formulasearch::set_target(const int& u, const int& c, const int& j, const int& e) {
    target.u = u; target.c = c; target.j = j; target.e = e;
    if ( verbose ) Rcpp::Rcout << "Setting target: " << to_string(target) << endl;
}

void formulasearch::set_options(const vector<int>& r) {

    trivial_id = false;
    format_do = true;
    index = 0;
    pool = 0;
    md_s = g->get_md_switches();
    md_p = g->get_md_proxies();
    md_t = md_s >> 1;
    md_u = 0;
    md = md_s > 0;

    tr = g->get_trnodes();
    sb = g->get_sbnodes();
    trsb = tr | sb;

    if ( r.size() > 0 ) rules = r;
    else {
        // Additional rules for missing data problems
        if ( md ) {
            rules = {4, 5, 8, 9, -1, -2, -3, 1, 2, 3, 6, -6, 7, -7};
        } else {
            if ( improve && (trsb == 0) ) rules = {4, 5, -2, -3, 2, 3, 6, -6};
            else rules = {4, 5, -1, -2, -3, 1, 2, 3, 6, -6};
        }
    }

    // Disable counts for now
    // rule_counts = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    // rule_counts_deriv = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    if ( improve ) {
        z_sets = get_subsets(n);
    } else {
        for (int i = 1; i < full_set(n); i++) {
            z_sets.push_back(i);
        }
    }

}

void formulasearch::add_known(const int& u, const int& c, const int& j, const int& e) {
    index++;
    p pp;
    distr iquery;
    pp.u = u; pp.c = c; pp.j = j; pp.e = e;
    iquery.name = to_string(pp);
    iquery.rule_name = "";
    iquery.rule_num = 0;
    iquery.formula = iquery.name;
    iquery.pp = pp;
    iquery.pa1 = 0;
    iquery.pa2 = 0;
    iquery.primitive = true;
    iquery.index = index;
    iquery.score = compute_score(pp);
    L[index] = iquery;
    ps[iquery.name] = index;
    if ( equal_p(pp, target) ) {
        trivial_id = true;
        target_dist.push_back(L[index]);
    }
    if ( heuristic ) Q.push(&L[index]);
    if ( md ) pool = (pool | u) | ((u & md_p) >> 2);
    else pool = pool | u;
    if ( verbose ) Rcpp::Rcout << "Adding known distribution: " << to_string(pp) << endl;
}

void formulasearch::set_labels(const Rcpp::StringVector& lab) {
    labels = vector<string>(2*n);
    for ( int i = 0; i < n; i++ ) {
        labels[i] = lab(i);
        labels[n+i] = "I(" + lab(i) + ")";
    }
}

void formulasearch::set_graph(dcongraph* g_) {
    g = g_;
}

void formulasearch::set_derivation(derivation* d_) {
    deriv = d_;
}

// Heuristic for search order
int formulasearch::compute_score(const p& pp) const {
    int score = 0;
    int common_y = pp.u & target.u;
    int common_x = pp.j & target.j;
    int pp_w = pp.c - pp.j;
    int target_w = target.c - target.j;
    int common_z = pp_w & target_w;

    score += 10 * set_size(common_y);
    score -= 2 * set_size(target.u - common_y);
    score += 5 * set_size(common_x);
    score -= 2 * set_size(pp.j - common_x);
    score -= 2 * set_size(target.j - common_x);
    score += 3 * set_size(common_z);
    score -= 1 * set_size(pp_w - common_z);
    score -= 1 * set_size(target_w - common_z);

    return(score);
}

// Heuristic for search order that takes proxy variables into account
int formulasearch::compute_score_md(const p& pp) const {
    int score = 0;
    int pp_w = pp.c - pp.j;

    int proxy_u = pp.u & md_p;
    int proxy_w = pp_w & md_p;

    int proxy_total = proxy_u | proxy_w;
    int switch_total = (pp.u | pp_w) & md_s;
    int proxy_needed = switch_total << 1;
    int switch_needed = proxy_total >> 1;

    int proxy_match = proxy_total & proxy_needed;
    int proxy_mismatch = proxy_total - proxy_needed;
    int switch_match = switch_total & switch_needed;
    int switch_mismatch = switch_total - switch_needed;
    int common_y = ((pp.u - proxy_u) | (proxy_u >> 2)) & target.u;
    int common_x = pp.j & target.j;
    int target_w = target.c - target.j;
    int common_z = ((pp_w - proxy_w) | (proxy_w >> 2)) & target_w;

    score += 10 * set_size(common_y);
    score += 6 * set_size(proxy_match);
    score += 6 * set_size(switch_match);
    score -= 2 * set_size(proxy_mismatch);
    score -= 2 * set_size(switch_mismatch);
    score -= 2 * set_size(target.u - common_y);
    score += 6 * set_size(common_x);
    score -= 5 * set_size(pp.j - common_x);
    score -= 2 * set_size(target.j - common_x);
    score += 4 * set_size(common_z);
    score -= 2 * set_size(pp_w - common_z);
    score -= 2 * set_size(target_w - common_z);
    score += 10 * set_size(pp.e);

    return(score);

}

bool formulasearch::required_exists(distr &required) {
    if ( info.rp.u > 0 ) {
        int req = ps[to_string(info.rp)];
        if ( req > 0 ) {
            distr& reqd = L[req];
            required.name = reqd.name;
            required.pp = reqd.pp;
            required.primitive = reqd.primitive;
            required.pa1 = reqd.pa1;
            required.pa2 = reqd.pa2;
            required.rule_num = reqd.rule_num;
            required.rule_name = reqd.rule_name;
            required.index = reqd.index;
            required.score = reqd.score;
            return true;
        }
        return false;
    }
    return true;
}

Rcpp::List formulasearch::search_init() {

    info.to.u = 0; info.to.c = 0; info.to.j = 0; info.to.e = 0;
    info.from.u = 0; info.from.c = 0; info.from.j = 0; info.from.e = 0;
    info.rp.u = 0; info.rp.c = 0; info.rp.j = 0; info.rp.e = 0;
    info.ri.xset = 0; info.ri.yset = 0; info.ri.c = 0; info.ri.j = 0;
    info.valid = false;
    info.rule_name = "";

    unsigned ntarget = 0;
    vector<string> formulas;
    vector<string> derivations;
    string full_derivation = "";
    Rcpp::Timer timer;

    if ( improve ) {
        bool trivial = true;
        if ( (pool & target.u) == target.u ) {
            trivial = FALSE;
        }
        if ( trivial ) {
            return Rcpp::List::create(
                Rcpp::Named("identifiable") = false,
                Rcpp::Named("formula") = formulas,
                Rcpp::Named("derivation") = derivations,
                Rcpp::Named("full_derivation") = full_derivation,
                Rcpp::Named("time") = Rcpp::NumericVector::create(
                    Rcpp::Named("start") = 0.0,
                    Rcpp::Named("end") = 0.0
                )
                // Rcpp::Named("distributions") = L.size(),
                // Rcpp::Named("rule_count1") = rule_counts,
                // Rcpp::Named("rule_count2") = rule_counts_deriv
            );
        }
    }

    if ( !trivial_id || derive_all ) {
        timer.step("start");
        search_bfs();
        timer.step("end");
    } else {
        timer.step("start");
        timer.step("end");
    }

    ntarget = target_dist.size();
    if ( !derive_all && ntarget > 1) ntarget = 1;
    formulas = std::vector<string>(ntarget);
    derivations = std::vector<string>(ntarget);

    for ( unsigned i = 0; i < ntarget; i++ ) {
        if ( formula ) {
            derive_formula(target_dist[i]);
            formulas[i] = target_dist[i].formula;
        }
        if ( draw_derivation && !draw_all ) {
            derivation temp_deriv;
            temp_deriv.init();
            draw(target_dist[i], TRUE, temp_deriv);
            temp_deriv.finish();
            derivations[i] = temp_deriv.get();
        }
    }

    if (draw_derivation && draw_all ) {
        deriv->init();
        for ( auto it = L.begin(); it != L.end(); ++it ) {
            distr& d = it->second;
            draw(d, FALSE, *deriv);
        }
        for ( unsigned i = 0; i < ntarget; i++ ) {
            draw(target_dist[i], FALSE, *deriv);
        }
        deriv->finish();
        full_derivation = deriv->get();
    }

    Rcpp::NumericVector times(timer);
    int nano = 1000000;
    for (int i = 0; i < times.size(); i++) {
        times[i] = times[i] / nano;
    }

    return Rcpp::List::create(
        Rcpp::Named("identifiable") = ntarget > 0,
        Rcpp::Named("formula") = formulas,
        Rcpp::Named("derivation") = derivations,
        Rcpp::Named("full_derivation") = full_derivation,
        Rcpp::Named("time") = times
        // Rcpp::Named("distributions") = L.size(),
        // Rcpp::Named("rule_count1") = rule_counts,
        // Rcpp::Named("rule_count2") = rule_counts_deriv
    );

}

void formulasearch::search_bfs() {

    distr required;
    distr * iptr;
    bool found = false;
    bool primi = true;
    unsigned int i = 0;
    int u, c, j, e, z, ruleid, exist;
    int remaining = L.size();

    while ( remaining > 0 && (!found || derive_all) ) {

        if ( heuristic ) {
            iptr = Q.top();
            Q.pop();
        } else {
            iptr = &L[i+1];
        }
        distr& iquery = *iptr;
        remaining--;

        u = iquery.pp.u;
        c = iquery.pp.c;
        j = iquery.pp.j;
        e = iquery.pp.e;
        primi = iquery.primitive;

        if ( verbose ) Rcpp::Rcout << "Expanding: " << iquery.name << endl;

        for ( unsigned int r = 0; r < rules.size(); r++ ) {

            ruleid = rules[r];

            if ( improve && !valid_do_rule(ruleid, u, c, j, primi) ) continue;

            for ( unsigned int z_ind = 0; z_ind < z_sets.size(); z_ind++ ) {

                required.primitive = TRUE;
                z = z_sets[z_ind];

                apply_do_rule(ruleid, u, c, j, e, z);

                if ( !info.valid ) continue;
                // rule_counts[r]++;

                exist = ps[to_string(info.to)];
                if ( exist > 0 ) {

                    if ( replace ) {

                        if ( L[exist].primitive ) continue;
                        if ( !required_exists(required) ) continue;
                        if ( !is_primitive(primi, required.primitive, ruleid) ) continue;
                        L[exist].primitive = true;
                        L[exist].pa1 = iquery.index;
                        L[exist].pa2 = 0;
                        L[exist].rule_num = ruleid;
                        L[exist].rule_name = info.rule_name;
                        if ( ruleid == 6 || ruleid == -6 || ruleid == 7 || ruleid == -7 ) {
                            L[exist].pa2 = required.index;
                        }

                    } else continue;

                } else {

                    if ( !required_exists(required) ) continue;
                    if ( info.ri.xset > 0 ) {
                        if ( !g->dsep_set(info.ri.xset, info.ri.yset, info.ri.c, info.ri.j) ) continue;
                    }

                    distr nquery;
                    nquery.name = to_string(info.to);
                    nquery.pp = info.to;
                    nquery.primitive = is_primitive(iquery.primitive, required.primitive, ruleid);
                    nquery.pa1 = iquery.index;
                    nquery.pa2 = 0;
                    nquery.rule_num = ruleid;
                    nquery.rule_name = info.rule_name;

                    if ( ruleid == 6 || ruleid == -6 || ruleid == 7 || ruleid == -7 ) {
                        nquery.pa2 = required.index;
                    }

                    if ( equal_p(info.to, target) ) {

                        if ( verbose ) {
                            Rcpp::Rcout << "!!!! Managed to hit the target !!!!" << endl;
                            Rcpp::Rcout << "index = " << index << endl;
                        }
                        found = true;
                        target_dist.push_back(nquery);

                    } else {

                        index++;
                        remaining++;
                        nquery.index = index;
                        // rule_counts_deriv[r]++;
                        if ( heuristic ) {
                            if ( md ) nquery.score = compute_score_md(nquery.pp);
                            else nquery.score = compute_score(nquery.pp);
                            L[index] = nquery;
                            ps[nquery.name] = index;
                            Q.push(&L[index]);
                        } else {
                            L[index] = nquery;
                            ps[nquery.name] = index;
                        }

                    }

                }

                if ( found && !derive_all ) break;

            } // for z

            if ( found && !derive_all ) break;

        } // for ruleid

        i++;

        if ( verbose )  Rcpp::Rcout << "Marking here " << iquery.name << " as expanded (" << i << "/" << index << ")" << endl;

        /*
        if (i > imax) {
            found = false;
            Rcpp::Rcout << "Breaking the infinite loop!" << endl;
            break;
        } */

    } // while

}

bool formulasearch::is_primitive(const bool& pa1_primitive, const bool& pa2_primitive, const int& ruleid) {
    if ( pa1_primitive && pa2_primitive ) {
        if ( ruleid * ruleid < 16 ) return false;
        return true;
    }
    return false;
}

void formulasearch::derive_formula(distr& dist) {
    if ( dist.pa1 > 0 ) {
        distr& pa1 = L[dist.pa1];
        if (pa1.formula == "") derive_formula(pa1);
        if ( dist.pa2 > 0 ) {
            distr& pa2 = L[dist.pa2];
            if ( pa2.formula == "" ) derive_formula(pa2);
            if ( dist.rule_num == 6 || dist.rule_num == -6) {
                if ( dist.primitive ) dist.formula = dist.name;
                else {
                    if ( pa1.formula.length() < pa2.formula.length() ) dist.formula = pa1.formula + "*" + pa2.formula;
                    else dist.formula = pa2.formula + "*" + pa1.formula;
                }
            } else if ( dist.rule_num == 7) {
                if ( dist.primitive ) dist.formula = dist.name;
                else dist.formula = "[[" + pa1.formula + "]/[" + pa2.formula + "]]";
            } else if ( dist.rule_num == -7 ) {
                if ( dist.primitive ) dist.formula = dist.name;
                else dist.formula = "[[" + pa2.formula + "]/[" + pa1.formula + "]]";
            }
        } else {
            if ( dist.rule_num * dist.rule_num < 16 ) {
                dist.formula = pa1.formula;
            } else if ( dist.rule_num == 5 ) {
                if ( dist.primitive ) dist.formula = dist.name;
                else dist.formula = "[[" + pa1.formula + "]/[sum_{" + dec_to_text(dist.pp.u, 0) + "} " + pa1.formula + "]]";
            } else if ( dist.rule_num == 4 ) {
                if ( dist.primitive ) dist.formula = dist.name;
                else dist.formula =  "[sum_{" + dec_to_text(pa1.pp.u - dist.pp.u, 0) + "} [" + pa1.formula + "]]";
            } else if ( dist.rule_num >= 8 ) {
                if ( dist.primitive ) dist.formula = dist.name;
                else dist.formula = pa1.formula;
            }
        }
    } else {
        dist.formula = dist.name;
    }
}

void formulasearch::draw(const distr& dist, const bool& recursive, derivation& d) {
    if ( dist.pa1 > 0 ) {
        distr& pa1 = L[dist.pa1];
        d.add_edge(pa1.name, dist.name, dist.rule_name);
        if ( recursive ) draw(pa1, recursive, d);
        if ( dist.pa2 > 0 ) {
            distr& pa2 = L[dist.pa2];
            d.add_edge(pa2.name, dist.name, dist.rule_name);
            if ( recursive ) draw(pa2, recursive, d);
        }
    }
}

string formulasearch::dec_to_text(const int& dec, const int& enabled) const {
    if ( dec == 0 ) return("");
    string s = "";
    int elems = set_size(dec);
    int el = get_next_element(dec, 0);
    if ( in_set(el, enabled) ) s += labels[el-1] + " = " + md_sym;
    else s += labels[el-1];
    for ( int i = 1; i < elems; i++ ) {
        el = get_next_element(dec, el);
        if ( in_set(el, enabled) ) s += "," + labels[el-1] + " = " + md_sym;
        else s += "," + labels[el-1];
    }
    return s;
}

string formulasearch::to_string(const p& pp) const {
    int u = pp.u;
    int c = pp.c;
    int j = pp.j;
    int e = pp.e;
    string s = "";

    if ( format_do ) {
        s += "p(" + dec_to_text(u, e);
        if ( j != 0 || c != 0 ) s += "|";
        if ( j != 0 ) s += "do(" + dec_to_text(j, e) + ")";
        c = c & (~j);
        if ( c != 0 ) {
            if ( j != 0 ) s += ",";
            s += dec_to_text(c, e);
        }
        s += ")";
    } else {
        s += "p(" + dec_to_text(u, e);
        if ( c != 0 ) s += "|" + dec_to_text(c, e);
        if ( j != 0 ) s += "||" + dec_to_text(j, e);
        s += ")";
    }

    return s;
}

bool formulasearch::equal_p(const p& pp1, const p& pp2) const {
    return ( (pp1.u == pp2.u) && (pp1.c == pp2.c) && (pp1.j == pp2.j) && (pp1.e == pp2.e) );
}

bool formulasearch::valid_do_rule(const int &ruleid, const int &u, const int &c, const int &j, const bool &primi) const {

    switch ( ruleid ) {

        case -1 : {
            // there must be observations to delete
            if ( (c - j) == 0 ) return false;
            else return true;
        }

        case -2 : {
             // there must be actions to exchange
            if ( j == 0 ) return false;
            else return true;
        }

        case -3 : {
            // there must be actions to delete
            if ( j == 0 ) return false;
            else return true;
        }

        case 2 : {
            // there must be observations to exchange
            if ( (c - j) == 0 ) return false;
            else return true;
        }

        case 4 : {
            // there must be other variables
            if ( set_size(u) == 1 ) return false;
            else return true;
        }

        case 5 : {
            // there must be other variables
            if ( set_size(u) == 1 ) return false;
            else return true;
        }

        case 6 : {
            // there must be conditioning variables
            if ( (c - j) == 0 ) return false;
            else return true;
        }

        case -7 : {
            // there must be conditioning variables
            if ( (c - j) == 0 ) return false;
            else return true;
        }

        case 7 : {
            // there must be other variables
            if ( set_size(u) == 1 ) return false;
            else return true;
        }

        /*
        case 8: {
            if ( primi ) return true;
            else return false;
        } */

        default : {
            return true;
        }
    }

    return true;
}

void formulasearch::apply_do_rule(const int &ruleid, const int &u, const int &c, const int &j, const int &e, const int &z) {

    int a, b, x, y, w, v, k;
    int d = 0;

    info.valid = false;

    switch ( ruleid ) {

        // Rule 1: Insertion of observations
        case 1 : {

            y = u;
            x = j;
            w = c;
            // z cannot be in any of the sets
            if ( (z & y) != 0 ) return; // z intersection y = 0
            if ( (z & x) != 0 ) return; // z intersection x = 0
            if ( (z & c) != 0 ) return; // z cannot have elements in c

            if ( md ) {
                d = e;
                a = ((y | c) & md_p) >> 2;
                if ( (a & z) != 0 ) return; // cannot add x if x* exists
                b = ((y | c) & md_t) << 2;
                if ( (b & z) != 0 ) return; // cannot add x* if x exists

                a = (z & md_p) >> 2;
                if ( (a & z) != 0) return; // cannot add both x and x*
                b = (z & md_t) << 2;
                if ( (b & z) != 0) return;// cannot add both x and x*
                // d = e | (z & md_s); // Always set inpedendent switches to zero
            }

            break;

        }

        // Rule 1: Deletion of observations
        case -1 : {

            y = u;
            x = j;
            if ( (z & y) != 0 ) return; // z intersection y = 0
            if ( (z & x) != 0 ) return; // z intersection x = 0
            // w OR z OR j = c => w = c - j - z
            if ( (z & (c - j)) != z ) return; // z must be a genuine subset of c
            w = c - j - z;

            if ( md ) d = e - (e & z); // Disable those z that were enabled

            break;

        }

        // Rule 2: Exchange actions to observations
        case -2 : {

            y = u;
            w = c - j;
            if ( (z & y) != 0 ) return; // z intersection y = 0
            if ( (z & w) != 0 ) return; // z intersection z = 0
            if ( (z & j) != z ) return; // z must be a subset of j
            // if ( (z & trsb) != 0 ) return; // z intersection s = 0
            // if ( (z & zu) != z ) return; // z must precede u;
            // if ( improve && ((z & zu) != z) ) return; // z must precede u;
            x = j - z;

            if ( md ) {
                if ( (z & md_p) != 0 ) return; // z intersection m = 0
                d = e;
            }

            break;

        }

        // Rule 2: Exchange observations to actions
        case 2 : {

            y = u;
            x = j;
            if ( (z & y) != 0 ) return; // z intersection y = 0
            if ( (z & x) != 0 ) return; // z intersection x = 0
            if ( (z & trsb) != 0 ) return; // z intersection s = 0
            if ( (z & c) != z ) return; // z must be a subset of c
            w = c - j - z;

            if ( md ) {
                if ( (z & md_p) != 0 ) return; // z intersection m = 0
                d = e | (md_s & z); // An action on a missing data switch will always set it to zero
            }

            break;

        }

        // Rule 3: Deletion of actions
        case -3 : {

            y = u;
            w = c - j;
            if ( (z & y) != 0 ) return; // z intersection y = 0
            if ( (z & w) != 0 ) return; // z intersection w = 0
            if ( (z & j) != z ) return; // z must be a subset of j
            x = j - z;

            if ( md ) {
                if ( (z & md_p) != 0 ) return; // z intersection m = 0
                d = e - (e & z); // Disable those z that were enabled
            }

            break;

        }

        // Rule 3: Insertion of actions
        case 3 : {

            y = u;
            w = c - j;
            x = j;
            if ( (z & y) != 0 ) return; // z intersection y = 0
            if ( (z & x) != 0 ) return; // z intersection x = 0
            if ( (z & w) != 0 ) return; // z intersection w = 0
            if ( (z & trsb) != 0 ) return; // z intersection s = 0

            if ( md ) {
                if ( (z & md_p) != 0 ) return; // z intersection m = 0

                a = ((y | c) & md_p) >> 2;
                if ( (a & z) != 0 ) return; // cannot add x if x* exists
                b = ((y | c) & md_t) << 2;
                if ( (b & z) != 0 ) return; // cannot add x* if x exists

                a = (z & md_p) >> 2;
                if ( (a & z) != 0) return; // cannot add both x and x*
                b = (z & md_t) << 2;
                if ( (b & z) != 0) return; // cannot add both x and x*
                d = e | (z & md_s); // An action on a missing data switch will always set it to zero
            }

            break;

        }

        // Rule 4: Marginalisation
        case 4 : {

            if ( (z & u) != z ) return; // z has to be in u
            if ( z == u ) return; // there has to be something else in u
            if ( (z & trsb) != 0 ) return; // z intersection s = 0
            x = j;
            y = u - z;
            w = c - j;

            if ( md ) {
                a = u & md_s;
                if ( (z & a) != 0 && (z & e) != 0 ) return; // Cannot marginalize over enabled switch
                d = e;
            }

            break;

        }

        // Rule 5: Conditioning
        case 5 : {

            if ( (z & u) != z ) return; // z has to be in u
            if ( z == u ) return; // there has to be something else in u
            if ( (z & trsb) != 0 ) return; // z intersection s = 0
            x = j;
            y = u - z;
            w = c - j;

            if ( md ) {
                a = u & md_s;
                if ( (y & a) != 0 && (y & e) != 0 ) return; // Cannot condition if switch is enabled
                d = e;
            }

            break;

        }

        // Rule 6: Product rule
        case 6 : {

            x = j;
            y = u;
            if ( (z & c) != z ) return; // z has to be in c
            if ( (z & j) != 0 ) return; // z cannot be in j
            if ( (z & trsb) != 0 ) return; // z cannot be in s
            w = c - x - z;
            d = e;

            break;

        }

        // Rule 6: Product rule
        case -6 : {

            x = j;
            y = u;
            w = c - j;
            // z cannot be in any of the sets
            if ( (z & y) != 0 ) return; // z intersection y = 0
            if ( (z & x) != 0 ) return; // z intersection x = 0
            if ( (z & w) != 0 ) return; // z intersection z = 0
            if ( (z & trsb) != 0 ) return; // z intersection s = 0
            d = e | (z & md_s);

            if ( md ) {
                a = ((y | c) & md_p) >> 2;
                if ( (a & z) != 0 ) return; // cannot add x if x* exists
                b = ((y | c) & md_t) << 2;
                if ( (b & z) != 0 ) return; // cannot add x* if x exists

                a = (z & md_p) >> 2;
                if ( (a & z) != 0) return; // cannot add both x and x*
                b = (z & md_t) << 2;
                if ( (b & z) != 0) return;// cannot add both x and x*
            }

            break;

        }

        // Rule 7: Conditioning via product rule (for missing data problems)
        case 7 : {

            if ( e == 0 ) return;
            if ( (z & u) != z ) return; // z has to be in u
            if ( z == u ) return; // there has to be something else in u
            if ( (z & trsb) != 0 ) return; // z intersection s = 0
            x = j;
            y = u - z;
            w = c - j;
            d = e - (e & z);

            break;

        }

        // Rule 7: Conditioning via product rule (for missing data problems)
        case -7 : {

            if ( e == 0 ) return;
            x = j;
            y = u;
            if ( (z & c) != z ) return; // z has to be in c
            if ( (z & j) != 0 ) return; // z cannot be in j
            if ( (z & trsb) != 0 ) return; // z cannot be in s
            w = c - j - z;
            d = e - (e & y);

            if ( md ) {
                a = ((y | c) & md_p) >> 2;
                if ( (a & z) != 0 ) return; // cannot add x if x* exists
                b = ((y | c) & md_t) << 2;
                if ( (b & z) != 0 ) return; // cannot add x* if x exists

                a = (z & md_p) >> 2;
                if ( (a & z) != 0) return; // cannot add both x and x*
                b = (z & md_t) << 2;
                if ( (b & z) != 0) return;// cannot add both x and x*
            }

            break;

        }

        // Rule 8: Enable missing data mechanism (for missing data problems)
        case 8 : {

            // if ( !iquery.primitive ) return;
            w = c - j;
            a = (u & md_s) | (w & md_s);
            if ( (z & a) != z ) return;
            if ( (z & e) != 0 ) return;
            y = u;
            x = j;
            d = e | z;

            break;

        }

        // Rule 9: Exchange proxy with true variable (for missing data problems)
        case 9 : {

            // if ( !iquery.primitive ) return;
            v = u & md_p; // Get proxy variables p(v*|.)
            k = c & md_p; // Get proxy variables p(.|k*) or p(.|do(k*))
            if ( (z & (v | k)) != z ) return; // z has to contain only proxy variables
            v = v & z;
            k = k & z;
            a = u & e; // Get missing data mechanisms that have been enabled p(.|r_a = 0)
            b = c & e; // Get missing data mechanisms that have been enabled p(r_b = 0|.)
            if ( (((a | b) << 1) & v) == v && ((b << 1) & k) == k ) {
                /* Check whether all proxy variables in z can be replaced
                 * Proxy variables p(v*|.) can be replaced by true variables either if p(r_v = 0,v*|.) or p(v*|r_v* = 0)
                 * However, proxies p(.|k*) can only be replaced by true variables if p(.|k*,r_k* = 0)
                 */
                y = (u - v) + (v >> 2); // Replace proxies with their respective true variables
                w = (c - k) + (k >> 2); // Replace proxies with their respective true variables
                x = j;
                w = c - j;
                d = e;
            } else return;

            break;

        }

    }

    info.valid = true;
    get_ruleinfo(ruleid, y, x, z, w, d, e);

}


void formulasearch::get_ruleinfo(const int& ruleid, const int& y, const int& x, const int& z, const int& w, const int& d, const int& e) {

    string r_str = "";

    int xw = x | w;
    int xz = x | z;
    int xzw = (x | z) | w;
    int yz = y | z;
    int xyw = (x | y) | w;

    switch ( ruleid ) {

        case 1 : {

            info.rule_name = "R1";

            info.to.u = y; info.to.c = xzw; info.to.j = x; info.to.e = d;
            info.from.u = y; info.from.c = xw; info.from.j = x; info.from.e = e;
            info.ri.xset = y; info.ri.yset = z; info.ri.c = xw; info.ri.j = x;
            info.rp.u = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " << to_string(info.from) << " - RULE 1: insertion of observations " << r_str << endl;
            }

            return;

        }

        case -1 : {

            info.rule_name = "R1";

            info.to.u = y; info.to.c = xw; info.to.j = x; info.to.e = d;
            info.from.u = y; info.from.c = xzw; info.from.j = x; info.from.e = e;
            info.ri.xset = y; info.ri.yset = z; info.ri.c = xw; info.ri.j = x;
            info.rp.u = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " << to_string(info.from) << " - RULE 1: deletion of observations " << r_str << endl;
            }

            return;

        }

        case 2 : {

            info.rule_name = "R2";

            info.to.u = y; info.to.c = xzw; info.to.j = xz; info.to.e = d;
            info.from.u = y; info.from.c = xzw; info.from.j = x; info.from.e = e;
            info.ri.xset = y; info.ri.yset = z << n; info.ri.c = xzw; info.ri.j = x;
            info.rp.u = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " << to_string(info.from) << " - RULE 2: observation to action exchange " << r_str << endl;
            }

            return;

        }

        case -2 : {

            info.rule_name = "R2";

            info.to.u = y; info.to.c = xzw; info.to.j = x; info.to.e = d;
            info.from.u = y; info.from.c = xzw; info.from.j = xz; info.from.e = e;
            info.ri.xset = y; info.ri.yset = z << n; info.ri.c = xzw; info.ri.j = x;
            info.rp.u = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " << to_string(info.from) << " - RULE 2: action to observation exchange " << r_str << endl;
            }

            return;

        }

        case 3 : {

            info.rule_name = "R3";

            info.to.u = y; info.to.c = xzw; info.to.j = xz; info.to.e = d;
            info.from.u = y; info.from.c = xw; info.from.j = x; info.from.e = e;
            info.ri.xset = y; info.ri.yset = z << n; info.ri.c = xw; info.ri.j = x;
            info.rp.u = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " << to_string(info.from) << " - RULE 3: insertion of actions " << r_str << endl;
            }

            return;

        }

        case -3 : {

            info.rule_name = "R3";

            info.to.u = y; info.to.c = xw; info.to.j = x; info.to.e = d;
            info.from.u = y; info.from.c = xzw; info.from.j = xz; info.from.e = e;
            info.ri.xset = y; info.ri.yset = z << n; info.ri.c = xw; info.ri.j = x;
            info.rp.u = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " << to_string(info.from) << " - RULE 3: deletion of actions " << r_str << endl;
            }

            return;

        }

        case 4 : {

            info.rule_name = "M";

            info.to.u = y; info.to.c = xw; info.to.j = x; info.to.e = d;
            info.from.u = yz; info.from.c = xw; info.from.j = x; info.from.e = e;
            info.rp.u = 0;
            info.ri.xset = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " << to_string(info.from) << " - RULE 4: marginalisation " << r_str << endl;
            }

            return;

        }

        case 5 : {

            info.rule_name = "C";

            info.to.u = y; info.to.c = xzw; info.to.j = x; info.to.e = d;
            info.from.u = yz; info.from.c = xw; info.from.j = x; info.from.e = e;
            info.rp.u = 0;
            info.ri.xset = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " << to_string(info.from) << " - RULE 5: conditioning " << r_str << endl;
            }

            return;

        }

        case 6 : {

            info.rule_name = "P";

            info.to.u = yz; info.to.c = xw; info.to.j = x; info.to.e = d;
            info.from.u = y; info.from.c = xzw; info.from.j = x; info.from.e = e;
            info.rp.u = z; info.rp.c = xw; info.rp.j = x; info.rp.e = e - (e & y);
            info.ri.xset = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " <<  to_string(info.from) << " AND [" << to_string(info.rp)  << " req.] - RULE 6: product rule (a) " << r_str << endl;
            }

            return;

        }

        case -6 : {

            info.rule_name = "P";

            info.to.u = yz; info.to.c = xw; info.to.j = x; info.to.e = d;
            info.from.u = y; info.from.c = xw; info.from.j = x; info.from.e = e;
            info.rp.u = z; info.rp.c = xyw; info.rp.j = x; info.rp.e = d;
            info.ri.xset = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " <<  to_string(info.from) << " AND [" << to_string(info.rp)  << " req.] - RULE 6: product rule (b) " << r_str << endl;
            }

            return;

        }

        case 7 : {

            info.rule_name = "D";

            info.to.u = y; info.to.c = xw; info.to.j = x; info.to.e = d;
            info.from.u = yz; info.from.c = xw; info.from.j = x; info.from.e = e;
            info.rp.u = z; info.rp.c = xyw; info.rp.j = x; info.rp.e = e;
            info.ri.xset = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " <<  to_string(info.from) << " AND [" << to_string(info.rp)  << " req.] - RULE 7: Bayes' rule (a)" << r_str << endl;
            }

            return;

        }

        case -7 : {

            info.rule_name = "D";

            info.to.u = z; info.to.c = xw; info.to.j = x; info.to.e = d;
            info.from.u = y; info.from.c = xzw; info.from.j = x; info.from.e = e;
            info.rp.u = yz; info.rp.c = xw; info.rp.j = x; info.rp.e = e;
            info.ri.xset = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " <<  to_string(info.from) << " AND [" << to_string(info.rp)  << " req.] - RULE 7: Bayes' rule (b)" << r_str << endl;
            }

            return;

        }

        case 8 : {

            info.rule_name = "0";

            info.to.u = y; info.to.c = xw; info.to.j = x; info.to.e = d;
            info.from.u = y; info.from.c = xw; info.from.j = x; info.from.e = e;
            info.rp.u = 0;
            info.ri.xset = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " << to_string(info.from) << " - RULE 8: Enable missing data mechanisms" << r_str << endl;
            }

            return;

        }

        case 9 : {

            info.rule_name = "EX";

            int v = (z >> 2) & y;
            int k = (z >> 2) & xw;
            int u = (y - v) + (v << 2);
            int c = (xw - k) + (k << 2);

            info.to.u = y; info.to.c = xw; info.to.j = x; info.to.e = d;
            info.from.u = u; info.from.c = c; info.from.j = x; info.from.e = e;
            info.rp.u = 0;
            info.ri.xset = 0;

            if ( verbose ) {
                Rcpp::Rcout << to_string(info.to) << " <= " << to_string(info.from) << " - RULE 9: Exchange proxy with true variable" << r_str << endl;
            }

            return;

        }

    }
}
