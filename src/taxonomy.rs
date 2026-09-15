//! Taxonomic assignment by lowest common ancestor (LCA), mirroring MMseqs2's
//! `taxonomy` / `lca` modules.
//!
//! Inputs are the standard NCBI taxdump pieces plus a sequence→taxon map:
//! * `nodes.dmp` — `taxid | parent_taxid | rank | ...` (the tree).
//! * `names.dmp` — `taxid | name | | name class | ...` (scientific names).
//! * a `seqid<TAB>taxid` table mapping each target accession to a taxon.
//!
//! For each query we take its hits, map every target to its taxon, and assign the
//! deepest node that is an ancestor of them all. A "weighted" variant keeps only
//! the hits whose bit score is within a fraction of the best (MMseqs2's
//! approximate-2bLCA idea), which resists a single spurious hit dragging the
//! assignment up to the root. We consume a taxdump rather than bundling one.

use std::collections::{HashMap, HashSet};

/// A parsed NCBI taxonomy tree.
pub struct TaxTree {
    parent: HashMap<u32, u32>,
    rank: HashMap<u32, String>,
    name: HashMap<u32, String>,
}

/// Split a `.dmp` line on the `\t|\t` (and trailing `\t|`) field separator.
fn dmp_fields(line: &str) -> Vec<&str> {
    let trimmed = line.trim_end_matches(|c| c == '|' || c == '\t' || c == '\n' || c == '\r');
    trimmed.split("\t|\t").map(|s| s.trim()).collect()
}

impl TaxTree {
    /// Parse `nodes.dmp` content: fields are taxid, parent, rank.
    pub fn parse_nodes(text: &str) -> Self {
        let mut parent = HashMap::new();
        let mut rank = HashMap::new();
        for line in text.lines() {
            let f = dmp_fields(line);
            if f.len() < 3 {
                continue;
            }
            if let (Ok(t), Ok(p)) = (f[0].parse::<u32>(), f[1].parse::<u32>()) {
                parent.insert(t, p);
                rank.insert(t, f[2].to_string());
            }
        }
        TaxTree {
            parent,
            rank,
            name: HashMap::new(),
        }
    }

    /// Merge scientific names from `names.dmp` content (name class == "scientific
    /// name"). Fields are taxid, name, unique-name, name-class.
    pub fn load_names(&mut self, text: &str) {
        for line in text.lines() {
            let f = dmp_fields(line);
            if f.len() < 4 {
                continue;
            }
            if f[3] == "scientific name" {
                if let Ok(t) = f[0].parse::<u32>() {
                    self.name.insert(t, f[1].to_string());
                }
            }
        }
    }

    pub fn rank_of(&self, taxid: u32) -> &str {
        self.rank.get(&taxid).map(|s| s.as_str()).unwrap_or("no rank")
    }

    pub fn name_of(&self, taxid: u32) -> &str {
        self.name.get(&taxid).map(|s| s.as_str()).unwrap_or("unclassified")
    }

    pub fn contains(&self, taxid: u32) -> bool {
        self.parent.contains_key(&taxid)
    }

    /// Ancestor chain `[taxid, parent, …, root]` (root is where `parent == self`
    /// or the parent is missing). Cycle-safe.
    pub fn ancestors(&self, taxid: u32) -> Vec<u32> {
        let mut path = vec![taxid];
        let mut seen = HashSet::new();
        seen.insert(taxid);
        let mut cur = taxid;
        while let Some(&p) = self.parent.get(&cur) {
            if p == cur || !seen.insert(p) {
                break;
            }
            path.push(p);
            cur = p;
        }
        path
    }

    /// Lowest common ancestor of a set of taxa. Unknown taxa are ignored; returns
    /// `None` only if none of the inputs are in the tree.
    pub fn lca(&self, taxids: &[u32]) -> Option<u32> {
        let known: Vec<u32> = taxids.iter().copied().filter(|&t| self.contains(t)).collect();
        let (first, rest) = known.split_first()?;
        // ancestors(first) is ordered deepest→root; keep only nodes ancestral to
        // every other taxon, then the first survivor is the deepest = the LCA.
        let mut common = self.ancestors(*first);
        for &t in rest {
            let anc: HashSet<u32> = self.ancestors(t).into_iter().collect();
            common.retain(|x| anc.contains(x));
            if common.is_empty() {
                break;
            }
        }
        common.first().copied()
    }
}

/// Parse a `seqid<TAB>taxid` mapping table.
pub fn parse_seqmap(text: &str) -> HashMap<String, u32> {
    let mut m = HashMap::new();
    for line in text.lines() {
        let mut it = line.split(|c| c == '\t' || c == ',');
        let id = it.next().map(|s| s.trim());
        let taxid = it.next().and_then(|s| s.trim().parse::<u32>().ok());
        if let (Some(id), Some(taxid)) = (id, taxid) {
            if !id.is_empty() {
                m.insert(id.to_string(), taxid);
            }
        }
    }
    m
}

/// Keep the taxa used for LCA. With `frac <= 0` **all** mapped hits are kept
/// (plain LCA of every hit). With `frac > 0` only hits whose score is within that
/// fraction of the best are kept (MMseqs2's approximate-2bLCA idea), which resists
/// a single spurious hit dragging the assignment up toward the root. `scored` is
/// `(taxid, score)`.
pub fn weighted_taxa(scored: &[(u32, f64)], frac: f64) -> Vec<u32> {
    if scored.is_empty() {
        return Vec::new();
    }
    if frac <= 0.0 {
        return scored.iter().map(|&(t, _)| t).collect();
    }
    let best = scored.iter().map(|&(_, s)| s).fold(f64::MIN, f64::max);
    let cutoff = best * (1.0 - frac);
    scored
        .iter()
        .filter(|&&(_, s)| s >= cutoff)
        .map(|&(t, _)| t)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // A tiny tree:
    //   root(1) ─ Bacteria(2) ─ Bacillus genus(3) ─ B.subtilis(4)
    //                                              └ B.cereus(5)
    //                         └ E.coli genus(6)   ─ E.coli(7)
    const NODES: &str = "\
1\t|\t1\t|\tno rank\t|\n\
2\t|\t1\t|\tsuperkingdom\t|\n\
3\t|\t2\t|\tgenus\t|\n\
4\t|\t3\t|\tspecies\t|\n\
5\t|\t3\t|\tspecies\t|\n\
6\t|\t2\t|\tgenus\t|\n\
7\t|\t6\t|\tspecies\t|\n";

    const NAMES: &str = "\
1\t|\troot\t|\t\t|\tscientific name\t|\n\
2\t|\tBacteria\t|\t\t|\tscientific name\t|\n\
3\t|\tBacillus\t|\t\t|\tscientific name\t|\n\
4\t|\tBacillus subtilis\t|\t\t|\tscientific name\t|\n\
5\t|\tBacillus cereus\t|\t\t|\tscientific name\t|\n\
6\t|\tEscherichia\t|\t\t|\tscientific name\t|\n\
7\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n";

    fn tree() -> TaxTree {
        let mut t = TaxTree::parse_nodes(NODES);
        t.load_names(NAMES);
        t
    }

    #[test]
    fn parses_tree_and_names() {
        let t = tree();
        assert_eq!(t.rank_of(3), "genus");
        assert_eq!(t.name_of(4), "Bacillus subtilis");
        assert_eq!(t.ancestors(4), vec![4, 3, 2, 1]);
    }

    #[test]
    fn lca_of_sister_species_is_their_genus() {
        let t = tree();
        // B. subtilis (4) + B. cereus (5) -> Bacillus (3)
        assert_eq!(t.lca(&[4, 5]), Some(3));
        assert_eq!(t.rank_of(t.lca(&[4, 5]).unwrap()), "genus");
    }

    #[test]
    fn lca_across_genera_is_the_superkingdom() {
        let t = tree();
        // B. subtilis (4) + E. coli (7) -> Bacteria (2)
        assert_eq!(t.lca(&[4, 7]), Some(2));
    }

    #[test]
    fn lca_of_single_taxon_is_itself() {
        let t = tree();
        assert_eq!(t.lca(&[7]), Some(7));
    }

    #[test]
    fn unknown_taxa_are_ignored() {
        let t = tree();
        // 999 not in tree; 4 and 5 -> genus 3
        assert_eq!(t.lca(&[4, 999, 5]), Some(3));
        assert_eq!(t.lca(&[999]), None);
    }

    #[test]
    fn weighted_filter_drops_low_scores() {
        // best 100; frac 0.1 keeps >=90. The spurious taxon at 50 is dropped, so
        // a weighted LCA of {sister species} stays at the genus instead of root.
        let scored = vec![(4u32, 100.0), (5u32, 95.0), (7u32, 50.0)];
        let kept = weighted_taxa(&scored, 0.1);
        assert!(kept.contains(&4) && kept.contains(&5));
        assert!(!kept.contains(&7));
        let t = tree();
        assert_eq!(t.lca(&kept), Some(3)); // genus, not superkingdom
    }

    #[test]
    fn frac_zero_keeps_all_hits() {
        // frac 0 = plain LCA of every hit: sisters 4 & 5 -> genus 3
        let scored = vec![(4u32, 100.0), (5u32, 80.0)];
        let kept = weighted_taxa(&scored, 0.0);
        assert_eq!(kept.len(), 2);
        assert_eq!(tree().lca(&kept), Some(3));
    }

    #[test]
    fn seqmap_parses() {
        let m = parse_seqmap("seqA\t4\nseqB\t7\nseqC,5\n");
        assert_eq!(m.get("seqA"), Some(&4));
        assert_eq!(m.get("seqB"), Some(&7));
        assert_eq!(m.get("seqC"), Some(&5));
    }
}
