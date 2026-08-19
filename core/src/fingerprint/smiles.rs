use std::collections::HashMap;
use anyhow::{anyhow, Result};
use crate::fingerprint::MolecularFingerprint;

/// Bond order representation in a molecular graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BondType {
    Single = 1,
    Double = 2,
    Triple = 3,
    Aromatic = 4,
}

/// Atom attributes perceived from a parsed SMILES representation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedAtom {
    pub atomic_num: u8,
    pub is_aromatic: bool,
    pub formal_charge: i8,
    pub explicit_h: u8,
    pub implicit_h: u8,
    pub degree: u8,
}

/// A parsed molecular chemical graph.
#[derive(Debug, Clone, Default)]
pub struct MolecularGraph {
    pub atoms: Vec<ParsedAtom>,
    pub adjacency: Vec<Vec<(usize, BondType)>>,
}

impl MolecularGraph {
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            adjacency: Vec::new(),
        }
    }

    pub fn add_atom(&mut self, atomic_num: u8, is_aromatic: bool, formal_charge: i8, explicit_h: u8) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(ParsedAtom {
            atomic_num,
            is_aromatic,
            formal_charge,
            explicit_h,
            implicit_h: 0,
            degree: 0,
        });
        self.adjacency.push(Vec::new());
        idx
    }

    pub fn add_bond(&mut self, u: usize, v: usize, bond_type: BondType) {
        if u < self.atoms.len() && v < self.atoms.len() && u != v {
            self.adjacency[u].push((v, bond_type));
            self.adjacency[v].push((u, bond_type));
            self.atoms[u].degree += 1;
            self.atoms[v].degree += 1;
        }
    }

    /// Computes typical organic implicit hydrogen counts based on valence.
    pub fn perceive_hydrogens(&mut self) {
        for i in 0..self.atoms.len() {
            let atom = &mut self.atoms[i];
            if atom.explicit_h > 0 {
                atom.implicit_h = atom.explicit_h;
                continue;
            }

            // Standard organic default valences
            let default_val = match atom.atomic_num {
                6 => if atom.is_aromatic { 3 } else { 4 }, // Carbon
                7 => if atom.is_aromatic { 3 } else { 3 }, // Nitrogen
                8 => if atom.is_aromatic { 2 } else { 2 }, // Oxygen
                9 => 1,                                    // Fluorine
                15 => 3,                                   // Phosphorus
                16 => if atom.is_aromatic { 2 } else { 2 },// Sulfur
                17 => 1,                                   // Chlorine
                35 => 1,                                   // Bromine
                53 => 1,                                   // Iodine
                5 => 3,                                    // Boron
                _ => 0,
            };

            let mut bond_sum: u8 = 0;
            for &(_, b_type) in &self.adjacency[i] {
                bond_sum += match b_type {
                    BondType::Single => 1,
                    BondType::Double => 2,
                    BondType::Triple => 3,
                    BondType::Aromatic => 1,
                };
            }

            if default_val > bond_sum {
                atom.implicit_h = default_val - bond_sum;
            } else {
                atom.implicit_h = 0;
            }
        }
    }

    /// Generates a 2048-bit Morgan / ECFP4 fingerprint (radius = 2).
    pub fn to_ecfp4(&mut self) -> MolecularFingerprint {
        self.perceive_hydrogens();
        let n_atoms = self.atoms.len();
        if n_atoms == 0 {
            return MolecularFingerprint::zeros();
        }

        let mut fp = MolecularFingerprint::zeros();

        // 1. Initial atom invariants (radius 0)
        let mut current_invariants = Vec::with_capacity(n_atoms);
        for atom in &self.atoms {
            let inv = hash_initial_atom(atom);
            fp.set_bit((inv % 2048) as usize);
            current_invariants.push(inv);
        }

        // 2. Circular neighborhood iteration (radius 1 and 2)
        for radius in 1..=2u32 {
            let mut next_invariants = Vec::with_capacity(n_atoms);

            for i in 0..n_atoms {
                let mut neighbor_tuples = Vec::new();
                for &(neighbor_idx, bond_type) in &self.adjacency[i] {
                    neighbor_tuples.push((bond_type as u32, current_invariants[neighbor_idx]));
                }

                // Lexicographical sorting guarantees canonical invariance
                neighbor_tuples.sort_unstable();

                let mut h = fnv1a_init();
                h = fnv1a_update(h, radius);
                h = fnv1a_update(h, current_invariants[i]);
                for (b_code, n_inv) in neighbor_tuples {
                    h = fnv1a_update(h, b_code);
                    h = fnv1a_update(h, n_inv);
                }

                let new_inv = h as u32;
                fp.set_bit((new_inv % 2048) as usize);
                next_invariants.push(new_inv);
            }

            current_invariants = next_invariants;
        }

        fp
    }
}

// 64-bit FNV-1a Hash functions for fast deterministic hashing
const FNV_OFFSET_BASIS: u64 = 0xcbf29ce484222325;
const FNV_PRIME: u64 = 0x100000001b3;

#[inline]
fn fnv1a_init() -> u64 {
    FNV_OFFSET_BASIS
}

#[inline]
fn fnv1a_update(mut hash: u64, value: u32) -> u64 {
    for byte in value.to_le_bytes() {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

#[inline]
fn hash_initial_atom(atom: &ParsedAtom) -> u32 {
    let mut h = fnv1a_init();
    h = fnv1a_update(h, atom.atomic_num as u32);
    h = fnv1a_update(h, atom.degree as u32);
    h = fnv1a_update(h, (atom.implicit_h + atom.explicit_h) as u32);
    h = fnv1a_update(h, (atom.formal_charge as i32) as u32);
    h = fnv1a_update(h, if atom.is_aromatic { 1 } else { 0 });
    h as u32
}

/// Parses a standard chemical SMILES string into a MolecularGraph.
pub fn parse_smiles(smiles: &str) -> Result<MolecularGraph> {
    let smiles = smiles.trim();
    if smiles.is_empty() {
        return Err(anyhow!("Empty SMILES string provided"));
    }

    let mut graph = MolecularGraph::new();
    let chars: Vec<char> = smiles.chars().collect();
    let mut pos = 0;

    let mut prev_atom: Option<usize> = None;
    let mut pending_bond: Option<BondType> = None;
    let mut branch_stack: Vec<Option<usize>> = Vec::new();
    let mut ring_closures: HashMap<u32, (usize, Option<BondType>)> = HashMap::new();

    while pos < chars.len() {
        let ch = chars[pos];

        match ch {
            // Branches
            '(' => {
                branch_stack.push(prev_atom);
                pos += 1;
            }
            ')' => {
                if let Some(parent) = branch_stack.pop() {
                    prev_atom = parent;
                } else {
                    return Err(anyhow!("Unmatched closing parenthesis ')' at pos {}", pos));
                }
                pos += 1;
            }

            // Bonds
            '-' => {
                pending_bond = Some(BondType::Single);
                pos += 1;
            }
            '=' => {
                pending_bond = Some(BondType::Double);
                pos += 1;
            }
            '#' => {
                pending_bond = Some(BondType::Triple);
                pos += 1;
            }
            ':' => {
                pending_bond = Some(BondType::Aromatic);
                pos += 1;
            }
            '/' | '\\' => {
                // Directional stereobond; treat as single bond for ECFP4
                pending_bond = Some(BondType::Single);
                pos += 1;
            }
            '.' => {
                // Disconnected structures (salts/mixtures)
                prev_atom = None;
                pending_bond = None;
                pos += 1;
            }

            // Ring closures (1..9 and %10..%99)
            '0'..='9' => {
                let ring_num = ch.to_digit(10).unwrap() as u32;
                handle_ring_closure(&mut graph, ring_num, prev_atom, &mut pending_bond, &mut ring_closures)?;
                pos += 1;
            }
            '%' => {
                pos += 1;
                if pos + 1 < chars.len() && chars[pos].is_ascii_digit() && chars[pos + 1].is_ascii_digit() {
                    let d1 = chars[pos].to_digit(10).unwrap();
                    let d2 = chars[pos + 1].to_digit(10).unwrap();
                    let ring_num = d1 * 10 + d2;
                    handle_ring_closure(&mut graph, ring_num, prev_atom, &mut pending_bond, &mut ring_closures)?;
                    pos += 2;
                } else {
                    return Err(anyhow!("Invalid two-digit ring closure '%' syntax at pos {}", pos));
                }
            }

            // Bracketed atom: [NH4+], [OH-], [13C], [nH], etc.
            '[' => {
                let end_bracket = chars[pos..].iter().position(|&c| c == ']');
                if let Some(offset) = end_bracket {
                    let bracket_content: String = chars[pos + 1..pos + offset].iter().collect();
                    let (atom_num, is_arom, charge, expl_h) = parse_bracket_atom(&bracket_content)?;
                    let current = graph.add_atom(atom_num, is_arom, charge, expl_h);

                    attach_atom(&mut graph, prev_atom, current, &mut pending_bond, is_arom);
                    prev_atom = Some(current);
                    pos += offset + 1;
                } else {
                    return Err(anyhow!("Unclosed bracket '[' at pos {}", pos));
                }
            }

            // Organic subset atoms
            'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'I' | 'B' | 'c' | 'n' | 'o' | 's' | 'p' => {
                let (atom_num, is_arom, consumed) = parse_organic_atom(&chars[pos..])?;
                let current = graph.add_atom(atom_num, is_arom, 0, 0);

                attach_atom(&mut graph, prev_atom, current, &mut pending_bond, is_arom);
                prev_atom = Some(current);
                pos += consumed;
            }

            // Cl and Br special 2-letter organic cases
            _ => {
                return Err(anyhow!("Unexpected character '{}' at pos {}", ch, pos));
            }
        }
    }

    if !ring_closures.is_empty() {
        return Err(anyhow!("Unclosed ring closure numbers: {:?}", ring_closures.keys()));
    }

    Ok(graph)
}

fn attach_atom(
    graph: &mut MolecularGraph,
    prev_atom: Option<usize>,
    current: usize,
    pending_bond: &mut Option<BondType>,
    is_aromatic: bool,
) {
    if let Some(prev) = prev_atom {
        let b_type = pending_bond.take().unwrap_or_else(|| {
            if is_aromatic && graph.atoms[prev].is_aromatic {
                BondType::Aromatic
            } else {
                BondType::Single
            }
        });
        graph.add_bond(prev, current, b_type);
    }
}

fn handle_ring_closure(
    graph: &mut MolecularGraph,
    ring_num: u32,
    current_atom: Option<usize>,
    pending_bond: &mut Option<BondType>,
    ring_closures: &mut HashMap<u32, (usize, Option<BondType>)>,
) -> Result<()> {
    let current = current_atom.ok_or_else(|| anyhow!("Ring number {} without a preceding atom", ring_num))?;

    if let Some((first_atom, first_bond)) = ring_closures.remove(&ring_num) {
        let b_type = pending_bond
            .take()
            .or(first_bond)
            .unwrap_or_else(|| {
                if graph.atoms[current].is_aromatic && graph.atoms[first_atom].is_aromatic {
                    BondType::Aromatic
                } else {
                    BondType::Single
                }
            });
        graph.add_bond(first_atom, current, b_type);
    } else {
        ring_closures.insert(ring_num, (current, pending_bond.take()));
    }
    Ok(())
}

fn parse_organic_atom(slice: &[char]) -> Result<(u8, bool, usize)> {
    if slice.is_empty() {
        return Err(anyhow!("End of input reached while parsing organic atom"));
    }

    // 2-letter elements: Cl, Br
    if slice.len() >= 2 {
        if slice[0] == 'C' && slice[1] == 'l' {
            return Ok((17, false, 2));
        }
        if slice[0] == 'B' && slice[1] == 'r' {
            return Ok((35, false, 2));
        }
    }

    let (atom_num, is_arom) = match slice[0] {
        'C' => (6, false),
        'N' => (7, false),
        'O' => (8, false),
        'S' => (16, false),
        'P' => (15, false),
        'F' => (9, false),
        'I' => (53, false),
        'B' => (5, false),
        'c' => (6, true),
        'n' => (7, true),
        'o' => (8, true),
        's' => (16, true),
        'p' => (15, true),
        ch => return Err(anyhow!("Unknown organic atom '{}'", ch)),
    };

    Ok((atom_num, is_arom, 1))
}

fn parse_bracket_atom(content: &str) -> Result<(u8, bool, i8, u8)> {
    let mut atom_str = String::new();
    let mut charge: i8 = 0;
    let mut explicit_h: u8 = 0;
    let chars: Vec<char> = content.chars().collect();
    let mut i = 0;

    // Skip optional leading isotope digits (e.g. 13C)
    while i < chars.len() && chars[i].is_ascii_digit() {
        i += 1;
    }

    // Extract chemical symbol
    if i < chars.len() {
        atom_str.push(chars[i]);
        i += 1;
        if i < chars.len() && chars[i].is_ascii_lowercase() && !['h', '@'].contains(&chars[i]) {
            atom_str.push(chars[i]);
            i += 1;
        }
    }

    // Process remainder (chirality, hydrogens, charges)
    while i < chars.len() {
        match chars[i] {
            '@' => {
                // Ignore stereochemistry for 2D topological ECFP4
                i += 1;
            }
            'H' => {
                i += 1;
                if i < chars.len() && chars[i].is_ascii_digit() {
                    explicit_h = chars[i].to_digit(10).unwrap() as u8;
                    i += 1;
                } else {
                    explicit_h = 1;
                }
            }
            '+' => {
                i += 1;
                if i < chars.len() && chars[i].is_ascii_digit() {
                    charge = chars[i].to_digit(10).unwrap() as i8;
                    i += 1;
                } else if i < chars.len() && chars[i] == '+' {
                    charge = 2;
                    i += 1;
                } else {
                    charge = 1;
                }
            }
            '-' => {
                i += 1;
                if i < chars.len() && chars[i].is_ascii_digit() {
                    charge = -(chars[i].to_digit(10).unwrap() as i8);
                    i += 1;
                } else if i < chars.len() && chars[i] == '-' {
                    charge = -2;
                    i += 1;
                } else {
                    charge = -1;
                }
            }
            _ => {
                i += 1;
            }
        }
    }

    let (atomic_num, is_aromatic) = match atom_str.as_str() {
        "H" => (1, false),
        "C" => (6, false),
        "c" => (6, true),
        "N" => (7, false),
        "n" => (7, true),
        "O" => (8, false),
        "o" => (8, true),
        "F" => (9, false),
        "P" => (15, false),
        "p" => (15, true),
        "S" => (16, false),
        "s" => (16, true),
        "Cl" => (17, false),
        "Br" => (35, false),
        "I" => (53, false),
        "Na" => (11, false),
        "K" => (19, false),
        "Ca" => (20, false),
        "Fe" => (26, false),
        "Zn" => (30, false),
        "Se" => (34, false),
        "se" => (34, true),
        "Si" => (14, false),
        other => return Err(anyhow!("Unsupported element symbol in bracket: [{}]", other)),
    };

    Ok((atomic_num, is_aromatic, charge, explicit_h))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smiles_parsing_simple() {
        let mut g = parse_smiles("CC(=O)O").expect("Valid Acetic Acid");
        assert_eq!(g.atoms.len(), 4);
        let fp = g.to_ecfp4();
        assert!(fp.popcount() > 0);
    }

    #[test]
    fn test_smiles_aromatic_rings() {
        let mut g_benzene = parse_smiles("c1ccccc1").expect("Benzene");
        assert_eq!(g_benzene.atoms.len(), 6);
        let fp_benzene = g_benzene.to_ecfp4();
        assert!(fp_benzene.popcount() > 0);

        let mut g_pyridine = parse_smiles("c1ccncc1").expect("Pyridine");
        let fp_pyridine = g_pyridine.to_ecfp4();
        
        let sim = fp_benzene.tanimoto(&fp_pyridine);
        assert!(sim > 0.3 && sim < 1.0);
    }

    #[test]
    fn test_smiles_aspirin_identity() {
        let mut g1 = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").expect("Aspirin");
        let mut g2 = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").expect("Aspirin");
        let fp1 = g1.to_ecfp4();
        let fp2 = g2.to_ecfp4();
        assert!((fp1.tanimoto(&fp2) - 1.0).abs() < 1e-6);
    }
}
