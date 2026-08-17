/// C-ABI Bridge types and interface for native graph isomorphism matchers (VF3 / RDKit C++).

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct FfiCandidateMolecule {
    pub cid: u32,
    pub mol_block: *const u8,
    pub mol_block_len: usize,
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct FfiMatchResult {
    pub cid: u32,
    pub is_match: bool,
    pub graph_score: f32,
}

extern "C" {
    // Stub definition for future Phase 6 native library binding
    // fn molir_verify_substructure(
    //     query_smarts: *const std::ffi::c_char,
    //     candidates: *const FfiCandidateMolecule,
    //     candidate_count: usize,
    //     out_results: *mut FfiMatchResult,
    // ) -> i32;
}
