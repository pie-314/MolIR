use serde::{Deserialize, Serialize};

/// Physicochemical properties used for multi-parameter scoring and filtering.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MolecularProperties {
    pub molecular_weight: f32,
    pub logp: f32,
    pub tpsa: f32,
    pub hbd: u16,
    pub hba: u16,
    pub rotatable_bonds: u16,
}

impl Default for MolecularProperties {
    fn default() -> Self {
        Self {
            molecular_weight: 0.0,
            logp: 0.0,
            tpsa: 0.0,
            hbd: 0,
            hba: 0,
            rotatable_bonds: 0,
        }
    }
}

/// Property constraint filter definition.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PropertyFilter {
    pub min_mw: Option<f32>,
    pub max_mw: Option<f32>,
    pub min_logp: Option<f32>,
    pub max_logp: Option<f32>,
    pub max_tpsa: Option<f32>,
}

impl PropertyFilter {
    pub fn matches(&self, props: &MolecularProperties) -> bool {
        if let Some(min) = self.min_mw {
            if props.molecular_weight < min {
                return false;
            }
        }
        if let Some(max) = self.max_mw {
            if props.molecular_weight > max {
                return false;
            }
        }
        if let Some(min) = self.min_logp {
            if props.logp < min {
                return false;
            }
        }
        if let Some(max) = self.max_logp {
            if props.logp > max {
                return false;
            }
        }
        if let Some(max) = self.max_tpsa {
            if props.tpsa > max {
                return false;
            }
        }
        true
    }
}
