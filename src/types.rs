// src/types.rs

use std::fmt;

use crate::{body::Body, condensation::RegionThermalProfile};

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum MassType {
    Star,
    Planet,
    Moon,
    GasGiant,
}

impl fmt::Display for MassType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            MassType::Star => write!(f, "Star"),
            MassType::Planet => write!(f, "Planet"),
            MassType::Moon => write!(f, "Moon"),
            MassType::GasGiant => write!(f, "Gas Giant"),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum StellarTopology {
    Single { primary_index: usize },
    Binary(BinaryTopology),
    HigherMultiplicity { primary_index: usize, star_count: usize },
}

#[derive(Debug, Clone, PartialEq)]
pub struct BinaryTopology {
    pub primary_index: usize,
    pub secondary_index: usize,
    pub semi_major_axis_au: f64,
    pub eccentricity: f64,
    pub total_mass_solar: f64,
    pub mass_ratio: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct StellarHierarchy {
    pub root_index: usize,
    pub nodes: Vec<StellarHierarchyNode>,
}

impl StellarHierarchy {
    pub fn single(star_index: usize) -> Self {
        Self {
            root_index: 0,
            nodes: vec![StellarHierarchyNode {
                parent_index: None,
                kind: StellarHierarchyNodeKind::Star { star_index },
            }],
        }
    }

    pub fn unresolved(star_indices: Vec<usize>) -> Self {
        Self {
            root_index: 0,
            nodes: vec![StellarHierarchyNode {
                parent_index: None,
                kind: StellarHierarchyNodeKind::Unresolved { star_indices },
            }],
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StellarHierarchyNode {
    pub parent_index: Option<usize>,
    pub kind: StellarHierarchyNodeKind,
}

#[derive(Debug, Clone, PartialEq)]
pub enum StellarHierarchyNodeKind {
    Star { star_index: usize },
    Binary(StellarHierarchyBinary),
    Unresolved { star_indices: Vec<usize> },
}

#[derive(Debug, Clone, PartialEq)]
pub struct StellarHierarchyBinary {
    pub primary_child_index: usize,
    pub secondary_child_index: usize,
    pub semi_major_axis_au: f64,
    pub eccentricity: f64,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct HostedSystemRef {
    pub hosted_system_index: usize,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct StarRef {
    pub hosted_system_index: usize,
    pub star_index: usize,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct PlanetRef {
    pub hosted_system_index: usize,
    pub region_index: usize,
    pub body_index: usize,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum HierarchyHostRef {
    Star {
        node_index: usize,
        global_star_index: usize,
    },
    Binary {
        node_index: usize,
        primary_child_index: usize,
        secondary_child_index: usize,
    },
}

impl HierarchyHostRef {
    pub fn node_index(&self) -> usize {
        match self {
            HierarchyHostRef::Star { node_index, .. } | HierarchyHostRef::Binary { node_index, .. } => *node_index,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum PlanetHost {
    Circumstellar {
        star_ref: StarRef,
    },
    Circumbinary {
        hosted_system_index: usize,
        primary_star_index: usize,
        secondary_star_index: usize,
    },
}

impl PlanetHost {
    pub fn hosted_system_index(&self) -> usize {
        match self {
            PlanetHost::Circumstellar { star_ref } => star_ref.hosted_system_index,
            PlanetHost::Circumbinary {
                hosted_system_index, ..
            } => *hosted_system_index,
        }
    }

    pub fn display_order_key(&self) -> (usize, u8, usize) {
        match self {
            PlanetHost::Circumstellar {
                star_ref:
                    StarRef {
                        hosted_system_index,
                        star_index: 0,
                    },
            } => (*hosted_system_index, 0, 0),
            PlanetHost::Circumstellar { star_ref } => (star_ref.hosted_system_index, 1, star_ref.star_index),
            PlanetHost::Circumbinary {
                hosted_system_index, ..
            } => (*hosted_system_index, 2, 0),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StableRegion {
    pub host: PlanetHost,
    pub hierarchy_host: Option<HierarchyHostRef>,
    pub inner_au: f64,
    pub outer_au: f64,
    pub stability_inner_au: f64,
    pub stability_outer_au: f64,
    pub truncation_inner_au: f64,
    pub truncation_outer_au: f64,
}

impl StableRegion {
    pub fn hierarchy_node_index(&self) -> Option<usize> {
        self.hierarchy_host.map(|host| host.node_index())
    }
}

#[derive(Debug, Clone)]
pub struct PlanetaryRegion {
    pub stable_region: StableRegion,
    pub host_mass_solar: f64,
    pub host_luminosity_solar: f64,
    pub thermal_profile: RegionThermalProfile,
    pub bodies: Vec<Body>,
}

impl PlanetaryRegion {
    pub fn new(
        stable_region: StableRegion,
        host_mass_solar: f64,
        host_luminosity_solar: f64,
        bodies: Vec<Body>,
    ) -> Self {
        let thermal_profile = Self::derive_thermal_profile(&stable_region, host_mass_solar, host_luminosity_solar);

        Self {
            stable_region,
            host_mass_solar,
            host_luminosity_solar,
            thermal_profile,
            bodies,
        }
    }

    pub fn derive_thermal_profile(
        _stable_region: &StableRegion,
        _host_mass_solar: f64,
        host_luminosity_solar: f64,
    ) -> RegionThermalProfile {
        RegionThermalProfile::from_host_luminosity(host_luminosity_solar)
    }
}

#[derive(Debug, Clone)]
pub struct HostedSystem {
    pub designation_suffix: Option<String>,
    pub stars: Vec<crate::star::Star>,
    pub topology: StellarTopology,
    pub planetary_regions: Vec<PlanetaryRegion>,
}
