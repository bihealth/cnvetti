use std::fmt;


/// Enum for selecting preset.
#[derive(Clone, Debug, PartialEq)]
pub enum OptionsPreset {
    Wgs,
    WesOffTarget,
    WesOnTarget,
}


impl OptionsPreset {
    /// Parse `CountKind` from `&str`.
    pub fn from_str(s: &str) -> Option<OptionsPreset> {
        match s {
            "Wgs" => Some(OptionsPreset::Wgs),
            "WesOffTarget" => Some(OptionsPreset::WesOffTarget),
            "WesOnTarget" => Some(OptionsPreset::WesOnTarget),
            _ => None,

        }
    }
}


/// Enum for selecting count type.
#[derive(Clone, Debug, PartialEq)]
pub enum CountKind {
    Coverage,
    Alignments,
}


impl CountKind {
    /// Parse `CountKind` from `&str`.
    pub fn from_str(s: &str) -> Option<CountKind> {
        match s {
            "Coverage" => Some(CountKind::Coverage),
            "Alignments" => Some(CountKind::Alignments),
            _ => None,
        }
    }
}


impl fmt::Display for CountKind {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}
