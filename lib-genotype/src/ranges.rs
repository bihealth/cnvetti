//! Simple interval arithmetics on ranges.

use std::ops::{Range, Sub};

/// Trait for computing overlaps.
pub trait Overlaps<Rhs = Self> {
    /// Return `true` if `self` overlaps `other`.
    fn overlaps(&self, other: &Rhs) -> bool;
}

impl<Idx: Ord> Overlaps<Range<Idx>> for Range<Idx> {
    fn overlaps(&self, other: &Range<Idx>) -> bool {
        return (self.end <= other.start) && (other.end <= self.start);
    }
}

/// Implementation of a simple set of `Range` objects.
#[derive(Debug)]
pub struct RangeSet<Idx> {
    pub entries: Vec<Range<Idx>>,
}

impl<Idx> RangeSet<Idx> {
    pub fn new(entries: Vec<Range<Idx>>) -> Self {
        RangeSet { entries: entries }
    }

    pub fn from_range(range: Range<Idx>) -> Self {
        RangeSet {
            entries: vec![range],
        }
    }
}

fn sub<Idx: Ord + Clone>(lhs: &Range<Idx>, rhs: &Range<Idx>) -> Vec<Range<Idx>> {
    let lhs = lhs.clone();
    let rhs = rhs.clone();
    if lhs.start >= rhs.end || lhs.end <= rhs.start {
        Vec::new()
    } else if rhs.start > lhs.start && rhs.end < lhs.end {
        vec![lhs.start..rhs.start, rhs.end..lhs.end]
    } else if rhs.start < lhs.start {
        vec![rhs.end..lhs.end]
    } else {
        vec![lhs.start..rhs.start]
    }
}

impl<Idx: Ord + Sub + Clone> RangeSet<Idx> {
    pub fn sub(&self, other: &Range<Idx>) -> RangeSet<Idx> {
        let mut result = Vec::new();

        for entry in &self.entries {
            if entry.overlaps(other) {
                result.append(&mut sub(entry, other));
            }
        }

        RangeSet::new(result)
    }
}
