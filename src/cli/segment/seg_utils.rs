/// Various utility code for segmentation.
use cli::segment::math::median;

// TODO: implement outlier detection
// TODO: implement segmentation refinement

/// Replace `vals` with segment medians, given by `breakpoints`.
pub fn replace_with_segment_medians(values: &[f64], breakpoints: &[usize]) -> Vec<f64> {
    let mut result: Vec<f64> = Vec::with_capacity(values.len());
    if values.is_empty() {
        return result;
    }

    assert!(
        !breakpoints.is_empty(),
        "At least begin and end pos must be given"
    );

    for i in 0..(breakpoints.len() - 1) {
        let x = median(&values[breakpoints[i]..breakpoints[i + 1]]);
        for _pos in breakpoints[i]..breakpoints[i + 1] {
            result.push(x);
        }
    }

    result
}
